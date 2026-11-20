"""Run the packaged Snakemake workflow for a project directory."""

from __future__ import annotations

import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

from . import project as proj


def snakemake_command(project_dir: Path, configfile: Path, profile_dir: Path,
                      targets: list[str], extra_args: list[str]) -> list[str]:
    return [
        "snakemake",
        *targets,
        "--snakefile", str(proj.snakefile()),
        "--directory", str(project_dir),
        "--configfile", str(configfile),
        "--profile", str(profile_dir),
        *extra_args,
    ]


def run_project(project_file: Path, extra_args: list[str], targets: list[str] | None = None,
                overrides: dict[str, Any] | None = None) -> int:
    """Derive the run files from the project file and stream a Snakemake run."""
    if shutil.which("snakemake") is None:
        print("snakemake was not found on PATH. Activate the environment that has Snakemake "
              "(see the README's Install section) and try again.", file=sys.stderr)
        return 1
    if not project_file.exists():
        print(f"No project file at {project_file}. Run `vine-mirus config` here first "
              "(or `vine-mirus test` for a mocked end-to-end run).", file=sys.stderr)
        return 1

    project_dir = project_file.parent
    data = proj.load_yaml(project_file)
    try:
        execution = proj.execution_settings(data)
    except ValueError as exc:
        print(f"{project_file.name}: {exc}", file=sys.stderr)
        return 1

    n_samples, problems = proj.inspect_samples(data.get("samples") or "samples.tsv", project_dir)
    if problems:
        for problem in problems:
            print(f"Sample sheet: {problem}", file=sys.stderr)
        return 1

    configfile, profile_dir = proj.write_run_files(project_file, data, overrides)
    db_status = proj.ensure_databases_link(project_dir, execution["databases_dir"])
    test_mode = str((overrides or {}).get("test_mode", data.get("test_mode"))).strip().lower() in proj.TRUE_STRINGS

    log_dir = project_dir / "logs" / "vine-mirus"
    log_dir.mkdir(parents=True, exist_ok=True)
    driver_log = log_dir / f"driver.{datetime.now().strftime('%Y%m%d-%H%M%S')}.log"
    command = snakemake_command(project_dir, configfile, profile_dir, targets or [], extra_args)

    print("  ><(((°>  vine-mirus")
    print(f"  Project:   {project_dir}")
    print(f"  Samples:   {n_samples}  ({data.get('samples')})")
    print(f"  Executor:  {execution['executor']}")
    print(f"  Mode:      {'TEST (mocked tools)' if test_mode else 'production'}")
    print(f"  {db_status}")
    print(f"  Conda:     {execution['conda_prefix']}")
    print(f"  Log:       {driver_log}")
    print("")

    with driver_log.open("w", encoding="utf-8") as handle:
        handle.write(" ".join(command) + "\n\n")
        process = subprocess.Popen(command, cwd=project_dir, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True, bufsize=1)
        assert process.stdout is not None
        try:
            for line in process.stdout:
                sys.stdout.write(line)
                handle.write(line)
        except KeyboardInterrupt:
            process.terminate()
        process.wait()
    if process.returncode == 0 and test_mode:
        print("\nTest mode: the tools were mocked. A pass shows the joins and logic work, "
              "nothing about biology.")
    return process.returncode
