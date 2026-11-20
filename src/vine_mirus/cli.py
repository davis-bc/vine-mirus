"""vine-mirus command-line interface: config, run, test, databases."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from . import __version__
from . import project as proj

THEMES = ("superior", "dunes")
BANNER = "  ~~~ ><(((°>  vine-mirus · Great Lakes virome survey  ~~~"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="vine-mirus",
        description=(
            f"{BANNER}\n\n"
            "Viral discovery in metagenomes, with Snakemake running under the hood.\n"
            "Settings live in ./vine-mirus.yaml in each project directory.\n\n"
            "  vine-mirus test        mocked end-to-end run, no databases (start here)\n"
            "  vine-mirus config      edit this project's settings in the terminal UI\n"
            "  vine-mirus databases   download/build the reference databases\n"
            "  vine-mirus run         run the workflow"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", metavar="command")

    p = sub.add_parser("config", help="Edit this project's settings in the terminal UI.",
                       formatter_class=argparse.RawTextHelpFormatter,
                       description="Create (if missing) and edit ./vine-mirus.yaml.\n\n"
                                   "Examples:\n  vine-mirus config\n  vine-mirus config --theme dunes\n"
                                   "  vine-mirus config --file /data/study1/vine-mirus.yaml")
    p.add_argument("--file", type=Path, help="Project file to edit (default: ./vine-mirus.yaml).")
    p.add_argument("--theme", choices=THEMES, default="superior", help="TUI theme (default: superior).")
    p.add_argument("--print-path", action="store_true", help="Print the project file path and exit.")
    p.set_defaults(handler=handle_config)

    run_help = ("Unrecognised arguments go straight to Snakemake.\n\nExamples:\n"
                "  vine-mirus run -n                 dry run (do this first)\n"
                "  vine-mirus run\n"
                "  vine-mirus run results/evidence_table.tsv --forcerun evidence_score")
    p = sub.add_parser("run", help="Run the workflow with this project's settings.",
                       formatter_class=argparse.RawTextHelpFormatter,
                       description="Run the packaged workflow for the project in the current directory.\n" + run_help)
    p.add_argument("--file", type=Path, help="Project file (default: ./vine-mirus.yaml).")
    p.set_defaults(handler=handle_run)

    p = sub.add_parser("test", help="Mocked end-to-end run on synthetic reads (no databases).",
                       formatter_class=argparse.RawTextHelpFormatter,
                       description="Generate synthetic reads and run the whole workflow in test mode.\n"
                                   "Every tool is mocked; a pass shows the joins and logic work.\n\n"
                                   "Examples:\n  vine-mirus test\n  vine-mirus test --dir /scratch/vm-test --cores 8")
    p.add_argument("--dir", type=Path, default=Path("vine-mirus-test"),
                   help="Test project directory (default: ./vine-mirus-test).")
    p.add_argument("--cores", type=int, default=4, help="Cores to use (default: 4).")
    p.set_defaults(handler=handle_test)

    p = sub.add_parser("databases", help="Download/build the reference databases (hours).",
                       formatter_class=argparse.RawTextHelpFormatter,
                       description="Build every database the enabled tools need, into the shared store\n"
                                   "(execution.databases_dir). nr is only fetched with download_nr: true.\n" + run_help.split("\n\n")[0])
    p.add_argument("--file", type=Path, help="Project file (default: ./vine-mirus.yaml).")
    p.set_defaults(handler=handle_databases)
    return parser


def handle_config(args: argparse.Namespace, _extra: list[str]) -> int:
    path = proj.ensure_project(proj.project_path(args.file))
    if args.print_path:
        print(path)
        return 0
    from .tui import launch_config_tui

    return launch_config_tui(path, theme_name=args.theme)


def handle_run(args: argparse.Namespace, extra: list[str]) -> int:
    from .runtime import run_project

    return run_project(proj.project_path(args.file), extra)


def handle_databases(args: argparse.Namespace, extra: list[str]) -> int:
    from .runtime import run_project

    return run_project(proj.project_path(args.file), extra, targets=["databases"],
                       overrides={"test_mode": False})


def handle_test(args: argparse.Namespace, extra: list[str]) -> int:
    from .runtime import run_project

    test_dir = args.dir.expanduser().resolve()
    data_dir = test_dir / "test_data"
    test_dir.mkdir(parents=True, exist_ok=True)
    if not all((data_dir / f"sample{i}_R{m}.fastq.gz").exists() for i in (1, 2) for m in (1, 2)):
        generator = proj.asset_dir("test") / "generate_test_data.py"
        rc = subprocess.run([sys.executable, "-I", str(generator), "--outdir", str(data_dir)]).returncode
        if rc != 0:
            return rc
    sheet = test_dir / "samples.tsv"
    sheet.write_text("sample\tr1\tr2\n" + "".join(
        f"sample{i}\ttest_data/sample{i}_R1.fastq.gz\ttest_data/sample{i}_R2.fastq.gz\n" for i in (1, 2)))
    project_file = test_dir / proj.PROJECT_FILE
    if not project_file.exists():
        proj.seed_project(project_file, samples="samples.tsv", overrides={
            "test_mode": True,
            ("execution", "executor"): "local",
            ("execution", "cores"): args.cores,
        })
    return run_project(project_file, ["--cores", str(args.cores), *extra], overrides={"test_mode": True})


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args, extra = parser.parse_known_args(sys.argv[1:] if argv is None else argv)
    if not hasattr(args, "handler"):
        parser.print_help()
        defaults = proj.execution_defaults()
        print(f"\nProject file here:  {proj.project_path()}"
              f"{'' if proj.project_path().exists() else '  (not created yet)'}")
        print(f"Shared conda envs:  {defaults['conda_prefix']}")
        print(f"Shared databases:   {defaults['databases_dir']}")
        print(f"Packaged workflow:  {proj.snakefile()}")
        return 0
    if args.command == "config" and extra:
        parser.error(f"unrecognized arguments: {' '.join(extra)}")
    return args.handler(args, extra)
