"""The per-project settings file and the run files derived from it.

A study directory holds one user-edited file, ``vine-mirus.yaml``: the full
pipeline config (seeded from the packaged ``config/config.yaml``, comments kept)
plus a ``samples:`` path and an ``execution:`` block that only the CLI reads.
Each run derives two files under ``.vine-mirus/``:

* ``config.yaml``  - the pipeline config without ``execution:``, passed to
  Snakemake as ``--configfile`` (so nested keys need no ``--config`` override);
* ``profile/config.yaml`` - the packaged local/slurm profile template merged
  with the ``execution:`` values.
"""

from __future__ import annotations

import copy
import csv
import os
import shutil
from datetime import datetime
from pathlib import Path
from typing import Any

from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap

PROJECT_FILE = "vine-mirus.yaml"
STATE_DIR = ".vine-mirus"
EXECUTORS = ("local", "slurm")
TRUE_STRINGS = {"true", "yes", "on", "1"}  # same spellings as cfg_bool() in the Snakefile


# ── Packaged assets ───────────────────────────────────────────────────────────

def asset_dir(name: str) -> Path:
    """Directory of a packaged asset tree (workflow, config, envs, profiles, test).

    Installed from a wheel the trees sit beside this module (vine_mirus/<name>).
    From a source checkout (editable install) they are the repo-root
    directories, two levels above src/vine_mirus/.
    """
    here = Path(__file__).resolve().parent
    for candidate in (here / name, here.parents[1] / name):
        if candidate.is_dir():
            return candidate
    raise FileNotFoundError(f"packaged '{name}' directory not found next to {here}")


def snakefile() -> Path:
    return asset_dir("workflow") / "Snakefile"


def default_config_path() -> Path:
    return asset_dir("config") / "config.yaml"


def _xdg(var: str, fallback: str) -> Path:
    configured = os.environ.get(var)
    return Path(configured).expanduser() if configured else Path.home() / fallback


def app_data_dir() -> Path:
    return _xdg("XDG_DATA_HOME", ".local/share") / "vine-mirus"


# ── YAML ──────────────────────────────────────────────────────────────────────

def make_yaml() -> YAML:
    yaml = YAML()
    yaml.preserve_quotes = True
    yaml.indent(mapping=2, sequence=4, offset=2)
    yaml.width = 4096
    return yaml


def load_yaml(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return make_yaml().load(handle)


def dump_yaml(data: Any, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as handle:
        make_yaml().dump(data, handle)
    tmp.replace(path)


def get_path(data: Any, keys: tuple[str, ...], default: Any = None) -> Any:
    node = data
    for key in keys:
        if not isinstance(node, dict) or key not in node:
            return default
        node = node[key]
    return node


def set_path(data: Any, keys: tuple[str, ...], value: Any) -> None:
    node = data
    for key in keys[:-1]:
        if key not in node or not isinstance(node[key], dict):
            node[key] = CommentedMap()
        node = node[key]
    node[keys[-1]] = value


# ── Project file ──────────────────────────────────────────────────────────────

def execution_defaults() -> dict[str, Any]:
    data = app_data_dir()
    return {
        "executor": "local",
        "cores": 8,
        "mem_mb": 32000,
        "slurm_account": "",
        "slurm_partition": "main",
        "slurm_mem_mb": 8000,
        "slurm_runtime": 60,
        "jobs": 100,
        "latency_wait": 60,
        "keep_going": True,
        "rerun_incomplete": True,
        "conda_prefix": str(data / "conda-envs"),
        "databases_dir": str(data / "databases"),
    }


_EXECUTION_COMMENTS = {
    "executor": "local (this machine) | slurm (submit each job to the cluster)",
    "cores": "local: total cores shared by all running jobs",
    "mem_mb": "local: total memory budget (~90% of RAM) shared by all running jobs",
    "slurm_account": "slurm: account to charge (blank = cluster default)",
    "slurm_partition": "slurm: default partition",
    "slurm_mem_mb": "slurm: memory for rules that declare none",
    "slurm_runtime": "slurm: walltime in minutes for rules that declare none",
    "jobs": "slurm: maximum jobs queued at once",
    "latency_wait": "seconds to wait for output files on a shared filesystem",
    "keep_going": "keep running independent jobs after one fails",
    "rerun_incomplete": "re-run jobs left incomplete by an interrupted run",
    "conda_prefix": "shared conda env store, so envs are built once, not per project",
    "databases_dir": "shared database store; <project>/databases links here",
}


def project_path(explicit: str | os.PathLike[str] | None = None) -> Path:
    return Path(explicit).expanduser().resolve() if explicit else Path.cwd() / PROJECT_FILE


def _scalar(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return str(value)
    return "''" if value == "" else str(value)


def seed_project(path: Path, samples: str = "samples.tsv", overrides: dict[Any, Any] | None = None) -> Path:
    """Create a project file: header, samples, execution block, then the
    packaged config.yaml verbatim (so every default keeps its comments)."""
    width = max(len(k) for k in _EXECUTION_COMMENTS) + 2
    lines = [
        "# vine-mirus project settings.",
        "#   edit:  vine-mirus config        run:  vine-mirus run -n, then vine-mirus run",
        "",
        f"samples: {samples}    # sample sheet (TSV): sample + r1/r2 and/or contigs",
        "",
        "# ── Execution ──────────────────────────────────────────────────────────────",
        "# Read by the vine-mirus CLI to build the Snakemake profile; not passed to the",
        "# workflow itself.",
        "execution:",
    ]
    for key, value in execution_defaults().items():
        entry = f"  {key}: {_scalar(value)}"
        lines.append(f"{entry:<{width + 14}}  # {_EXECUTION_COMMENTS[key]}")
    lines += ["", "# ══ Pipeline config (packaged defaults) ════════════════════════════════════", ""]
    text = "\n".join(lines) + "\n" + default_config_path().read_text(encoding="utf-8")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    if overrides:
        data = load_yaml(path)
        for keys, value in overrides.items():
            set_path(data, keys if isinstance(keys, tuple) else (keys,), value)
        dump_yaml(data, path)
    return path


def ensure_project(path: Path) -> Path:
    if not path.exists():
        seed_project(path)
    return path


def execution_settings(data: Any) -> dict[str, Any]:
    """Execution block with defaults filled in for any key a user deleted."""
    merged = execution_defaults()
    merged.update({k: v for k, v in (data.get("execution") or {}).items() if v is not None})
    if merged["executor"] not in EXECUTORS:
        raise ValueError(f"execution.executor must be one of {EXECUTORS}, got {merged['executor']!r}")
    return merged


def backup(path: Path) -> Path:
    target_dir = path.parent / STATE_DIR / "backups"
    target_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S-%f")
    target = target_dir / f"{path.stem}.{stamp}{path.suffix}.bak"
    shutil.copy2(path, target)
    return target


# ── Sample sheet ──────────────────────────────────────────────────────────────

def inspect_samples(sheet: str | os.PathLike[str], project_dir: Path) -> tuple[int, list[str]]:
    """(sample count, problems). Paths resolve against the project directory,
    which is Snakemake's working directory for a `vine-mirus run`."""
    path = Path(sheet).expanduser()
    if not path.is_absolute():
        path = project_dir / path
    if not path.is_file():
        return 0, [f"sample sheet not found: {path}"]
    problems: list[str] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        cols = set(reader.fieldnames or [])
        rows = list(reader)
    if "sample" not in cols:
        return 0, [f"{path.name}: no 'sample' column (found {sorted(cols)})"]
    has_reads, has_contigs = {"r1", "r2"} <= cols, "contigs" in cols
    if not (has_reads or has_contigs):
        problems.append(f"{path.name}: needs r1+r2 and/or contigs columns")
    names = [(row.get("sample") or "").strip() for row in rows]
    if len(set(names)) != len(names):
        problems.append("duplicate sample names")
    missing = 0
    for row in rows:
        for col in ("r1", "r2", "contigs"):
            value = (row.get(col) or "").strip()
            if value:
                fp = Path(value).expanduser()
                if not (fp if fp.is_absolute() else project_dir / fp).exists():
                    missing += 1
    if missing:
        problems.append(f"{missing} file path(s) in the sheet do not exist")
    return len(rows), problems


# ── Derived run files ─────────────────────────────────────────────────────────

def _profile_template(executor: str) -> dict[str, Any]:
    safe = YAML(typ="safe")
    with (asset_dir("profiles") / executor / "config.yaml").open(encoding="utf-8") as handle:
        return safe.load(handle) or {}


def build_profile(execution: dict[str, Any]) -> dict[str, Any]:
    profile = _profile_template(execution["executor"])
    profile["use-conda"] = True
    profile["conda-prefix"] = str(Path(execution["conda_prefix"]).expanduser().resolve())
    profile["latency-wait"] = int(execution["latency_wait"])
    profile["keep-going"] = bool(execution["keep_going"])
    profile["rerun-incomplete"] = bool(execution["rerun_incomplete"])
    if execution["executor"] == "local":
        profile["cores"] = int(execution["cores"])
        profile["resources"] = [f"mem_mb={int(execution['mem_mb'])}"]
    else:
        defaults = dict(profile.get("default-resources") or {})
        defaults["slurm_partition"] = str(execution["slurm_partition"])
        defaults["mem_mb"] = int(execution["slurm_mem_mb"])
        defaults["runtime"] = int(execution["slurm_runtime"])
        if str(execution.get("slurm_account") or "").strip():
            defaults["slurm_account"] = str(execution["slurm_account"]).strip()
        profile["default-resources"] = defaults
        profile["jobs"] = int(execution["jobs"])
    return profile


def write_run_files(project_file: Path, data: Any, overrides: dict[str, Any] | None = None) -> tuple[Path, Path]:
    """Write .vine-mirus/config.yaml and .vine-mirus/profile/config.yaml."""
    state = project_file.parent / STATE_DIR
    execution = execution_settings(data)
    pipeline = copy.deepcopy(data)
    pipeline.pop("execution", None)
    for key, value in (overrides or {}).items():
        pipeline[key] = value
    configfile = state / "config.yaml"
    dump_yaml(pipeline, configfile)
    profile_dir = state / "profile"
    profile_dir.mkdir(parents=True, exist_ok=True)
    plain = YAML()
    plain.default_flow_style = False
    with (profile_dir / "config.yaml").open("w", encoding="utf-8") as handle:
        handle.write("# Generated by vine-mirus from the execution: block of vine-mirus.yaml.\n"
                     "# Do not edit; it is rewritten on every run.\n")
        plain.dump(build_profile(execution), handle)
    return configfile, profile_dir


def ensure_databases_link(project_dir: Path, databases_dir: str) -> str:
    """Point <project>/databases at the shared store. Returns a status line.

    Every database path in the workflow is relative to the working directory
    (databases/...), so one symlink shares a single store across projects.
    """
    target = Path(databases_dir).expanduser().resolve()
    link = project_dir / "databases"
    if link.is_symlink():
        if link.resolve() == target:
            return f"databases -> {target}"
        link.unlink()
    elif link.exists():
        return f"WARNING: {link} is a real directory; using it instead of {target}"
    target.mkdir(parents=True, exist_ok=True)
    link.symlink_to(target, target_is_directory=True)
    return f"databases -> {target} (linked)"
