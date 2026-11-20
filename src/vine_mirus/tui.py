"""`vine-mirus config`: a Textual editor for a project's vine-mirus.yaml.

Layout: a brand bar with status buoys, a wave strip, the settings viewport
(one tab at a time) beside a side panel of tabs, and a "Field notes" log.
Widgets are generated from schema.FIELDS. Only changed values are written back,
through ruamel's round-trip mode, so the file's comments survive a save.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from . import project as proj
from .schema import FIELDS, TABS, TOOLS, Field

try:
    from ruamel.yaml.comments import CommentedSeq
    from textual.app import App, ComposeResult
    from textual.binding import Binding
    from textual.containers import Horizontal, Vertical, VerticalScroll
    from textual.theme import Theme
    from textual.widgets import Button, ContentSwitcher, Footer, Input, Label, RichLog, Select, Static, Switch
except ModuleNotFoundError as exc:  # pragma: no cover - depends on the user's environment
    MISSING_TUI_DEPENDENCY: str | None = exc.name or "a required dependency"
else:
    MISSING_TUI_DEPENDENCY = None


# ── Themes ────────────────────────────────────────────────────────────────────
# Superior: deep cold lake water, teal shoals, walleye gold, Pictured Rocks
# sandstone, pine and birch. Dunes: Sleeping Bear sand and Lake Michigan blue.

THEME_SPECS = {
    "superior": dict(
        name="vine-superior", dark=True,
        primary="#2fa4a9", secondary="#1d6f8c", accent="#d9a520",
        foreground="#e8efe9", background="#0b2233", surface="#123247", panel="#081824",
        success="#5aa36b", warning="#c9875a", error="#d0503f",
        variables={
            "lake-deep": "#081824", "lake": "#123247", "shoal": "#1f4d63", "foam": "#7fc8c4",
            "walleye": "#d9a520", "sandstone": "#c9875a", "birch": "#e8efe9", "reed": "#8fa89a",
            "notes-bg": "#0e2a3b", "notes-ink": "#d8e4dc",
        },
    ),
    "dunes": dict(
        name="vine-dunes", dark=False,
        primary="#1f6f9a", secondary="#2f8f9d", accent="#b5651d",
        foreground="#1d2b33", background="#f3ead7", surface="#fbf6ea", panel="#e6d8b8",
        success="#3f7d4e", warning="#b5651d", error="#b03a2e",
        variables={
            "lake-deep": "#e6d8b8", "lake": "#fbf6ea", "shoal": "#d9c9a3", "foam": "#1f6f9a",
            "walleye": "#9a6a12", "sandstone": "#b5651d", "birch": "#1d2b33", "reed": "#5d6b5f",
            "notes-bg": "#fffaf0", "notes-ink": "#2a3a33",
        },
    ),
}

NOTE_COLORS = {
    "superior": {"info": "#d8e4dc", "ok": "#7fd18f", "warn": "#e0a46e", "error": "#ff7a66", "fish": "#7fc8c4"},
    "dunes": {"info": "#2a3a33", "ok": "#2f6b3d", "warn": "#9a5214", "error": "#b03a2e", "fish": "#1f6f9a"},
}

WAVES = "  ~  ≈  ~~   ><>    ≈ ~  ~~   ≈   <><  ~   ≈≈  ~    ~~ ≈   ><(((°>   ~  ≈  ~~    ≈ ~   <><   ~~  ≈  ~ " * 3


# ── Value conversion ──────────────────────────────────────────────────────────

def to_text(fld: Field, value: Any) -> str:
    if value is None:
        return ""
    if fld.kind == "list":
        return ", ".join(str(v) for v in (value or []))
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def _number(text: str) -> int | float:
    return int(text) if text.lstrip("+-").isdigit() else float(text)


def parse_value(fld: Field, raw: Any) -> Any:
    """Widget value -> YAML value, or raise ValueError with a readable reason."""
    if fld.kind == "bool":
        return bool(raw)
    if fld.kind == "choice":
        if raw not in fld.choices:
            raise ValueError(f"choose one of {', '.join(fld.choices)}")
        return raw
    text = str(raw).strip()
    if fld.kind in ("int", "float", "evalue"):
        if not text:
            raise ValueError("a number is required")
        try:
            number = int(text) if fld.kind == "int" else _number(text) if fld.kind == "float" else float(text)
        except ValueError:
            raise ValueError(f"'{text}' is not {'a whole number' if fld.kind == 'int' else 'a number'}") from None
        if fld.lo is not None and number < fld.lo:
            raise ValueError(f"must be ≥ {fld.lo:g}")
        if fld.hi is not None and number > fld.hi:
            raise ValueError(f"must be ≤ {fld.hi:g}")
        return text if fld.kind == "evalue" else number
    if fld.kind == "list":
        items = [item.strip() for item in text.split(",") if item.strip()]
        bad = [i for i in items if fld.choices and i not in fld.choices]
        if bad:
            raise ValueError(f"unknown item(s) {', '.join(bad)}; allowed: {', '.join(fld.choices)}")
        return items
    if fld.kind == "words":
        words = text.split()
        bad = [w for w in words if fld.choices and w not in fld.choices]
        if bad or not words:
            raise ValueError(f"space-separated, from: {' '.join(fld.choices)}")
        return " ".join(words)
    return text  # str / path


def same(a: Any, b: Any) -> bool:
    if isinstance(a, (list, tuple)) or isinstance(b, (list, tuple)):
        return list(a or []) == list(b or [])
    if isinstance(a, bool) or isinstance(b, bool):
        return a is b or a == b
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return float(a) == float(b)
    return str(a if a is not None else "") == str(b if b is not None else "")


def replace_seq(old: Any, new: list) -> None:
    """Replace a ruamel sequence's items in place, keeping its style and the
    comment block that trails its last item (ruamel files that comment under
    the last index, so a plain slice assignment would drop it)."""
    # ruamel does not mark an empty `[]` as flow style, so treat empty as flow:
    # writing block items under it would land them after any following comments.
    flow = not len(old) or bool(old.fa.flow_style())
    ca = getattr(old, "ca", None)
    trailing = ca.items.pop(len(old) - 1, None) if ca is not None and len(old) else None
    if ca is not None:
        ca.items.clear()
    old[:] = new
    if trailing is not None and len(new):
        ca.items[len(new) - 1] = trailing
    if flow or not new:
        old.fa.set_flow_style()


def require_tui_dependencies() -> None:
    if MISSING_TUI_DEPENDENCY:
        raise SystemExit(f"`vine-mirus config` needs {MISSING_TUI_DEPENDENCY}. "
                         "Install with: pip install 'textual>=6,<7' 'ruamel.yaml>=0.18,<0.19'")


if MISSING_TUI_DEPENDENCY is None:

    class FieldBlock(Vertical):
        """Label + control on one row, the description underneath."""

        def __init__(self, fld: Field, control: Any) -> None:
            classes = "field" + (f" show-{fld.show}" if fld.show else "")
            super().__init__(classes=classes)
            self.fld = fld
            self.control = control

        def compose(self) -> ComposeResult:
            with Horizontal(classes="field-row"):
                yield Label(self.fld.label, classes="field-label")
                yield self.control
            hint = self.fld.help
            if self.fld.kind == "list" and self.fld.choices:
                hint += f"  [{', '.join(self.fld.choices)}]"
            yield Static(hint, classes="field-help", markup=False)

    class Section(Vertical):
        def __init__(self, title: str, blocks: list[FieldBlock], show: str | None) -> None:
            super().__init__(classes="section" + (f" show-{show}" if show else ""))
            self.border_title = title
            self.blocks = blocks

        def compose(self) -> ComposeResult:
            yield from self.blocks

    class ConfigEditorApp(App[int]):
        CSS_PATH = "tui.tcss"
        TITLE = "vine-mirus config"
        BINDINGS = [
            Binding("ctrl+s", "save", "Save", priority=True),
            Binding("ctrl+r", "reload", "Reload", priority=True),
            Binding("ctrl+q", "quit_app", "Quit", priority=True),
            Binding("escape", "quit_app", "Quit", show=False),
            *[Binding(str(i + 1), f"show_tab('{key}')", title, show=False) for i, (key, _g, title) in enumerate(TABS)],
        ]

        def __init__(self, project_file: Path, theme_name: str = "superior") -> None:
            super().__init__()
            self.project_file = project_file
            self.project_dir = project_file.parent
            self.theme_key = theme_name
            self.data: Any = None
            self.dirty = False
            self.quit_armed = False
            self._loading = True
            # Themes before the stylesheet is parsed: tui.tcss uses their variables.
            for spec in THEME_SPECS.values():
                self.register_theme(Theme(**spec))
            self.theme = THEME_SPECS[theme_name]["name"]
            self._load()

        # ── data ──────────────────────────────────────────────────────────────

        def _load(self) -> None:
            self.data = proj.load_yaml(self.project_file)

        def current(self, fld: Field) -> Any:
            value = proj.get_path(self.data, fld.path)
            if value is None and fld.path[0] == "execution":
                value = proj.execution_defaults().get(fld.path[1])
            return value

        def make_control(self, fld: Field) -> Any:
            value = self.current(fld)
            wid = fld.widget_id
            if fld.kind == "bool":
                return Switch(value=bool(value), id=wid)
            if fld.kind == "choice":
                options = list(fld.choices)
                if value not in options and value is not None:
                    options.append(str(value))  # surfaced, then rejected on save
                return Select([(o, o) for o in options], value=value if value in options else Select.BLANK,
                              allow_blank=value not in options, id=wid)
            placeholder = {"path": "path", "list": "comma, separated", "words": "space separated",
                           "evalue": "e.g. 1e-5"}.get(fld.kind, "")
            return Input(value=to_text(fld, value), placeholder=placeholder, id=wid)

        # ── layout ────────────────────────────────────────────────────────────

        def compose(self) -> ComposeResult:
            with Horizontal(id="topbar"):
                yield Static("><(((°>  [b]vine-mirus[/b]\n[i]Great Lakes virome survey[/i]", id="brand")
                yield Static(id="buoy-samples", classes="buoy")
                yield Static(id="buoy-tools", classes="buoy")
                yield Static(id="buoy-mode", classes="buoy")
                yield Static(id="buoy-executor", classes="buoy")
            yield Static(WAVES, id="waves", markup=False)
            with Horizontal(id="client"):
                with ContentSwitcher(initial=f"tab-{TABS[0][0]}", id="viewport"):
                    for key, glyph, title in TABS:
                        with VerticalScroll(id=f"tab-{key}", classes="tab-body"):
                            yield from self.tab_sections(key)
                with Vertical(id="sidepanel"):
                    yield Static("Chart room", id="side-title")
                    for i, (key, glyph, title) in enumerate(TABS):
                        yield Button(f"{i + 1} {glyph} {title}", id=f"side-{key}", classes="side-tab")
                    yield Static(self.side_info(), id="side-info")
                    yield Button("Save  ^s", id="save", classes="action", variant="primary")
                    yield Button("Reload  ^r", id="reload", classes="action")
            notes = RichLog(id="notes", markup=True, wrap=True, max_lines=400)
            notes.border_title = "Field notes"
            yield notes
            yield Footer()

        def tab_sections(self, tab: str) -> list[Any]:
            groups: dict[str, list[Field]] = {}
            for fld in FIELDS:
                if fld.tab == tab:
                    groups.setdefault(fld.group, []).append(fld)
            sections = []
            for title, flds in groups.items():
                shows = {f.show for f in flds}
                section_show = shows.pop() if len(shows) == 1 else None
                blocks = [FieldBlock(f, self.make_control(f)) for f in flds]
                sections.append(Section(title, blocks, section_show))
            return sections

        def side_info(self) -> str:
            return (f"[b]Project file[/b]\n{self.project_file}\n\n"
                    "[b]Then run[/b]\nvine-mirus run -n\nvine-mirus run")

        # ── lifecycle ─────────────────────────────────────────────────────────

        def on_mount(self) -> None:
            self.query_one("#viewport").border_title = f"{TABS[0][1]}  {TABS[0][2]}"
            self.mark_tab(TABS[0][0])
            self.apply_executor(self.widget_value(self.field("execution.executor")))
            self.refresh_status()
            self.note(f"><> Casting off: {self.project_file}", "fish")
            self.note("Tabs: side panel or keys 1-7 · Save ^s · Reload ^r · Quit ^q/esc", "info")
            self.call_after_refresh(self._done_loading)

        def _done_loading(self) -> None:
            self._loading = False

        # ── helpers ───────────────────────────────────────────────────────────

        @staticmethod
        def field(dotted: str) -> Field:
            return next(f for f in FIELDS if f.dotted == dotted)

        def widget_value(self, fld: Field) -> Any:
            widget = self.query_one(f"#{fld.widget_id}")
            value = widget.value
            return None if value is Select.BLANK else value

        def note(self, message: str, level: str = "info") -> None:
            color = NOTE_COLORS[self.theme_key].get(level, NOTE_COLORS[self.theme_key]["info"])
            self.query_one("#notes", RichLog).write(f"[{color}]{message}[/]")

        def mark_tab(self, tab: str) -> None:
            for key, _glyph, _title in TABS:
                self.query_one(f"#side-{key}", Button).set_class(key == tab, "active")

        def apply_executor(self, executor: Any) -> None:
            for other in proj.EXECUTORS:
                for widget in self.query(f".show-{other}"):
                    widget.display = other == executor

        def refresh_status(self) -> None:
            sheet = self.widget_value(self.field("samples")) or ""
            n, problems = proj.inspect_samples(sheet, self.project_dir) if sheet else (0, ["no sample sheet set"])
            samples = self.query_one("#buoy-samples", Static)
            samples.update(f"Samples\n{n}" + (" ⚠" if problems else ""))
            samples.set_class(bool(problems), "buoy-warn")
            on = sum(bool(self.widget_value(self.field(f"virus_identification.tools.{t}"))) for t in TOOLS)
            tools = self.query_one("#buoy-tools", Static)
            tools.update(f"Tools\n{on}/{len(TOOLS)}")
            tools.set_class(on == 0, "buoy-warn")
            test = bool(self.widget_value(self.field("test_mode")))
            mode = self.query_one("#buoy-mode", Static)
            mode.update("Mode\nTEST" if test else "Mode\nPRODUCTION")
            mode.set_class(test, "buoy-test")
            self.query_one("#buoy-executor", Static).update(
                f"Executor\n{self.widget_value(self.field('execution.executor')) or '?'}")

        def form_differs(self) -> bool:
            for fld in FIELDS:
                try:
                    value = parse_value(fld, self.widget_value(fld))
                except ValueError:
                    return True
                if not same(value, self.current(fld)):
                    return True
            return False

        def _changed(self) -> None:
            if self._loading:
                return
            self.dirty = self.form_differs()
            self.quit_armed = False
            self.sub_title = "● unsaved changes" if self.dirty else ""
            self.refresh_status()

        # ── events ────────────────────────────────────────────────────────────

        def on_button_pressed(self, event: Button.Pressed) -> None:
            bid = event.button.id or ""
            if bid.startswith("side-"):
                self.action_show_tab(bid.removeprefix("side-"))
            elif bid == "save":
                self.action_save()
            elif bid == "reload":
                self.action_reload()

        def on_switch_changed(self, _event: Switch.Changed) -> None:
            self._changed()

        def on_input_changed(self, _event: Input.Changed) -> None:
            self._changed()

        def on_select_changed(self, event: Select.Changed) -> None:
            if event.select.id == self.field("execution.executor").widget_id:
                self.apply_executor(event.value)
            self._changed()

        # ── actions ───────────────────────────────────────────────────────────

        def action_show_tab(self, tab: str) -> None:
            self.query_one("#viewport", ContentSwitcher).current = f"tab-{tab}"
            glyph, title = next((g, t) for k, g, t in TABS if k == tab)
            self.query_one("#viewport").border_title = f"{glyph}  {title}"
            self.mark_tab(tab)

        def action_reload(self) -> None:
            self._loading = True
            self._load()
            for fld in FIELDS:
                widget = self.query_one(f"#{fld.widget_id}")
                value = self.current(fld)
                if isinstance(widget, Switch):
                    widget.value = bool(value)
                elif isinstance(widget, Select):
                    widget.value = value if value in fld.choices else Select.BLANK
                else:
                    widget.value = to_text(fld, value)
            self.apply_executor(self.current(self.field("execution.executor")))
            self.dirty = False
            self.sub_title = ""
            self.refresh_status()
            self.note("↻ Reloaded from disk.", "info")
            self.call_after_refresh(self._done_loading)

        def action_quit_app(self) -> None:
            if self.dirty and not self.quit_armed:
                self.quit_armed = True
                self.note("Unsaved changes. ^s to save, or quit again to drop them.", "warn")
                return
            self.exit(0)

        def collect(self) -> tuple[dict[tuple[str, ...], Any], list[str]]:
            values, errors = {}, []
            for fld in FIELDS:
                try:
                    values[fld.path] = parse_value(fld, self.widget_value(fld))
                except ValueError as exc:
                    errors.append(f"{fld.label} ({fld.dotted}): {exc}")
            return values, errors

        def cross_check(self, v: dict[tuple[str, ...], Any]) -> tuple[list[str], list[str]]:
            """(errors that block the save, warnings that do not)."""
            g = lambda dotted: v.get(tuple(dotted.split(".")))  # noqa: E731
            errors, warnings = [], []
            tools_on = [t for t in TOOLS if g(f"virus_identification.tools.{t}")]
            if not tools_on:
                errors.append("Enable at least one detection tool.")
            t1, t2, t3 = (g(f"collation.confidence_tiers.tier_{i}_min_tools") for i in (1, 2, 3))
            if None not in (t1, t2, t3) and not (t1 >= t2 >= t3):
                errors.append(f"Tier thresholds must run tier 1 ≥ tier 2 ≥ tier 3 (got {t1}, {t2}, {t3}).")
            if g("preprocessing.host_removal.enabled") and not g("preprocessing.host_removal.reference_genomes"):
                errors.append("Host removal is on but no host reference FASTA is set.")
            if g("execution.executor") == "slurm" and not g("execution.slurm_partition"):
                errors.append("Slurm needs a default partition.")
            if "vibrant" not in tools_on:
                warnings.append("VIBRANT is off: the phage-exclusion rules depend on its calls.")
            scored = [t for t in ("genomad", "virsorter2", "vibrant", "pprmeta") if t in tools_on]
            if g("collation.candidate_rule") != "union" and sum(g(f"tool_weights.{t}") or 0 for t in scored) <= 0:
                errors.append("The score/both rule needs a positive weight on at least one enabled scoring tool.")
            sheet = g("samples")
            _n, problems = proj.inspect_samples(sheet, self.project_dir) if sheet else (0, ["no sample sheet set"])
            warnings += [f"Sample sheet: {p}" for p in problems]
            if not g("test_mode"):
                if g("virus_identification.tools.diamond_nr") and not g("download_nr"):
                    nr = Path(str(g("diamond_nr_db"))).expanduser()
                    nr = nr if nr.is_absolute() else self.project_dir / nr
                    if not nr.exists():
                        warnings.append(f"DIAMOND vs nr is on, download_nr is off and {nr} does not exist: "
                                        "the run will stop and explain how to provide nr.")
                db_dir = Path(str(g("execution.databases_dir"))).expanduser()
                if not db_dir.is_dir():
                    warnings.append(f"Database store {db_dir} does not exist yet: run `vine-mirus databases` first.")
            return errors, warnings

        def action_save(self) -> None:
            values, errors = self.collect()
            if not errors:
                more_errors, warnings = self.cross_check(values)
                errors += more_errors
            if errors:
                self.note("✗ Not saved - fix these first:", "error")
                for err in errors:
                    self.note(f"   • {err}", "error")
                return
            changes = []
            for fld in FIELDS:
                new, old = values[fld.path], proj.get_path(self.data, fld.path)
                if old is None and fld.path[0] == "execution":
                    old = proj.execution_defaults().get(fld.path[1])
                if same(new, old) and proj.get_path(self.data, fld.path) is not None:
                    continue
                old_text = to_text(fld, old)
                if isinstance(old, list) and isinstance(new, list):
                    replace_seq(old, new)
                elif isinstance(new, list):
                    proj.set_path(self.data, fld.path, CommentedSeq(new))
                else:
                    proj.set_path(self.data, fld.path, new)
                changes.append(f"{fld.dotted}: {old_text or '(blank)'} → {to_text(fld, new) or '(blank)'}")
            if not changes:
                self.note("Nothing changed; file left as it was.", "info")
                self.dirty = False
                self.sub_title = ""
                return
            bak = proj.backup(self.project_file)
            proj.dump_yaml(self.data, self.project_file)
            self.dirty = False
            self.quit_armed = False
            self.sub_title = ""
            self.note(f"✓ Saved {len(changes)} change(s) to {self.project_file.name} (backup: {bak.name})", "ok")
            for change in changes:
                self.note(f"   {change}", "ok")
            for warn in warnings:
                self.note(f"⚠ {warn}", "warn")


def launch_config_tui(project_file: Path, theme_name: str = "superior") -> int:
    require_tui_dependencies()
    result = ConfigEditorApp(project_file, theme_name=theme_name).run()
    return int(result or 0)
