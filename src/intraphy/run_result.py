"""Explicit, portable ownership of the artifacts produced by a completed run."""
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple


@dataclass(frozen=True)
class RunResult:
    model: str
    analysis_scope: str
    tree_file: str
    artifacts: Tuple[str, ...]
    status: str = "completed"

    def write(self, directory):
        path = Path(directory) / "run_result.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        data = {"schema_version": 1, "status": self.status, "model": self.model,
                "analysis_scope": self.analysis_scope, "tree_file": self.tree_file,
                "artifacts": list(self.artifacts)}
        temporary = path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        temporary.replace(path)


def record_run_result(output_dir, input_dir, model, analysis_scope):
    names = ("node_structural_states.tsv", "branch_structural_events.tsv", "structural_site_summary.tsv", "gene_change_summary.tsv", "minimum_change_history.tsv", "character_catalogue.tsv") if model == "parsimony" else (
        "model_fits.tsv", "model_tests.tsv", "node_state_posteriors.tsv",
        "branch_transition_posteriors.tsv", "structural_changes.tsv")
    from intraphy.storage.tabular import read_tsv
    tree_file = "species_tree.tsv"
    scope_rows = read_tsv(Path(output_dir) / "phylogeny_scope.tsv", optional=True)
    if scope_rows:
        tree_file = scope_rows[0].get("tree_file") or tree_file
    if analysis_scope == "experimental-multicopy":
        names = ("candidate_structural_events.tsv", "event_support_summary.tsv",
                 "model_fit_summary.tsv", "phylogeny_scope.tsv")
    result = RunResult(model, analysis_scope, tree_file,
                       tuple(name for name in names if (Path(output_dir) / name).exists()))
    result.write(output_dir)
    return result


def result_model(result_dir):
    """Read explicit run metadata. Legacy artifacts are accepted only if unambiguous."""
    directory = Path(result_dir)
    manifest = directory / "run_result.json"
    if manifest.exists():
        data = json.loads(manifest.read_text(encoding="utf-8"))
        if data.get("status") != "completed":
            raise ValueError("cannot render an incomplete run")
        model = data.get("model")
        if model not in {"parsimony", "er-ard", "foreground", "experimental-multicopy", "exon-parsimony", "exon-ctmc", "exon-locus-ctmc"}:
            raise ValueError(f"invalid run result model: {model!r}")
        return model
    parameters = directory / "run_parameters.json"
    if parameters.exists():
        data = json.loads(parameters.read_text(encoding="utf-8"))
        model = data.get("model", data.get("model_test"))
        if model in {"parsimony", "er-ard", "foreground"}:
            return model
    possible = [model for model, name in (
        ("parsimony", "branch_structural_events.tsv"),
        ("er-ard", "structural_changes.tsv"),
        ("experimental-multicopy", "event_support_summary.tsv"),
    ) if (directory / name).exists()]
    formal = [model for model in possible if model != "experimental-multicopy"]
    if len(formal) == 1:
        return formal[0]
    if len(possible) > 1:
        raise ValueError("ambiguous legacy results: select or regenerate run_result.json; old artifacts will not be guessed")
    return possible[0] if possible else None


def begin_run(output_dir, model, analysis_scope):
    """Invalidate an older completed result before attempting a new inference."""
    result = RunResult(model, analysis_scope, "species_tree.tsv", (), status="running")
    result.write(output_dir)
    return result
