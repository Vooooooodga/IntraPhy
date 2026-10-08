"""Command execution, explicit output ownership, logs and failure provenance."""
from contextlib import contextmanager
from datetime import datetime, timezone
import json
import logging
from pathlib import Path
import sys

from .environment import environment_report
from .preflight import preflight

GUARDED_COMMANDS = {"build-case", "run", "infer-phylogeny", "visualize", "import-orthofinder",
                    "extract-loci", "normalize-annotation", "explain", "analyze", "fit-exon-rates", "realign-exons",
                    "compare-exon-foreground"}
GUARDED_COMMANDS.update({"infer-exon-repertoires", "fit-exon-repertoire-rates"})
GUARDED_COMMANDS.update({"locus-statistics", "prepare-locus-evidence", "prepare-genomic-evidence",
                         "calibrate-exons"})
OWNERS = {"execution.json", "run_result.json", "case_build_report.tsv", "case_provenance.tsv",
          "visualization_manifest.tsv", "target_manifest.tsv", "figure_manifest.json",
          "locus_history.json", "locus_fit.json", "locus_statistics.json", "locus_model.json",
          "genomic_sites.tsv", "genomic_members.tsv", "dna_observations.tsv", "dna_fit.json", "dna_history.json",
          "intron_families.tsv", "intron_positions.tsv", "intron_observations.tsv", "intron_alignment.tsv",
          "intron_fit.json", "intron_history.json"}
OWNERS.update({"exon_structure_fit.json", "exon_history.json", "exon_configurations.jsonl"})
OWNERS.add("calibration_metadata.json")
OWNERS.update({"exon_foreground_comparison.json", "foreground_profile.tsv",
               "foreground_family_fits.tsv", "foreground_unit_diagnostics.tsv",
               "foreground_branches.tsv", "state_space_diagnostics.jsonl"})


def _write_json(path, data):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")
    temporary.replace(path)


def reserve_output(args):
    path = Path(args.output_dir).expanduser()
    if path.is_symlink():
        raise ValueError("Output directory must not be a symbolic link")
    target = path.resolve()
    if target in {Path('/'), Path.home(), Path.cwd()}:
        raise ValueError("Use a dedicated output directory, not the filesystem root, home or working directory")
    alignment_source = getattr(args, "alignment_evidence_dir", None)
    if alignment_source and target.is_relative_to(Path(alignment_source).resolve()):
        raise ValueError("Output directory must not be inside --alignment-evidence-dir.")
    for field in ("input_dir", "result_dir", "genomic_evidence_dir", "alignment_evidence_dir", "structural_site_matrix", "manifest", "species_tree", "locus_model", "evidence_json",
                  "fasta", "gff", "orthologs", "config", "exon_configurations", "exon_rates", "repertoire_model", "profile_dir", "codon_matrix"):
        value = getattr(args, field, None)
        paths = value if isinstance(value, (list, tuple)) else [value]
        if any(item and Path(item).resolve().is_relative_to(target) for item in paths):
            raise ValueError(f"Output directory contains --{field.replace('_', '-')}; inputs cannot be overwritten")
    if path.exists() and not path.is_dir():
        raise ValueError("Output path exists and is not a directory")
    if path.exists() and any(path.iterdir()):
        if (path / '.intraphy.lock').exists():
            raise ValueError("Output has an active/stale .intraphy.lock; inspect the previous run before removing it")
        if not getattr(args, "force", False):
            raise ValueError("Output directory is nonempty; choose a new directory or use --force to preserve it as a backup")
        if not any((path / name).exists() for name in OWNERS):
            raise ValueError("--force refuses to move a directory without an IntraPhy output manifest")
        stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        backup = path.with_name(path.name + '.previous-' + stamp)
        path.rename(backup)
        print(f"Previous output preserved: {backup}", file=sys.stderr)
    path.mkdir(parents=True, exist_ok=True)
    return path


@contextmanager
def command_session(args):
    preflight(args)
    if args.command not in GUARDED_COMMANDS:
        yield
        return
    directory = reserve_output(args)
    lock = directory / '.intraphy.lock'
    with lock.open('x', encoding='utf-8') as handle:
        handle.write('IntraPhy command in progress\n')
    logger = logging.getLogger('intraphy')
    logger.setLevel(logging.INFO)
    handler = logging.FileHandler(directory / 'intraphy.log', encoding='utf-8')
    handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(message)s'))
    logger.addHandler(handler)
    state = {"command": args.command, "status": "running", "arguments": {k: v for k, v in vars(args).items() if not k.startswith("_")},
             "started_at_utc": datetime.now(timezone.utc).isoformat()}
    try:
        _write_json(directory / 'environment.json', environment_report())
        _write_json(directory / 'execution.json', state)
        logger.info("Starting %s", args.command)
        if not getattr(args, 'quiet', False):
            print(f"IntraPhy: {args.command} -> {directory}", file=sys.stderr)
        if getattr(args, 'analysis_range', 'all') != 'all':
            logger.warning("Coverage subset requested; the threshold has no calibrated biological interpretation")
        if getattr(args, 'command', None) == 'compare-exon-foreground':
            logger.warning("Foreground comparison is conditional composite-likelihood inference; its p-value and confidence-interval calibration are not available")
        elif getattr(args, 'model', 'parsimony') == 'exon-structure-ctmc':
            logger.warning("Exon-structure estimates are conditional on the mapped finite catalogue and use a composite likelihood for linked local units")
        elif getattr(args, 'model', 'parsimony') in {'er-ard', 'foreground', 'exon-ctmc', 'dna-presence-ctmc', 'intron-position-ctmc'}:
            logger.warning("CTMC results are conditional; finite-sample LRT calibration remains unassessed")
        yield
        state['status'] = 'completed'
        logger.info("Completed %s", args.command)
        if not getattr(args, 'quiet', False):
            print(f"Completed. Results: {directory}", file=sys.stderr)
    except BaseException as exc:
        state.update(status='failed', error_type=type(exc).__name__, error=str(exc))
        logger.error("Failed: %s", exc)
        result = directory / 'run_result.json'
        if result.exists():
            data = json.loads(result.read_text())
            data['status'] = 'failed'
            _write_json(result, data)
        raise
    finally:
        state['finished_at_utc'] = datetime.now(timezone.utc).isoformat()
        _write_json(directory / 'execution.json', state)
        logger.removeHandler(handler)
        handler.close()
        lock.unlink(missing_ok=True)
