"""Fail before expensive analysis when arguments or declared inputs are invalid."""
import math
from pathlib import Path
import shutil

from ..storage.tabular import read_tsv
from ..topology import SpeciesTree


def required_tools(args):
    if args.command == "prepare-genomic-evidence":
        return []
    if getattr(args, "model", None) == "exon-locus-ctmc":
        return []
    names = set()
    raw_analyze = (args.command == "analyze" and getattr(args, "model", None) == "dna-presence-ctmc"
                   and not getattr(args, "input_dir", None))
    if args.command in {"build-case", "derive-tables"} or raw_analyze:
        names.add("mafft")  # Family protein alignment and exon-pair alignment.
        names.add(getattr(args, "context_aligner", "minimap2"))
        names.add(getattr(args, "aligner", "mafft"))
    new_model = str(getattr(args, "model", "")).startswith("exon-")
    if (args.command == "run" and not (new_model and getattr(args, "exon_configurations", None))) or (args.command == "infer-phylogeny" and new_model and not getattr(args, "exon_configurations", None)):
        names.update({"mafft", getattr(args, "evidence_aligner", "minimap2")})
        if getattr(args, "evidence_aligner", "minimap2") == "miniprot":
            names.add("minimap2")  # Nucleotide evidence is a separate channel.
    if args.command == "realign-exons":
        names.add(args.cesar)
    if args.command == "normalize-annotation":
        names.add("agat_convert_sp_gxf2gxf.pl")
    return sorted(names - {"internal", "auto"})


def validate_arguments(args):
    if getattr(args, "command", None) == "analyze" and getattr(args, "model", None) is None:
        args.model = "exon-locus-ctmc" if getattr(args, "locus_model", None) else "dna-presence-ctmc"
    if getattr(args, "command", None) == "analyze" and getattr(args, "model", None) == "dna-presence-ctmc":
        if getattr(args, "root_frequency", "stationary") not in {"stationary", "fixed"}:
            raise ValueError("dna-presence-ctmc supports --root-frequency stationary or fixed.")
        if getattr(args, "locus_model", None):
            raise ValueError("--locus-model requires --model exon-locus-ctmc")
        if getattr(args, "exon_configurations", None) or getattr(args, "exon_rates", None):
            raise ValueError("Exon-configuration inputs require an explicit exon-parsimony or exon-ctmc model.")
        legacy_controls = (
            ("observation_view", "evidence"), ("max_locus_bases", 100000),
            ("alignment_timeout", 600), ("exon_identity", .7),
            ("anchor_bases", 12), ("anchor_identity", .8),
        )
        if any(getattr(args, name, default) != default for name, default in legacy_controls) or any(
                getattr(args, name, None) not in (None, [], ()) for name in
                ("max_states", "max_origin_scenarios", "max_observation_scenarios", "origin_root_sensitivity")):
            raise ValueError("Exon-configuration controls do not apply to dna-presence-ctmc.")
        if getattr(args, "input_dir", None) and any(getattr(args, name, None) for name in
                                                      ("fasta", "manifest", "gff", "orthologs")):
            raise ValueError("Use either --input-dir prepared inputs or raw FASTA/GFF inputs.")
        if getattr(args, "branch_length_mode", "supplied") != "supplied":
            raise ValueError("dna-presence-ctmc uses supplied tree branch lengths.")
        if getattr(args, "parameter_mode", "fit") != "fit":
            raise ValueError("--parameter-mode applies only to exon-locus-ctmc.")
        if (getattr(args, "dna_gain_rate", None) is None) != (getattr(args, "dna_loss_rate", None) is None):
            raise ValueError("--dna-gain-rate and --dna-loss-rate must be supplied together.")
        if getattr(args, "root_frequency", "stationary") == "stationary" and getattr(args, "root_presence", .5) != .5:
            raise ValueError("--root-presence is used only with --root-frequency fixed.")
        if getattr(args, "survey_max_dp_cells", 250000) < 1:
            raise ValueError("--survey-max-dp-cells must be positive.")
        for field in ("survey_min_identity", "survey_min_coverage", "dna_gain_rate", "dna_loss_rate"):
            value = getattr(args, field, None)
            if value is not None and (not math.isfinite(value) or value < 0 or
                                      (field.startswith("survey_") and value > 1)):
                raise ValueError(f"--{field.replace('_', '-')} is outside its valid range")
        if getattr(args, "root_frequency", "stationary") == "fixed" and not 0 <= float(getattr(args, "root_presence", .5)) <= 1:
            raise ValueError("--root-presence must lie in [0, 1] for a fixed root distribution.")
    if getattr(args, "command", None) == "prepare-genomic-evidence":
        if getattr(args, "threads", 1) < 1:
            raise ValueError("--threads must be at least 1")
        if getattr(args, "survey_max_dp_cells", 250000) < 1:
            raise ValueError("--survey-max-dp-cells must be positive.")
        for field in ("survey_min_identity", "survey_min_coverage"):
            value = getattr(args, field, None)
            if value is not None and (not math.isfinite(value) or not 0 <= value <= 1):
                raise ValueError(f"--{field.replace('_', '-')} must be finite and in [0, 1]")
        return
    if args.command == "locus-statistics":
        from .locus_statistics import validate_locus_statistics_arguments
        validate_locus_statistics_arguments(args)
        return
    if args.command == "prepare-locus-evidence":
        if not all(getattr(args, field, None) for field in
                   ("locus_model", "species_tree", "evidence_json", "output_dir")):
            raise ValueError("prepare-locus-evidence requires model, tree, evidence JSON, and output directory")
        return
    if getattr(args, "model", None) == "exon-locus-ctmc":
        if not getattr(args, "locus_model", None):
            raise ValueError("The default exon-locus-ctmc model requires --locus-model with a qualified intraphy.exon-locus-model/2 DNA record.")
        if not getattr(args, "species_tree", None):
            raise ValueError("The exon-locus-ctmc model requires --species-tree (rooted Newick or TSV).")
        if getattr(args, "exon_configurations", None) or getattr(args, "exon_rates", None):
            raise ValueError("--exon-configurations and --exon-rates belong to the separate exon-configuration model.")
        if getattr(args, "repertoire_model", None):
            raise ValueError("--repertoire-model belongs to infer-exon-repertoires.")
        if getattr(args, "analysis_scope", "single-copy") != "single-copy":
            raise ValueError("The exon-locus-ctmc model uses labelled copies in one supplied locus and does not fit copy genealogy.")
        if getattr(args, "branch_length_mode", "supplied") != "supplied":
            raise ValueError("The exon-locus-ctmc model uses the supplied tree branch lengths and does not support --branch-length-mode unit.")
        if getattr(args, "bootstrap_replicates", 0) or getattr(args, "stochastic_maps", 0):
            raise ValueError("Bootstrap and stochastic-map options are not implemented for exon-locus-ctmc.")
        if getattr(args, "foreground_branches", None) or getattr(args, "structural_site_matrix", None):
            raise ValueError("Foreground branches and structural-site matrices are not part of exon-locus-ctmc.")
        if getattr(args, "analysis_range", "all") != "all" or getattr(args, "min_callable_fraction", .70) != .70:
            raise ValueError("Coverage filtering options are not implemented for exon-locus-ctmc.")
        if getattr(args, "root_frequency", "estimated") != "estimated" or getattr(args, "root_presence", .5) != .5:
            raise ValueError("Binary root-frequency options are not used; provide the full root distribution in --locus-model.")
        if (getattr(args, "dna_gain_rate", None) is not None or getattr(args, "dna_loss_rate", None) is not None
                or getattr(args, "survey_min_identity", .7) != .7
                or getattr(args, "survey_min_coverage", .8) != .8
                or getattr(args, "survey_max_dp_cells", 250000) != 250000):
            raise ValueError("DNA-presence rates and survey thresholds apply only to --model dna-presence-ctmc.")
        if getattr(args, "ascertainment", "observed-at-least-one") != "observed-at-least-one":
            raise ValueError("Legacy ascertainment options are not used; discovery is declared in --locus-model.")
        if getattr(args, "command", None) in {"run", "infer-phylogeny"} and getattr(args, "input_dir", None):
            raise ValueError("Use either --locus-model/--species-tree or --input-dir prepared inputs; these routes cannot be combined.")
        if getattr(args, "command", None) == "analyze":
            if getattr(args, "fasta", None) or getattr(args, "manifest", None) or getattr(args, "gff", None) or getattr(args, "orthologs", None):
                raise ValueError("--locus-model cannot be combined with raw FASTA/GFF/manifest inputs.")
        for name in ("max_states", "max_origin_scenarios", "max_observation_scenarios"):
            if getattr(args, name, None) is not None:
                raise ValueError(f"--{name.replace('_', '-')} belongs to the exon-configuration model and does not limit locus closure.")
        if getattr(args, "observation_view", "evidence") != "evidence" or getattr(args, "annotation_view", "repertoire") != "repertoire":
            raise ValueError("Legacy annotation-view selectors are not used by exon-locus-ctmc.")
        if getattr(args, "origin_root_sensitivity", []):
            raise ValueError("--origin-root-sensitivity is not implemented for exon-locus-ctmc; supply a root distribution in the model.")
        if getattr(args, "evidence_aligner", "minimap2") != "minimap2":
            raise ValueError("Alignment backend options do not apply to an explicit locus model.")
        if getattr(args, "command", None) == "analyze":
            if getattr(args, "family_id", None):
                raise ValueError("--family-id belongs to raw genomic input selection and is not used with --locus-model.")
            if getattr(args, "flank", 1000) != 1000 or getattr(args, "max_extension", 10000) != 10000:
                raise ValueError("--flank and --max-extension belong to raw genomic preparation.")
        for name, default in (("exon_identity", .7), ("anchor_bases", 12), ("anchor_identity", .8),
                              ("max_locus_bases", 100000), ("alignment_timeout", 600)):
            if getattr(args, name, default) != default:
                raise ValueError(f"--{name.replace('_', '-')} belongs to raw genomic/configuration inference.")
    elif getattr(args, "locus_model", None) or (getattr(args, "parameter_mode", "fit") != "fit"):
        raise ValueError("--locus-model and --parameter-mode apply only to --model exon-locus-ctmc.")
    if args.command == "fit-exon-repertoire-rates":
        if len(args.log_scale_bounds) != 2 or any(not math.isfinite(x) for x in args.log_scale_bounds) or args.log_scale_bounds[0] >= args.log_scale_bounds[1]:
            raise ValueError("--log-scale-bounds must be two ordered values")
        try:
            endpoints = [math.exp(float(x)) for x in args.log_scale_bounds]
        except (OverflowError, ValueError):
            raise ValueError("--log-scale-bounds must produce finite positive scales")
        if any(not math.isfinite(x) or x <= 0 for x in endpoints):
            raise ValueError("--log-scale-bounds must produce finite positive scales")
        if not args.collection_provenance.strip():
            raise ValueError("--collection-provenance is required")
    if args.command == "infer-exon-repertoires" and not isinstance(args.expected_edits, bool):
        raise ValueError("--expected-edits must be Boolean")
    if args.command == "realign-exons" and args.timeout < 1:
        raise ValueError("--timeout must be positive")
    if args.command == "fit-exon-rates":
        if (args.max_states is not None and args.max_states < 1) or args.gene_bootstrap < 0 or args.parametric_bootstrap < 0:
            raise ValueError("State budget must be positive and replicate counts nonnegative")
    if getattr(args, "threads", 1) < 1:
        raise ValueError("--threads must be at least 1")
    for field in ("min_callable_fraction", "root_presence", "identity_threshold", "min_identity", "min_coverage"):
        value = getattr(args, field, None)
        if value is not None and (not math.isfinite(value) or not 0 <= value <= 1):
            raise ValueError(f"--{field.replace('_', '-')} must be finite and in [0, 1]")
    for field in ("flank", "max_extension", "short_context_max_length"):
        if getattr(args, field, 0) < 0:
            raise ValueError(f"--{field.replace('_', '-')} must be nonnegative")
    if str(getattr(args, "model", "")).startswith("exon-") and getattr(args, "model", "") != "exon-locus-ctmc":
        for name in ("max_states", "max_origin_scenarios", "max_observation_scenarios", "max_locus_bases", "alignment_timeout", "anchor_bases"):
            value = getattr(args, name)
            if value is None and name in {"max_observation_scenarios", "max_states", "max_origin_scenarios"}:
                continue
            if value < 1:
                raise ValueError(f"--{name.replace('_', '-')} must be positive")
        for name in ("exon_identity", "anchor_identity"):
            value = getattr(args, name)
            if not math.isfinite(value) or not 0 <= value <= 1:
                raise ValueError(f"--{name.replace('_', '-')} must lie in [0,1]")
        if getattr(args, "analysis_scope", "single-copy") != "single-copy":
            raise ValueError("The exon-configuration model does not implement multicopy genealogy")
        if getattr(args, "structural_site_matrix", None):
            raise ValueError("Structural-site matrices are not part of the exon-configuration model")
        if getattr(args, "analysis_range", "all") != "all":
            raise ValueError("The exon-configuration model does not use a percentage trim rule")
        if getattr(args, "annotation_view", "repertoire") != "repertoire":
            raise ValueError("The exon-configuration model uses explicit local annotation scenarios; --annotation-view is unsupported")
        if getattr(args, "root_frequency", "estimated") != "estimated" or getattr(args, "root_presence", .5) != .5:
            raise ValueError("Binary root-frequency flags do not apply to the exon-configuration model")
        if getattr(args, "ascertainment", "observed-at-least-one") != "observed-at-least-one":
            raise ValueError("Discovery is declared in the exon-configuration catalogue; --ascertainment is unsupported")
        if getattr(args, "evidence_aligner", "minimap2") != "minimap2":
            raise ValueError("The exon-configuration model uses the MAFFT/minimap2 evidence route")
        if getattr(args, "expected_edits", False) and args.model != "exon-ctmc":
            raise ValueError("--expected-edits requires exon-ctmc")
        if args.command == "analyze" and args.exon_configurations:
            raise ValueError("Use infer-phylogeny to read an existing exon-configuration catalogue without rebuilding raw inputs")
        if (args.model == "exon-ctmc") != bool(args.exon_rates):
            raise ValueError("exon-ctmc requires --exon-rates; parsimony must omit it")
    if args.command in {"run", "infer-phylogeny"} and args.model != "exon-locus-ctmc" and args.analysis_scope == "single-copy":
        if args.bootstrap_replicates or args.stochastic_maps:
            raise ValueError("Bootstrap and stochastic-map flags belong to the experimental model; formal values must be 0")
        if (args.model == "foreground") != bool(args.foreground_branches):
            raise ValueError("--model foreground requires --foreground-branches; other models must omit it")


def validate_input_paths(args):
    if args.command == "prepare-genomic-evidence":
        directory = Path(args.input_dir)
        if not directory.is_dir():
            raise FileNotFoundError(f"--input-dir: directory does not exist: {args.input_dir}")
        tree_path = directory / "species_tree.tsv"
        if not tree_path.is_file():
            raise FileNotFoundError("Prepared genomic input lacks species_tree.tsv")
        from ..inference.locus_io import _tree
        tree, _ = _tree(tree_path)
        from .genomic_preflight import validate_prepared_input
        validate_prepared_input(directory, tree)
        return
    if args.command in {"locus-statistics", "prepare-locus-evidence"}:
        fields = (("locus_model", args.locus_model), ("species_tree", args.species_tree))
        if args.command == "prepare-locus-evidence":
            fields += (("evidence_json", args.evidence_json),)
        for field, value in fields:
            if not Path(value).is_file():
                raise FileNotFoundError(f"--{field.replace('_', '-')}: file does not exist: {value}")
        from ..inference.locus_io import load_locus_model, locus_model_from_record
        if args.command == "locus-statistics":
            bundle = load_locus_model(args.locus_model, args.species_tree)
            from .locus_statistics import validate_loaded_locus_statistics
            validate_loaded_locus_statistics(args, bundle)
            args._locus_bundle = bundle
        else:
            from ..inference.locus_codec import strict_json
            from ..inference.locus_evidence import apply_material_evidence
            template = strict_json(args.locus_model)
            rows = strict_json(args.evidence_json)
            if not isinstance(rows, list):
                raise ValueError("--evidence-json must contain a JSON array of evidence rows")
            updated = apply_material_evidence(template, rows)
            args._prepared_locus_bundle = locus_model_from_record(updated, args.species_tree)
            args._prepared_model_record = updated
        return
    for field in ("manifest", "genome", "annotation", "species_tree", "locus_model", "structural_site_matrix", "foreground_branches", "exon_configurations", "exon_rates", "repertoire_model"):
        value = getattr(args, field, None)
        values = value if isinstance(value, (tuple, list)) else [value]
        for item in values:
            if item and not Path(item).is_file():
                raise FileNotFoundError(f"--{field.replace('_', '-')}: file does not exist: {item}")
    if args.command in {"infer-exon-repertoires", "fit-exon-repertoire-rates"}:
        from ..structure.serialization import read_catalogues
        from ..inference.repertoire_run import load_repertoire_bundle
        catalogues = read_catalogues(args.exon_configurations)
        args._repertoire_bundle = load_repertoire_bundle(catalogues, args.species_tree, args.repertoire_model)
    if args.command == "normalize-annotation":
        from ..inputs.resources import expand_files, GFF_SUFFIXES
        expand_files(args.gff, GFF_SUFFIXES)
        if args.config and not Path(args.config).is_file():
            raise FileNotFoundError(f"AGAT configuration does not exist: {args.config}")
        if args.timeout < 1:
            raise ValueError("--timeout must be positive")
    if args.command == "analyze" and getattr(args, "model", None) == "dna-presence-ctmc":
        from ..inference.locus_io import _tree
        if getattr(args, "input_dir", None):
            if getattr(args, "species_tree", None):
                raise ValueError("--species-tree is already fixed in --input-dir prepared inputs; omit the override.")
            directory = Path(args.input_dir)
            if not directory.is_dir():
                raise FileNotFoundError(f"--input-dir: directory does not exist: {args.input_dir}")
            tree_path = directory / "species_tree.tsv"
            if not tree_path.is_file():
                raise FileNotFoundError("Prepared genomic input lacks species_tree.tsv")
            args._prepared_input_dir = directory
            args._species_tree, _ = _tree(tree_path)
            from .genomic_preflight import validate_prepared_input
            validate_prepared_input(directory, args._species_tree)
        else:
            if not getattr(args, "species_tree", None):
                raise ValueError("Raw genomic analyze requires --species-tree.")
            args._species_tree, _ = _tree(args.species_tree)
            if not (getattr(args, "fasta", None) or getattr(args, "manifest", None)):
                raise ValueError("Raw genomic analyze requires --fasta or --manifest.")
            from ..inputs.selection import resolve_inputs
            args._input_selection = resolve_inputs(args)
            if set(args._input_selection.tree_species) != set(args._species_tree.leaf_by_label):
                raise ValueError("Resolved raw input species do not match the supplied tree tips.")
        if getattr(args, "genomic_evidence_dir", None):
            if not getattr(args, "input_dir", None):
                raise ValueError("--genomic-evidence-dir reuse requires --input-dir prepared inputs.")
            evidence = Path(args.genomic_evidence_dir) / "dna_observations.tsv"
            if not evidence.is_file():
                raise FileNotFoundError(f"--genomic-evidence-dir lacks dna_observations.tsv: {evidence}")
            from .genomic_preflight import validate_reusable_evidence
            validate_reusable_evidence(
                args.genomic_evidence_dir, args._prepared_input_dir, args._species_tree,
                min_identity=args.survey_min_identity,
                min_coverage=args.survey_min_coverage,
                max_dp_cells=args.survey_max_dp_cells,
            )
        return
    if getattr(args, "model", None) == "exon-locus-ctmc":
        if not getattr(args, "locus_model", None):
            return
        from ..inference.locus_io import load_locus_model
        args._locus_bundle = load_locus_model(args.locus_model, args.species_tree)
        return
    if args.command in {"build-case", "check", "extract-loci", "analyze"}:
        if args.command == "analyze" and not (getattr(args, "fasta", None) or getattr(args, "manifest", None)):
            raise ValueError("Raw-input analyze requires --fasta or --manifest; the default model requires --locus-model.")
        if args.command == "analyze" and not getattr(args, "species_tree", None):
            raise ValueError("Raw-input analyze requires --species-tree.")
        from ..inputs.selection import resolve_inputs
        args._input_selection = resolve_inputs(args)
    if args.command in {"run", "infer-phylogeny"}:
        if not args.input_dir:
            raise ValueError("Prepared-input run/infer-phylogeny requires --input-dir; the default model requires --locus-model and --species-tree.")
        directory = Path(args.input_dir)
        tree_path = directory / "species_tree.tsv"
        if not tree_path.is_file():
            raise FileNotFoundError(f"Prepared input lacks {tree_path.name}; run build-case with --species-tree")
        SpeciesTree(read_tsv(tree_path))
        if not getattr(args, "exon_configurations", None) and (args.command == "run" or not args.structural_site_matrix):
            required = ("gene_loci.tsv", "gene_loci.fasta", "transcript_paths.tsv") if str(args.model).startswith("exon-") else ("segment_occurrences.tsv", "segment_homology.tsv", "segment_sequences.fasta")
            for name in required:
                if not (directory / name).is_file():
                    raise FileNotFoundError(f"Prepared input lacks {name}")


def preflight(args):
    validate_arguments(args)
    validate_input_paths(args)
    missing = [name for name in required_tools(args) if not shutil.which(name)]
    if missing:
        raise RuntimeError("Required external tools are unavailable on PATH: " + ", ".join(missing)
                           + ". See docs/installation.md and run intraphy inspect-aligners.")
