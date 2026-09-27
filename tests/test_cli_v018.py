"""Command-line validation and output ownership."""
from contextlib import redirect_stderr, redirect_stdout
import io
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from intraphy.cli import main
from intraphy.commands.parser import build_parser
from intraphy.commands.preflight import validate_arguments, validate_input_paths
from intraphy.commands.session import command_session, reserve_output
from intraphy.preparation.manifests import load_manifest
from intraphy.verification.native_cases import build_native_example


class CommandTests(unittest.TestCase):
    def test_example_is_raw_and_portable(self):
        with tempfile.TemporaryDirectory() as tmp, redirect_stdout(io.StringIO()):
            root=Path(tmp)/'native'
            self.assertEqual(main(['example','--output-dir',str(root)]),0)
            manifest=load_manifest(root/'manifest.tsv')
            self.assertEqual(len(manifest),4)
            self.assertTrue(all(Path(r['genome_fasta']).is_absolute() for r in manifest))
            self.assertFalse((root/'segment_homology.tsv').exists())
            self.assertFalse((root/'structural_site_matrix.tsv').exists())

    def test_example_refuses_nonempty_directory(self):
        with tempfile.TemporaryDirectory() as tmp, redirect_stderr(io.StringIO()):
            Path(tmp,'keep.txt').write_text('keep')
            self.assertEqual(main(['example','--output-dir',tmp]),2)
            self.assertEqual(Path(tmp,'keep.txt').read_text(),'keep')

    def test_nonempty_output_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'out'; root.mkdir(); (root/'keep').write_text('x')
            with self.assertRaises(ValueError):
                reserve_output(SimpleNamespace(output_dir=str(root),force=False))

    def test_force_preserves_previous_output(self):
        with tempfile.TemporaryDirectory() as tmp, redirect_stderr(io.StringIO()):
            root=Path(tmp)/'out'; root.mkdir()
            (root/'execution.json').write_text('{}'); (root/'keep').write_text('original')
            reserve_output(SimpleNamespace(output_dir=str(root),force=True))
            backups=list(Path(tmp).glob('out.previous-*'))
            self.assertEqual(len(backups),1)
            self.assertEqual((backups[0]/'keep').read_text(),'original')
            self.assertFalse((root/'keep').exists())

    def test_force_rejects_unowned_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'out'; root.mkdir(); (root/'other_program').write_text('x')
            with self.assertRaises(ValueError): reserve_output(SimpleNamespace(output_dir=str(root),force=True))

    def test_output_cannot_contain_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'work'; root.mkdir(); source=root/'input'; source.mkdir()
            with self.assertRaises(ValueError):
                reserve_output(SimpleNamespace(output_dir=str(root),input_dir=str(source),force=True))

    def test_output_symlink_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); (root/'real').mkdir(); (root/'link').symlink_to(root/'real',target_is_directory=True)
            with self.assertRaises(ValueError): reserve_output(SimpleNamespace(output_dir=str(root/'link'),force=True))

    def test_existing_lock_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'out'; root.mkdir(); (root/'.intraphy.lock').write_text('locked')
            with self.assertRaises(ValueError): reserve_output(SimpleNamespace(output_dir=str(root),force=True))

    def test_failure_records_status_and_releases_lock(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/'out'
            args=SimpleNamespace(command='run',output_dir=str(out),force=False,quiet=True)
            with patch('intraphy.commands.session.preflight'), patch('intraphy.commands.session.environment_report',return_value={}):
                with self.assertRaisesRegex(RuntimeError,'test failure'):
                    with command_session(args): raise RuntimeError('test failure')
            self.assertEqual(json.loads((out/'execution.json').read_text())['status'],'failed')
            self.assertFalse((out/'.intraphy.lock').exists())
            self.assertIn('test failure',(out/'intraphy.log').read_text())

    def test_finite_fraction_validation(self):
        parser=build_parser()
        for value in ('nan','inf','-0.1','1.1'):
            with self.subTest(value=value):
                args=parser.parse_args(['run','--model','exon-parsimony','--input-dir','in','--output-dir','out','--min-callable-fraction',value])
                with self.assertRaises(ValueError): validate_arguments(args)

    def test_default_scope_is_all(self):
        args=build_parser().parse_args(['run','--input-dir','in','--output-dir','out'])
        self.assertEqual(args.analysis_range,'all')
        self.assertEqual(args.model,'exon-locus-ctmc')

    def test_formal_bootstrap_not_silently_accepted(self):
        args=build_parser().parse_args(['run','--model','exon-parsimony','--input-dir','in','--output-dir','out','--bootstrap-replicates','1'])
        with self.assertRaises(ValueError): validate_arguments(args)
