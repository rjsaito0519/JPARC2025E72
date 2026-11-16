"""Decode must read the same DCGEO that residual calibration updates."""
import importlib.util
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('dcgeo_create_runlist', REPO/'calibration/scripts/create_runlist.py')
create = importlib.util.module_from_spec(spec)
spec.loader.exec_module(create)


class DcgeoRunlistTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        paths = {name: root/name for name in ('PARAM_DIR','OUTPUT_DIR','SCRATCH_DIR','DECODE_DIR','DATA_DIR','ANALYZER_DIR')}
        paths["PARAM_DIR"] = paths["ANALYZER_DIR"] / "param"
        for directory in paths.values():
            directory.mkdir(parents=True, exist_ok=True)
        self.addCleanup(patch.stopall)
        patch.multiple(create.config, **paths, SUB_DIR='e72').start()
        self.template = paths['PARAM_DIR']/'DCGEO/example/DCGeomParam_e72_example'
        self.template.parent.mkdir(parents=True)
        self.template.write_text('1001 BLC2a-U1 0 0 0 0 0 0 0 .2 16 -5 .35\n')
        conf = paths['PARAM_DIR']/'conf/example/analyzer_e72_example1.conf'
        conf.parent.mkdir(parents=True)
        conf.write_text('DCGEO: param/DCGEO/example/DCGeomParam_e72_example\nHDPHC: param/HDPHC/example/HodoPHCParam_e72_example\n')
        self.phc_template=paths['PARAM_DIR']/'HDPHC/example/HodoPHCParam_e72_example'
        self.phc_template.parent.mkdir(parents=True)
        self.phc_template.write_text('initial correction\n')
        for mode in ('bcin','bcout'):
            runlist = paths['ANALYZER_DIR']/f"runmanager/runlist/{create.MODES[mode]['tpl_yml']}"
            runlist.parent.mkdir(parents=True, exist_ok=True)
            runlist.write_text('DEFAULT: {}\nRUN:\n')
        self.target = paths['PARAM_DIR']/'DCGEO/e72/DCGeomParam_run00012_K'

    def setup_run(self, mode='bcin', stage='resi', ref=None):
        create.setup_calibration_run(12, create.MODES[mode], ['K_'+stage], ref, None)
        conf = create.config.PARAM_DIR/f'conf/e72/analyzer_run00012_{mode}_K.conf'
        return create.config.PARAM_DIR / Path(conf.read_text().split()[1]).relative_to('param')

    def test_missing_run_geometry_is_seeded_before_decode(self):
        source = self.setup_run()
        self.assertEqual(self.target, source)
        self.assertEqual(self.template.read_bytes(), source.read_bytes())
        source.write_text('updated by residual calibration\n')
        self.assertEqual('updated by residual calibration\n', self.setup_run('bcout').read_text())
        self.assertEqual('1001', self.template.read_text().split()[0])

    def test_existing_geometry_is_preserved_across_stages(self):
        self.target.parent.mkdir(parents=True)
        self.target.write_text('manual geometry\n')
        for stage in ('tdc','drift','resi'):
            self.assertEqual('manual geometry\n', self.setup_run(stage=stage).read_text())

    def test_explicit_reference_does_not_create_run_geometry(self):
        reference = self.target.with_name('DCGeomParam_run00011_K')
        reference.parent.mkdir(parents=True)
        reference.write_text('reference geometry\n')
        self.assertEqual(reference, self.setup_run(ref=11))
        self.assertFalse(self.target.exists())

    def test_missing_phc_is_seeded_before_first_hodo_decode(self):
        self.setup_run(mode='hodo',stage='hdphc')
        conf=create.config.PARAM_DIR/'conf/e72/analyzer_run00012_hodo_K.conf'
        line=next(line for line in conf.read_text().splitlines() if line.startswith('HDPHC:'))
        target=create.config.PARAM_DIR/Path(line.split()[1]).relative_to('param')
        self.assertEqual('HodoPHCParam_run00012_K',target.name)
        self.assertEqual(self.phc_template.read_bytes(),target.read_bytes())
        target.write_text('updated PHC correction\n')
        self.assertEqual('updated PHC correction\n',target.read_text())
        self.assertEqual('initial correction\n',self.phc_template.read_text())

    def test_existing_phc_is_not_replaced_by_example(self):
        target=create.config.PARAM_DIR/'HDPHC/e72/HodoPHCParam_run00012_K'
        target.parent.mkdir(parents=True)
        target.write_text('manual correction\n')
        self.setup_run(mode='hodo',stage='hdphc')
        self.assertEqual('manual correction\n',target.read_text())

    def test_failed_publication_does_not_leave_partial_geometry(self):
        self.target.parent.mkdir(parents=True)
        with patch.object(create.os, 'replace', side_effect=OSError('publication failed')):
            with self.assertRaises(OSError):
                create.seed_run_parameter(self.target, self.template)
        self.assertFalse(self.target.exists())
        self.assertEqual([], list(self.target.parent.glob(self.target.name+'.*.*')))


if __name__ == '__main__':
    unittest.main()
