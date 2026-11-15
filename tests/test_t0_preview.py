"""ROOT Gaussian T0 preview and curve replay, without calibration writes."""
import json
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np
import uproot

REPO=Path(__file__).resolve().parents[1]
BINARY=REPO/'bin/T0_OFFSET_PREVIEW'


@unittest.skipUnless(BINARY.is_file(),'build T0_OFFSET_PREVIEW first')
class T0PreviewTests(unittest.TestCase):
    def setUp(self):
        self.directory=tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.base=Path(self.directory.name)
        edges=np.linspace(-10.,10.,2001)
        x=(edges[:-1]+edges[1:])/2
        counts=np.round(400.*np.exp(-.5*((x+.7)/.15)**2))
        with uproot.recreate(self.base/'input.root') as root:
            for segment in range(15):root[f'T0_seg{segment}_TimeOffset_Pi']=counts,edges

    def preview(self,mode='fixed',bounds=(-1.4,0.),function='gaussian',automatic=False):
        settings=self.base/'settings.txt'
        settings.write_text(f'bh2-3-t0:Pi {function} {mode} {bounds[0]} {bounds[1]} 400 -.6 .2 0\n')
        command=[str(BINARY),str(self.base/'input.root'),'Pi','3','T0',str(self.base/'result.json')]
        if not automatic:command.append(str(settings))
        result=subprocess.run(command,capture_output=True,text=True,timeout=60)
        return result,json.loads((self.base/'result.json').read_text()) if (self.base/'result.json').exists() else None

    def test_fixed_fit_preserves_bounds_and_negative_offset(self):
        process,result=self.preview()
        self.assertEqual(0,process.returncode,process.stderr[-1000:])
        self.assertTrue(result['success'])
        self.assertEqual(3,len(result['parameters']))
        self.assertAlmostEqual(-.7,result['parameters'][1],places=3)
        self.assertAlmostEqual(.15,result['parameters'][2],places=3)
        self.assertEqual([-1.4,0.],result['range'])

    def test_auto_range_is_inside_requested_search_window(self):
        process,result=self.preview(mode='auto')
        self.assertEqual(0,process.returncode)
        self.assertTrue(result['success'])
        self.assertGreater(result['range'][0],-1.4)
        self.assertLess(result['range'][1],0.)
        self.assertAlmostEqual(-.7,result['parameters'][1],places=3)

    def test_empty_fit_is_not_successful(self):
        process,result=self.preview(bounds=(5.,6.))
        self.assertEqual(0,process.returncode)
        self.assertFalse(result['success'])

    def test_non_gaussian_and_outside_range_are_rejected(self):
        for arguments in ({'function':'landau'},{'bounds':(-11.,0.)}):
            process,_=self.preview(**arguments)
            self.assertNotEqual(0,process.returncode)

    def test_curve_only_replays_gaussian_without_fitting(self):
        settings=self.base/'curve.txt'
        settings.write_text('bh2-3-t0:Pi gaussian fixed -1. 1. 10 0 .2 0\n')
        output=self.base/'curve.json'
        result=subprocess.run([str(BINARY),'--curve',str(settings),'bh2-3-t0:Pi',str(output)],capture_output=True,text=True,timeout=60)
        self.assertEqual(0,result.returncode,result.stderr[-1000:])
        curve=json.loads(output.read_text())['curve']
        self.assertAlmostEqual(10./(.2*np.sqrt(2*np.pi)),max(point[1] for point in curve),places=8)

    def test_default_fit_exposes_saved_gaussian_parameters(self):
        process,result=self.preview(automatic=True)
        self.assertEqual(0,process.returncode,process.stderr[-1000:])
        self.assertAlmostEqual(-.7,result['parameters'][1],places=3)
        self.assertEqual('gaussian',result['function'])
        self.assertTrue(result['curve'])


if __name__=='__main__':unittest.main()
