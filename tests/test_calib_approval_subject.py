"""Subject composition for final approval (#37)."""
import json
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from common.lib import calib_approval, config, qa_session  # noqa: E402


class BuildSubjectTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.out = self.root / "out"
        self.param = self.root / "param"
        self.data = self.root / "data"
        self.decode = self.root / "decode"
        self.scratch = self.root / "scratch"
        for path in (self.out, self.param, self.data, self.decode, self.scratch):
            path.mkdir(parents=True, exist_ok=True)
        self._old = {
            "OUTPUT_DIR": config.OUTPUT_DIR,
            "PARAM_DIR": config.PARAM_DIR,
            "DATA_DIR": config.DATA_DIR,
            "DECODE_DIR": config.DECODE_DIR,
            "SCRATCH_DIR": config.SCRATCH_DIR,
            "SUB_DIR": config.SUB_DIR,
        }
        config.OUTPUT_DIR = self.out
        config.PARAM_DIR = self.param
        config.DATA_DIR = self.data
        config.DECODE_DIR = self.decode
        config.SCRATCH_DIR = self.scratch
        config.SUB_DIR = "e72"

    def tearDown(self):
        for key, value in self._old.items():
            setattr(config, key, value)

    def _write_param(self, cat, prefix, suffix="Pi", ext=""):
        pdir = self.param / cat / "e72"
        pdir.mkdir(parents=True, exist_ok=True)
        path = pdir / f"{prefix}_run00001_{suffix}{ext}"
        path.write_bytes(b"v1")
        return path

    def _manifest(self, exec_id, family, mode, particle="pion"):
        man_dir = self.out / "calib/run00001/qa_manifests"
        man_dir.mkdir(parents=True, exist_ok=True)
        path = man_dir / f"{exec_id}.json"
        path.write_text(json.dumps({
            "exec_id": exec_id,
            "run": 1,
            "meta": {"family": family, "mode": mode, "particle": particle},
            "finished_at": "2026-01-01T00:00:00Z",
            "outputs": [{"name": f"{exec_id}.pdf", "sha": "dead"}],
            "params": [],
        }) + "\n", encoding="utf-8")
        return path

    def test_hodo_only_does_not_require_dc(self):
        st = {
            "plan": {
                "calib_particle": "pion",
                "steps": [{"name": "hodo", "opts": {}}],
            },
            "steps": [],
        }
        subject = calib_approval.build_subject(1, st)
        self.assertEqual(subject["families"], ["hodo"])
        self.assertTrue(any(m.startswith("param:HDPRM") for m in subject["missing"]))
        self.assertFalse(any("DCTDC" in m or "DCGEO" in m for m in subject["missing"]))

    def test_dcgeo_change_changes_subject_sha(self):
        self._write_param("HDPRM", "HodoParam")
        self._write_param("HDPHC", "HodoPHCParam")
        self._write_param("DCTDC", "DCTdcParam")
        self._write_param("DCDRFT", "DCDriftParam", ext=".root")
        geo = self._write_param("DCGEO", "DCGeomParam")
        (self.data / "run00001_Hodo.root").write_bytes(b"hodo")
        (self.data / "run00001_BcOut.root").write_bytes(b"bc")
        (self.out / "calib/run00001/hdphc_rounds").mkdir(parents=True)
        (self.out / "calib/run00001/hdphc_rounds/summary.txt").write_text("ok\n")
        self._manifest("h1", "hodo", "hdprm")
        self._manifest("d1", "dc", "t0")
        st = {
            "plan": {
                "calib_particle": "pion",
                "steps": [
                    {"name": "hodo", "opts": {}},
                    {"name": "dc_chamber", "opts": {}},
                ],
            },
            "steps": [],
        }
        s1 = calib_approval.build_subject(1, st)
        self.assertEqual(s1["missing"], [])
        sha1 = calib_approval.subject_sha(s1)
        geo.write_bytes(b"v2-changed")
        s2 = calib_approval.build_subject(1, st)
        sha2 = calib_approval.subject_sha(s2)
        self.assertNotEqual(sha1, sha2)

    def test_particle_filters_manifest(self):
        self._manifest("k1", "hodo", "hdprm", particle="kaon")
        self._manifest("p1", "hodo", "hdprm", particle="pion")
        st = {
            "plan": {
                "calib_particle": "kaon",
                "kaon": True,
                "steps": [{"name": "hodo", "opts": {}}],
            },
            "steps": [],
            "kaon": True,
        }
        subject = calib_approval.build_subject(1, st)
        exec_ids = [m["exec_id"] for m in subject["qa_manifests"]]
        self.assertIn("k1", exec_ids)
        self.assertNotIn("p1", exec_ids)


class FinalQaIndexFilterTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.img = self.root / "img/run00001"
        self.img.mkdir(parents=True)
        self._old = config.OUTPUT_DIR
        config.OUTPUT_DIR = self.root
        from common.lib import final_qa_index
        self.idx = final_qa_index

    def tearDown(self):
        config.OUTPUT_DIR = self._old

    def test_filters_particle_and_family(self):
        (self.img / "run00001_hodo_hdprm_Pi.pdf").write_bytes(b"a")
        (self.img / "run00001_hodo_hdprm_K.pdf").write_bytes(b"b")
        (self.img / "run00001_BLC2_TDC_Pi.pdf").write_bytes(b"c")
        (self.img / "d5_pos_review_run00001_Pi.pdf").write_bytes(b"d")
        rows = self.idx.list_qa_pdfs(1, particle="pion", families=["hodo"])
        names = [r["name"] for r in rows]
        self.assertEqual(names, ["run00001_hodo_hdprm_Pi.pdf"])


if __name__ == "__main__":
    unittest.main()
