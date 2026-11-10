"""Run with python3 -m unittest discover -s tests -v -p 'test_qa_session.py'."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

REPO = Path(__file__).resolve().parents[1]


def _make_pdf(path, pages=1, stamp=b""):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # Distinct content without ghostscript: minimal PDF-like bytes + stamp.
    path.write_bytes(b"%PDF-1.4\n" + stamp + b"\n%%EOF\n" + bytes([pages]))


class QaSessionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.out = Path(self.tmp.name) / "out"
        self.out.mkdir()
        if str(REPO) not in sys.path:
            sys.path.insert(0, str(REPO))
        from common.lib import config as cfg_mod
        from common.lib import hodo_qa_pdf, qa_session

        self.cfg_mod = cfg_mod
        self.qa = qa_session
        self.hodo = hodo_qa_pdf
        self._old_out = cfg_mod.OUTPUT_DIR
        cfg_mod.OUTPUT_DIR = self.out

    def tearDown(self):
        self.cfg_mod.OUTPUT_DIR = self._old_out

    def test_verify_rejects_stale_partial_update(self):
        img = self.out / "img/run00001"
        a = img / "run00001_BHT_HDPRM_Pi.pdf"
        b = img / "run00001_BH2_HDPRM_Pi.pdf"
        _make_pdf(a, stamp=b"old-a")
        _make_pdf(b, stamp=b"old-b")
        session = self.qa.begin(1, "hodo_hdprm", [a, b], meta={"family": "hodo"})
        _make_pdf(a, stamp=b"new-a")  # only one refreshed
        ok, msg, hashes = self.qa.verify_sources(session)
        self.assertFalse(ok)
        self.assertIn("stale PDF", msg)
        self.assertIn("run00001_BH2_HDPRM_Pi.pdf", msg)
        self.assertIn("run00001_BHT_HDPRM_Pi.pdf", hashes)

    def test_merge_session_fails_on_stale(self):
        img = self.out / "img/run00002"
        a = img / "run00002_BHT_HDPRM_Pi.pdf"
        b = img / "run00002_BH2_HDPRM_Pi.pdf"
        _make_pdf(a, stamp=b"old-a")
        _make_pdf(b, stamp=b"old-b")
        session = self.qa.begin(2, "hodo_hdprm", [a, b])
        _make_pdf(a, stamp=b"new-a")
        # BH2 left stale — merge must refuse before pdfunite.
        ok, msg, out = self.hodo.merge(2, "hdprm", session=session)
        self.assertFalse(ok)
        self.assertIn("stale PDF", msg)
        merged = self.out / "img/run00002/run00002_hodo_hdprm_Pi.pdf"
        self.assertFalse(merged.is_file())

    def test_verify_ok_writes_manifest(self):
        img = self.out / "img/run00003"
        a = img / "run00003_BHT_HDPRM_Pi.pdf"
        b = img / "run00003_BH2_HDPRM_Pi.pdf"
        _make_pdf(a, stamp=b"old-a")
        _make_pdf(b, stamp=b"old-b")
        root = self.out / "run00003_Hodo.root"
        root.write_bytes(b"root")
        session = self.qa.begin(
            3, "hodo_hdprm", [a, b], input_path=root,
            meta={"family": "hodo", "mode": "hdprm"},
        )
        _make_pdf(a, stamp=b"new-a")
        _make_pdf(b, stamp=b"new-b")
        ok_v, msg_v, hashes = self.qa.verify_sources(session)
        self.assertTrue(ok_v, msg_v)
        path = self.qa.write_manifest(
            session, source_hashes=hashes, output_pdfs=[a],
        )
        data = json.loads(path.read_text())
        self.assertEqual(data["exec_id"], session["exec_id"])
        self.assertEqual(data["run"], 3)
        self.assertEqual(len(data["sources"]), 2)
        self.assertTrue(data["input_token"].startswith("run00003_Hodo.root:"))
        self.assertTrue(path.name.endswith(".json"))

    def test_dc_try_post_skips_stale(self):
        from common.lib import dc_qa_slack

        pdf = self.out / "img/run00004/run00004_BLC2_TDC_Pi.pdf"
        _make_pdf(pdf, stamp=b"old")
        session = self.qa.begin(4, "dc_t0", [pdf], meta={"family": "dc"})
        # unchanged -> stale
        with mock.patch.object(dc_qa_slack, "load_slack_config",
                               return_value=({"token": "x", "channel_id": "C"}, "ok")):
            ok, msg = dc_qa_slack.try_post_after_dcprm(4, "t0", session=session)
        self.assertFalse(ok)
        self.assertIn("stale PDF", msg)

    def test_post_hodo_cli_requires_pdf_or_manifest(self):
        script = REPO / "calibration/scripts/post_hodo_qa_slack.py"
        result = subprocess.run(
            [sys.executable, str(script), "1", "hdprm"],
            capture_output=True, text=True, timeout=15,
            env={**os.environ, "PYTHONPATH": str(REPO)},
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("--pdf", result.stdout + result.stderr)

    def test_seal_and_resolve_rejects_overwritten_live_path(self):
        """Manifest points at sealed copy; mutating live path must not change post sha."""
        img = self.out / "img/run00005"
        live = img / "run00005_hodo_hdprm_Pi.pdf"
        _make_pdf(live, stamp=b"old-qa")
        session = self.qa.begin(
            5, "hodo_hdprm", [live],
            meta={"family": "hodo", "mode": "hdprm", "particle": "pion"},
        )
        param = self.out / "HodoParam_run00005_Pi"
        param.write_bytes(b"fit-input-param")
        self.qa.attach_params(session, [param], role="fit_input")
        _make_pdf(live, stamp=b"new-qa")
        ok_v, msg_v, hashes = self.qa.verify_sources(session)
        self.assertTrue(ok_v, msg_v)
        sealed = self.qa.seal_outputs(session, [live])
        self.qa.attach_params(session, [param], role="written")
        man_path = self.qa.write_manifest(
            session, source_hashes=hashes, output_pdfs=sealed,
        )
        data = json.loads(man_path.read_text())
        self.assertEqual(data["params"][0]["role"], "fit_input")
        self.assertEqual(data["params"][0]["sha"], self.qa.file_sha256(param))
        self.assertTrue(str(sealed[0]).endswith(f"/{session['exec_id']}/{live.name}"))
        sealed_sha = data["outputs"][0]["sha"]

        # Overwrite the usual img path after sealing — resend must still use seal.
        _make_pdf(live, stamp=b"even-newer")
        ok_r, msg_r, resolved = self.qa.resolve_sealed_output(
            data, run=5, mode="hdprm", particle="pion",
        )
        self.assertTrue(ok_r, msg_r)
        self.assertEqual(self.qa.file_sha256(resolved), sealed_sha)
        self.assertNotEqual(self.qa.file_sha256(live), sealed_sha)

    def test_resolve_rejects_sha_mismatch_on_manifest_path(self):
        img = self.out / "img/run00006"
        live = img / "run00006_hodo_hdprm_Pi.pdf"
        _make_pdf(live, stamp=b"old")
        session = self.qa.begin(
            6, "hodo_hdprm", [live],
            meta={"family": "hodo", "mode": "hdprm", "particle": "pion"},
        )
        _make_pdf(live, stamp=b"new")
        ok_v, _, hashes = self.qa.verify_sources(session)
        self.assertTrue(ok_v)
        # Legacy-style manifest: outputs.path is the mutable live path.
        man_path = self.qa.write_manifest(
            session, source_hashes=hashes, output_pdfs=[live],
        )
        data = json.loads(man_path.read_text())
        old_sha = data["outputs"][0]["sha"]
        _make_pdf(live, stamp=b"tampered")
        data["outputs"][0]["sha"] = old_sha  # keep recorded sha, disk changed
        ok_r, msg_r, _resolved = self.qa.resolve_sealed_output(
            data, run=6, mode="hdprm", particle="pion",
        )
        self.assertFalse(ok_r)
        self.assertIn("sha mismatch", msg_r)

    def test_resolve_rejects_mode_mismatch(self):
        img = self.out / "img/run00007"
        live = img / "run00007_hodo_hdprm_Pi.pdf"
        _make_pdf(live, stamp=b"a")
        session = self.qa.begin(
            7, "hodo_hdprm", [live],
            meta={"family": "hodo", "mode": "hdprm", "particle": "pion"},
        )
        _make_pdf(live, stamp=b"b")
        ok_v, _, hashes = self.qa.verify_sources(session)
        self.assertTrue(ok_v)
        sealed = self.qa.seal_outputs(session, [live])
        man_path = self.qa.write_manifest(
            session, source_hashes=hashes, output_pdfs=sealed,
        )
        data = json.loads(man_path.read_text())
        ok_r, msg_r, _ = self.qa.resolve_sealed_output(
            data, run=7, mode="hdphc", particle="pion",
        )
        self.assertFalse(ok_r)
        self.assertIn("mode mismatch", msg_r)


if __name__ == "__main__":
    unittest.main()
