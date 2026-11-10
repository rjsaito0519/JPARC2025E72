"""Subject composition / correspondence for final approval (#37)."""
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

    def _write_param(self, cat, prefix, suffix="Pi", ext="", content=b"v1"):
        pdir = self.param / cat / "e72"
        pdir.mkdir(parents=True, exist_ok=True)
        path = pdir / f"{prefix}_run00001_{suffix}{ext}"
        path.write_bytes(content)
        return path

    def _sealed_manifest(
        self,
        exec_id,
        family,
        mode,
        *,
        particle="pion",
        bcin=False,
        all_chambers=False,
        decode_path=None,
        param_paths=None,
        finished_at="2026-01-01T00:00:00Z",
        seal_bytes=b"%PDF sealed",
    ):
        man_dir = self.out / "calib/run00001/qa_manifests"
        man_dir.mkdir(parents=True, exist_ok=True)
        seal_dir = man_dir / exec_id
        seal_dir.mkdir(parents=True, exist_ok=True)
        pdf = seal_dir / f"{exec_id}.pdf"
        pdf.write_bytes(seal_bytes)
        sha = qa_session.file_sha256(pdf)
        meta = {"family": family, "mode": mode, "particle": particle}
        if family == "dc":
            meta["bcin"] = bool(bcin)
            meta["all"] = bool(all_chambers)
        written = []
        for path in param_paths or []:
            path = Path(path)
            written.append({
                "name": path.name,
                "path": str(path),
                "sha": qa_session.file_sha256(path),
                "role": "written",
            })
        tok = None
        if decode_path is not None and Path(decode_path).is_file():
            tok = qa_session.input_token(decode_path)
        path = man_dir / f"{exec_id}.json"
        path.write_text(json.dumps({
            "exec_id": exec_id,
            "run": 1,
            "meta": meta,
            "finished_at": finished_at,
            "input_token": tok,
            "outputs": [{"name": pdf.name, "path": str(pdf), "sha": sha}],
            "params": [],
            "params_written": written,
        }) + "\n", encoding="utf-8")
        return path

    def _seed_hodo_complete(self):
        """Time-ordered hodo: hdprm → decode → t0 → decode → hdphc (final)."""
        hdprm = self._write_param("HDPRM", "HodoParam", content=b"hdprm-v1")
        hdphc = self._write_param("HDPHC", "HodoPHCParam", content=b"hdphc-v1")
        hodo = self.data / "run00001_Hodo.root"
        (self.out / "calib/run00001/hdphc_rounds").mkdir(parents=True)
        (self.out / "calib/run00001/hdphc_rounds/summary.txt").write_text("ok\n")
        hodo.write_bytes(b"hodo-after-hdprm")
        self._sealed_manifest(
            "h_hdprm", "hodo", "hdprm",
            decode_path=hodo, param_paths=[hdprm],
            finished_at="2026-01-01T00:00:00Z",
        )
        hdprm.write_bytes(b"hdprm-after-t0")
        hodo.write_bytes(b"hodo-after-t0")
        self._sealed_manifest(
            "h_t0", "hodo", "t0",
            decode_path=hodo, param_paths=[hdprm],
            finished_at="2026-01-01T00:01:00Z",
        )
        hdphc.write_bytes(b"hdphc-final")
        hodo.write_bytes(b"hodo-final")
        self._sealed_manifest(
            "h_hdphc", "hodo", "hdphc",
            decode_path=hodo, param_paths=[hdphc],
            finished_at="2026-01-01T00:02:00Z",
        )
        return hdprm, hdphc, hodo

    def _seed_dc_chamber(self, *, bcin=False, modes=None):
        """Time-ordered dc chamber: rewrite decode between modes; final=reso."""
        dctdc = self._write_param("DCTDC", "DCTdcParam", content=b"tdc0")
        dcdrft = self._write_param(
            "DCDRFT", "DCDriftParam", ext=".root", content=b"drft0",
        )
        geo = self._write_param("DCGEO", "DCGeomParam", content=b"geo0")
        chamber = "BcIn" if bcin else "BcOut"
        root = self.data / f"run00001_{chamber}.root"
        tag = "in" if bcin else "out"
        mode_params = {
            "t0": [dctdc],
            "drift": [dcdrft],
            "resi": [geo],
            "reso": [geo],
        }
        contents = {
            "t0": (b"tdc-final", b"root-t0"),
            "drift": (b"drft-final", b"root-drift"),
            "resi": (b"geo-resi", b"root-resi"),
            "reso": (b"geo-final", b"root-final"),
        }
        for i, mode in enumerate(modes or ("t0", "drift", "resi", "reso")):
            param_bytes, root_bytes = contents[mode]
            for path in mode_params[mode]:
                path.write_bytes(param_bytes)
            root.write_bytes(root_bytes)
            self._sealed_manifest(
                f"dc_{tag}_{mode}",
                "dc",
                mode,
                bcin=bcin,
                decode_path=root,
                param_paths=mode_params[mode],
                finished_at=f"2026-01-01T00:{i:02d}:00Z",
            )
        return dctdc, dcdrft, geo, root

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

    def test_required_modes_and_matching_set_ok(self):
        self._seed_hodo_complete()
        self._seed_dc_chamber(bcin=False)
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
        subject = calib_approval.build_subject(1, st)
        self.assertEqual(subject["missing"], [])
        req = {(r["family"], r["mode"], r["chamber"]) for r in subject["qa_required"]}
        self.assertIn(("hodo", "hdprm", ""), req)
        self.assertIn(("dc", "reso", "bcout"), req)
        self.assertTrue(all(
            m.get("checks_ok") for m in subject["qa_manifests"] if m.get("required")
        ))

    def test_dcgeo_change_blocks_and_changes_sha(self):
        self._seed_hodo_complete()
        _tdc, _drft, geo, _root = self._seed_dc_chamber(bcin=False)
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
        self.assertNotEqual(sha1, calib_approval.subject_sha(s2))
        self.assertTrue(any(m.startswith("qa_param:") for m in s2["missing"]))

    def test_sealed_pdf_missing_blocks(self):
        self._seed_hodo_complete()
        st = {
            "plan": {
                "calib_particle": "pion",
                "steps": [{"name": "hodo", "opts": {"skip_t0": True}}],
            },
            "steps": [],
        }
        # Corrupt sealed file for hdprm.
        pdf = self.out / "calib/run00001/qa_manifests/h_hdprm/h_hdprm.pdf"
        pdf.unlink()
        subject = calib_approval.build_subject(1, st)
        self.assertIn("qa_sealed:hodo:hdprm", subject["missing"])

    def test_stale_final_input_token_blocks(self):
        self._seed_hodo_complete()
        st = {
            "plan": {
                "calib_particle": "pion",
                "steps": [{"name": "hodo", "opts": {}}],
            },
            "steps": [],
        }
        (self.data / "run00001_Hodo.root").write_bytes(b"hodo-rewritten")
        subject = calib_approval.build_subject(1, st)
        self.assertIn("qa_input:hodo:hdphc", subject["missing"])
        self.assertNotIn("qa_input:hodo:hdprm", subject["missing"])
        self.assertNotIn("qa_input:hodo:t0", subject["missing"])

    def test_time_ordered_hodo_approvable(self):
        self._seed_hodo_complete()
        st = {
            "plan": {
                "calib_particle": "pion",
                "steps": [{"name": "hodo", "opts": {}}],
            },
            "steps": [],
        }
        subject = calib_approval.build_subject(1, st)
        self.assertEqual(subject["missing"], [])
        finals = {(r["family"], r["mode"], r["chamber"]) for r in subject["qa_final"]}
        self.assertEqual(finals, {("hodo", "hdphc", "")})
        by_mode = {m["mode"]: m for m in subject["qa_manifests"]}
        self.assertTrue(by_mode["hdprm"]["checks_ok"])
        self.assertFalse(by_mode["hdprm"]["final"])
        self.assertTrue(by_mode["hdphc"]["final"])
        from common.lib import final_qa_index
        idx = final_qa_index.index_from_subject(1, subject)
        names = {p["name"] for p in idx["pdfs"]}
        self.assertIn("h_hdprm.pdf", names)
        self.assertIn("h_hdphc.pdf", names)

    def test_final_hdphc_param_change_blocks(self):
        _hdprm, hdphc, _hodo = self._seed_hodo_complete()
        st = {
            "plan": {
                "calib_particle": "pion",
                "steps": [{"name": "hodo", "opts": {}}],
            },
            "steps": [],
        }
        self.assertEqual(calib_approval.build_subject(1, st)["missing"], [])
        hdphc.write_bytes(b"hdphc-tampered")
        subject = calib_approval.build_subject(1, st)
        self.assertTrue(any(
            m.startswith("qa_param:hodo:hdphc:") for m in subject["missing"]
        ))

    def test_time_ordered_dc_approvable(self):
        self._seed_dc_chamber(bcin=False)
        st = {
            "plan": {
                "calib_particle": "pion",
                "steps": [{"name": "dc_chamber", "opts": {}}],
            },
            "steps": [],
        }
        subject = calib_approval.build_subject(1, st)
        self.assertEqual(subject["missing"], [])
        finals = {(r["family"], r["mode"], r["chamber"]) for r in subject["qa_final"]}
        self.assertEqual(finals, {("dc", "reso", "bcout")})
        by_mode = {m["mode"]: m for m in subject["qa_manifests"] if m["chamber"] == "bcout"}
        self.assertTrue(by_mode["t0"]["checks_ok"])
        self.assertFalse(by_mode["t0"]["final"])
        self.assertTrue(by_mode["reso"]["final"])

    def test_skip_hdphc_final_is_t0(self):
        hdprm = self._write_param("HDPRM", "HodoParam", content=b"hdprm-v1")
        hodo = self.data / "run00001_Hodo.root"
        hodo.write_bytes(b"hodo-after-hdprm")
        self._sealed_manifest(
            "h_hdprm", "hodo", "hdprm",
            decode_path=hodo, param_paths=[hdprm],
            finished_at="2026-01-01T00:00:00Z",
        )
        hdprm.write_bytes(b"hdprm-after-t0")
        hodo.write_bytes(b"hodo-after-t0")
        self._sealed_manifest(
            "h_t0", "hodo", "t0",
            decode_path=hodo, param_paths=[hdprm],
            finished_at="2026-01-01T00:01:00Z",
        )
        st = {
            "plan": {
                "calib_particle": "pion",
                "steps": [{"name": "hodo", "opts": {"skip_hdphc": True}}],
            },
            "steps": [],
        }
        subject = calib_approval.build_subject(1, st)
        self.assertEqual(subject["missing"], [])
        finals = {(r["family"], r["mode"], r["chamber"]) for r in subject["qa_final"]}
        self.assertEqual(finals, {("hodo", "t0", "")})
        self.assertFalse(any(m.startswith("param:HDPHC") for m in subject["missing"]))
        self.assertFalse(any(p.get("key") == "HDPHC" for p in subject["params"]))

    def test_final_missing_input_token_blocks(self):
        self._seed_hodo_complete()
        man = self.out / "calib/run00001/qa_manifests/h_hdphc.json"
        data = json.loads(man.read_text(encoding="utf-8"))
        data["input_token"] = None
        man.write_text(json.dumps(data) + "\n", encoding="utf-8")
        st = {
            "plan": {
                "calib_particle": "pion",
                "steps": [{"name": "hodo", "opts": {}}],
            },
            "steps": [],
        }
        subject = calib_approval.build_subject(1, st)
        self.assertIn("qa_input:hodo:hdphc", subject["missing"])

    def test_particle_filters_manifest(self):
        self._sealed_manifest("k1", "hodo", "hdprm", particle="kaon")
        self._sealed_manifest("p1", "hodo", "hdprm", particle="pion")
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

    def _both_chamber_plan(self):
        return {
            "plan": {
                "calib_particle": "pion",
                "steps": [
                    {"name": "dc_chamber", "opts": {}},
                    {"name": "dc_chamber", "opts": {"bcin": True}},
                ],
            },
            "steps": [],
        }

    def test_keeps_bcout_and_bcin_same_mode(self):
        self._sealed_manifest("out_t0", "dc", "t0", bcin=False)
        self._sealed_manifest(
            "in_t0", "dc", "t0", bcin=True, finished_at="2026-01-02T00:00:00Z",
        )
        subject = calib_approval.build_subject(1, self._both_chamber_plan())
        by_chamber = {
            (m["mode"], m["chamber"]): m["exec_id"]
            for m in subject["qa_manifests"]
        }
        self.assertEqual(by_chamber.get(("t0", "bcout")), "out_t0")
        self.assertEqual(by_chamber.get(("t0", "bcin")), "in_t0")

    def test_missing_bcin_mode_when_plan_needs_both(self):
        self._sealed_manifest("out_t0", "dc", "t0", bcin=False)
        subject = calib_approval.build_subject(1, self._both_chamber_plan())
        self.assertIn("qa_manifest:dc:t0:bcin", subject["missing"])
        self.assertNotIn("qa_manifest:dc:t0:bcout", subject["missing"])

    def test_bcin_does_not_drop_bcout(self):
        self._sealed_manifest("out_t0", "dc", "t0", bcin=False)
        st = self._both_chamber_plan()
        s1 = calib_approval.build_subject(1, st)
        self.assertEqual(
            [m["exec_id"] for m in s1["qa_manifests"] if m["chamber"] == "bcout"],
            ["out_t0"],
        )
        self._sealed_manifest(
            "in_t0", "dc", "t0", bcin=True, finished_at="2026-01-03T00:00:00Z",
        )
        s2 = calib_approval.build_subject(1, st)
        by_chamber = {
            (m["mode"], m["chamber"]): m["exec_id"]
            for m in s2["qa_manifests"]
        }
        self.assertEqual(by_chamber.get(("t0", "bcout")), "out_t0")
        self.assertEqual(by_chamber.get(("t0", "bcin")), "in_t0")

    def test_reso_all_is_separate_chamber_key(self):
        self._sealed_manifest("out_t0", "dc", "t0", bcin=False)
        self._sealed_manifest("all_reso", "dc", "reso", all_chambers=True)
        subject = calib_approval.build_subject(1, self._both_chamber_plan())
        by_key = {
            (m["mode"], m["chamber"]): m["exec_id"]
            for m in subject["qa_manifests"]
        }
        self.assertEqual(by_key.get(("t0", "bcout")), "out_t0")
        self.assertEqual(by_key.get(("reso", "all")), "all_reso")
        self.assertIn("qa_manifest:dc:t0:bcin", subject["missing"])
        # --all does not satisfy per-chamber reso requirement.
        self.assertIn("qa_manifest:dc:reso:bcout", subject["missing"])

    def test_reconcile_stale_approval(self):
        self._seed_hodo_complete()
        st = {
            "status": "approved",
            "finished": True,
            "plan": {
                "calib_particle": "pion",
                "steps": [{"name": "hodo", "opts": {"skip_t0": True}}],
            },
            "steps": [],
            "approval": {"subject_sha": "deadbeefdeadbeef", "via": "test"},
        }
        state_path = self.out / "calib/run00001/calib_dispatch_state.json"
        state_path.parent.mkdir(parents=True, exist_ok=True)
        state_path.write_text(json.dumps(st) + "\n", encoding="utf-8")
        changed = calib_approval.reconcile_stale_approval(1)
        self.assertTrue(changed)
        new_st = json.loads(state_path.read_text(encoding="utf-8"))
        self.assertEqual(new_st["status"], "awaiting_final_approval")
        self.assertNotIn("approval", new_st)


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

    def test_index_from_subject_uses_sealed(self):
        subject = {
            "qa_manifests": [
                {
                    "required": True,
                    "checks_ok": True,
                    "family": "hodo",
                    "mode": "hdprm",
                    "chamber": "",
                    "exec_id": "e1",
                    "sealed": [{
                        "name": "a.pdf",
                        "path": "/x/a.pdf",
                        "sha": "aa",
                        "exec_id": "e1",
                    }],
                },
                {
                    "required": True,
                    "checks_ok": False,
                    "family": "hodo",
                    "mode": "hdphc",
                    "chamber": "",
                    "exec_id": "e2",
                    "sealed": [],
                },
            ],
        }
        idx = self.idx.index_from_subject(1, subject)
        self.assertEqual(idx["n_pdf"], 1)
        self.assertEqual(idx["pdfs"][0]["exec_id"], "e1")
        self.assertEqual(idx["source"], "subject_sealed")


if __name__ == "__main__":
    unittest.main()
