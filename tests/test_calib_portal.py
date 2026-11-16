"""Run-independent calibration portal, routing and token isolation."""
import importlib.util
import json
from pathlib import Path
import tempfile
import threading
import unittest
import urllib.error
import urllib.request
from unittest.mock import patch

REPO = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("portal_web", REPO / "calibration/mip_web/mip_web.py")
web = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(web)


class PortalTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.folder = Path(self.directory.name)
        for name in ("OUTPUT_DIR", "DATA_DIR", "WORK_DIR", "PARAM_DIR"):
            patcher = patch.object(web.config, name, self.folder/name)
            patcher.start()
            self.addCleanup(patcher.stop)
        self.portal = web.Portal()
        for run in (12, 13):
            path = Path(web.config.DATA_DIR)/f"run{run:05d}_Hodo.root"
            path.parent.mkdir(parents=True, exist_ok=True)
            path.touch()
            self.fit(run, "Pi", "hdprm")
            self.fit(run, "Pi", "t0")

    def fit(self, run, particle, mode):
        kind = "T0_Offset" if mode == "t0" else "BH2_HDPRM"
        path = Path(web.config.OUTPUT_DIR)/"root"/f"run{run:05d}"/f"run{run:05d}_{kind}_{particle}.root"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()

    def test_available_only_saved_results_with_matching_decode(self):
        self.fit(14, "Pi", "hdprm")
        self.assertEqual(4, len(self.portal.available()))
        self.assertEqual({12, 13}, {row["run"] for row in self.portal.available()})
        self.assertFalse(self.portal.applications)
        self.assertTrue(all(not row["ready"] for row in self.portal.available()))

    def test_workspaces_are_independent_and_reused(self):
        first = self.portal.resolve({"run": ["12"]})
        self.assertIs(first, self.portal.resolve({"run": ["12"]}))
        others = [self.portal.resolve({"run": ["13"]}), self.portal.resolve({"run": ["12"], "mode": ["t0"]}),
                  self.portal.resolve({"run": ["12"], "particle": ["K"]})]
        self.assertEqual(4, len({app.token for app in [first, *others]}))
        self.assertEqual(4, len({app.staging_path for app in [first, *others]}))
        for app in others:
            self.assertIs(first.preview_lock, app.preview_lock)

    def test_particle_ambiguity_requires_selection_without_a_session(self):
        self.fit(12, "K", "hdprm")
        with self.assertRaisesRegex(ValueError, "PiとK"):
            self.portal.resolve({"run": ["12"]})
        path = Path(web.config.OUTPUT_DIR)/"calib/run00012/lanes/hodo/calib_dispatch_state.json"
        web.write_json(path, {"status": "held", "calib_particle": "kaon"})
        self.assertEqual("K", self.portal.resolve({"run": ["12"]}).particle)

    def test_invalid_context_never_creates_a_workspace(self):
        for query in ({}, {"run": ["../12"]}, {"run": ["12", "13"]}, {"run": ["99999"]},
                      {"run": ["12"], "mode": ["unknown"]}, {"run": ["12"], "particle": ["Proton"]}):
            with self.assertRaises(ValueError):
                self.portal.resolve(query)
        self.assertFalse(self.portal.applications)

    def test_http_root_links_and_cross_workspace_token_rejection(self):
        first = self.portal.resolve({"run": ["12"]})
        second = self.portal.resolve({"run": ["13"]})
        for app in (first, second):
            app.state = lambda app=app: {"run": app.run, "particle": app.particle, "mode": app.mode}
            app.edit = lambda request: {"edited": True}
        server = web.ThreadingHTTPServer(("127.0.0.1", 0), web.handler(portal=self.portal))
        worker = threading.Thread(target=server.serve_forever, daemon=True)
        worker.start()
        self.addCleanup(worker.join)
        self.addCleanup(server.server_close)
        self.addCleanup(server.shutdown)
        base = f"http://127.0.0.1:{server.server_port}"
        with urllib.request.urlopen(base+"/") as response:
            self.assertIn("較正レビュー", response.read().decode())
        with urllib.request.urlopen(base+"/runs") as response:
            self.assertEqual(4, len(json.load(response)["runs"]))
        for run in (12, 13):
            with urllib.request.urlopen(base+f"/state?run={run}") as response:
                self.assertEqual(run, json.load(response)["run"])
        request = urllib.request.Request(base+"/edit?run=13", data=b"{}", headers={"X-MIP-Token": first.token})
        with self.assertRaises(urllib.error.HTTPError) as error:
            urllib.request.urlopen(request)
        self.assertEqual(403, error.exception.code)
        with urllib.request.urlopen(base+"/?run=13&mode=hdprm") as response:
            page = response.read().decode()
            self.assertIn(second.token, page)
            self.assertNotIn(first.token, page)
            self.assertNotIn("__PARTICLE__", page)
        with urllib.request.urlopen(base+"/?run=99999") as response:
            self.assertIn("Hodo decode", response.read().decode())


if __name__ == "__main__":
    unittest.main()
