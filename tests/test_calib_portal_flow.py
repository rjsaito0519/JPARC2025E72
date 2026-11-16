"""Exercise portal HTTP editing and replay boundaries using isolated calibration records."""
import copy
import importlib.util
import json
from pathlib import Path
import re
import tempfile
import threading
import unittest
import urllib.error
import urllib.parse
import urllib.request
from unittest.mock import patch

import numpy as np
import uproot

REPO = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("portal_flow_web", REPO/"calibration/mip_web/mip_web.py")
web = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(web)


class PortalFlowTests(unittest.TestCase):
    def setUp(self):
        work = Path(web.config.WORK_DIR)
        self.directory = tempfile.TemporaryDirectory(dir=work if work.is_dir() else None)
        self.addCleanup(self.directory.cleanup)
        self.folder = Path(self.directory.name)
        for name in ("OUTPUT_DIR", "WORK_DIR", "DATA_DIR", "PARAM_DIR"):
            patcher = patch.object(web.config, name, self.folder/name)
            patcher.start()
            self.addCleanup(patcher.stop)
        root = Path(web.config.DATA_DIR)/"run00012_Hodo.root"
        root.parent.mkdir(parents=True)
        root.touch()
        self.portal = web.Portal()
        self.server = web.ThreadingHTTPServer(("127.0.0.1", 0), web.handler(portal=self.portal))
        worker = threading.Thread(target=self.server.serve_forever, daemon=True)
        worker.start()
        self.addCleanup(worker.join)
        self.addCleanup(self.server.server_close)
        self.addCleanup(self.server.shutdown)
        self.url = f"http://127.0.0.1:{self.server.server_port}"

    def workspace(self, mode):
        app = self.portal.resolve({"run": ["12"], "mode": [mode], "particle": ["Pi"]})
        app.fingerprint = lambda: "source-"+mode
        settings = ({"parameters": [20., -.3, .1], "range": [-.7, .1], "function": "gaussian", "range_mode": "fixed"}
                    if mode == "t0" else
                    {"parameters": [12000., 250., 30., 4.], "range": [180., 340.], "function": "auto", "range_mode": "auto"})
        result = ({"parameters": [21., -.31, .09], "range": [-.7, .1], "function": "gaussian"}
                  if mode == "t0" else
                  {"parameters": [13000., 255., 28., 3.], "range": [200., 320.], "function": "gaussian"})
        result.update(success=True, status=0, covariance=3, ndf=20, curve=[[0., 1.]], origin="isolated fixture")
        app.execute_preview = lambda *args, **kwargs: copy.deepcopy(result)
        app.automatic = lambda *args: dict(copy.deepcopy(result), parameters=list(settings["parameters"]), range=list(settings["range"]))
        if mode == "t0":
            output = app.fit_output()
            output.parent.mkdir(parents=True, exist_ok=True)
            with uproot.recreate(output) as saved:
                for segment in range(web.config.NUM_OF_CH["BH2"]):
                    saved[f"T0_seg{segment}_TimeOffset_Pi"] = np.ones(60), np.linspace(-10., 10., 61)
        state_path = app.base/"lanes/hodo/calib_dispatch_state.json"
        hold = {"id": "old-"+mode, "step": app.expected_step, "rc": 0}
        state = {"status": "held", "hold": hold, "calib_particle": "pion"}
        web.write_json(state_path, state)
        (app.base/"calib_session.yml").write_text("auto: false\n")
        query = urllib.parse.urlencode({"run": 12, "mode": mode, "particle": "Pi"})
        with urllib.request.urlopen(self.url+"/?"+query) as response:
            token = re.search(r"token='([^']+)'", response.read().decode())[1]
        self.assertEqual(app.token, token)
        return app, settings, result, state_path, state, hold, query

    def request(self, path, query, body=None, token=None):
        headers = {"X-MIP-Token": token} if token else {}
        payload = json.dumps(body).encode() if body is not None else None
        if body is not None:
            headers["Content-Type"] = "application/json"
        request = urllib.request.Request(self.url+path+"?"+query, data=payload, headers=headers)
        with urllib.request.urlopen(request, timeout=10) as response:
            return json.load(response)

    def test_adc_and_t0_http_write_replay_and_new_review(self):
        for mode in ("hdprm", "t0"):
            with self.subTest(mode=mode):
                app, settings, result, state_path, state, hold, query = self.workspace(mode)
                key = "bh2-3-t0" if mode == "t0" else "bh2-3-u"
                # Other channels, particle settings and locks must survive the HTTP Write.
                original = {"fit": {"T1": {"hdprm": {"t1-0-u": {"mip_left": 123.}}},
                                    "BH2": {mode: {key: {app.override_key: {"K": settings}}}}},
                            "calib_locks": {"BH2": {"hdprm": ["existing-lock"]}}}
                web.fit_overrides.save(12, original)
                before = self.request("/state", query)
                fit_request = {"fingerprint": before["fingerprint"], "key": key, "settings": settings}
                fitted = self.request("/fit", query, fit_request, app.token)
                self.assertIn("receipt", fitted)
                self.assertFalse(self.request("/state", query)["staged"])
                staged = self.request("/stage", query, dict(fit_request, receipt=fitted["receipt"]), app.token)
                review = self.request("/review", query, {"fingerprint": staged["fingerprint"]}, app.token)
                self.assertEqual([key], [row["key"] for row in review["rows"]])
                write_request = {"fingerprint": review["fingerprint"], "revision": review["revision"],
                                 "hold_id": review["gate"]["hold_id"]}

                def retry(operation, run, hold_id, background_resume):
                    self.assertEqual(("retry", 12, hold["id"], True), (operation, run, hold_id, background_resume))
                    web.write_json(state_path, dict(state, status="running"))
                    return True, "isolated replay accepted"

                with patch("common.lib.calib_ops.run_dispatch", side_effect=retry) as replay:
                    written = self.request("/write", query, write_request, app.token)
                    self.assertEqual("running", written["job"]["status"])
                    with self.assertRaises(urllib.error.HTTPError) as error:
                        self.request("/write", query, write_request, app.token)
                    self.assertEqual(409, error.exception.code)
                    replay.assert_called_once()
                actual = web.fit_overrides.load(12)
                self.assertEqual(original["fit"]["T1"], actual["fit"]["T1"])
                self.assertEqual(original["calib_locks"], actual["calib_locks"])
                values = actual["fit"]["BH2"][mode][key][app.override_key]
                self.assertEqual(settings, values["K"])
                self.assertEqual(result["parameters"], values["Pi"]["parameters"])
                # Model only the dispatcher's completed QA gate; no ROOT job or Slack API is invoked.
                new_hold = {"id": "new-"+mode, "step": app.expected_step, "rc": 0}
                web.write_json(state_path, dict(state, hold=new_hold))
                completed = self.request("/state", query)
                self.assertEqual("done", completed["job"]["status"])
                self.assertFalse(completed["staged"])
                self.assertEqual(new_hold["id"], completed["gate"]["hold_id"])
                self.assertIn("Slack", completed["job"]["message"])

    def test_stale_review_gate_and_other_particle_block_http_write(self):
        app, settings, result, state_path, state, hold, query = self.workspace("hdprm")
        request = {"fingerprint": app.state()["fingerprint"], "key": "bh2-3-u", "settings": settings}
        fitted = self.request("/fit", query, request, app.token)
        staged = self.request("/stage", query, dict(request, receipt=fitted["receipt"]), app.token)
        write_request = {"fingerprint": staged["fingerprint"], "revision": staged["revision"], "hold_id": hold["id"]}
        with patch("common.lib.calib_ops.run_dispatch") as replay:
            for changed in (dict(state, hold=dict(hold, id="changed")), dict(state, calib_particle="kaon")):
                web.write_json(state_path, changed)
                with self.assertRaises(urllib.error.HTTPError) as error:
                    self.request("/write", query, write_request, app.token)
                self.assertEqual(409, error.exception.code)
            replay.assert_not_called()
        self.assertFalse(web.fit_overrides.overrides_path(12).exists())


if __name__ == "__main__":
    unittest.main()
