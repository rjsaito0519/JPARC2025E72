"""Generated schema directories must not invalidate otherwise fixed PHC inputs."""
import sys
import tempfile
import unittest
from pathlib import Path
import uproot
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.lib import hodo_btof_monitor as monitor

XML = ('<unpacker xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
       'xsi:noNamespaceSchemaLocation="{schema}">\n'
       '<control><skip>0</skip><max_loop>50000</max_loop>'
       '<error_check><format check="false"/></error_check></control>\n</unpacker>')

class UnpackFingerprintTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)
        self.schema = self.directory/"source.xsd"
        self.schema.write_text("<schema/>")
        self.unpack = self.directory/"unpacker.xml"
        self.unpack.write_text(XML.format(schema="source.xsd"))
        self.files = {"UNPACK": self.unpack}

    def fingerprint(self, schema="/project/runmanager/tmp/tmpabc/hodo_run00001.xsd", **changes):
        root = self.directory/"input.root"
        content = XML.format(schema=schema)
        for before, after in changes.items():
            content=content.replace(before,after)
        with uproot.recreate(root) as output:
            output["param/UNPACK"]=content
            output["param/HDPRM"]="fixed gain and offset"
            output["param/USER"]="fixed window"
        return monitor._fixed_parameters(root,self.files)

    def test_generated_schema_directory_is_stable(self):
        self.assertEqual(self.fingerprint(),self.fingerprint(
            "/project/runmanager/tmp/tmp_def/hodo_run00001.xsd"))

    def test_actual_control_changes_remain_visible(self):
        initial=self.fingerprint()
        for before,after in (("<skip>0","<skip>1"),("50000","25000"),('check="false"','check="true"')):
            with self.subTest(before=before):
                self.assertNotEqual(initial,self.fingerprint(**{before:after}))

    def test_non_generated_schema_path_or_name_is_not_ignored(self):
        initial=self.fingerprint()
        for schema in ("/other/runmanager/tmp/tmpabc/hodo_run00001.xsd",
                       "/project/runmanager/tmp/tmpabc/different.xsd",
                       "/project/runmanager/tmp/permanent/hodo_run00001.xsd"):
            self.assertNotEqual(initial,self.fingerprint(schema))

    def test_schema_contents_remain_visible(self):
        initial=self.fingerprint()
        self.schema.write_text("<schema changed='yes'/>")
        self.assertNotEqual(initial,self.fingerprint())

    def test_missing_schema_does_not_disable_validation(self):
        self.schema.unlink()
        with self.assertRaises(OSError):self.fingerprint()

    def test_root_macro_and_external_fallback_normalize_consistently(self):
        generated=[]
        for directory in ("tmp_a","tmp_b"):
            schema=self.directory/"runmanager/tmp"/directory/"hodo_run00001.xsd"
            schema.parent.mkdir(parents=True)
            schema.write_text(self.schema.read_text())
            generated.append(schema)
        root=self.directory/"with_macro.root"
        empty=self.directory/"empty.root"
        with uproot.recreate(root) as output:
            output["param/UNPACK"]=XML.format(schema=generated[0])
        with uproot.recreate(empty):pass
        self.unpack.write_text(XML.format(schema=generated[1]))
        recorded=monitor._fixed_parameters(root,self.files)
        fallback=monitor._fixed_parameters(empty,self.files)
        self.assertEqual(recorded["UNPACK"],fallback["UNPACK"])
        self.assertEqual(recorded["UNPACK_SCHEMA"],fallback["UNPACK_SCHEMA"])

if __name__=="__main__":
    unittest.main()
