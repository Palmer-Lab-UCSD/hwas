"""Test hwas/_config

"""


from unittest import TestCase, main
import io


from hwas import _config

MOCK_CONFIG_FILE = """
[common]
option_a
option_b = b
option_c = c
option_d = ${option_a}
option_e = ${option_b}


[query]
q = value
qa = ${common:option_a}
qb = ${common:option_b}
qc = ${common:option_d}
qd = ${common:option_e}
qf = $$
qg = $${a:b}
qh = {a:b}
qi = {common:option_c}
qj = common:option_c
"""



class TestConfigParser(TestCase):
    def setUp(self):
        self.cfg = _config.ConfigParser()
        self.cfg._cfg.read_string(MOCK_CONFIG_FILE)

    def test_is_interpolation(self):
        self.assertTrue(self.cfg.is_interpolation("common", "option_d"))
        self.assertTrue(self.cfg.is_interpolation("query", "qa"))


        self.assertFalse(self.cfg.is_interpolation("query", "q"))
        self.assertFalse(self.cfg.is_interpolation("common", "option_a"))
        self.assertFalse(self.cfg.is_interpolation("query", "qf"))
        self.assertFalse(self.cfg.is_interpolation("query", "qg"))
        self.assertFalse(self.cfg.is_interpolation("query", "qh"))
        self.assertFalse(self.cfg.is_interpolation("query", "qi"))
        self.assertFalse(self.cfg.is_interpolation("query", "qj"))


    def test_get_key_val_pair(self):
        self.assertEqual(self.cfg.get("common", "option_b"), "b")
        self.assertEqual(self.cfg.get("query", "q"), "value")
        self.assertEqual(self.cfg.get("query", "qd"), "b")

        self.assertIsNone(self.cfg.get("common", "option_a"))
        self.assertIsNone(self.cfg.get("common", "option_d"))
        self.assertIsNone(self.cfg.get("query", "qa"))
        self.assertIsNone(self.cfg.get("query", "qc"))


    def test_get_option_interpolator(self):
        
        self.assertIsNone(self.cfg.get_option_interpolator("common", "option_a"))
        self.assertIsNone(self.cfg.get_option_interpolator("common", "option_b"))
        self.assertIsNone(self.cfg.get_option_interpolator("query", "q"))
        self.assertIsNone(self.cfg.get_option_interpolator("query", "qg"))
        self.assertIsNone(self.cfg.get_option_interpolator("query", "qh"))
        self.assertIsNone(self.cfg.get_option_interpolator("query", "qi"))
        self.assertIsNone(self.cfg.get_option_interpolator("query", "qj"))


        self.assertTupleEqual(
                self.cfg.get_option_interpolator("common", "option_d"),
                ("common", "option_a")
                )
        self.assertTupleEqual(
                self.cfg.get_option_interpolator("query", "qa"),
                ("common", "option_a")
                )
        self.assertTupleEqual(
                self.cfg.get_option_interpolator("query", "qb"),
                ("common", "option_b")
                )
        self.assertTupleEqual(
                self.cfg.get_option_interpolator("query", "qd"),
                ("common", "option_e")
                )

    def test_set(self):
        self.cfg.set("common", "option_b", "c")
        self.assertEqual(self.cfg.get("common", "option_b"), "c")

        self.cfg.set("common", "option_a", "a")
        self.assertEqual(self.cfg.get("common", "option_a"), "a")

        self.cfg.set("common", "option_d", "a")
        self.assertEqual(self.cfg.get("common", "option_d"), "a")

        self.cfg.set("common", "option_c", None)
        self.assertIsNone(self.cfg.get("common", "option_c"))

        self.cfg.set("query", "qd", "interpolating")
        self.assertEqual(self.cfg.get("common", "option_b"), "interpolating")
        self.assertTrue(self.cfg.is_interpolation("query", "qd"))
        self.assertTupleEqual(
                self.cfg.get_option_interpolator("query", "qd"),
                ("common", "option_e")
                )
        self.assertTrue(self.cfg.is_interpolation("common", "option_e"))
        self.assertTupleEqual(
                self.cfg.get_option_interpolator("common", "option_e"),
                ("common", "option_b")
                )




if __name__ == "__main__":
    main()
