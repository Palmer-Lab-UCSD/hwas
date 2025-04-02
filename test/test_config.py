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
option_f = ${option_g}
option_g = ${option_f}


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


    def test_init(self):
        tmp = _config.ConfigParser(max_recursion_depth = 15)
        self.assertEqual(tmp.max_recursion_depth, 15)

        with self.assertRaises(ValueError):
            _config.ConfigParser(max_recursion_depth = -3)

        with self.assertRaises(ValueError):
            _config.ConfigParser(max_recursion_depth = 3.5)

        with self.assertRaises(ValueError):
            _config.ConfigParser(max_recursion_depth = "4")
            
        with self.assertRaises(ValueError):
            _config.ConfigParser(max_recursion_depth = [4])


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


    def test_get(self):
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


    def test_traverse_interpolators(self):
        with self.assertRaises(RecursionError):
            self.cfg._traverse_interpolators("common", "option_g")

        with self.assertRaises(RecursionError):
            self.cfg._traverse_interpolators("common", "option_f")


class TestDynamicConfigSection(TestCase):
    def test_init(self):
        tmp = _config.DynamicConfigSection("common")
        self.assertEqual(tmp.name, "common")

        self.assertListEqual(tmp._dynamic_option_names, [])

    def test_set_property(self):
        tmp = _config.DynamicConfigSection("common")
        tmp._x = 0
        tmp._set_property("x")
        self.assertEqual(tmp.x, tmp._x)

        tmp.x = 25
        self.assertEqual(tmp.x, 25)
        self.assertEqual(tmp._x, 25)

        with self.assertRaises(TypeError):
            tmp._set_property(None)


        # because I use a classmethod to setup the property, the
        # property will be defined at the class level.  Meaning
        # that it will be applied to all existing and subsequent
        # instances of the class.  
        tmp_2 = _config.DynamicConfigSection("query")
        self.assertTrue("x" in tmp_2.__dir__())
        tmp_2._x = 45
        self.assertEqual(tmp_2.x, 45)

        tmp_2._y = "a"
        tmp_2._set_property("y")
        self.assertEqual(tmp_2.y, tmp_2._y)

        # Here we see that the `tmp` class instance now has property
        # y without having explicity made the property in the 
        # instance `tmp`

        self.assertTrue("y" in tmp.__dir__())
        with self.assertRaises(AttributeError):
            a = tmp.y

        tmp._y = "Weird behavior"
        self.assertEqual(tmp.y, "Weird behavior")


    def test_opt_to_attr(self):
        tmp = _config.DynamicConfigSection("common")
        self.assertEqual(tmp._opt_to_attr("x"), "_x")
        self.assertEqual(tmp._opt_to_attr("_x"), "__x")

        self.assertEqual(tmp._opt_to_attr(None), "_None")
        self.assertEqual(tmp._opt_to_attr(False), "_False")

        self.assertEqual(tmp._opt_to_attr(1235.3422), "_1235.3422")


# TODO
# class TestGetConfigSection(TestCase):
# class TestLoadDefaultConfig(TestCase):
# class TestInit(TestCase):



if __name__ == "__main__":
    main()
