import unittest
import numpy as np
import lssps

class TestConfig(unittest.TestCase):
    #def setUp(self):

    def test_sizeof_float(self):
        """Test that sizeof(Float) is consistent between the library and the
        Python extension. If this fails check CPPFLAG -DDOUBLEPRECISION.
        """

        s = lssps._lssps._config_sizeof_float()
        self.assertEqual(s[0], s[1])

if __name__ == '__main__':
    unittest.main()
