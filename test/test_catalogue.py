import unittest
import numpy as np
import lssps

class TestCatalogue(unittest.TestCase):
    def setUp(self):
        self.cat = lssps.Catalogue()
        self.cat.loadtxt('../data/wizcola_realspace.txt')

    def test_len(self):
        self.assertEqual(len(self.cat), 798698)

    def test_particles(self):
        """Test conversion of the catalogue to an array"""
        a = self.cat[:]
        
        # first particle
        self.assertAlmostEqual(a[0, 0], 8.490359e-01)
        self.assertAlmostEqual(a[0, 1], 4.980907e+00)
        self.assertAlmostEqual(a[0, 2], 1.014142e+01)
        self.assertAlmostEqual(a[0, 3], 1.0)

        # second particle
        self.assertAlmostEqual(a[1, 0], 5.946966e-01)
        self.assertAlmostEqual(a[1, 1], 2.829677e-01)
        self.assertAlmostEqual(a[1, 2], 2.601586e+01)

        # last particle
        n = a.shape[0]
        self.assertAlmostEqual(a[n - 1, 0], 5.968730e+02)
        self.assertAlmostEqual(a[n - 1, 1], 5.272345e+02)
        self.assertAlmostEqual(a[n - 1, 2], 5.882148e+02)
        
if __name__ == '__main__':
    unittest.main()
