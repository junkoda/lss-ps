import unittest
import numpy as np
import lssps

omega_m = 0.31
z_max= 1.2

class TestCatalogue(unittest.TestCase):
    def setUp(self):
        lssps.cosmology.init(omega_m, z_max, 2001)

    def test_distance_redshift(self):
        d = []
        for z in [0.6, 0.9, 1.1]:
            d.append(lssps.cosmology.compute_comoving_distance(z))
        z= lssps.cosmology.redshift(np.array(d))

        print(z)
        
        self.assertAlmostEqual(z[0], 0.6, places=4)
        self.assertAlmostEqual(z[1], 0.9, places=4)
        self.assertAlmostEqual(z[2], 1.1, places=4)

if __name__ == '__main__':
    unittest.main()
