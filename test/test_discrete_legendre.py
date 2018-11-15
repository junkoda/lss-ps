import unittest
import numpy as np
import lssps

nc = 64
boxsize = 128

class TestDiscreteLegendre(unittest.TestCase):
    #def setUp(self):
    #    lssps.cosmology.init(omega_m, z_max, 2001)

    def test_discrete_legendre0(self):
        """
        (2l + 1)/2 int_-1^1 ~P_l(mu) dmu~
        """
        mu0 = lssps.grid.mu2(nc, boxsize, axis=2)
        mu0[:] = 1.0

        ps = lssps.power_spectrum.compute_discrete_power_multipoles(mu0,
                                            k_min=0.0, k_max=1.0, dk=0.1,
                                            line_of_sight=2)
        nbin = len(ps)

        for i in range(nbin):
            if ps.nmodes[i] > 0:
                self.assertAlmostEqual(ps.P0[i], 1.0, places=12)
                self.assertAlmostEqual(ps.P2[i], 0.0, places=12)
                self.assertAlmostEqual(ps.P4[i], 0.0, places=12)

    def test_discrete_legendre2(self):
        """
        (2l + 1)/2 int_-1^1 mu^2 ~P_l(mu) dmu~
        """
        mu2 = lssps.grid.mu2(nc, boxsize, axis=2)

        ps = lssps.power_spectrum.compute_discrete_power_multipoles(mu2,
                                            k_min=0.0, k_max=1.0, dk=0.1,
                                            line_of_sight=2)
        nbin = len(ps)

        for i in range(nbin):
            if ps.nmodes[i] > 0:
                self.assertAlmostEqual(ps.P0[i], 1.0/3.0, places=12)
                self.assertAlmostEqual(ps.P2[i], 2.0/3.0, places=12)
                self.assertAlmostEqual(ps.P4[i], 0.0, places=12)

    def test_discrete_legendre4(self):
        """
        (2l + 1)/2 int_-1^1 mu^4 ~P_l(mu) dmu~
        """
        mu4 = lssps.grid.mu2(nc, boxsize, axis=2)
        mu4[:] = mu4[:]**2

        ps = lssps.power_spectrum.compute_discrete_power_multipoles(mu4,
                                            k_min=0.0, k_max=1.0, dk=0.1,
                                            line_of_sight=2)
        nbin = len(ps)

        for i in range(nbin):
            if ps.nmodes[i] > 0:
                self.assertAlmostEqual(ps.P4[i], 8.0/35.0, places=12)


if __name__ == '__main__':
    unittest.main()
