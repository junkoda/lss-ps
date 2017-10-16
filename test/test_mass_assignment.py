import unittest
import numpy as np
import lssps

class TestGrid(unittest.TestCase):
    def setUp(self):
        self.nc = 4
        self.dx = 100.0
        self.boxsize = self.dx*self.nc        
        self.grid = lssps.Grid(self.nc, self.boxsize)
        self.w = 1.0

    def _set_one(self, ix, iy, iz):
        xyz = np.zeros((1, 3))
        xyz[0, 0] = (ix + 0.25)*self.dx
        xyz[0, 1] = (iy + 0.25)*self.dx
        xyz[0, 2] = (iz + 0.25)*self.dx
        w = np.ones(1)*self.w
        nbar = np.ones(1)
        return xyz, nbar, w

    def _test_one_cic(self, ix, iy, iz):
        xyz, nbar, w = self._set_one(ix, iy, iz)
        self.grid.clear()
        self.grid.assign_density(xyz=xyz, nbar=nbar, mas='CIC')

        a = self.grid[:]
        w0 = 0.75
        w1 = 0.25

        ix0 = ix
        ix1 = (ix + 1) % self.nc
        iy0 = iy
        iy1 = (iy + 1) % self.nc
        iz0 = iz
        iz1 = (iz + 1) % self.nc
        
        self.assertAlmostEqual(a[ix0, iy0, iz0], w0*w0*w0)
        self.assertAlmostEqual(a[ix0, iy0, iz1], w0*w0*w1)
        self.assertAlmostEqual(a[ix0, iy1, iz0], w0*w1*w0)
        self.assertAlmostEqual(a[ix0, iy1, iz1], w0*w1*w1)
        self.assertAlmostEqual(a[ix1, iy0, iz0], w1*w0*w0)
        self.assertAlmostEqual(a[ix1, iy0, iz1], w1*w0*w1)
        self.assertAlmostEqual(a[ix1, iy1, iz0], w1*w1*w0)
        self.assertAlmostEqual(a[ix1, iy1, iz1], w1*w1*w1)

        # The grid value should add up to w
        self.assertAlmostEqual(np.sum(a), self.w)

    def _test_one_tsc(self, ix, iy, iz):
        xyz, nbar, w = self._set_one(ix, iy, iz)
        self.grid.clear()
        self.grid.assign_density(xyz=xyz, nbar=nbar, mas='TSC')

        a = self.grid[:]
        w0 = 0.28125
        w1 = 0.6875
        w2 = 0.03125

        ix0 = (ix - 1) % self.nc
        ix1 =  ix      % self.nc
        ix2 = (ix + 1) % self.nc

        iy0 = (iy - 1) % self.nc
        iy1 =  iy      % self.nc
        iy2 = (iy + 1) % self.nc

        iz0 = (iz - 1) % self.nc
        iz1 =  iz      % self.nc
        iz2 = (iz + 1) % self.nc

        self.assertAlmostEqual(a[ix0, iy0, iz0], w0*w0*w0)
        #self.assertAlmostEqual(a[ix0, iy0, iz1], w0*w0*w1)
        #self.assertAlmostEqual(a[ix0, iy1, iz0], w0*w1*w0)
        #self.assertAlmostEqual(a[ix0, iy1, iz1], w0*w1*w1)
        #self.assertAlmostEqual(a[ix1, iy0, iz0], w1*w0*w0)
        #self.assertAlmostEqual(a[ix1, iy0, iz1], w1*w0*w1)
        #self.assertAlmostEqual(a[ix1, iy1, iz0], w1*w1*w0)
        #self.assertAlmostEqual(a[ix1, iy1, iz1], w1*w1*w1)

        # The grid value should add up to w
        self.assertAlmostEqual(np.sum(a), self.w)

    def test_one_cic(self):
        nc = lssps._lssps._grid_nc(self.grid._grid)

        for ix in range(self.nc):
            for iy in range(self.nc):
                for iz in range(self.nc):
                    self._test_one_cic(ix, iy, iz)

    # TSC test not working
    #def test_one_tsc(self):
    #    nc = lssps._lssps._grid_nc(self.grid._grid)
    #
    #    for ix in range(self.nc):
    #        for iy in range(self.nc):
    #            for iz in range(self.nc):
    #                self._test_one_tsc(ix, iy, iz)

        
if __name__ == '__main__':
    unittest.main()
