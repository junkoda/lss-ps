import unittest
import numpy as np
import lssps

class TestGrid(unittest.TestCase):
    def setUp(self):
        self.nc = 4
        self.grid = lssps.Grid(self.nc)

    def test_nc(self):
        nc = lssps._lssps._grid_nc(self.grid._grid)
        self.assertEqual(nc, self.nc)

    def test_mode(self):
        self.assertEqual(self.grid.mode, 'real-space')

    def test_fx(self):
        a = self.grid[:]
        self.assertEqual(a.shape, (self.nc, self.nc, self.nc))
        self.assertTrue(np.all(a == 0.0))
        
if __name__ == '__main__':
    unittest.main()
