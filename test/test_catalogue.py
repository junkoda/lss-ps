import unittest
import numpy as np
import lssps

class TestCatalogue(unittest.TestCase):
    def setUp(self):
        self.cat = lssps.Catalogue()
        self.cat.loadtxt('../data/wizcola_realspace.txt')
        print('len = ', len(self.cat))

    def test_len(self):
        self.assertEqual(len(self.cat), 798698)

if __name__ == '__main__':
    unittest.main()
