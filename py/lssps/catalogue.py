"""Catalogue

Catalogue is a collection of particles

cat = lssps.Catalogue()
"""

import lssps._lssps as c

class Catalogue:
    def __init__(self, filename = None):
        self._cat = c._catalogue_alloc()

    def __len__(self):
        return c._catalogue_len(self._cat)
    
    def loadtxt(self, filename):
        c._catalogue_read_text(self._cat, filename)
