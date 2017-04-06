"""Catalogue

Catalogue object contains an array of particles,
or information to read particles from a catalogue file

cat = lssps.open("catalogue.txt")
"""

import lssps._lssps as c

class Catalogue:
    def __init__(self):
        self._cat = c._catalogue_alloc()

    def __len__(self):
        return c._catalogue_len(self._cat)

    def __getitem__(self, index):
        return c._catalogue_asarray(self._cat)[index]

class CatalogueFile:
    def __init__(self, fname):
        self._filename = fname

    def __repr__(self):
        return 'CatalogueFile(\'%s\')' % self._filename

    @property
    def filename():
        return self._filename
    
    pass
        

def ascii(filename, *,
          xyz=(0, 1, 2), radec=None, r=None,
          weights=(), nbar=-1, Pest=0.0):
    f = CatalogueFile(filename)
    f._f = c._catalogue_file_ascii_alloc(filename, xyz, weights, nbar, Pest)
    
    return f

def gadget(filename, *,
           particle_types=(1), rsd=False):
    pass

def loadtxt(filename, *, xyz=[0, 1, 2], radec=None, r=None,
             weights=None, nbar=None):
    cat = Catalogue()
    c._catalogue_read_text(cat._cat, filename)
    
    return cat
