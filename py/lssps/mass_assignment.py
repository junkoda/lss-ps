import lssps._lssps as c


def cic(cat, x0, boxsize, grid):
    c._mass_assignment_cic(cat._cat, x0, boxsize, grid._grid)

