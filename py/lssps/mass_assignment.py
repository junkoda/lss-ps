from lssps.grid import Grid
from lssps.catalogue import CatalogueFile, loadtxt
import lssps._lssps as c


_mass_assignment_scheme = {'NGP':1, 'CIC':2, 'TSC': 3}

def compute_density(catalogue_file, mas, nc, *,
                    x0=(0.0, 0.0, 0.0), boxsize=None,
                    interlacing=False):
    """
    Args:
        catalogue_file: CatalogueFile object
        mas: mass assignment scheme 'NGP', 'CIC', or 'TSC'
        nc: number of grids per dimension
        boxsize: length of the cubic box on a side
        centre: corrdinate of the centre of the box
    Returns:
        density with nc**3 grids if interlacing=False,
        (grid, grid_shifted) if interlacing=True
    """

    if not isinstance(catalogue_file, CatalogueFile):
        raise TypeError('CatalogueFile is not given to compute_density')
    
    grid = Grid(nc)
    n_mas = _mass_assignment_scheme[mas.upper()]

    if interlacing:
        grid_shifted = Grid(nc)
        c._mass_assignment(catalogue_file._f, x0, boxsize,
                           n_mas, (grid._grid, grid_shifted._grid))

        return (grid, grid_shifted)

    #if isinstance(cat, str):
    # catalogue is the filename
    #    cat = loadtxt(cat)

    c._mass_assignment(catalogue_file._f, x0, boxsize,
                       n_mas, grid._grid)

    return grid

