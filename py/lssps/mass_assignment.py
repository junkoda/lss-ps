import lssps.grid
from lssps.catalogue import CatalogueFile, loadtxt
import lssps._lssps as c


_mass_assignment_scheme = {'NGP':1, 'CIC':2, 'TSC': 3}

def compute_density(catalogue_file, nc, *,
                    mas='CIC',
                    x0=(0.0, 0.0, 0.0), boxsize=None,
                    interlacing=False):
    """
    Compute density of particles in `catalogue_file` on a grid using a 
    mass assignment scheme `grid`.

    Args:
        catalogue_file: CatalogueFile object
        nc: number of grids per dimension
        boxsize: length of the cubic box on a side
        centre: corrdinate of the centre of the box
        mas: mass assignment scheme 'NGP', 'CIC', or 'TSC'


    Returns:
        A tuple of density grids
        (grid,) if interlacing=False,
        (grid, grid_shifted) if interlacing=True

    Example:
        grid = compute_density(lssps.catalogue.ascii('catalogue.txt'), 64)[0]
        a = grid[:] # 3D array of density field
    """

    if not isinstance(catalogue_file, CatalogueFile):
        raise TypeError('CatalogueFile is not given to compute_density')
    
    grid = lssps.grid.zeros(nc)
    n_mas = _mass_assignment_scheme[mas.upper()]

    if interlacing:
        grid_shifted = lssps.grid.zeros(nc)
        c._mass_assignment(catalogue_file._f, x0, boxsize,
                           n_mas, grid._grid, grid_shifted._grid)

        return (grid, grid_shifted)

    #if isinstance(cat, str):
    # catalogue is the filename
    #    cat = loadtxt(cat)

    c._mass_assignment(catalogue_file._f, x0, boxsize,
                       n_mas, grid._grid, None)

    return (grid,)


def assign_density(grid, *,
                   cat=None,
                   xyz=None, weight=None, nbar=None,
                   mas='CIC'):
    """Assign density to an existing grid.
    Args:
        grid: A Grid object
        mas:  The mass assignment scheme 'NGP', 'CIC', or 'TSC'
    """

    if grid.boxsize <= 0.0:
        raise AssertionError('grid.boxsize is not set.')
        
    
    
    if xyz is not None:
        c._mass_assignment_from_array(xyz, weight, nbar,
                                      _mass_assignment_scheme[mas],
                                      grid._grid)
    else:
        RuntimeError('xyz not provided')
        # TODO cat
