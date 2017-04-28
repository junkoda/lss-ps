"""
High-level functions
"""
import lssps._lssps as c
import lssps.power_spectrum
import lssps.grid

    
def compute_power_spectrum(data, rand, nc, estimator='plane-parallel', *,
                           mas='CIC',
                           x0=(0,0,0), boxsize=None,
                           k_min=0.0, k_max=1.0, dk=0.01, nmu=0,
                           subtract_shotnoise=True,
                           correct_mas=True,
                           interlacing=False):
    """compute_power_spectrum(data, rand, nc, 
                              mas='CIC', x0=None, boxsize=None)
    Args:
        data:    data CatalogueFile
        rand:    random CatalogueFile
        nc:      number of grids per dimension
        x0:      left bottom coordinate of the box
        boxsize: length of the cubic box on a side
        subtract_shotnoise: subtract constant shotnoise from power spectrum
        correct_mas: correct for mass assignment window function
        interlacing: reduce aliasing using interlacing
    """
    
    assert(nc > 0)

    if boxsize is None:
        boxsize = data.boxsize
        if boxsize is None:
            raise TypeError('boxsize is not given')

    # Read catalogue and compute density on grids
    grid_data = lssps.compute_density(data, nc, mas=mas,
                             x0=x0, boxsize=boxsize, interlacing=interlacing)

    if rand is None:
        grid_rand = None
    else:
        grid_rand = lssps.compte_density(rand, nc, mas=mas,
                             x0=x0, boxsize=boxsize, interlacing=interlacing)
 
    # compute fluctuation
    grid_delta = lssps.compute_fluctuation(grid_data, grid_rand)
        
    # FFT
    for d in grid_delta:
        d.fft()

    # interlacing
    if interlacing:
        assert(len(grid_delta) == 2)
        c._interlacing(grid_delta[0]._grid, grid_delta[1]._grid)

    # power spectrum
    _ps = c._power_spectrum_compute_plane_parallel(k_min, k_max, dk, nmu,
                grid_delta[0]._grid,
                int(subtract_shotnoise), int(correct_mas))
    

    return lssps.PowerSpectrum(_ps)
