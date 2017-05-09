import lssps._lssps as c


def interlacing(grid):
    """interlacing(grid)
    Arg:
       grid (lssps.grid.Grid) A grid with shifted grid, grid.shifted.
    """

    if grid.shifted is None:
        raise TypeError('grid.shifted is required for interlacing')

    c._interlacing(grid._grid, grid.shifted._grid)

    return grid
