import lssps._lssps as c

def interlacing(grid):
    """interlacing((grid, grid_shifted))
    Arg:
       grid, a pair 
    
    """
    assert(len(grid) == 2)
    
    c._interlacing(grid[0]._grid, grid[1]._grid)

    return grid[0]
