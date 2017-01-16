import lssps._lssps as c

class Grid:
    """Grid is a 3-dimensional cubic grid with nc points per dimension"""
    
    def __init__(self, nc):
        self._grid = c._grid_alloc(nc)


    def __getitem__(self, index):
        if self.mode == 'real-space':
            return c._grid_fx_asarray(self._grid)[index]
        elif self.mode == 'fourier-space':
            return c._grid_fk_asarray(self._grid)[index]
        return None
        

    def fft(self):
        c._grid_fft(self._grid)

    @property
    def mode(self):
        """Current mode of the grid:
        'real-space' or 'fourier-space'
        """
        return c._grid_mode(self._grid)

