import lssps._lssps as c
import numpy as np

def compute_multipoles(xyz, nbar, r_max, *, weight=None, dr=1.0,
                       kind='end-point'):
    if kind == 'end-point':
        kind = -1
    elif isinstance(kind, int) and (0 <= kind < 3):
        pass
    else:
        raise ValueError('Unknown kind of rr multipoles: {}'.format(kind))
    
    _rr= c._rr_compute_multipoles(xyz, weight, nbar, r_max, dr, kind)

    n = round(r_max/dr) + 1

    d = {}
    d['r'] = (np.arange(n) + 0.5)*dr
    d['n'] = n

    for n in [0, 1, 2, 3, 4]:
        for l in [0, 1, 2, 3, 4]:
            d['rr%d%d' % (n, l)] = np.copy(c._rr_asarray(_rr, n, l))

    return d
