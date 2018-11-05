import lssps._lssps as c
import numpy as np

def compute_multipoles(xyz, nbar, r_max, *, weight=None, dr=1.0):
    _rr= c._rr_compute_multipoles(xyz, weight, nbar, r_max, dr)

    n = round(r_max/dr) + 1

    d = {}
    d['r'] = (np.arange(n) + 0.5)*dr
    d['n'] = n

    for i in [0, 1, 2, 3, 4]:
        d['rr%d' % i] = np.copy(c._rr_asarray(_rr, i))

    return d
