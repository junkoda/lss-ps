from lssps.catalogue import Catalogue
from lssps.grid import Grid
from lssps.mass_assignment import compute_density, assign_density
from lssps.compute import compute_power_spectrum

from lssps.power_spectrum import PowerSpectrum
import lssps.msg

loglevel = lssps.msg.LogLevel()

#import lssps.power_spectrum
_mass_assignment_scheme = {'NGP':1, 'CIC':2, 'TSC': 3}    


