#
# $ make
#

from distutils.core import setup, Extension
import numpy as np
import os


# Remove the "-Wstrict-prototypes" compiler option, which isn't valid for C++.
# https://stackoverflow.com/questions/8106258
import distutils.sysconfig
cfg_vars = distutils.sysconfig.get_config_vars()
for key, value in cfg_vars.items():
    if type(value) == str:
        cfg_vars[key] = value.replace("-Wstrict-prototypes", "")

# directories for include -I(idir)
idirs = os.environ["IDIRS"]
if idirs:
    idirs = idirs.split()
else:
    idirs = []

idirs = ['../lib', np.get_include()] + idirs

# directories for libraries -L(dir)
ldirs = os.environ["LDIRS"]
if ldirs:
    ldirs = ldirs.split()
else:
    ldirs = []

libs = os.environ['LIBS'].split()

setup(name='lssps',
      version='0.0.1',
      author='Jun Koda',
      py_modules=['lssps.msg', 'lssps.comology', 
                  'lssps.grid', 
      ],
      ext_modules=[
          Extension('lssps._lssps',
                    ['../lib/msg.cpp', '../lib/config.cpp',
                     '../lib/cosmology.cpp', '../lib/grid.cpp',
                     '../lib/mass_assignment.cpp', '../lib/interlacing.cpp',
                     '../lib/multipole.cpp', '../lib/power_spectrum.cpp',
                     '../lib/mass_assignment_serial.cpp',
                     '../lib/mass_assignment_atomic.cpp',
                     'py_package.cpp', 'py_config.cpp', 'py_msg.cpp',
                     'py_cosmology.cpp',
                     'py_grid.cpp', 'py_mass_assignment.cpp',
                     'py_power_spectrum.cpp', 'py_interlacing.cpp',
                     'py_performance.cpp',
                     'py_mass_assignment_variations.cpp',
                    ],
                    include_dirs = idirs,
                    extra_compile_args = [os.environ['OPT']],
                    library_dirs = ldirs,
                    libraries = libs,
                    undef_macros = ['NDEBUG'],
          )
      ],
      packages=['lssps'],
)
