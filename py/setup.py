from distutils.core import setup, Extension
import numpy as np
import os

#dir = os.getcwd()
#print(dir)
#print("np.get_include()")
#print(np.get_include())

idirs = ['../lib', np.get_include()] + os.environ["IDIRS"].split(' ')
ldirs = os.environ["LDIRS"].split(' ')
libs = os.environ['LIBS'].split(' ')

setup(name='lssps',
      version='0.0.1',
      author='Jun Koda',
      py_modules=['lssps.catalogue', 'lssps.grid',
      ],
      ext_modules=[
          Extension('lssps._lssps',
                    ['../lib/msg.cpp', '../lib/config.cpp', 
                     '../lib/catalogue.cpp', '../lib/grid.cpp',
                     '../lib/mass_assignment.cpp', '../lib/interlacing.cpp',
                     '../lib/multipole.cpp', '../lib/power_spectrum.cpp',
                     'py_package.cpp', 'py_config.cpp', 'py_catalogue.cpp',
                     'py_grid.cpp', 'py_mass_assignment.cpp',
                     'py_power_spectrum.cpp', 'py_interlacing.cpp',
                    ],
                    include_dirs = idirs,
                    extra_compile_args = [os.environ["OPT"]],
                    library_dirs = ldirs,
                    libraries = libs,
          )
      ],
      packages=['lssps'],
)
