from distutils.core import setup, Extension
import numpy as np
import os

print("np.get_include()")
print(np.get_include())

setup(name='lssps',
      version='0.0.1',
      author='Jun Koda',
      py_modules=['lssps.catalogue',
      ],
      ext_modules=[
          Extension('lssps._lssps',
                    ['py_package.cpp', 'py_config.cpp', 'py_catalogue.cpp',
                     'py_grid.cpp',
                    ],
                    include_dirs = ['../lib', np.get_include()],
                    extra_compile_args = [os.environ["OPT"]],
                    library_dirs =  ['../lib'],
                    libraries = ['lssps'],
          )
      ],
      packages=['lssps'],
)
