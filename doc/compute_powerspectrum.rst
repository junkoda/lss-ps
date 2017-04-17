Compute power spectrum
======================

.. py:function:: compute_power_spectrum(data, rand, nc, \
		 esimator='plane-parallel',\
		 mas='CIC',\
                 x0=(0,0,0), boxsize=None,\
		 los=2,
                 k_min=0.0, k_max=1.0, dk=0.01, nmu=0,\
                 subtract_shotnoise=True,\
                 correct_mas=True,\
                 interlacing=False)

   :param CatalogueFile data: The data catalogue
   :param rand: The random catalogue
   :type rand: CatalogueFile or None
   :param int nc: The number of grid pointsThe body of the message
   :param str estimator: The power spectrum estimator: 'plane-parallel', 'Yamamoto-Scoccimarro', or 'Yamamoto-Bianchi'		  
   :param str mas: The Mass assignment scheme: 'NGP', 'CIC' or 'TSC'
   :param sequience x0: The coordinate of the box corner
   :type x0: tuple(float, float, float)
   :param float boxsize: The length of the box on a side
   :param int los: The line-of-sight direction 0, 1, or 2 for plane-parallel estimator		 
   :param float k_min: The lower range of the power spectrum k bins
   :param float k_max: The upper range of the power spectrum k bins
   :param float dk:    The bin width of the power spectrum k bins
   :param int nmu:     The number of mu bins for 2D power spectrum
   :param bool coorect: Correct for the mass assignment smoothing effect
   :param bool interlacing: Reduce aliasing by using interlacing
   :return: PowerSpectrum object
   :rtype: PowerSpectrum
   :raises FileNotFoundError: if the catalogue files are not found


Examples
========

Partiles in a periodic box
--------------------------

.. code-block:: python

   import lssps

   data = lssps.catalogue.ascii('data/wizcola_realspace.txt', xyz=(0, 1, 2))
   rand = None
   nc = 256

   ps = lssps.compute_power_spectrum(data, rand, nc, 'plane-parallel',
                                     boxsize=600)

   for k, nmodes, P0, P2, P4 in ps:
       if nmodes > 0.0:
           print('%e %e %e' % (k, nmodes, P0))
