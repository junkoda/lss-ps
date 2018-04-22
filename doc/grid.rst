Grid
======================

.. py:class:: Grid

   .. py:attribute:: nc

   Number of grid points per dimension

   .. py:attribute:: boxsize

   Length of the cubic box on a side

   .. py:attribute:: x0
   
   Coordinate of the corner of the box
		     
   .. py:attribute:: offset (float)

   Position of the grid points in units of grid spacing.
   Grid points are located at x_i = (n + offset) \Delta x
   where \Delta x = boxsize / nc, n = 0, 1, 2, .., nc - 1

   .. py:attribute:: mode

   'real-space' or 'fourier-space'

   .. py:attribute:: pk_normalisation

   P(k) = pk_normalisation <F(k) F(k)^*>

   Set to
   1/grid_rand.nw2_sum
   in compute_fluctuation()

   .. py::atribute:: shot_noise
   
		     
   .. py:method:: clear()

   Clear the grid value to 0, and reset mode to 'real-space'.

   .. py:method:: compute_fluctuation(grid_rand = None)

   Compute the fluctuation

   grid <- grid - alpha n_{rand}
   
   Uses the homogeneous mean density $$ \\bar{n} = N_{data}/boxsize^3 $$ if
   grid_rand = None.
   grid <- grid - \bar{n}
		     
   .. py:method:: interlace()

   .. py:method:: correct_mas()


.. py:function:: lssps.grid.zeros(nc, boxsize, x0, offset=0.0, *, interlacing=True)

.. py:function:: lssps.grid.empty(nc, boxsize, x0, offset=0.0, *, interlacing=True)
   
Examples
