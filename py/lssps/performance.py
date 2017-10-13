import lssps._lssps as c
import numpy as np

def xyz_sum(xyz):
    """
    Args:
      array with 3 columns
    Returns:
      sums (Tuple) of xyz
      None if the length of the array is zero
    """

    return c._performance_xyz_sum(xyz)
