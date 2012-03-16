
__all__ = ['NetCDFFile']
__doc__ = """
.. _netcdf
:mod:`netcdf` -- netcdf import point
====================================

.. module:: netcdf
   :platform: Unix, Windows
   :synopsis: Povides a single import point for a package.  If
              a user has one of many netcdf interfaces, this module
              selects it and provides it.
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""
from netCDF4 import Dataset
class NetCDFFile(Dataset):
    def __new__(cls, *args, **kwds):
        kwds.setdefault('format', 'NETCDF3_CLASSIC')
        return Dataset(*args, **kwds)

