__all__ = ['NetCDFFile']
from warnings import warn
netcdf_pkgs = [('pynetcdf', 'NetCDFFile'), ('netCDF3', 'Dataset'), \
               ('netCDF4', 'Dataset'), ('Scientific.IO.NetCDF', 'NetCDFFile'), \
               ('pupynere', 'NetCDFFile')]
for pkg, reader in netcdf_pkgs:
    try:
        NetCDFFile = getattr(__import__(pkg, fromlist = [reader]),reader)
        break
    except ImportError, e:
        warn(e.message)
        pass
else:
    raise ImportError, "Did not find a NetCDFFile object"
