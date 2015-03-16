# pyPA Installation #

Known Platforms:
  * Unix
  * Mac
  * Linux

Requirements:
  * [Python >=2.5](http://python.org)
  * [NumPy >=1.2](http://numpy.scipy.org)
  * [YAML](http://www.yaml.org)
  * One of the below netcdf readers (in descending order of preference):
    * [netCDF4](http://code.google.com/p/netcdf4-python)
    * [netCDF3](http://code.google.com/p/netcdf4-python)
    * [pynetcdf](http://pypi.python.org/pypi/pynetcdf/0.7)
    * [Scientific](http://dirac.cnrs-orleans.fr/plone/software/scientificpython/)
    * [pupynere](http://pypi.python.org/pypi/pupynere/)

Instructions
  1. [Download current source distribution](https://dawes.sph.unc.edu/trac/pyPA/changeset/HEAD/trunk?old_path=%2F&format=zip&filename=pyPA-src)
  1. unzip pyPA-`*`.zip -d pyPA-src
  1. cd pyPA-src/trunk
  1. python setup.py install