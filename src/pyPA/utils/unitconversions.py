toextrinsic = {'ppmV': ('micromoles', 'AIRMOLS'),
               'ppm': ('micromoles', 'AIRMOLS'),
               'ppm/h': ('micromoles/h', 'AIRMOLS'),
               'ppmV/h': ('micromoles/h', 'AIRMOLS'),
               'ppm/hr': ('micromoles/hr', 'AIRMOLS'),
               'ppmV/hr': ('micromoles/hr', 'AIRMOLS'),
               'ugram/m**3': ('micrograms', 'VOL'),
               'ugrams/m**3': ('micrograms', 'VOL'),
               'microgram/m**3': ('micrograms', 'VOL'),
               'micrograms/m**3': ('micrograms', 'VOL'),
               'm**2/m**3': ('m**2', 'VOL'),
               'number/m**3': ('number', 'VOL'),
               'umol/m**3': ('micromoles', 'VOL'),
               'umols/m**3': ('micromoles', 'VOL'),
               'micromol/m**3': ('micromoles', 'VOL'),
               'micromols/m**3': ('micromoles', 'VOL'),
               'micromole/m**3': ('micromoles', 'VOL'),
               'micromoles/m**3': ('micromoles', 'VOL'),
               'm**3/mol': ('m**3', 'AIRMOLS'),
               'm**3/mole': ('m**3', 'AIRMOLS'),
               'm**3/moles': ('m**3', 'AIRMOLS'),
               'm**3': ('m**3', '1'),
               'mol': ('moles', '1'),
               'mole': ('moles', '1'),
               'moles': ('moles', '1'),
               'None': ('None', '1')}

from PseudoNetCDF import PseudoNetCDFVariable
import numpy as np

def getextrinsic(nfile, varkey):
    """
    Given a variable with the property units,
    this function will return an extrinsic property
    """
    var = nfile.variables[varkey]
    unit = getattr(var, 'units', 'None')
    newunit, denominator = toextrinsic[unit.strip()]
    denom = eval(denominator, None, nfile.variables)
    if not isinstance(denom, (float, int)):
        denom = denom[:]
        if (denom.shape[0] + 1) == var[:].shape[0]:
            start = 1
        else:
            start = None
    if start is None:
        ext = var[:] * denom
    else:
        ext = np.zeros_like(var[:, :denom.shape[1], :denom.shape[2], :denom.shape[3]])
        ext[1:] = var[1:, :denom.shape[1]] * denom
        ext[0] = var[0] * denom[0]
    
    extvar = PseudoNetCDFVariable(nfile, varkey, var.dtype.char, var.dimensions, values = ext)
    for k in var.ncattrs():
        setattr(extvar, k, getattr(var, k))
    setattr(extvar, 'units', varkey)
    setattr(extvar, 'oldunits', varkey)
    setattr(extvar, 'expression', 'OLD * %s' % denominator)
    return ext
