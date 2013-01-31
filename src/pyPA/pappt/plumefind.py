from ..netcdf import NetCDFFile
from numpy import *
from PseudoNetCDF.camxfiles.Memmaps import uamiv
from ..utils.unitconversions import getextrinsic
import sys

def getlen(dk, dv):
    if dk == 'VAR':
        return 1
    elif dk == 'TSTEP':
        if isinstance(dv, int):
            return dv + 1
        else:
            return len(dv) + 1
    else:
        if isinstance(dv, int):
            return dv
        else:
            return len(dv)

def mass_threshold(nfile, propkey, massthreshhold = .9, tracer = False, volmax = inf):
    conc = eval(propkey, nfile.variables)
    firstkey = propkey.split('+')[0]
    mass = getextrinsic(nfile, propkey)
    if not tracer:
        factor = mass / conc
        conc2 = conc.copy()
        conc2.shape = conc2.shape[:2] + (prod(conc2.shape[2:]),)
        conc = conc[:] - percentile(conc2, [84.13], axis = 2)[0][:, :, None, None]
        conc = where(conc < 0, 0, conc)
        mass = factor * conc
    total_mass = mass.sum(1).sum(1).sum(1)
    vol = nfile.variables['VOL']
    tflag = nfile.variables['TFLAG']
    shape = ones(conc.shape , 'i')
    print >> sys.stdout, "[Day Hr]", "Mass Fraction", "Mass Captured", "Total Mass Available"
    for hri, (hrconc, hrmass, hrtot) in enumerate(zip(conc, mass, total_mass)):
        hrpeak = hrconc.max()
        factor = 0.
        thresh_mass = massthreshhold * hrtot
        capt_mass = thresh_mass * 2
        capt_mass_old = 0
        same_mass_count = 0
        n = 0
        while same_mass_count < 10:
            if capt_mass < thresh_mass:
                factor -= .5**n
            else:
                factor += .5**n
            n += 1
            hrshape = hrconc >= factor*hrpeak
            if capt_mass_old == capt_mass:
                same_mass_count += 1
            else:
                same_mass_count = 0
            capt_mass_old = capt_mass
            capt_mass = hrmass[hrshape].sum()
        
        print >> sys.stdout, tflag[hri,0,:], capt_mass/hrtot, capt_mass, hrtot
        shape[hri] = hrshape
    f = NetCDFFile('massthresh.nc', 'w', format = 'NETCDF3_CLASSIC')
    for dk, dv in nfile.dimensions.iteritems():
        f.createDimension(dk, getlen(dk, dv))
    for pk in nfile.ncattrs():
        setattr(f, pk, getattr(nfile, pk))
    
    newtflag = f.createVariable('TFLAG', 'i', tflag.dimensions)
    newtflag[1:] = tflag[:, [0], :]
    newtflag[0] = tflag[0, [0], :]
    newtflag[0, 0, 1] -= (newtflag[3, 0, 1] - newtflag[2, 0, 1])
    for pk in tflag.ncattrs():
        setattr(newtflag, pk, getattr(tflag, pk))
    
    shapev = f.createVariable(('PLUME%.2f' % massthreshhold).replace('.', 'pt'), 'i', nfile.variables[firstkey].dimensions)
    shapev[1:] = shape[:]
    shapev[0] = shape[0]
    shapev.units = 'ON/OFF'
    shapev.var_desc = 'ON/OFF'
    shapev.long_name = ('PLUME%.2f' % massthreshhold).replace('.', 'pt')
    firstvar = nfile.variables[firstkey]
    propvar = f.createVariable(firstkey, firstvar.dtype.char, firstvar.dimensions)
    for pk in tflag.ncattrs():
        setattr(propvar, pk, getattr(tflag, pk))
    propvar[1:] = conc[:]
    propvar.expression = 'val - 90(%) in horizontal layer'
    f.close()
    

if __name__ == '__main__':
    from pyPA.utils.CMAQTransforms import cmaq_pa_master
    from pyPA.utils.CAMxTransforms import camx_pa_master
    from PseudoNetCDF.MetaNetCDF import file_master
    import yaml
    pafile_master = file_master
    try:
        input = yaml.load(file(sys.argv[1], 'r'))
        pa_master = eval(input.get('metawrapper', 'file_master'))(input['files'])
        mass_threshold(pa_master, sys.argv[2], massthreshhold = .9, volmax = inf, tracer = eval(sys.argv[3]))
    except Exception, e:
        print e
        print
        print
        print '======================================================='
        print 'python -m pyPA.pappt.plumefind <config.yaml> <VARIABLE> <TRACER>'
        print ' - config.yaml: path to a pyPA configuration file.'
        print ' - VARIABLE: string identifying the variable to use.'
        print ' - TRACER: 0 if this is a standard compound; 1 if it is'
        print '           emitted exclusively to trace the plume.'
        print '======================================================='
        print
