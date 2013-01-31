import sys
from warnings import warn
from os.path import exists
import re

from yaml import load

from PseudoNetCDF.MetaNetCDF import file_master
from PseudoNetCDF import PseudoNetCDFVariable

from pyPA.utils.CMAQTransforms import cmaq_pa_master
from pyPA.utils.CAMxTransforms import camx_pa_master
from pyPA.pappt.lagrangian import boxes, box_id
from pyPA.pappt.pseudo_procs import simple_pseudo_procs
from pyPA.netcdf import NetCDFFile


def ext_mrg(input):
    from numpy import ndarray, newaxis, where, fromfile, zeros, ones
    from numpy.ma import masked_where, masked_invalid
    from PseudoNetCDF import PseudoNetCDFFile, PseudoNetCDFVariable
    from PseudoNetCDF.pncgen import pncgen

    initial = input.get('initial', 'INIT')
    final = input.get('final', 'FCONC')
    
    # Open a single PseudoNetCDF file that can access
    # all variables necessary from input files
    pa_master = eval(input.get('metawrapper', 'file_master'))(input['files'])
    
    # Analysis volume shape variable name
    shape_name = input.get('shape', 'DEFAULT_SHAPE')

    # Variable dimension order
    #
    # Provide the order of spatiotemporal dimensions for variables
    # *assuming* consistent for all variables
    #
    # For convenience dimensions are also ordered for later use
    dimensions = input.get('dimensions', dict(TSTEP=0, LAY = 1, ROW = 2, COL = 3))    
    spatial_dimensions_ordinals = [v for k,v in dimensions.iteritems() if k != 'TSTEP']
    spatial_dimensions_ordinals.sort(reverse = True)
    dimensions_ordered = [(v, k) for k,v in dimensions.iteritems()]
    dimensions_ordered.sort(reverse = False)
    dimensions_ordered = [k for v,k in dimensions.iteritems()]


    # Unit conversion dictionary
    #
    # each element of the dictionary has a key that represents the input unit
    # and a dictionary value.  The dictionary value should have two keys: expression
    # and new_unit.  expression is a textual expression to perform on the values which
    # should be entered as <values>
    #
    # *Providing default conversion of umol/m**3 values to ppm for CAMx
    # **<values> is replaced by var[:] which evaluates in the global environment
    #   to the current values
    unitconversions = {'umol/m**3': dict(expression = "<values> * VOL / AIRMOLS", new_unit = 'ppm')}
    unitconversions.update(input.get('unitconversions', {}))
    for unit, unit_dict in unitconversions.iteritems():
        unitconversions[unit]['expression'] = unitconversions[unit]['expression'].replace('<value>', 'var[:]')

    # intrinsic to extrinsic conversion dictionary
    #
    # each element has a key for the input unit that it converts
    # the value is the name of a variable that can be multiplied to
    # an intrinic unit to obtain an extrinsic unit
    contributions = {'ppmV': 'AIRMOLS',
                     'ppm': 'AIRMOLS',
                     'ppm/h': 'AIRMOLS',
                     'ppmV/h': 'AIRMOLS',
                     'ppm/hr': 'AIRMOLS',
                     'ppmV/hr': 'AIRMOLS',
                     'ugram/m**3': 'VOL',
                     'ugrams/m**3': 'VOL',
                     'microgram/m**3': 'VOL',
                     'micrograms/m**3': 'VOL',
                     'm**2/m**3': 'VOL',
                     'number/m**3': 'VOL',
                     'umol/m**3': 'VOL',
                     'umols/m**3': 'VOL',
                     'micromol/m**3': 'VOL',
                     'micromols/m**3': 'VOL',
                     'micromole/m**3': 'VOL',
                     'micromoles/m**3': 'VOL',
                     'm**3/mol': 'AIRMOLS',
                     'm**3/mole': 'AIRMOLS',
                     'm**3/moles': 'AIRMOLS',
                     'm**3': '1',
                     'mol': '1',
                     'mole': '1',
                     'moles': '1',
                     'None': '1'}
                     
    contributions.update(input.get('contributions', {}))
    
    # extrinsic to intrinsic conversion dictionary
    #
    # each element has a key for the extrinsic unit that it converts
    # the value is the name of a variable that can divide a extrinsic
    # unit to return an intrinsic unit
    # *by default all contribution variables are also normalizers
    normalizers = contributions.copy()
    normalizers.update(input.get('normalizers', {}))


    # analysis volume shape
    #
    # shape has ones (on) and zeros (off) to indicate
    # which cells in the process analysis volume will
    # used in the analysis
    #
    # if an ascii mask is provided, it further restricts
    # the cells used following the same on/off convection
    # the ascii mask is an ascii map of the domain cells
    # delimmited by spaces.
    try:
        shape = eval(shape_name, locals(), pa_master.variables)[:]
    except:
        raise UserWarning, '%s is not valid in your files as processed by the metawrapper' % shape_name

    if input.has_key('ascii_mask'):
        ascii_mask = input['ascii_mask']
        if exists(ascii_mask):
            cols = len(file(ascii_mask, 'r').readline().split(" "))
            rows = len([l for l in file(ascii_mask, 'r').readlines() if l != '\n'])
            
            ascii_mask = fromfile(ascii_mask, sep = " ", dtype = 'bool').reshape(1, 1, rows, cols)[:, :, ::-1, :]
            shape = shape * ascii_mask
            
        else:
            print >> file(ascii_mask, 'w'), '\n'.join([' '.join(['1']*shape.shape[3])]*shape.shape[2])
            raise ValueError, """File %s was not found; instead, a template was created.
    Edit the template to include only those cells of interest"""


    def reduce_space(var):
        """
        Convenience function that removes the spatial dimensions
        by summing across them.  This should be used for extrinsic
        variables to get the total volume value.  Dimension order is
        taken from global environment
        
        example: mass_o3 might have 4 dimensions (TIME, LAYERS, ROWS, COLS)
                 reduce(mass_o3) returns total mass_o3 with only a time dimension
        """
        out = var.copy()
        for d in spatial_dimensions_ordinals:
            out = out.sum(d)
        return out
                
    mask = shape.astype('bool') == False
    
    species = input.get('species', None)
    species = species or [key[len(initial) + 1:] for key in pa_master.variables.keys() if key[:len(initial) + 1] == initial + '_']
    species = list(set(species))
    
    temp_spc = species[0]
    processes = input.get('processes', None)
    processes = processes or [key[:-len(temp_spc)-1] for key in pa_master.variables.keys() if '_' + temp_spc == key[-len(temp_spc)-1:]]
    processes = list(set(processes))

    reactions = input.get('reactions', None)
    reactions = reactions or '(RXN|IRR)_.*'
    if isinstance(reactions, str):
        reactions = re.compile(reactions)
        reactions = [key for key in pa_master.variables.keys() if reactions.match(key)]
        ordre = re.compile('\d+')
        try:
            reactions = [(int(ordre.search(key).group()), key) for key in reactions]
            reactions.sort()
            reactions = [key for ord, key in reactions]
        except:
            pass
    else:
        reactions = list(set(reactions))
        
    
    agg_keys = [(shape_name, 's')]
    agg_keys.extend([(k, 'a') for k in list(set([v for k, v in contributions.iteritems()]+[v for k, v in normalizers.iteritems()]))])
    agg_keys.extend([(rxn, 'r') for rxn in reactions])
    
    for spc in species:
        for prc in processes:
            ptype = {initial: 'i', final: 'f'}.get(prc, 'p')
            agg_keys.append(('_'.join([prc, spc]), ptype))
            
    for k, v in contributions.iteritems():
        try:
            unit_contribution = eval(v, locals(), pa_master.variables)
        except:
            warn('pyPA could not evaluate %s with the files you provided; if you need %s, make sure you provided all necessary inputs' % (v, v))
            continue
            
        if isinstance(unit_contribution, ndarray):
            unit_contribution = unit_contribution[:]
            
        contributions[k] = unit_contribution
    
    bxs=boxes(shape,kaxis = dimensions['LAY'])
    norm_bxs = {}
    for unit, v in normalizers.iteritems():
        try:
            unit_normalizer = eval(v, locals(), pa_master.variables)
        except:
            warn('pyPA could not evaluate %s with the files you provided; if you need %s, make sure you provided all necessary inputs' % (v, v))
            continue
            
        if isinstance(unit_normalizer, ndarray):
            unit_normalizer = unit_normalizer[:]
            norm_bxs[unit] = reduce_space(bxs * unit_normalizer[..., newaxis])
        else:
            norm_bxs[unit] = ones(bxs[:, 0, 0, 0, :].shape) * unit_normalizer
            
        normalizers[unit] = unit_normalizer
    
    denominators = {}
    init_denominators = {}

    for k, v in normalizers.iteritems():
        if isinstance(v, ndarray):
            denominators[k] = reduce_space(masked_where(mask[1:], v))
            init_denominators[k] = reduce_space(masked_where(mask[:-1], v))
        else:
            denominators[k] = 1
            init_denominators[k] = 1

        
    outputfile = PseudoNetCDFFile()
    
    def getdimlen(dim):
        try:
            return len(dim)
        except TypeError:
            return dim

    outputfile.createDimension('TSTEP', getdimlen(pa_master.dimensions['TSTEP']))
    outputfile.createDimension('TSTEP_STAG', getdimlen(outputfile.dimensions['TSTEP'])+1)
        
    outputfile.createDimension('LAY', mask.shape[dimensions['LAY']])
    outputfile.createDimension('ROW', mask.shape[dimensions['ROW']])
    outputfile.createDimension('COL', mask.shape[dimensions['COL']])
    outputfile.createDimension('VAR', len(agg_keys)+1)
    outputfile.createDimension('DATE-TIME', 2)
    
    var = pa_master.variables['TFLAG']
    outputfile.variables['TFLAG'] = PseudoNetCDFVariable(outputfile, 'TFLAG', 'i', ('TSTEP', 'VAR', 'DATE-TIME'), units = var.units, long_name = var.long_name, var_desc = var.var_desc, values = var[:][:, [0], :].repeat(len(agg_keys)+1, 1))
    
    for key, ktype in agg_keys:
        print >> sys.stdout, key, 
        dimensions = ('TSTEP',)
        print key
        try:
            var = pa_master.variables[key]
        except KeyError, (e):
            if ktype in ('p', 'i', 'f'):
                warn("No %s process variable" % key)
                var = PseudoNetCDFVariable(pa_master, 'temp', 'f', dimensions_ordered, units = 'None', long_name = key, var_desc = "Dummy values (0) for %s" % key, values = zeros(mask.shape, 'f')[1:])
            elif ktype == 'a':
                continue
            else:
                raise KeyError, "No %s variable" % key

        unit = var.units.strip()
        if unit not in contributions and ktype not in ('a', 's'):
            warn("Ignoring %s; cannot process unit (%s)" % (key, unit))
            warn("To add processing for a unit, update the contributions and/or normalizations dictionary with appropriate variable")
            continue
        mask_slice = slice(1, None)
        
        if ktype in ('a',):
            print >> sys.stdout, key
            numerator = reduce_space(masked_where(mask[mask_slice], var[:]))
            denominator = 1
        elif ktype in ('s',):
            numerator = shape
            dimensions = ('TSTEP_STAG', 'LAY', 'ROW', 'COL')
            denominator = 1
        else:
            if ktype in ('i',):
                mask_slice = slice(0, -1)
                denominator = init_denominators[unit]
            else:
                denominator = denominators[unit]
            numerator = reduce_space(masked_where(mask[mask_slice], var[:]*contributions[unit]))
        
        values = numerator / denominator
        outputfile.variables[key] = PseudoNetCDFVariable(outputfile, key, 'f', dimensions, values = values, units = unit, long_name = var.long_name, var_desc = var.var_desc)
    print >> sys.stdout
    
    if reduce_space(bxs).sum(0)[[box_id.HENT, box_id.HDET, box_id.VENT, box_id.VDET]].astype('bool').any():
        simple_pseudo_procs(pa_master = pa_master, outputfile = outputfile, spcs = species, initial = initial, bxs = bxs, norm_bxs = norm_bxs, contributions = contributions, reduce_space = reduce_space)
        processes.extend("VDET VENT HDET HENT EDHDIL EDVDIL".split())
    
    # For each variable, change the unit until there
    # is no change left.
    for name, var in outputfile.variables.iteritems():
        in_unit = var.units.strip()
        while in_unit in unitconversions:
            out_unit = unitconversions[in_unit]['new_unit'].strip()
            expression_str = unitconversions[in_unit]['expression'].replace('<values>', 'var[:]')
            var[:] = eval(expression_str, locals(), outputfile.variables)
            var.units = out_unit.ljust(16)
            in_unit = out_unit
    
    outputfile.Processes = '\t'.join([p.ljust(16) for p in processes])
    outputfile.Species = '\t'.join([p.ljust(16) for p in species])
    outputfile.Reactions = '\t'.join([p.ljust(16) for p in reactions])
    outputfile.PYPAVERSION = '1'
    from pyPA.netcdf import NetCDFFile
    outputfile = pncgen(outputfile, NetCDFFile(input['outfile'], mode = 'w', format = 'NETCDF3_CLASSIC'))
    outputfile.sync()
    return outputfile

if __name__ == '__main__':
    input = load("""
outfile: test.mrg.nc
metawrapper: cmaq_pa_master
files:
  - [/Users/barronh/Development/CMAQ/v4.7/data/cctm/CCTM_e1a_Darwin9_i386_PA_1.ncf, NetCDFFile]
  - [/Users/barronh/Development/CMAQ/v4.7/data/cctm/CCTM_e1a_Darwin9_i386_PA_2.ncf, NetCDFFile]
  - [/Users/barronh/Development/CMAQ/v4.7/data/cctm/CCTM_e1a_Darwin9_i386_PA_3.ncf, NetCDFFile]
  - [/Users/barronh/Development/CMAQ/v4.7/data/cctm/CCTM_e1a_Darwin9_i386_IRR_1.ncf, NetCDFFile]
  - [/Users/barronh/Development/CMAQ/v4.7/data/cctm/CCTM_e1a_Darwin9_i386_CONC.ncf, NetCDFFile]
  - [/Users/barronh/Development/CMAQ/v4.7/data/cctm/../mcip3/M_36_2001/METCRO2D_010722, NetCDFFile]
  - [/Users/barronh/Development/CMAQ/v4.7/data/cctm/../mcip3/M_36_2001/METCRO3D_010722, NetCDFFile]
initial: 'INIT'
final: 'FCONC'
species: [NH3, ANH4I] #, ANH4J, ANH4K] #O3, 'NO', NO2, ALD2]
processes: [AERO, CHEM, CLDS, DDEP, EMIS, HADV, HDIF, VDIF, ZADV, INIT, FCONC]
reactions: ['IRR_1']
shape: DEFAULT_SHAPE
#ascii_mask: ascii_mask.txt
unitconversions:
    ppmV:
        expression: <value>*1000.
        new_unit: ppb
    ppm:
        expression: <value>*1000.
        new_unit: ppb
""")
    ext_mrg(input)
    
