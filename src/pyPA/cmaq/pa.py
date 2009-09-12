from pyPA.netcdf import NetCDFFile
from PseudoNetCDF.sci_var import PseudoNetCDFFile
from PseudoNetCDF.MetaNetCDF import file_master
from warnings import warn

def mrgidx(ipr_paths, irr_paths, idx, conc_paths = None):
    """
    ipr_paths - list of strings or a string for IPR (aka PA) paths
    irr_paths - list of strings or a string for IRR paths
    conc_paths - list of strings or a string for CONC paths
    idx - 4-D index usually a time slice with a 3-D location 
          (e.g. [slice(None), 3, 4, 5])
    """
    if isinstance(irr_paths,str):
        irrf = NetCDFFile(irr_paths)
    else:
        irrf = file_master([NetCDFFile(irr_path) for irr_path in irr_paths])
    
    if isinstance(ipr_paths,str):
        iprf = NetCDFFile(ipr_paths)
    else:
        iprf = file_master([NetCDFFile(ipr_path) for ipr_path in ipr_paths])
        
    if conc_paths is None:
        concf = None
        concidx = None
    else:
        if isinstance(conc_paths,str):
            concf = NetCDFFile(conc_paths)
        else:
            concf = file_master([NetCDFFile(conc_path) for conc_path in conc_paths])
        xoffset = (concf.XORIG -  iprf.XORIG) / iprf.XCELL
        assert(float(xoffset) == float(int(xoffset)))
        xoffset = int(xoffset)

        yoffset = (concf.YORIG -  iprf.YORIG) / iprf.YCELL
        assert(float(yoffset) == float(int(yoffset)))
        yoffset = int(yoffset)
        concidx = [i for i in idx]
        concidx[2] += yoffset
        concidx[3] += xoffset
    
        
    
    # Process and Reaction keys should exclude TFLAG
    pr_keys = [pr for pr in iprf.variables.keys() if pr not in ('TFLAG',)]
    rr_keys = [rr for rr in irrf.variables.keys() if rr not in ('TFLAG',)]
    try:
        conc_keys = [c for c in concf.variables.keys() if c not in ('TFLAG',)]
    except:
        conc_keys = []
    
    # Attempt to order reactions by number
    # this is not necessary, but is nice and clean
    try:
        rr_keys = [(int(rr.split('_')[1]), rr) for rr in rr_keys]
        rr_keys.sort()
        rr_keys = [rr[1] for rr in rr_keys]
    except:
        warn("Cannot sort reaction keys")
    
    # Processes are predicated by a delimiter
    prcs = ['INIT']+list(set([pr.split('_')[0] for pr in pr_keys]))+['FCONC']
    # Species are preceded by a delimiter
    spcs = list(set(['_'.join(pr.split('_')[1:]).replace(' ', '_') for pr in pr_keys]))
    
    # Select a dummy variable for extracting properties
    pr_tmp = iprf.variables[pr_keys[0]]
    
    # Create an empty file and decorate
    # it as necessary
    outf = PseudoNetCDFFile()
    outf.Species = "".join([spc.ljust(16) for spc in spcs])
    outf.Process = "".join([prc.ljust(16) for prc in prcs])
    outf.Reactions = "".join([rr_key.ljust(16) for rr_key in rr_keys])
    outf.createDimension("PROCESS", len(prcs))
    outf.createDimension("SPECIES", len(spcs))
    outf.createDimension("RXN", len(rr_keys))
    outf.createDimension("TSTEP", pr_tmp[:,0,0,0].shape[0])
    outf.createDimension("TSTEP_STAG", outf.dimensions["TSTEP"]+1)
    outf.createDimension("ROW", 1)
    outf.createDimension("LAY", 1)
    outf.createDimension("COL", 1)
    outf.createDimension("VAR", 3)
    outf.createDimension("DATE-TIME", 2)
    tflag = outf.createVariable("TFLAG", "i", ('TSTEP', 'VAR', 'DATE-TIME'))
    tflag.__dict__.update(dict(units = "<YYYYJJJ,HHDDMM>", var_desc = 'TFLAG'.ljust(16), long_name = 'TFLAG'.ljust(16)))
    tflag[:,:,:] = iprf.variables['TFLAG'][:][:,[0],:]
    shape = outf.createVariable("SHAPE", "i", ("TSTEP", "LAY", "ROW", "COL"))
    shape.__dict__.update(dict(units = "ON/OFF", var_desc = "SHAPE".ljust(16), long_name = "SHAPE".ljust(16)))
    shape[:] = 1
    irr = outf.createVariable("IRR", "f", ("TSTEP", "RXN"))
    irr.__dict__.update(dict(units = pr_tmp.units, var_desc = "IRR".ljust(16), long_name = "IRR".ljust(16)))
    ipr = outf.createVariable("IPR", "f", ("TSTEP", "SPECIES", "PROCESS"))
    ipr.__dict__.update(dict(units = pr_tmp.units, var_desc = "IPR".ljust(16), long_name = "IPR".ljust(16)))

    for rr, var in zip(rr_keys,irr.swapaxes(0,1)):
        var[:] = irrf.variables[rr][:][idx]
        
    for prc, prcvar in zip(prcs,ipr.swapaxes(0,2)[1:-1]):
        if not prc in ('INIT', 'FCONC'):
            for spc, spcvar in zip(spcs,prcvar):
                try:
                    spcvar[:] = iprf.variables['_'.join([prc,spc])][:][idx]
                except KeyError, es:
                    warn(str(es))

    if concf is not None:
        for prc_slice, prcvar in zip([slice(0,-1), slice(1,None)], ipr.swapaxes(0,2)[[1,-1]]):
            for spc, spcvar in zip(spcs,prcvar):
                if concf.variables.has_key(spc):
                    spcvar[:] = concf.variables[spc][:][concidx][prc_slice]
                else:
                    warn("No concentration given for %s" % spc)
            
    return outf
    
if __name__ == '__main__':
    from PseudoNetCDF.pncdump import pncdump
    from pyPA.cmaq.pa import mrgidx
    import sys
    
    x = ', '.join(sys.argv[1:])
    try:
        pncdump(eval('mrgidx(%s)' % (x,)), name = mrgidx)
    except:
        print >> sys.stderr, "\nUsage: python -m pyPA.cmaq.pa \"ipr_paths, irr_paths, idx, conc_paths\"\n  ipr_paths - a single string or list of strings\n  irr_paths - a single string or list of strings\n  idx - 4D numpy slice indexes a time series\n  conc_paths - a single string or list of strings\n\n\nExample:\n  $ python -m pyPA.cmaqfiles.pa \"['CMAQ_IPR1.nc','CMAQ_IPR2.nc'], 'CMAQ_IRR.nc', (slice(None),0,1,1)\""
        raise
