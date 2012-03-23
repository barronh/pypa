from netcdf import NetCDFFile
import numpy as np
def combine(paths):
    files = map(NetCDFFile, paths)
    temp_file = files[0]
    outputfile = NetCDFFile('CAMx_PA_MASTER.nc', 'w', format = 'NETCDF4_CLASSIC')
    nstep = sum([len(f.dimensions['TSTEP']) for f in files])
    nstep_stag = nstep + 1
    outputfile.createDimension('TSTEP', nstep)
    outputfile.createDimension('TSTEP_STAG', nstep_stag)
    outputfile.createDimension('LAY', len(temp_file.dimensions['LAY']))
    outputfile.createDimension('ROW', len(temp_file.dimensions['ROW']))
    outputfile.createDimension('COL', len(temp_file.dimensions['COL']))
    outputfile.createDimension('VAR', len(temp_file.dimensions['VAR']))
    outputfile.createDimension('DATE-TIME', 2)

    for k in temp_file.ncattrs():
        setattr(outputfile, k, getattr(temp_file, k))

    for key, temp_var in temp_file.variables.iteritems():
        print "Adding", key
        new_var = outputfile.createVariable(key, temp_var.dtype.char, temp_var.dimensions, contiguous = True)
        for k in temp_var.ncattrs():
            setattr(new_var, k, getattr(temp_var, k))
    
        print "Populating", key
        if 'TSTEP_STAG' in temp_var.dimensions:
            new_var[:] = np.concatenate([f.variables[key][:-1] for f in files] + [files[-1].variables[key][-1:]], axis = 0)
        else:
            new_var[:] = np.concatenate([f.variables[key][:] for f in files], axis = 0)
        outputfile.sync()
    outputfile.close()

if __name__ == '__main__':
    import sys
    paths = sys.argv[1:]
    combine(paths)