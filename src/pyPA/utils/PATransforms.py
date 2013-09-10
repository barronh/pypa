from pyPA.netcdf import NetCDFFile
from PseudoNetCDF.MetaNetCDF import file_master

def pafile_master(paths_and_readers):
    for i, (p, r) in enumerate(paths_and_readers):
        if r[:5] == 'from ' and 'import' in r:
            exec r in globals(), locals()
            paths_and_readers[i][1] = r.split(' ')[-1]
    # Create a list of opened files
    files=[eval(r)(p) for p,r in paths_and_readers]
    
    # Create a file master object from files
    return file_master(files)
