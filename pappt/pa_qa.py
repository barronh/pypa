HeadURL="$HeadURL$"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy$"
__version__ = RevisionNum

import sys
from numpy import array, where, isnan, zeros, logical_and
from pynetcdf import NetCDFFile as ncf
from pyPA.utils.CAMxFiles import ipr,ipr_mem

# Used to delimit species and process names
delim="_"

def PAdC(pa_file,spc_list=None,prc_list=None,model=None,mechanism=None,init='INIT',final='FCONC',verbose=True):
    """PAdC
    Process Analysis tallies the change in concentration (conc.) due to each 
    process in the model.  The sum of these process changes should be
    equal to the change in concentration....
    
    If the sum of process changes is dC, Initial conc. + dC should equal Final conc.
    This function checks that for a list of species and processes, which can optionally
    be supplied by the program for a given model/mechanism.
    
    pa_file - NetCDF like object (i.e. has a variables dictionary with array like values)
    spc_list - optional list of species to check
    prc_list - optional list of processes to check
    model - optional name of model cmaq or camx
    mechanism - optional mechanism cbivmech3
    init - name of initial concentration "process"
    final - name of final concentration "process"
    verbose - display to screen errors and values
    
    Either spc_list/prc_list or model/mechanism must be provided
    """
    
    # if model or mechanism is not provided, spc_list and prc_list must be
    if model==None or mechanism==None:
        if spc_list==None or prc_list==None:
            raise ValueError, "spc_list and prc_list are only optional if model and mechanism are provided"
    
    # if not provided get model and mechanism defaults
    if spc_list==None:
        spc_list=yaml.load(file(os.path.join(os.path.dirname(__file__),model.lower()+'_'+mechanism.lower()+'.yaml'),'r')).species
        
    # if not provided get model defaults
    if prc_list==None:
        prc_list=yaml.load(file(os.path.join(os.path.dirname(__file__),model.lower()+'_'+mechanism.lower()+'.yaml'),'r')).processes
    
    prc_list=[prc for prc in prc_list if prc not in [init,final,'FCONC','INIT']]
        
    
    # Accumulate stats for output
    spc_min_mean_nmean_max=[]
    
    # Initialize the change accumulator with the shape from a variable
    PAdC=zeros(pa_file.variables[prc_list[0]+delim+spc_list[0]].shape,'f')

    # Iterate over species before processes; this order is best for CAMx and for CMAQ it doesn't matter
    for spc in spc_list:
        # Initialize PAdC for each variable
        PAdC[...]=0.
        
        # Iterate over processes
        for prc in prc_list:
            # Combine process and species into a key
            key=prc+delim+spc
            # If the key is not available (e.g. emissions of O3), skip it
            if key in pa_file.variables.keys():
                # Add individual process to change accumulator
                PAdC+=pa_file.variables[key]
                if verbose:
	                print >> sys.stderr,"Found %s for %s" % (prc,spc)
            else:
                # If verbose, report skip to user
                if verbose:
                    print >> sys.stderr,prc, "does not exist for", spc
        if PAdC.sum()==0:
        	pritn >> sys.stderr,prc,"has no PAdC"
        # Calculate the change based on initial and final values
        ConcdC=array(pa_file.variables[final+delim+spc])-array(pa_file.variables[init+delim+spc])
        
        # Calculate the fraction of change that PA captured
        dc_frac=PAdC/where(ConcdC==0,PAdC/2.,ConcdC)
        
        # When Final-Init and dC equal 0, nans show up and must be removed
        dc_frac=where(isnan(dc_frac),1,dc_frac)

        count_within=logical_and(dc_frac>.99,dc_frac<1.01).sum()
        pct_within=count_within/float(dc_frac.size)
        
        # Calculate a mean fraction by weighting the actual change in concentration
        ndc_frac_mean=(ConcdC*dc_frac).sum()/ConcdC.sum()
        
        # Accumulate min, mean, normalized mean, and max for returning and reporting
        spc_min_mean_nmean_max.append((spc,dc_frac.min(),dc_frac.mean(),ndc_frac_mean,dc_frac.max(),count_within,pct_within))
        
    # If verbose, report stats to user
    if verbose:
        print >> sys.stdout, "(Spc, Min, Mean, Norm Mean, Max, Count within .99-1.01, %Total)"
        for l in spc_min_mean_nmean_max:
            print >> sys.stdout, "%s, %8.5f, %8.5f, %8.5f, %8.5f, %d, %5.2f" % l
    
    # Return result to user
    return spc_min_mean_nmean_max,PAdC,ConcdC

if __name__=='__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-m", "--model", dest="model",default="camx",
                      help="Model can either be camx or cmaq", metavar="MODEL")

    parser.add_option("-c", "--mechanism", dest="mechanism", default="cbivmech3",
                      help="Chemical mechanisms: cbivmech3", metavar="MECHANISM")
                      
    parser.add_option("-s", "--species", dest="species",
                      help="List of species", metavar="SPECIES")
                      
    parser.add_option("-P", "--process", dest="process",
                      help="List of process", metavar="PROCESS")
                      
    parser.add_option("-i", "--initial", dest="initial", default="INIT")
    parser.add_option("-f", "--final", dest="final", default="FCONC")
    
    (options, args) = parser.parse_args()
    if len(args)!=1 or \
        options.model.lower() not in ['camx','cmaq'] or \
        options.mechanism.lower() not in ['cbivmech3']:
        parser.print_help()
        sys.exit()
        
    iprfile=args[0]
    
    try:
        pa_file=ncf(iprfile,'r+')
    except IOError:
        try:
            pa_file=ipr(iprfile)
        except OverflowError:
            pa_file=ipr_mem(iprfile)

    PAdC(pa_file,spc_list=options.species,prc_list=options.process,model=options.model,mechanism=options.mechanism,verbose=True,init=options.initial,final=options.final)