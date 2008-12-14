#!/usr/bin/env python
from pyPA.utils.MetaNetCDF import add_derived,file_master,window,newresolution
from pyPA.utils.sci_var import PseudoNetCDFFile, PseudoNetCDFVariable, PseudoNetCDFVariables, Pseudo2NetCDF
from pyPA.utils.util import AttrDict

from pyPA.pappt.pappt import ext_mrg,MergedWriter
from pyPA.utils.CAMxFiles import *
from pyPA.utils.CMAQTransforms import *
from pyPA.utils.CAMxTransforms import *
from pyPA.pappt.kvextract import tops2shape,vertcamx
from pyPA.pappt.legacy import LegacyMerged,LegacyMergedCMAQ
from pyPA.utils.ArrayTransforms import CenterTime
from numpy import ones,zeros,array
try:
    from Scientific.IO.NetCDF import NetCDFFile as ncf
except:
    from pynetcdf import NetCDFFile as ncf
import os,sys,yaml

__all__ = ['camxshapemaker', 'LoadPAQAFromYAML', 'LoadPyPAFromYAML', 'template']

def camxshapemaker(iprfile,hpfile=None,kvfile=None,outpath='shape.nc',pagrid=0):
    if not ((hpfile!=None and kvfile!=None) or pblfile!=None):
        raise ValueError, "User must provide height_pressure/vertical_diffusivity or pbl"
    try:
        iprfile=ipr(iprfile)
    except OverflowError:
        iprfile=ipr_mem(iprfile)
    padomain=iprfile.padomain[pagrid]
    grid=iprfile.grid[padomain['grid']]
    nrows=grid['nrow']
    ncols=grid['ncol']
    
    metfiles=[]
    if hpfile!=None: 
        metfiles.append(CenterTime(height_pressure(hpfile,nrows,ncols)))
    if kvfile!=None: 
        metfiles.append(CenterTime(vertical_diffusivity(kvfile,nrows,ncols)))

    metfiles=file_master(metfiles)
    nlays=metfiles.dimensions['LAY']
    ntimes=metfiles.dimensions['TSTEP']
    iprtflag=iprfile.variables['TFLAG'][:,0,:].tolist()
    iprtflagstart=iprtflag[0]
    mettflag=metfiles.variables['TFLAG'][:,0,:].tolist()
    tstart=mettflag.index(iprtflagstart)
    tend=len(iprtflag)
    
    tslice=slice(tstart,tend)
    islice=slice(padomain['istart']-1,padomain['iend'])
    jslice=slice(padomain['jstart']-1,padomain['jend'])
    kslice=slice(padomain['blay']-1,padomain['tlay'])
    metfiles=window(file_master(metfiles),tslice=tslice,kslice=kslice,jslice=jslice,islice=islice)
    shape=tops2shape(vertcamx(metfiles.variables['KV'],metfiles.variables['HGHT']),metfiles.variables['HGHT'].shape)
    outfile=ncf(outpath,'r+')
    for d in ['LAY','ROW','COL']:
        outfile.createDimension(d,ipr.dimensions[d])
    outfile.createDimension('TSTEP',ipr.dimensions['TSTEP']+1)
    new_shape=[i for i in shape.shape]
    new_shape[0]+=1
    new_shape=zeros(new_shape,'b')
    new_shape[1:,:,:,:]=shape
    new_shape[0,:,:,:]=shape[0,:,:,:]
    v=outfile.createVariable('SHAPE','b',('TSTEP','LAY','ROW','COL'))
    v.assignValue(new_shape)
    v.units='ON/OFF'
    v.long_name=v.var_desc='SHAPE'.ljust(16)

def LoadPAQAFromYAML(yamlfile): 
    # Step 0: load YAML file as job, which is an attribute dictionary
    job=AttrDict(yaml.load(file(yamlfile)))
    
    AddDefaults(job)

    pr_rr=eval(job.metawrapper)(job.files)
    spc_iter=job.species
    prc_iter=job.processes
    from pyPA.pappt.pa_qa import PAdC
    PAdC(pr_rr,spc_list=job.species,prc_list=job.processes,model=job.model,mechanism=job.mechanism,verbose=True,init=job.init_conc,final=job.final_conc)

def LoadPyPAFromYAML(yamlfile):
    # Step 0: load YAML file as job, which is an attribute dictionary
    job=AttrDict(yaml.load(file(yamlfile)))
    
    AddDefaults(job)

    pr_rr=eval(job.metawrapper)(job.files)
    spc_iter=job.species
    prc_iter=job.processes
    rxn_iter=job.reactions
    shape=array(pr_rr.variables[job.shape])

    # May be scalar
    if type(job.irr_unitconversion)==str:
        irr_ucnv=pr_rr.variables[job.irr_unitconversion]
    else:
        irr_ucnv=job.irr_unitconversion

    # May be scalar
    if type(job.ipr_unitconversion)==str:
        ipr_ucnv=pr_rr.variables[job.ipr_unitconversion]
    else:
        ipr_ucnv=job.ipr_unitconversion

    # May be scalar
    if type(job.irr_contribution)==str:
        irr_contribution=pr_rr.variables[job.irr_contribution]
    else:
        irr_contribution=job.irr_contribution

    if job.ipr_contribution==job.irr_contribution:
        ipr_contribution=irr_contribution
    else:
        if type(job.ipr_contribution)==str:
            ipr_contribution=pr_rr.variables[job.ipr_contribution]
        else:
            ipr_contribution=job.ipr_contribution

    if job.normalizer==job.irr_contribution:
        normalizer=irr_contribution
    else:
        if type(job.normalizer)==str:
            normalizer=pr_rr.variables[job.normalizer]
        else:
            normalizer=job.normalizer
        
    
    kaxis=job.kaxis
    # Step 2: run ext_mrg and received 5 volume output
    output_5volumes=ext_mrg(pr_rr,spc_iter,prc_iter,rxn_iter,shape,ipr_ucnv,irr_ucnv,ipr_contribution,irr_contribution,normalizer,kaxis)
    
    # Step 3: merge 5 volume output creating pseudo-processes en(de)train and dilution
    ipr_irr=output_5volumes.merge()
    

    # BEGIN STEP 4: add meta data for storage
    ipr_irr.createDimension('TSTEP_STAG',shape.shape[0])
    ipr_irr.createDimension('LAY',shape.shape[1])
    ipr_irr.createDimension('ROW',shape.shape[2])
    ipr_irr.createDimension('COL',shape.shape[3])
    ipr_irr.SDATE=pr_rr.SDATE
    ipr_irr.STIME=pr_rr.STIME
    ipr_irr.TSTEP=pr_rr.TSTEP
    ipr_irr.irrfile='\n'.join([p for p,r in job.files])
    ipr_irr.iprfile='\n'.join([p for p,r in job.files])
    try:
        SDATE=ipr_irr.SDATE[0]
        STIME=ipr_irr.STIME[0]
        TSTEP=ipr_irr.TSTEP[0]
    except:
        SDATE=ipr_irr.SDATE
        STIME=ipr_irr.STIME
        TSTEP=ipr_irr.TSTEP
    
    # BEGIN STEP 5: write merged output to disk and sync buffer to disk
    outf=MergedWriter(job.outfile,ipr_irr,shape,pr_rr.variables['TFLAG'][:,0,:])
    outf.sync()
    # END STEP 5: write merged output to disk and sync buffer to disk
    
    # BEGIN STEP 6: send old merge text file to standard out
    try:
        {'camx': LegacyMerged, 'cmaq': LegacyMergedCMAQ}[job.model.lower()](sys.stdout,outf)
    except:
        pass
    # END STEP 6: send old merge text file to standard out

def AddDefaults(job):
    try:
        defaults=AttrDict(yaml.load(file(os.path.join(os.path.dirname(__file__),'defaults','_'.join([job.model.lower(),job.mechanism.lower()])+'.yaml'),'r')))
    except IOError:
        defaults=AttrDict(yaml.load(file(os.path.join(os.path.dirname(__file__),'defaults',job.model.lower()+'.yaml'),'r')))
    except:
        defaults={}
    for k in defaults:
        if not job.has_key(k):
            job[k]=defaults[k]
            
def template(model,mechanism,outfile=sys.stdout):
    if type(outfile)==str:
        outfile=file(outfile,'w')
    try:
        print >> outfile, file(os.path.join(os.path.dirname(__file__),'defaults',model.lower()+'_'+mechanism.lower()+'.yaml'),'r').read()
    except:
        print >> outfile, file(os.path.join(os.path.dirname(__file__),'defaults','new_template.yaml'),'r').read()
            
if __name__=='__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.set_usage("Usage: %prog [-tq] [-m <model>] [-c <mechanism>]\n"+(" "*16)+" [-i <init name>] [-f <final name>] <yamlfile>")
    parser.add_option("-q", "--qa", dest="qa",action="store_true",default=False,
                        help="Check PA data", metavar="QAPA")
    
    parser.add_option("-t", "--template", dest="template",action="store_true",default=False,
                        help="Output template on standard out (configurable with -m and -c", metavar="Template")
    
    parser.add_option("-m", "--model", dest="model",default="new",
                        help="Model can either be camx, cmaq or wrfchem (for use with -t)", metavar="MODEL")
    
    parser.add_option("-c", "--mechanism", dest="mechanism", default="template",
                        help="Chemical mechanisms: cbivmech3 (for use with -t)", metavar="MECHANISM")
    
    parser.add_option("-i", "--initial", dest="initial", default="INIT")
    parser.add_option("-f", "--final", dest="final", default="FCONC")
    
    (options, args) = parser.parse_args()
    if options.template:
        template(options.model,options.mechanism)
        parser.exit()
    if len(args)<1:
        parser.error(msg="Requires a yaml file as an argument.  For a template use the -t option.  The template will be output to the stdout.")
    else:
        yamlpath=args[0]
    
    if len(args)==2:
        netbalance=args[1]
    
    if options.qa:
        LoadPAQAFromYAML(yamlpath)
    else:
        LoadPyPAFromYAML(yamlpath)
