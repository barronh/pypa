__all__ = ['time_avg_new_unit', 'wind_center_time_cell', 'pypass_cmaq_met_master', 'pypass_cmaq_emiss_master', 'MetaMetPlusAirMols', 'cmaq_pa_master']

import os
from datetime import datetime
from warnings import warn
from numpy import array,ones
from PseudoNetCDF.MetaNetCDF import add_derived, \
                                  file_master
from ..netcdf import NetCDFFile
from PseudoNetCDF.ArrayTransforms import CenterCMAQWind, \
                                       CenterTime
from PseudoNetCDF.units import F2K
from PseudoNetCDF.sci_var import PseudoNetCDFFile, \
                    PseudoNetCDFVariable, \
                    PseudoIOAPIVariable, \
                    PseudoNetCDFVariables, \
                    PseudoNetCDFVariableConvertUnit

#==================================================================
#                                                             time_avg_new_unit 
class time_avg_new_unit(PseudoNetCDFFile):
#    __reader__=None
    def __init__(self,rffile,outunit={},endhour=False):
        self.__file=NetCDFFile(rffile)
        self.dimensions={}
        self.createDimension('TSTEP',self.__file.variables[self.__file.variables.keys()[0]].shape[0]-1)
        self.createDimension('LAY',self.__file.dimensions['LAY'])
        self.createDimension('ROW',self.__file.dimensions['ROW'])
        self.createDimension('COL',self.__file.dimensions['COL'])
        self.createDimension('VAR',self.__file.dimensions['VAR'])
        self.createDimension('DATE-TIME',self.__file.dimensions.get('DATE-TIME',2))
        self.__outunit=outunit
        self.variables=PseudoNetCDFVariables(self.__variables,self.__file.variables.keys())
        self.__timeslice={True:slice(1,None),False:slice(None,-1)}[endhour]
        v=self.createVariable('TFLAG','i',('TSTEP','VAR','DATE-TIME'),keep=True)
        v[:] = self.__file.variables['TFLAG'][self.__timeslice]
        v.long_name='Time flag'
        v.units='DATE-TIME'

    def __variables(self,k):
        outunit=self.__outunit.get(k,None)
        var=self.__file.variables[k]
        if outunit==None:
            outunit=var.units
        return PseudoNetCDFVariableConvertUnit(self.__decorator(var,PseudoNetCDFVariable(self,k,var.typecode(),var.dimensions,values=CenterTime(var))),outunit)
    
    def __decorator(self,ovar,nvar):
        for a,v in ovar.__dict__.iteritems():
            setattr(nvar,a,v)
        return nvar

#==================================================================
#                                                             NetCDFFile_center_time 
# class NetCDFFile_center_time(time_avg_new_unit):
#     __reader__=NetCDFFile

#==================================================================
#                                                             NetCDFFile_center_time 
class wind_center_time_cell(PseudoNetCDFFile):
    """
    CMAQ Files
    """
    def __init__(self,rffile,outunit='m/s',endhour=False):
        self.__windfile=NetCDFFile(rffile)
        self.dimensions={}
        self.createDimension('TSTEP',self.__windfile.variables['UWIND'].shape[0]-1)
        self.createDimension('LAY',self.__windfile.dimensions['LAY'])
        self.createDimension('ROW',self.__windfile.dimensions['ROW']-1)
        self.createDimension('COL',self.__windfile.dimensions['COL']-1)
        self.createDimension('VAR',self.__windfile.dimensions['VAR'])
        self.createDimension('DATE-TIME',self.__windfile.dimensions.get('DATE-TIME',2))
        self.__outunit=outunit
        self.variables=PseudoNetCDFVariables(self.__variables,['VWIND','UWIND','TFLAG'])
        self.__timeslice={True:slice(1,None),False:slice(None,-1)}[endhour]
    def __variables(self,k):
        self.__add_variables()
        return self.variables[k]
        
    def __add_variables(self):
        v=self.createVariable('TFLAG','i',('TSTEP','VAR','DATE-TIME'),keep=True)
        v[:] = self.__windfile.variables['TFLAG'][self.__timeslice]
        v.long_name='Time flag'
        v.units='DATE-TIME'
        for k in ['UWIND','VWIND']:
            preproc=CenterCMAQWind   
            var=self.__windfile.variables[k]
            v=PseudoNetCDFVariable(self,k,'f',('TSTEP','LAY','ROW','COL'),values=preproc(var))
            v.units=var.units
            v.long_name=k.ljust(16)
            v.var_desc=(k+' at center').ljust(16)
            self.variables[k]=PseudoNetCDFVariableConvertUnit(v,self.__outunit)
    
    

#==================================================================

def pypass_cmaq_met_master(metcro2d_path,metcro3d_path,metdot3d_path,rows,cols,endhour=True):
     windf=wind_center_time_cell(metdot3d_path,outunit='km/h',endhour=endhour)
     dim2f=time_avg_new_unit(metcro2d_path,outunit={'PBL':'km', 'TEMPG':'deg_F'},endhour=endhour)
     dim3f=time_avg_new_unit(metcro3d_path,outunit={'TA':'deg_F','ZF':'km'},endhour=endhour)
     return file_master([windf,dim2f,dim3f])
     
def pypass_cmaq_emiss_master(emiss_path,rows,cols,endhour=True):
     emiss = time_avg_new_unit(emiss_path,outunit={'CO':'moles/h', 'NO':'moles/h', 'NO2':'moles/h', 'ALD2':'moles/h', 'ETH':'moles/h','ETOH':'moles/h', 'MEOH':'moles/h','FORM':'moles/h', 'ISOP':'moles/h', 'OLE':'moles/h', 'PAR':'moles/h', 'TOL':'moles/h', 'XYL':'moles/h'},endhour=endhour)
     return file_master([emiss])

class MetaMetPlusAirMols(add_derived):
    __childclass__=lambda *args: NetCDFFile(*args[1:])
    __addvars__=['AIRMOLS','AREA']
    def __AREA__(self):
        #self.XSIZE*self.YSIZE
        return 4000**2
        
    def __AIRMOLS__(self):
        V=array(self.variables['DEPTH'])*array(self.variables['AREA']) #m**3
        T=F2K(array(self.variables['AIRTEMP'])) #K
        P=array(self.variables['PRES'])*100. #PA
        
        R=8.314472 #m**3 x Pa x K**-1 x mol**-1
        airmols=PseudoNetCDFVariable(self,'AIRMOLS','f',('TSTEP','LAY','ROW','COL'),values=P*V/R/T)
        airmols.units='moles'
        airmols.long_name='AIRMOLS'.ljust(16)
        airmols.var_desc='AIRMOLS'.ljust(16)
        return airmols

def cmaq_pa_master(paths_and_readers,tslice=None,kslice=None,jslice=None,islice=None):
    """
    CMAQ PA Master presents a single interface for CMAQ PA, IRR, and 
    Instantaneous concentration files.
    
    paths_and_readers - iterable of iterables (n x 2) where each element of the
                        primary iterable is an iterable containing a file path 
                        and a reader for that path.  The reader is expected to
                        present the Scientific.IO.NetCDF.NetCDFFile file inter-
                        face.
    
    optional:
       tslice - slice object used to window the time period from the 
                instantaneous concentration file to match the PA domain
       kslice - same as tslice, but for layers
       jslice - same as tslice, but for rows (defaults to an offset based on PA and conc difference)
       islice -  - same as tslice, but for columns (defaults to an offset based on PA and conc difference)
              
    """
    from ..pappt.kvextract import tops2shape,pblhghts2tops
    files=[]
    iprf = None
    concf = None
    for p,r in paths_and_readers:
        if not os.path.exists(p):
            raise ValueError, "File at %s does not exist" % p
        
        try:
            thisf = eval(r)(p)
        except:
            warn("Could not read %s with %s" % (p, r))
            raise
        files.append(thisf)
        if hasattr(thisf, 'FILEDESC'):
            if iprf is None and ("Integrated Process Rates Output File" in thisf.FILEDESC or "Integrated Reaction Rates Output File" in thisf.FILEDESC):
                iprf = thisf
            if "Concentration file output" in thisf.FILEDESC:
                concf = thisf

    master_file=file_master(files)
    
    if (iprf and concf) is not None:
        ipr_date = datetime.strptime('%dT%06d' % (iprf.SDATE, iprf.STIME), '%Y%jT%H%M%S')
        conc_date = datetime.strptime('%dT%06d' % (concf.SDATE, concf.STIME), '%Y%jT%H%M%S')
        step_sec = (iprf.TSTEP / 10000) * 3600 + iprf.TSTEP % 10000 / 100 * 60 + iprf.TSTEP % 100
        toffset = (ipr_date - conc_date).seconds / step_sec -1
        tend = toffset + len(iprf.dimensions['TSTEP'])
        tslice = tslice = slice(toffset, tend)
        tmslice = slice(tslice.start, tslice.stop+1)
        for vi, sigma in enumerate(concf.VGLVLS):
            if sigma == iprf.VGLVLS[0]:
                kslice = slice(vi, vi + iprf.NLAYS)
                break
        else:
            warn("""Could not align concentration and ipr vertical levels:
            conc: %s
            ipr: %s
            """ % (str(iprf.VGLVLS), str(concf.VGLVLS)))
            kslice = kslice or slice(None)
        
        
        xoffset = abs(concf.XORIG -  iprf.XORIG) / iprf.XCELL
        assert(float(xoffset) == float(int(xoffset)))
        xoffset = int(xoffset)
        xend = iprf.NCOLS + xoffset
        islice = islice or slice(xoffset, xend)
        
        yoffset = abs(concf.YORIG -  iprf.YORIG) / iprf.YCELL
        assert(float(yoffset) == float(int(yoffset)))
        yoffset = int(yoffset)
        yend = iprf.NROWS + yoffset
        jslice = jslice or slice(yoffset, yend)
        
        warn("""Automatically offseting based on subdomain.
\ttstart: %d
\ttend: %d
\tzstart: %d
\tzend: %d
\txstart: %d
\txend: %d
\tystart: %d
\tyend: %d""" % (tslice.start+1, tslice.stop, kslice.start+1, kslice.stop, islice.start+1, islice.stop, jslice.start+1, jslice.stop))
    else:
        warn("""Unable to identify CONC outputs and IPR outputs for automated slicing.
                Looking at FILEDESC properties for "Integrated Process Rates Output File" or
                "Integrated Reaction Rates Output File" and "Concentration file output"
                """)

    def InitLambda(x,tslice,kslice,jslice,islice):
        return lambda self: PseudoIOAPIVariable(self,x,'f',('TSTEP','LAY','ROW','COL'),values=self.variables[x][:-1,:,:,:][tslice,kslice,jslice,islice],units=self.variables[x].units)
    def FinalLambda(x,tslice,kslice,jslice,islice):
        return lambda self: PseudoIOAPIVariable(self,x,'f',('TSTEP','LAY','ROW','COL'),values=self.variables[x][1:,:,:,:][tslice,kslice,jslice,islice],units=self.variables[x].units)
    def MetLambda(x,tslice,kslice,jslice,islice):
        return lambda self: PseudoIOAPIVariable(self,x,'f',('TSTEP','LAY','ROW','COL'),values=CenterTime(self.variables[x])[tslice,kslice,jslice,islice],units=self.variables[x].units)
    for k in master_file.variables.keys():
        if '_' not in k and k!='TFLAG':
            master_file.addMetaVariable('INIT_'+k,InitLambda(k,tslice,kslice,jslice,islice))
            master_file.addMetaVariable('FCONC_'+k,FinalLambda(k,tslice,kslice,jslice,islice))
            master_file.addMetaVariable('INITIAL_'+k,InitLambda(k,tslice,kslice,jslice,islice))
            master_file.addMetaVariable('FINAL_'+k,FinalLambda(k,tslice,kslice,jslice,islice))
    master_file.addMetaVariable('VOL',lambda self: PseudoIOAPIVariable(self,'AIRMOLS','f',('TSTEP','LAY','ROW','COL'),values=self.XCELL*self.YCELL*2*CenterTime(self.variables['ZF'][:]-self.variables['ZH'][:])[tslice,kslice,jslice,islice],units='m**3'))
    master_file.addMetaVariable('CONC_AIRMOLS',lambda self: PseudoIOAPIVariable(self,'CONC_AIRMOLS','f',('TSTEP','LAY','ROW','COL'),values=CenterTime(self.variables['PRES'][:]/8.314472/self.variables['TA'][:])[tslice,kslice,jslice,islice],units='moles/m**3'))
    master_file.addMetaVariable('AIRMOLS',lambda self: PseudoIOAPIVariable(self,'AIRMOLS','f',('TSTEP','LAY','ROW','COL'),values=self.variables['CONC_AIRMOLS'][:]*self.XCELL*self.YCELL*2*CenterTime(self.variables['ZF'][:]-self.variables['ZH'][:])[tslice,kslice,jslice,islice],units='moles'))
    master_file.addMetaVariable('INVAIRMOLS',lambda self: PseudoIOAPIVariable(self,'INVAIRMOLS','f',('TSTEP','LAY','ROW','COL'),values=1/self.variables['AIRMOLS'][:],units='moles'))
    master_file.addMetaVariable('DEFAULT_SHAPE',lambda self: PseudoIOAPIVariable(self,'DEFAULT_SHAPE','f',('TSTEP','LAY','ROW','COL'),values=tops2shape(pblhghts2tops(self.variables['PBL'][:],self.variables['ZF'][:]),self.variables['ZF'][:].shape)[tmslice,kslice,jslice,islice],units='on/off'))
    
    return master_file
