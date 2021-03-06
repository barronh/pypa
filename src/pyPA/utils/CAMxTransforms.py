__all__ = ['point_source_newvarnames', 'pypass_camx_met_master', 'camx_pa_master', 'hght2dz', 'hght2zh', 'pypass_camx_emiss_master']


from warnings import warn
from ..netcdf import NetCDFFile
ncf = NetCDFFile
from numpy import zeros, ones
from PseudoNetCDF.ArrayTransforms import CenterTime
from PseudoNetCDF.camxfiles.wind.Transforms import wind_center_time_cell
from PseudoNetCDF.camxfiles.height_pressure.Transforms import height_pressure_center_time_plus, \
                                            height_pressure_center_time, \
                                            height_pressure_plus
from PseudoNetCDF.camxfiles.temperature.Transforms import temperature_center_time
from PseudoNetCDF.camxfiles.vertical_diffusivity.Transforms import vertical_diffusivity_center_time
from PseudoNetCDF.camxfiles.humidity.Transforms import humidity_center_time
from PseudoNetCDF.camxfiles.cloud_rain.Transforms import cloud_rain_center_time_plus, \
                                       cloud_rain_center_time, \
                                       cloud_rain_plus
from CAMxFiles import *
from PseudoNetCDF.MetaNetCDF import *
from CAMxFiles import point_source as reg_point_source, \
                      gridded_emissions as reg_gridded_emissions
from PseudoNetCDF.sci_var import PseudoNetCDFFile, \
                    PseudoNetCDFVariables, \
                    PseudoNetCDFVariable, \
                    PseudoIOAPIVariable


#==================================================================

#                                                             point_source 

class point_source_newvarnames(PseudoNetCDFFile):
    __childclass__=reg_point_source
    __newvars__=['ptCO','ptNO','ptNO2','ptALD2','ptALDX','ptETH','ptETHA','ptFORM','ptISOP','ptNR','ptOLE','ptIOLE','ptPAR','ptTOL','ptXYL', 'ptNH3', 'ptSO2', 'ptSULF', 'ptPEC','ptPNO3','ptPOA','ptPSO4','ptETOH','ptMEOH']
    __oldvars__=['XSTK','YSTK','HSTK']  
    def __init__(self,rffile):
        self.__child=self.__childclass__(rffile)
        self.dimensions=self.__child.dimensions
        self.variables=PseudoNetCDFVariables(self.__variables__,self.__newvars__ + self.__oldvars__)
    def __variables__(self,key):
        if key in self.__newvars__:
           val = self.__child.variables[key[2:]]
           var=PseudoNetCDFVariable(self,key,'f',('TSTEP','STK'),values=val)
        elif key in self.__oldvars__:
           val = self.__child.variables[key]
           var=PseudoNetCDFVariable(self,key,'f',('TSTEP','STK'),values=val)
#          return getattr(self,var)
        return var
                
#==================================================================

def pypass_camx_met_master(wind_path,hp_path,temp_path,kv_path,hum_path,cr_path,rows,cols,endhour=True):
     windf=wind_center_time_cell(wind_path,rows,cols,outunit='km/h',endhour=endhour,forcestaggered=False)
     hpf=height_pressure_center_time_plus(hp_path,rows,cols,outunit={'HGHT':'km', 'PRES':'hPA'},endhour=endhour)
     tempf=temperature_center_time(temp_path,rows,cols,outunit={'AIRTEMP':'deg_F','SURFTEMP':'deg_F'},endhour=endhour)
     kvf=vertical_diffusivity_center_time(kv_path,rows,cols,outunit={'KV':'m**2/s'},endhour=endhour)
     humf=humidity_center_time(hum_path,rows,cols,outunit={'HUM':'ppm'},endhour=endhour)
     crf=cloud_rain_center_time_plus(cr_path,rows,cols,outunit={'CLOUD':'g/m**3','RAIN':'g/m**3','SNOW':'g/m**3','GRAUPEL':'g/m**3','PRECIP':'g/m**3','PRECIPRATE':'mm/h','COD':'None'},endhour=endhour)
     return file_master([windf,hpf,tempf,kvf,humf,crf])

def camx_pa_master(paths_and_readers,tslice=None,kslice=None,jslice=None,islice=None, padom = 1):
    """
    CAMx PA Master presents a single interface for CAMx IPR, IRR, and shape definitions.
    If the shape is not defined, CAMx PA Master can provide a default shape horiztonally
    bound by the domain and vertically bound by the planetary boundary layer.  The
    planetary boundary layer is diagnosed by from the vertical diffusivity and vertical
    layer structure.
    
    paths_and_readers - iterable of iterables (n x 2) where each element of the
                        primary iterable is an iterable containing a file path 
                        and a reader for that path.  The reader is expected to
                        present the Scientific.IO.NetCDF.NetCDFFile file inter-
                        face.
    
    optional:
       tslice - slice object used to window the time period from the 
                instantaneous meteorological files to match the PA domain
       kslice - same as tslice, but for layers
       jslice - same as tslice, but for rows
       islice -  - same as tslice, but for columns
       padom - pa domain 1 is first pa domain
       
       if kslice, jslice, or islice are not provided, they will be calculated
       from pafiles
              
    """
    def defaultshape(self, kvname = 'KV', hghtname = 'HGHT'):
        from ..pappt.kvextract import tops2shape, vertcamx
        
        old_shape=[i for i in self.variables['UCNV_O3'].shape]
        new_shape=[i for i in self.variables['UCNV_O3'].shape]
        new_shape[0]+=1
        new_shape=zeros(tuple(new_shape),dtype='bool')
        new_shape[1:,:,:,:]=tops2shape(vertcamx(CenterTime(self.variables[kvname]),CenterTime(self.variables[hghtname]))[tslice,jslice,islice],old_shape)
        new_shape[0,:,:,:]=new_shape[1,:,:,:]
        return PseudoIOAPIVariable(self,'DEFAULT_SHAPE','i',('TSTEP', 'LAY', 'ROW', 'COL'),values=new_shape,units='on/off')

    
    for i, (p, r) in enumerate(paths_and_readers):
        if r[:5] == 'from ' and 'import' in r:
            exec r in globals(), locals()
            paths_and_readers[i][1] = r.split(' ')[-1]

    # Find irr or ipr and open as pafile; pafiles have grid
    # and padomain information that can be used to shape and
    # subset input files
    pafiles = ([(p, r) for p,r in paths_and_readers if r == 'ipr'] + [(p, r) for p,r in paths_and_readers if r == 'irr'])

    if pafiles != []:
        p, r = pafiles[0]
        pafile = eval(r)(p)
            
        padomain = pafile.padomains[0]
        grid = pafile.grids[padomain['grid']-1]
        
        tslice = tslice or slice(None)
        kslice = kslice or slice(padomain['blay']-1, padomain['tlay'])
        islice = islice or slice(padomain['istart']-1, padomain['iend'])
        jslice = jslice or slice(padomain['jstart']-1, padomain['jend'])
        nrows = grid['nrow']
        ncols = grid['ncol']
    else:
        tslice = tslice or slice(None)
        kslice = kslice or slice(None)
        islice = islice or slice(None)
        jslice = jslice or slice(None)
        for p, r in paths_and_readers:
            tempfile = eval(r)(p)
            if 'COL' in tempfile.dimensions and 'ROW' in tempfile.dimensions:
                nrows = tempfile.dimensions['ROW']
                ncols = tempfile.dimensions['COL']
                if not isinstance(nrows, int):
                    nrows = len(nrows)
                    ncols = len(ncols)
                break
        else:
            raise KeyError, "No files had COL and ROW dimensions"
            
    paths_and_readers = [(p, {'height_pressure': 'lambda path: height_pressure(path, %d, %d)' % (nrows, ncols), 'vertical_diffusivity': 'lambda path: vertical_diffusivity(path, %d, %d)' % (nrows, ncols)}.get(r, r)) for p,r in paths_and_readers]

    # Create a list of opened files
    files=[eval(r)(p) for p,r in paths_and_readers]
    
    # Create a file master object from files
    master_file=file_master(files)
    
    ## Add derived variables
    
    # Rename AVOL_O3 to VOL -- simplicity
    #   The choice of O3 was arbitrary.  All AVOL 
    #   values (e.g. AVOL_O3, AVOL_NO2, etc) are equal.
    master_file.addMetaVariable('VOL',lambda self: PseudoIOAPIVariable(self,'VOL','f',('TSTEP','LAY','ROW','COL'),values=self.variables['AVOL_O3'],units='m**3'))
    
    # Calculate AIRMOLS from IPR outputs
    #   UCNV [=] m**3/mol_{air}
    #   AVOL [=] m**3
    #   AIRMOLS = AVOL/UCNV [=] mol_air
    master_file.addMetaVariable('AIRMOLS',lambda self: PseudoIOAPIVariable(self,'AIRMOLS','f',('TSTEP','LAY','ROW','COL'),values=self.variables['AVOL_O3'][:]/self.variables['UCNV_O3'][:],units='moles'))

    # Calculate INVAIRMOLS from IPR outputs
    #   UCNV [=] m**3/mol_{air}
    #   AVOL [=] m**3
    #   AIRMOLS = UCNV/AVOL [=] 1/mol_air
    master_file.addMetaVariable('INVAIRMOLS',lambda self: PseudoIOAPIVariable(self,'INVAIRMOLS','f',('TSTEP','LAY','ROW','COL'),values=self.variables['UCNV_O3']/self.variables['AVOL_O3'],units='moles**-1'))

    # Calculate the well mixed region of the troposphere based
    # on vertical diffusivity and layer structure.
    mhas_key = master_file.variables.has_key
    if not mhas_key('DEFAULT_SHAPE'):
        if mhas_key('KV') and mhas_key('HGHT'):
            master_file.addMetaVariable('DEFAULT_SHAPE',lambda self: defaultshape(self, kvname = 'KV', hghtname = 'HGHT'))
        elif mhas_key('KV_M2pS') and mhas_key('ZGRID_M'):
            master_file.addMetaVariable('DEFAULT_SHAPE',lambda self: defaultshape(self, kvname = 'KV_M2pS', hghtname = 'ZGRID_M'))
        else:
            warn('Cannot diagnose 3-D DEFAULT_SHAPE\ncamx_pa_master needs KV from vertical_diffusivity (kv) file and HGHT from height/pressure (zp) file to diagnose the planetary boundary layer and provide DEFAULT_SHAPE')
            warn('Providing all cells.')
            all_vals = ones(map(lambda x: len(master_file.dimensions[x]) + 1 if x == 'TSTEP' else len(master_file.dimensions[x]), ('TSTEP', 'LAY', 'ROW', 'COL')), dtype = 'i')
            master_file.addMetaVariable('DEFAULT_SHAPE', lambda self: PseudoIOAPIVariable(self,'DEFAULT_SHAPE','i',('TSTEP', 'LAY', 'ROW', 'COL'),values=all_vals,units='on/off'))
    return master_file


def hght2dz(layer_tops):
    dz=layer_tops.copy()
    dz[:,1:,:,:]-=layer_tops[:,:-1,:,:]
    return dz

def hght2zh(layer_tops):
    return layer_tops-.5*dz

def pypass_camx_emiss_master(pointemiss_path,gridemiss_path,rows,cols,endhour=True):
    pointemiss=point_source_newvarnames(pointemiss_path)
    griddedemiss=reg_gridded_emissions(gridemiss_path)
    return file_master([pointemiss,griddedemiss])
    