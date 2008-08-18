__all__=['wind_center_time_cell']
from numpy import array

from pyPA.utils.MetaNetCDF import add_derived, time_avg_new_unit
from pyPA.utils.sci_var import PseudoNetCDFFile, PseudoNetCDFVariables, PseudoNetCDFVariable, PseudoNetCDFVariableConvertUnit
from pyPA.utils.CAMxFiles import wind as reg_wind
from pyPA.utils.ArrayTransforms import CenterTime,CenterCAMxU,CenterCAMxV
class wind_center_time_cell(PseudoNetCDFFile):
    """
    CAMx Files
    """
    def __init__(self,rffile,rows,cols,outunit='m/s',endhour=True,forcestaggered=False):
        self.__windfile=reg_wind(rffile,rows,cols)
        self.dimensions={}
        self.createDimension('TSTEP',self.__windfile.dimensions['TSTEP']-1)
        self.createDimension('LAY',self.__windfile.dimensions['LAY'])
        self.createDimension('ROW',self.__windfile.dimensions['ROW'])
        self.createDimension('COL',self.__windfile.dimensions['COL'])
        self.createDimension('VAR',self.__windfile.dimensions['VAR'])
        self.createDimension('DATE-TIME',self.__windfile.dimensions.get('DATE-TIME',2))
        self.__outunit=outunit
        self.__force_stagger=forcestaggered
        self.variables=PseudoNetCDFVariables(self.__variables,['V','U','TFLAG'])
        self.__timeslice={True:slice(1,None),False:slice(None,-1)}[endhour]
    def __variables(self,k):
        self.__add_variables()
        return self.variables[k]
        
    def __add_variables(self):
        v=self.createVariable('TFLAG','i',('TSTEP','VAR','DATE-TIME'),keep=True)
        v.assignValue(self.__windfile.variables['TFLAG'][self.__timeslice])
        v.long_name='Time flag'
        v.units='DATE-TIME'
        if self.__force_stagger and self.__windfile.LSTAGGER==0:
            warn('Cell centered values are being averaged as though staggered'+ \
                 'Could just be pre v4.3 file that was actually staggered')

        for k in ['U','V']:
            if self.__force_stagger or self.__windfile.LSTAGGER!=0:
                if k=='U':
                    preproc=CenterCAMxU
                elif k=='V':
                    preproc=CenterCAMxV
            else:
                preproc=CenterTime
                
            var=self.__windfile.variables[k]
            v=PseudoNetCDFVariable(self,k,'f',('TSTEP','LAY','ROW','COL'),preproc(var))
            v.units=var.units
            v.long_name=k.ljust(16)
            v.var_desc=(k+' at center').ljust(16)
            self.variables[k]=PseudoNetCDFVariableConvertUnit(v,self.__outunit)

