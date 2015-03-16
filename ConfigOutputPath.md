# Configuration Output Path #

This configuration key should be set to a valid path on the users system for the output file to  be created.  The output file is a NetCDF file and is described below using the NetCDF Common Data Language (CDL).  The output file has 4 variables and 11 dimensions.  The IPR variable contains all integrated process rates for each chemical species defined in the configuration file.  The IRR variable contains all integrated reaction rates for each reaction defined in the configuration file.  The shape variable is an array of 0s and 1s that define the Process Analysis Volume [wiki:ConfigShape (for more information, see the shape section)].  The TFLAG variable contains julian date and time (HHMMSS) data for each Process Analysis output time step.


```
netcdf outfile.mrg.nc {
dimensions:
        DATE-TIME = 2 ;
        RXN = 156 ;
        PROCESS = 31 ;
        COL = 54 ;
        LAY = 20 ;
        VAR = 3 ;
        TSTEP_STAG = 19 ;
        ROW = 46 ;
        SPECIES = 27 ;
        TSTEP = UNLIMITED ; // (18 currently)
variables:
        float IPR(TSTEP, SPECIES, PROCESS) ;
                IPR:units = "ppmV/time" ;
        float IRR(TSTEP, RXN) ;
                IRR:units = "ppmV/time" ;
        int SHAPE(TSTEP_STAG, LAY, ROW, COL) ;
                SHAPE:units = "onoff" ;
                SHAPE:var_desc = "SHAPE           " ;
                SHAPE:long_name = "SHAPE           " ;
        int TFLAG(TSTEP, VAR, DATE-TIME) ;
                TFLAG:units = "YYYYJJJ,HHDDMM" ;
                TFLAG:var_desc = "TFLAG           " ;
                TFLAG:long_name = "TFLAG           " ;

// global attributes:
                :Reactions = "RXN_01          RXN_02          RXN_03          ...          RXN_N";
                :SDATE = 2006233 ;
                :Process = "INIT             EMIS            CHEM          ...            FCONC";
                :iprfile = "camx451_cb05.20060821.bc06aqs1.reg8si_pscfv2.2006ep1_eta_dbemis_fddats_uhsst_utcsrlulc.0821_take2.ipr.nc\n",
                        "/work/yosuke/lpa/output/0821_take2/camx451_cb05.20060821.bc06aqs1.reg8si_pscfv2.2006ep1_eta_dbemis_fddats_uhsst_utcsrlulc.0821_take2.irr\n",
                        "shape.nc" ;
                :irrfile = "camx451_cb05.20060821.bc06aqs1.reg8si_pscfv2.2006ep1_eta_dbemis_fddats_uhsst_utcsrlulc.0821_take2.ipr.nc\n",
                           "/work/yosuke/lpa/output/0821_take2/camx451_cb05.20060821.bc06aqs1.reg8si_pscfv2.2006ep1_eta_dbemis_fddats_uhsst_utcsrlulc.0821_take2.irr\n",
                           "shape.nc" ;
                :STIME = 7.987401e-42f ;
                :Species = "O3              NO2             NO              ...";
                :TSTEP = 100. ;
}
```