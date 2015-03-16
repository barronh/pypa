[[PageOutline(2)]]
# Creating Model PA Output #
The methods and options for generating Process Analysis output are model specific.  We do our best to enable pyPA users to run PA on their host model, but the responsibility for documenting the process lies with the model user community.  The instructions below are provided as a ''best effort.'' Any discrepancies with the model's documentation should be assumed incorrect.  Occasionally, we find errors in model documentation and report it to the model community.  Anywhere our documentation should supersede the model documentation will be explicitly commented.

## CAMx ##
  1. Create a CAMx run and execute it without PA (see Air Quality Models),
  1. modify the namelist to enable PA according to the CAMx manual (see inset below), and
  1. rerun the CAMx simulation.

```
Probing_Tool = ’PA’, ! (None,OSAT,GOAT,APCA,DDM,PA,RTRAC) 
Number_of_Output_Species = 23, 
Output_Species_Names(1) = ’NO’, 
Output_Species_Names(2) = ’NO2’, 
Output_Species_Names(3) = ’O3’, 
Output_Species_Names(4) = ’OLE’, 
Output_Species_Names(5) = ’PAN’, 
Output_Species_Names(6) = ’NXOY’, 
Output_Species_Names(7) = ’PAR’, 
Output_Species_Names(8) = ’TOL’, 
Output_Species_Names(9) = ’XYL’, 
Output_Species_Names(10) = ’FORM’, 
Output_Species_Names(11) = ’ALD2’, 
Output_Species_Names(12) = ’ETH’, 
Output_Species_Names(13) = ’CRES’, 
Output_Species_Names(14) = ’MGLY’, 
Output_Species_Names(15) = ’OPEN’, 
Output_Species_Names(16) = ’PNA’, 
Output_Species_Names(17) = ’CO’, 
Output_Species_Names(18) = ’HONO’, 
Output_Species_Names(19) = ’H2O2’, 
Output_Species_Names(20) = ’HNO3’, 
Output_Species_Names(21) = ’ISOP’, 
Output_Species_Names(22) = ’MEOH’, 
Output_Species_Names(23) = ’ETOH’, 

Chemistry_Parameters = ’$INPUT/CAMx4.5.chemparam.3, 

&PA_Control 
PA_File_Root = ’CAMx.${RUN}.200206${CAL}.PA’, 

Number_of_PA_Domains = 1, 
Within_CAMx_Grid(1) = 1, 
PA_Beg_I_Index(1) = 41, 
PA_End_I_Index(1) = 50, 
PA_Beg_J_Index(1) = 51, 
PA_End_J_Index(1) = 60, 
PA_Beg_K_Index(1) = 1, 
PA_End_K_Index(1) = 16, 
&
```


## CMAQ ##
  1. Create a CMAQ run and execute it without PA,
  1. create a PACL input ﬁle,
  1. compile PROCAN with new PACL ﬁle,
  1. recompile CMAQ executable,
  1. modify the CMAQ execution script to enable PA,
  1. rerun the CMAQ simulation.
  * Details available from [CMAQ Operational Guidance](http://www.cmaq-model.org/op_guidance_4.6/html)

## WRF-Chem ##
  1. Create a WRF-Chem run and execute it without PA,
  1. modify the namelist to enable PA for the domain of interest, and
  1. rerun the WRF-chem simulation.
