# Configuration metawrapper #

The metawrapper configuration key tells pyPA how to provide meta-data for your analysis.  You must select one of the options delimited by the pipe ("|").  pyPA has model specific metawrappers (e.g. camx\_pa\_master, cmaq\_pa\_master, wrfchem\_pa\_master) or a generic (pafile\_master).  The model specific metawrappers use meteorological files to provide special metadata to help with pyPA analyses. There are two places that pyPA needs metadata: 1) defining the Process Analysis Volume shape and 2) aggregating Process Analysis data.

  1. Process Analysis Volume (PAV) provides the spatial/temporal definition for the pyPA data extraction (for more information see PAV).  It is often advantageous to have the PAV conform to the planetary boundary layer to account for the well-mixed region as a whole.  The height of the PBL is stored differently by each model (CAMx, CMAQ, WRF-Chem, etc), so the metawrapper for each model knows how to calculate the shape of the PAV with meteorology data inputs for a specific model.

  1. pyPA takes individual grid cell Process Analysis data and combines it to make a Process Analysis Volume (PAV).  To combine the data, intrinsic variables are converted to extrinsic values, aggregated, and converted back to intrinsic values.  The conversion is unit specific (ppb, micrograms/m3, etc) and often requires metadata (such as volume or air density) for conversion.  The metawrapper knows how to calculate extrinsic variables appropriate for conversion of Process Analysis data.

cmaq\_pa\_master:
> The CMAQ PA master needs the METCRO3D file for temperature and pressure and the METCRO2D file for the planetary boundary layer.

camx\_pa\_master:
> The CAMx PA master needs the vertical diffusivity to diagnose the planetary boundary layer and the height/pressure file to calculate air density.