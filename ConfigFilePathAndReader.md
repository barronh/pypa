# Configuration Paths and Readers #

pyPA reads Process Analysis and meteorology files in their model specific formats, but needs to be told the file paths and formats.  The paths can be any valid path on your system to a file.  The formats are generally model specific.  For some models (e.g. CMAQ and WRF-Chem), the inputs and outputs all have a single format (NetCDF).  For other models (e.g. CAMx), each file has its own unique format (arg...).  While pyPA attempts to provide file interfaces for all models, we cannot presume to know all models or future file formats.  If your model or data format is not explicitly handled, there is a way to build your own data interpreter.

CAMx:
> ipr: Integrated process rate file reader<br />
> irr: Integrated reaction rate file reader<br />
> vertical\_diffusivity: vertical diffusivity (aka kv) file reader<br />
> height\_pressure: height/pressure (aka zp) file reader<br />
> `*` vertical diffusivity and height pressure are only required when using DEFAULT\_SHAPE

CMAQ:
> NetCDFFile reads all required files: IRR, PA, METCRO2D and METCRO3D<br />
> `*` METCRO2D is only required when using DEFAULT\_SHAPE

WRF:
> NetCDFFile: reads all required files: IRR, PA, and met data

OTHER:
> Any 4-d (time/space) data can be used as long as you provide a reader
> function for the data.  The reader function must return a NetCDFFile
> object or a PseudoNetCDF object.  To specify your own reader put a
> "from package\_name import function\_name" statement in place of a reader.