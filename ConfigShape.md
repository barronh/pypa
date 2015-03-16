# Configuration Shape #

The shape variable is used to spatially/temporally subset model Process Analysis outputs before aggregating rate data to create the Process Analysis Volume.  This is useful when a Process Analysis domain cannot be defined given the limitations of the host model PA domain definitions or when a dynamically changing shape is required.  Host models typically allow for rectangular prizms, but pyPA allows for arbitrarily complex shapes in all 4 dimensions (space and time).  This allows for simple spatially static cubes, [planetary boundary layer conforming](stationary_pav.md), or complex [plume tracking analyses](pseudo_lagrangian_pav.md).  Any shape you can dream up, can be provided as a 4 dimensional variable of 1s and 0s.  When a variable element is 1, that space/time intersection is part of the analysis.  When a variable element is 0, that space/time intersection is excluded from the analysis.

pyPA provides several default shape variables or you can specify your own.  The DEFAULT\_SHAPE variable is the most common Process Analysis Volume (PAV).  The DEFAULT\_SHAPE variables uses the whole spatial domain of the Process Analysis output, but vertically subsets the data include only the well-mixed region of the troposphere within the planetary boundary layer (PBL).  The default shape can be horizontally subset using the [ascii\_mask](ConfigAsciiMask.md).  There are many alternative analyses that we cannot even begin to enumerate.  If your analysis requires a more complex shape, you can provide it by creating the variable and saving it as a file.  Include the path and reader in [files](ConfigFilePathAndReader.md) and set the shape configuration key to your new variable (as illustrated below).

```
files:
 - [/path/to/my/custom/shape_file, NetCDFFile]

shape: MY_SHAPE_NAME
```

The shape variable is always output in the pyPA output file.  This is a convenient starting point for defining your own custom shape following the steps below.  You use any NetCDF compatible program to edit the shape variable.

  1. Run pyPA with the shape: DEFAULT\_SHAPE
  1. Copy the output to a new file name
  1. Open the output file and edit the shape variable
  1. add the output file to the files paths and readers section
  1. set shape: shape
  1. Re-run pyPA to use your edited shape
