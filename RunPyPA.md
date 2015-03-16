# Post Processing PA Outputs #

Our group developed a PA post-processor, named pyPA, which uses Python, NumPy, and our own suite of file interface libraries to achieve the spatial ?exibility necessary to analyze a PA volume. pyPA uses a 4-D matrix to define the PA volume's 3-D spatial location at the beginning
and ending of each hour. This method allows us to quantify the changes in concentration due to PA volume motion and/or expansion. In a nutshell, pyPA calculates all the process in an individual PA cell and aggregates them for the chosen PA domain.


## pyPA Post-Processing Instructions ##
### Install pyPA ###

[Installation instructions](InstallInstructions.md)


### Create pyPA Configuration File ###

The configuration file allows pyPA correctly process your model PA data. By running the command below, a template configuration file is created for you to edit. If optional arguments are specified, pyPA can fill in almost all of the fields needed by pyPA.  The first optional argument is [model](model.md) (i.e. camx, cmaq, or wrfchem). The second argument is [Mechanism](Mechanism.md) (i.e. cbivmech3, cb05\_ae4\_aq, cb05cl\_ae5\_aq, cbmz). Note that the mechanism names follow the conventions of the host model.  A full list of models and mechanisms run pyPA without any options (i.e. `python -m pyPA`) <br />
```
$ python -m pyPA -t -m cmaq -c cb0cl_ae5_aq > configuration.yaml
```


### Configuration File ###

The configuration file has 5 sections.  Output path, model/mechanism identification, input file and readers, PA meta variables, PA IRR/IPR variables.  If you provided a model and mechanism parameter to generate your configuration file, you will only have to edit the output path and  input files options.  Other options can be edited to configure a more specific analysis.

#### OUTPUT PATH ####

> outfile::
> > pyPA output file path (directory and filename) <br />

#### MODEL/MECHANISM IDENTIFICATION ####


> model::
> > The model used to create PA output. This is already set when you make the default template. Model options are cmaq, camx, and wrfchem <br />


> mechanism::
> > The chemical mechanism used by the model run that generated the PA outputs. <br />


> kaxis::
> > Model outputs are arrays with dimensions.  The kaxis is the position of the vertical dimension of the array.  Most model outputs have dimensions (time, vertical, north-sout, east-west).  The position is a 0-based number, so the vertical axis is typically 1 and is pre-filled in the template. <br />


#### INPUT FILES AND READERS ####


> metawrapper::
> > The metawrapper's primary function is to take all inputs files and make them act as one. If any variables that the user specified are not in the raw data, a model specific wrapper (i.e. camx\_pa\_master or cmaq\_pa\_master) will attempt to use the available variables to calculate them. Parameters for this input are cmaq\_pa\_master, camx\_pa\_master, and file\_master (described below). <br />

  * camx\_pa\_master - Metawrapper that has special knowledge of CAMx. It knows how to process the meteorology files to create necessary IRR and IPR weighting and normalizing variables.  It can also process CAMx meteorology files to provide a default shape that is bounded vertically by the ground and the planetary boundary layer.  The planetary boundary layer is diagnosed by from the vertical diffusivity and vertical layer structure. The horizontal bounds are either the PA domain or, if provided, a text file mask (see shape)<br />
  * cmaq\_pa\_master - Metawrapper that has special knowledge of CAMx. It knows how to process the meteorology files to create necessary IRR and IPR weighting and normalizing variables.  It can also process CAMx meteorology files to provide a default shape that is bounded vertically by the ground and the planetary boundary layer.  The planetary boundary is used directly from the METCRO2D file and mapped to the model's vertical layer structure. The horizontal bounds are either the PA domain or, if provided, a text file mask (see shape)<br />
  * file\_master - A generic wrapper that has no knowledge of any model.  This wrapper should be used when the data contains all necessary variables and requires no in-line processing.
<br />


> files::
> > All necessary inputs files paths are included in this section. Each file should be on a line by itself following the pattern below.
```
 - [filepath, reader]
```
> > Required files for pyPA are the 'IPR' (Integrated Process Rate) and 'IPR' (Integrated Reaction Rate), which are outputs from the model. In addition to the IRR and IPR files, CMAQ needs instantaneous model concentration output files. The filepath is the path to the file (path syntax is machine dependent).  The reader to use depends on the model.  For CMAQ, the reader will always be NetCDFFile.  For CAMx, the reader will depend on the type of file. CAMx IRR files use irr as the reader; IPR files use ipr as the reader.


> In addition the required files, there are 2 optional file sets.  If the user wants to provide their own 'shape' variable, it should be provided as a NetCDFFile.  If the user wants pyPA to dynamically provide a shape that is vertically bounded by the planetary boundary layer, the user must provide meteorological inputs.  For CMAQ, the meteorology files required are METCRO2D and METCRO3D.  For CAMx, the meteorology files include the height pressure file (sometimes named camx\_zp...) and the vertical diffusivity file (sometimes named camx\_kv...).  For more information on the model mixing volume that is automatically generated see 'shape' in the PA META VARIABLES section.


#### PA META VARIABLES ####
PA needs the name of variables to be used for the following purposes:

> shape::
> > The shape values is used to generate a mask of 1's (on) and 0's (off) that tells pyPA which cells are part of the analysis volume.  If the shape value is a path to a file, the file should be an ascii file with 1's and 0's separated by a space.  Each line of the ascii file is a row in the model (top-to-bottom = north-to-south).  Each column (separated by spaces) corresponds to a column in the model (left-to-right = west-to-east).  The text below would be an ascii mask file for a 4 row and 5 column domain.

```
0 1 0 0 0
0 1 1 0 0
0 1 1 0 0
1 1 1 1 0
```


> If the shape is not the path to a file, it must be a variable in one of the specified files.  The variable must be a 4-D array of integers (1/0) or booleans (True/False) that includes or excludes cells in IRR/IPR.  The variable should have 1 more time than the data files.  The first location is the location before the analysis starts.

> init\_conc::
> > Name of the initial concentration variable in a specified file.


> final\_conc::
> > This include the variables name for the final concentrations in file.


> irr/ipr\_contribution::
> > Individual cell contribution to overall IRR/IPR weighted averages (RXN contribution and ipr contribution `[`defaults to 1`]`)


> normalizer::
> > Cell values to sum for normalization of IRR/IPR weighted averages (RXN normalizer and ipr normalizer `[`defaults to irr/ipr contribution`]`)


> irr/ipr\_unitconversion::
> > Conversion of units (RXN unitconversion and ipr unitconversion [to 1](defaults.md))



#### PA IRR/IPR VARIABLES ####
In these parameters you can specify which reactions, process and species are to be outputted.
The default setting are set to all species, process, and reactions. Keep in mind that these are
dependent of the atmospheric model used.

  * species
  * processes
  * reactions


### Edit pyPA template ###

The following are notes that should be taken into account when editing the pyPA template.


  * Edit paths for ipr, irr, vertical diffusivity and height pressure files. <br />
> > This workflow assumes absolute paths, nevertheless relative paths may also be applied.
  * Add readers for each of the files.
    * Files that use the IPR reader and are larger than 2G should use the mem reader. (i.e. instead of ipr reader use ipr mem)
    * Review what arguments are needed for each reader. (i.e vertical diffusivity takes 3 arguments [path, row, column])
  * Edit PA META VARIABLES
    * irr/ipr contribution : In the model, the PBL changes for every cell, and concentrations from the models are outputted in ppm for each cell. One can explain ppm as an extensive quantity (that is a physical quantity, whose value is proportional to the size of the system it describes ). The contribution variables change the units of these results to intensive quantities, which in this case are moles, so they can be accurately aggregated. CAMx outputs grid cell volumes for ipr cells so the ipr contributor is usually VOL. On the other hand, CMAQ doesn't have such variable, so a number of variables are used to calculate the ipr contributor (AIRMOLS). Generally, ipr contribution for both CAMx and CMAQ is AIRMOLS.
    * normalizer : These parameters convert the intensive quantities back to ppm and is usually AIRMOLS.


### Run pyPA ###
```
$ python -m pyPA configuration.yaml > pypa.log
```