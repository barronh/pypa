# Configuration Yaml Description #
## Overview ##
The configuration template file below has embedded links with more information on each component of the configuration file.

## Default Configuration Template ##
###### PA INPUT<br />
#### Process Analysis uses model generated Integrated Process and Reaction <br />
#### Rates (IPR/IRR) to create a conceputually well mixed volume that can<br />
#### be analyzed to conceptually represent physical and chemical evolution.<br />
####<br />
#### There are 5 sections.  Commented options are optional.<br />
<br />
###############################################<br />
###### OUTPUT PATH                       ######<br />
###############################################<br />
#### outfile is the path for results<br />
[outfile: <path for output file>](ConfigOutputPath.md)<br />
<br />
###############################################<br />
###### MODEL/MECHANISM IDENTIFICATION    ######<br />
###############################################<br />
#### models supported are cmaq, camx, wrfchem<br />
#### mechanisms supported cbivmech3<br />
#[model: cmaq|camx|wrfchem](ConfigModel.md)<br />
#[mechanism: cbivmech3|more-to-come](ConfigMechanism.md)<br />
<br />
###############################################<br />
###### PA META WRAPPERS                  ######<br />
###############################################<br />
#### Data that is needed by pyPA often has different<br />
#### names, forms, and or meta data conventions.  pyPA<br />
#### provides meta wrappers that will help you work with<br />
#### specific models.  Please choose the wrapper that goes<br />
#### with the model you are using.  If your model is not <br />
#### provided, try pafile\_master. pafile\_master is the most<br />
#### generic.<br />
[metawrapper: camx\_pa\_master|cmaq\_pa\_master|wrfchem\_pa\_master|pafile\_master](ConfigMetaWrapper.md)<br />
<br />
###############################################<br />
###### INPUT FILES AND READERS           ######<br />
###############################################<br />
#### pyPA uses integrated reaction rate, integrated<br />
#### process rate, instantaneous concentration data,<br />
#### and meteorological data from files.  The variables<br />
#### may be in any number of files with any number of<br />
#### formats that each have their own readers.  Provide<br />
#### the path for each file and the name of its reader<br />
#### function.<br />
####<br />
#### READER FUNCTIONS BY MODEL<br />
#### CAMx:<br />
####  ipr: Integrated process rate file reader<br />
####  irr: Integrated reaction rate file reader<br />
####  vertical\_diffusivity: vertical diffusivity (aka kv) file reader<br />
####  height\_pressure: height/pressure (aka zp) file reader<br />
#### CMAQ:<br />
####  NetCDFFile: reads all CMAQ IRR, PA, and met data<br />
#### WRF:<br />
####  NetCDFFile: reads all WRF IRR, PA, and met data<br />
#### OTHER:<br />
####  Any 4-d (time/space) data can be used as long as you provide a reader <br />
####  function for the data.  The reader function must return a NetCDFFile<br />
####  object or a PseudoNetCDF object.  To specify your own reader put a<br />
####  "from package\_name import function\_name" statement in place of a reader.<br />
[files:](ConfigFilePathAndReader.md) <br />
[- [<path to file>, &lt;reader&gt;](ConfigFilePathAndReader.md)]<br />
<br />
###############################################<br />
###### PA META VARIABLES                 ######<br />
###############################################<br />
#### PA needs the name of variables to be used<br />
#### for the following purposes:<br />
####  shape: name of PA volume shape variable<br />
####         that is a 4-D array of booleans that <br />
####         include or exclude cells in IRR/IPR<br />
####  init\_conc: Initial concentration<br />
####  final\_conc: Final concentration<br />
[shape: DEFAULT\_SHAPE](ConfigShape.md)<br />
[init\_conc: INIT](ConfigInitAndFinal.md)<br />
[final\_conc: FCONC](ConfigInitAndFinal.md)<br />
<br />
###############################################<br />
###### 2-D MASK FOR HORIZONTAL SUBSETTING######<br />
###############################################<br />
#### PA can subset your domain based on a text<br />
#### file with 1s (include) and 0s (excluded)<br />
#### <br />
#### ascii\_mask: path to text file.  If it does<br />
####             not exist, a template is created<br />
#[ascii\_mask: ascii\_mask.txt](ConfigAsciiMask.md)<br />
<br />
###############################################<br />
###### UNIT CONVERSION                   ######<br />
###############################################<br />
#### unitconversion: dictionary of units to <br />
####      convert.  Each unit to convert has<br />
####      an expression that is evaluated to<br />
####      convert the value to a new unit.<br />
unitconversions:<br />
> ppmV:<br />
> > expression: 

&lt;value&gt;

**1000.**<br />
> > new\_unit: ppb<br />

> ppm:<br />
> > expression: 

&lt;value&gt;

**1000.**<br />
> > new\_unit: ppb<br />
<br />
###############################################<br />
###### PA IRR/IPR VARIABES               ######<br />
###############################################<br />
#### species: list of species names<br />
#### processes: list of process names<br />
#### reactions: list of reaction variables or a<br />
####     regular expression pattern for reaction<br />
####     variable names<br />
[species: <list of species names>](ConfigSpecies.md)<br />
[processes: <list of process names>](ConfigProcesses.md)<br />
[reactions: '(IRR|RXN)\_.\*'](ConfigReactions.md)<br />