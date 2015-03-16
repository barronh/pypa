# Configuration ASCII Mask #

Process Analysis model outputs are often rectangular prizms, but pyPA allows for horizontal subsetting to any resolvable shape.  The ascii\_mask configuration key provides a path to a text file with 1s and 0s.  A 1 tells pyPA to include a cell and a 0 tells pyPA to exclude that cell.  If the ascii\_mask path does not exist, pyPA will create a template for you to edit.
