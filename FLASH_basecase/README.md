Example case for FLASH with input density/tele/tion from txt files
=====================================================================================

## Setup

Base case tested on FLASH 4.8 (Documentation[https://flash.rochester.edu/site/flashcode/user_support/flash4_ug_4p8.pdf])

## Changes to Simulation_initBlock.F90

The main change to this case is the updating of Simulation_initBlock.F90: 

#### readTxtData module

The creation of a new module readTxtData was made, including the following sunroutine:
- read_txt_grid: open the given txt file. the file must be a one dimension array
- read_txt_mesh: open the given txt file. the file must be a two dimension array
- linear_interpolate_clamped: linear interpolation function. Clamped between the maximum of the file and a given minimal value
- cubic_interpolate_clamped: cubic interpolation function. Clamped between the maximum of the file and a given minimal value
- init_geom_from_txt: create the interpolated values for the block being initialised.

#### Simulation_initBlock subroutine  

The subroutine has been changed to take into account a 'fromtxt' initialisation mode calling init_geom_from_txt at the beginning of the block treatment and defining rho, tele, trad (defined as tele) and tion from the txt file given.

#### Changes to flash.par/Config/Simulation_data.F90/Simulation_init.F90  

The Addition of a variable:
- sim_inputFilenameRoot: relative path to the txt files. Stop as the end of the common prefix.  

Added the 'fromtxt' option to sim_initGeom.  

Do not forget to add every single input files to Config.  