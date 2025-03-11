Conversion tool to generate FLASH initial state from FLASH outputs and openpmd_files
=====================================================================================

FLASH => FLASH  
PIConGPU  => FLASH

## QuickStart
To begin, the following dependencies are required:

#### Dependencies:

- argparse, os, numpy, scipy, yt
- openPMD-api

#### Running the conversion:

- go to the wished folder 

- Execute the python script in a shell.

#### FLASH to txt:

usage: flash2txt.py [-h] -run_directory RUN_DIRECTORY -filename FILENAME [-output_name OUTPUT_NAME]  
                    [-nx_output NX_OUTPUT] [-ny_output NY_OUTPUT] [-xmin XMIN] [-ymin YMIN] [-xmax XMAX]  
                    [-ymax YMAX] [-level LEVEL] [-debug]  

Convert FLASH output to FLASH-compatible format on a different grid.  

options:  
  -h, --help            show this help message and exit  
  -run_directory RUN_DIRECTORY  
                        Path to the directory containing FLASH input.  
  -filename FILENAME    FLASH input filename.  
  -output_name OUTPUT_NAME  
                        Base name for output files.  
  -nx_output NX_OUTPUT  Number of grid points in x direction for the output file.  
  -ny_output NY_OUTPUT  Number of grid points in y direction for the output file.  
  -xmin XMIN            Minimum x boundary in cm for the output file.  
  -ymin YMIN            Minimum y boundary in cm for the output file.  
  -xmax XMAX            Maximum x boundary in cm for the output file.  
  -ymax YMAX            Maximum y boundary in cm for the output file.  
  -level LEVEL          Refinement level in yt.  
  -debug                Enable debug output.  

#### OpenPMD to txt:

usage: openpmd2txt.py [-h] -input_file INPUT_FILE [-output_name OUTPUT_NAME] [-debug] [-nx_flash NX_FLASH]  
                      [-ny_flash NY_FLASH] [-xmin XMIN] [-ymin YMIN] [-xmax XMAX] [-ymax YMAX] [-xcenter XCENTER]  
                      [-ycenter YCENTER] [-delta_x_pic DELTA_X_PIC] [-min_temp MIN_TEMP] [-min_dens MIN_DENS]  
                      [-density_conversion_factor DENSITY_CONVERSION_FACTOR]  

Convert PIConGPU output to FLASH-compatible format.  

options:  
  -h, --help            show this help message and exit  
  -input_file INPUT_FILE  
                        Path to the PIConGPU output file.  
  -output_name OUTPUT_NAME  
                        Base name for output files.  
  -debug                Enable debug output.  
  -nx_flash NX_FLASH    Number of grid points in x direction for FLASH.  
  -ny_flash NY_FLASH    Number of grid points in y direction for FLASH.  
  -xmin XMIN            Minimum x boundary in the PIConGPU output.  
  -ymin YMIN            Minimum y boundary in the PIConGPU output.  
  -xmax XMAX            Maximum x boundary in the PIConGPU output.  
  -ymax YMAX            Maximum y boundary in the PIConGPU output.  
  -xcenter XCENTER      Center x coordinate in the PIConGPU output. Will become 0 in FLASH.  
  -ycenter YCENTER      Center y coordinate in the PIConGPU output. Will become 0 in FLASH.  
  -delta_x_pic DELTA_X_PIC  
                        Grid spacing in the PIConGPU simulation.  
  -min_temp MIN_TEMP    Minimum temperature in K for FLASH.  
  -min_dens MIN_DENS    Minimum density in g/cm³ for FLASH.  
  -density_conversion_factor DENSITY_CONVERSION_FACTOR  
                        Density conversion factor from 1/cm³ to g/cm³. Ion dependant  
