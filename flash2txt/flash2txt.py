import argparse
import numpy as np
import os
from scipy import interpolate
import yt

# Constants
kb = 1.380649e-23  # Boltzmann constant in J/K

def read_flash_data(run_directory, filename, nx_output, ny_output, xmin_output, ymin_output, xmax_output, ymax_output, level, debug):
    """
    Reads FLASH output from an HDF5 file and interpolates the density and temperature data.
    """
    path_to_sim = os.path.join(run_directory, filename)
    ds = yt.load(path_to_sim)
    
    if level == 0:
        scale = 1
    else:
        scale = np.array([level**2, level**2, 1])
    
    ds.force_periodicity()
    all_data = ds.covering_grid(level=level, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions * scale)
    
    density = all_data["flash", "dens"].v
    t_electron = all_data["flash", "tele"].v
    t_ion = all_data["flash", "tion"].v
    t_rad = all_data["flash", "trad"].v
    
    if debug:
        print("density range: [", np.min(density), ",", np.max(density), "]")
    
    y_max, x_max, _ = ds.domain_right_edge.v
    y_min, x_min, _ = ds.domain_left_edge.v
    
    nx, ny, _ = density.shape
    x = np.linspace(x_min, x_max, nx)
    y = np.linspace(y_min, y_max, ny)
    
    xnew = np.linspace(xmin_output, xmax_output, nx_output)
    ynew = np.linspace(ymin_output, ymax_output, ny_output)
    xx, yy = np.meshgrid(xnew, ynew)
    
    f_dens = interpolate.interpn((y, x), density[:, :, 0], (yy, xx), method="linear", bounds_error=False, fill_value=np.min(density))
    f_tele = interpolate.interpn((y, x), t_electron[:, :, 0], (yy, xx), method="linear", bounds_error=False, fill_value=np.min(t_electron))
    f_tion = interpolate.interpn((y, x), t_ion[:, :, 0], (yy, xx), method="linear", bounds_error=False, fill_value=np.min(t_ion))
    f_trad = interpolate.interpn((y, x), t_rad[:, :, 0], (yy, xx), method="linear", bounds_error=False, fill_value=np.min(t_rad))
    
    return f_dens, f_tele, f_tion, f_trad, xnew, ynew

def save_data(outputname, x, y, density, t_electron, t_ion, t_rad):
    """ Saves processed data to output files. """
    os.makedirs("output", exist_ok=True)
    np.savetxt(f"output/{outputname}_xfile.txt", x)
    np.savetxt(f"output/{outputname}_yfile.txt", y)
    np.savetxt(f"output/{outputname}_dens.txt", density)
    np.savetxt(f"output/{outputname}_tele.txt", t_electron)
    np.savetxt(f"output/{outputname}_tion.txt", t_ion)
    np.savetxt(f"output/{outputname}_trad.txt", t_rad)

def main():
    parser = argparse.ArgumentParser(description="Convert FLASH output to FLASH-compatible format on a different grid.")
    parser.add_argument("-run_directory", type=str, required=True, help="Path to the directory containing FLASH input.")
    parser.add_argument("-filename", type=str, required=True, help="FLASH input filename.")
    parser.add_argument("-output_name", type=str, default="flash_conversion", help="Base name for output files.")
    parser.add_argument("-nx_output", type=int, default=400, help="Number of grid points in x direction for the output file.")
    parser.add_argument("-ny_output", type=int, default=400, help="Number of grid points in y direction for the output file.")
    parser.add_argument("-xmin", type=float, default=-20e-4, help="Minimum x boundary in cm for the output file.")
    parser.add_argument("-ymin", type=float, default=-20e-4, help="Minimum y boundary in cm for the output file.")
    parser.add_argument("-xmax", type=float, default=20e-4, help="Maximum x boundary in cm for the output file.")
    parser.add_argument("-ymax", type=float, default=20e-4, help="Maximum y boundary in cm for the output file.")
    parser.add_argument("-level", type=int, default=4, help="Refinement level in yt.")
    parser.add_argument("-debug", action="store_true", help="Enable debug output.")
    
    args = parser.parse_args()
    
    density, t_electron, t_ion, t_rad, x, y = read_flash_data(
        args.run_directory, args.filename, args.nx_output, args.ny_output, args.xmin, args.ymin, args.xmax, args.ymax, 
        args.dens_cutoff, args.dens_limit, args.level, args.debug
    )
    
    save_data(args.output_name, x, y, density, t_electron, t_ion, t_rad)
    
    print("Conversion complete. Output files are stored in the 'output' directory.")

if __name__ == "__main__":
    main()
