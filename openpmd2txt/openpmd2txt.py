import argparse
import numpy as np
import scipy.constants as const
import os
from scipy import interpolate
import openpmd_api as io

# Constants
kb = const.Boltzmann

def read_flash_data(file_path, delta_x_pic, nx_flash, ny_flash, xmin_flash, ymin_flash, xmax_flash, ymax_flash, debug):
    """
    Reads FLASH output from an openPMD file and interpolates the density and temperature data.
    """
    ds = io.Series(file_path, access=io.Access_Type.read_only)

    i=0
    
    for it in ds.iterations: 
        i = ds.iterations[it]
    
    density_h = i.meshes["H_all_density"][io.Mesh_Record_Component.SCALAR]
    density = density_h[:,:]
    density_SI = density_h.unit_SI
    
    edensity_h = i.meshes["e_all_density"][io.Mesh_Record_Component.SCALAR]
    edensity = edensity_h[:,:]
    edensity_SI = edensity_h.unit_SI
    
    t_electron_h = i.meshes["e_all_energyDensity"][io.Mesh_Record_Component.SCALAR]
    t_electron = t_electron_h[:,:]
    t_electron_SI = t_electron_h.unit_SI
    
    t_ion_h = i.meshes["H_all_energyDensity"][io.Mesh_Record_Component.SCALAR]
    t_ion = t_ion_h[:,:]
    t_ion_SI = t_ion_h.unit_SI

    ds.flush()
    
    density *= density_SI
    edensity *= edensity_SI
    t_electron *= t_electron_SI
    t_ion *= t_ion_SI
    
    if debug:
        print("density range: [", np.min(density), ",", np.max(density), "]")
        print("edensity range: [", np.min(edensity), ",", np.max(edensity), "]")
        print("Ed_ele range: [", np.min(t_electron), ",", np.max(t_electron), "]")
        print("Ed_ion range: [", np.min(t_ion), ",", np.max(t_ion), "]")
    
    nx, ny = density.shape
    x = np.linspace(0, nx * delta_x_pic, nx)
    y = np.linspace(0, ny * delta_x_pic, ny)
    
    xnew = np.linspace(xmin_flash, xmax_flash, nx_flash)
    ynew = np.linspace(ymin_flash, ymax_flash, ny_flash)
    
    xx, yy = np.meshgrid(xnew, ynew)
    
    f_dens = interpolate.interpn((x, y), density, (xx, yy), method="nearest", bounds_error=False, fill_value=np.min(density))
    f_edens = interpolate.interpn((x, y), edensity, (xx, yy), method="nearest", bounds_error=False, fill_value=np.min(edensity))
    f_tele = interpolate.interpn((x, y), t_electron, (xx, yy), method="nearest", bounds_error=False, fill_value=np.min(t_electron))
    f_tion = interpolate.interpn((x, y), t_ion, (xx, yy), method="nearest", bounds_error=False, fill_value=np.min(t_ion))
    
    if debug:
        print("After interpolation: density range: [", np.min(f_dens), ",", np.max(f_dens), "]")
        print("After interpolation: edensity range: [", np.min(f_edens), ",", np.max(f_edens), "]")
        print("After interpolation: Ed_ele range: [", np.min(f_tele), ",", np.max(f_tele), "]")
        print("After interpolation: Ed_ion range: [", np.min(f_tion), ",", np.max(f_tion), "]")
    
    return f_dens, f_edens, f_tele, f_tion, xnew, ynew

def process_data(density, edensity, t_electron, t_ion, x, y, xcenter, ycenter, outputname, min_temp, min_dens, density_conversion_factor, debug):
    """
    Processes and converts the data to FLASH-compatible format.
    """
    x = (x - xcenter) * 1e2
    y = (y - ycenter) * 1e2
    
    with np.errstate(divide='ignore', invalid='ignore'):
        valid_both = (density != 0) & (edensity != 0)
        tele = np.where(valid_both, (t_electron / edensity) * 2 / (3 * kb), min_temp)
        tion = np.where(density != 0, (t_ion / density) * 2 / (3 * kb), min_temp)
    
    density = np.maximum(density * 1e-6 * density_conversion_factor, min_dens)
    
    if debug:
        print("After conversion: density range: [", np.min(density), ",", np.max(density), "]")
        print("tele range (K): [", np.min(tele), ",", np.max(tele), "]")
        print("tion range (K): [", np.min(tion), ",", np.max(tion), "]")
        print("tele range (eV): [", np.min(tele)/11606, ",", np.max(tele)/11606, "]")
        print("tion range (eV): [", np.min(tion)/11606, ",", np.max(tion)/11606, "]")
    
    os.makedirs("output", exist_ok=True)
    np.savetxt(f"output/{outputname}_xfile.txt", x)
    np.savetxt(f"output/{outputname}_yfile.txt", y)
    np.savetxt(f"output/{outputname}_dens.txt", density)
    np.savetxt(f"output/{outputname}_tion.txt", tion)
    np.savetxt(f"output/{outputname}_tele.txt", tele)

def main():
    parser = argparse.ArgumentParser(description="Convert PIConGPU output to FLASH-compatible format.")
    parser.add_argument("-input_file", type=str, required=True, help="Path to the PIConGPU output file.")
    parser.add_argument("-output_name", type=str, default="Jet_test", help="Base name for output files.")
    parser.add_argument("-debug", action="store_true", help="Enable debug output.")
    parser.add_argument("-nx_flash", type=int, default=800, help="Number of grid points in x direction for FLASH.")
    parser.add_argument("-ny_flash", type=int, default=800, help="Number of grid points in y direction for FLASH.")
    parser.add_argument("-xmin", type=float, default=5.0e-6, help="Minimum x boundary in the PIConGPU output.")
    parser.add_argument("-ymin", type=float, default=5.0e-6, help="Minimum y boundary in the PIConGPU output.")
    parser.add_argument("-xmax", type=float, default=15.0e-6, help="Maximum x boundary in the PIConGPU output.")
    parser.add_argument("-ymax", type=float, default=15.0e-6, help="Maximum y boundary in the PIConGPU output.")
    parser.add_argument("-xcenter", type=float, default=10e-6, help="Center x coordinate  in the PIConGPU output. Will become 0 in FLASH.")
    parser.add_argument("-ycenter", type=float, default=10e-6, help="Center y coordinate  in the PIConGPU output. Will become 0 in FLASH.")
    parser.add_argument("-delta_x_pic", type=float, default=0.8e-6 / 96., help="Grid spacing in the PIConGPU simulation.")
    parser.add_argument("-min_temp", type=float, default=400, help="Minimum temperature in K for FLASH.")
    parser.add_argument("-min_dens", type=float, default=1e-6, help="Minimum density in g/cm³ for FLASH.")
    parser.add_argument("-density_conversion_factor", type=float, default=1.00784 * 1.66e-24, help="Density conversion factor from 1/cm³ to g/cm³. Ion dependant")
    
    args = parser.parse_args()
    
    density, edensity, t_electron, t_ion, x, y = read_flash_data(
        args.input_file, args.delta_x_pic, args.nx_flash, args.ny_flash, args.xmin, args.ymin, args.xmax, args.ymax, args.debug
    )
    
    process_data(density, edensity, t_electron, t_ion, x, y, args.xcenter, args.ycenter, 
                 args.output_name, args.min_temp, args.min_dens, args.density_conversion_factor, args.debug
                 )
    
    print("Conversion complete. Output files are stored in the 'output' directory.")

if __name__ == "__main__":
    main()
