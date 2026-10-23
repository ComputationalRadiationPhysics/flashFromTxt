import argparse
import os
import yt
import numpy as np
from scipy import interpolate
import openpmd_api as io
from scipy import constants as sc


class Convert:
    def __init__(self, run_directory, filename):
        """Initialize the conversion process by specifying the FLASH output directory and filename."""
        if run_directory is None:
            raise ValueError('The run_directory parameter cannot be None!')
        self.path_to_sim = os.path.join(run_directory, filename)

    def read4Flash(
        self,
        nx_pic,
        ny_pic,
        xmin,
        ymin,
        xmax,
        ymax,
        dens_cutoff,
        dens_limit,
        level=4,
        add_cpl=False,
        cpl_cells=0,
        cpl_sides="all",
        cpl_zero_all=False,
        theta=0.0,
        rot_center_x=None,
        rot_center_y=None,
    ):
        """Read FLASH output and interpolate to the desired resolution.

        - Loads a FLASH checkpoint via yt.
        - Converts from CGS to SI units where appropriate.
        - Interpolates onto a uniform (nx_pic, ny_pic) grid in the x-y plane at z=0,
          optionally rotating the map by an angle theta (in degrees, counterclockwise)
          about a rotation center (defaults to the center of the new grid).
          The rotation is performed directly at interpolation time: the target grid
          coordinates are inverse-rotated and the FLASH data is sampled at these
          rotated positions. This avoids a second interpolation pass and its
          associated smoothing artifacts. Points that map outside the FLASH domain
          are filled with 0.
        - Optionally applies a density cutoff and a CPL boundary layer.
        """
        ds = yt.load(self.path_to_sim)

        if level == 0:
            scale = 1
        else:
            # AMR refinement scaling (same convention as in the original script)
            scale = np.array([level**2, level**2, 1])

        ds.force_periodicity()
        all_data = ds.covering_grid(
            level=level,
            left_edge=ds.domain_left_edge,
            dims=ds.domain_dimensions * scale,
        )

        # FLASH -> SI units
        density = all_data["flash", "dens"] * 1e3  # g/cm^3 -> kg/m^3
        t_electron = all_data["flash", "tele"]
        t_ion = all_data["flash", "tion"]
        t_rad = all_data["flash", "trad"]
        ye = all_data["flash", "ye"]

        # Some basic diagnostics on the raw slice at z = 0
        print("Dens range: [" + str(np.min(density[:, :, 0].v)) + ", " + str(np.max(density[:, :, 0].v)) + "]")
        print("Tele range: [" + str(np.min(t_electron[:, :, 0].v)) + ", " + str(np.max(t_electron[:, :, 0].v)) + "]")
        print("Tion range: [" + str(np.min(t_ion[:, :, 0].v)) + ", " + str(np.max(t_ion[:, :, 0].v)) + "]")
        print("Ye range: [" + str(np.min(ye[:, :, 0].v)) + ", " + str(np.max(ye[:, :, 0].v)) + "]")

        # FLASH domain coordinates: assume cm, convert to m
        x = np.linspace(
            ds.domain_left_edge[0].v * 1e-2,
            ds.domain_right_edge[0].v * 1e-2,
            density.shape[0],
        )
        y = np.linspace(
            ds.domain_left_edge[1].v * 1e-2,
            ds.domain_right_edge[1].v * 1e-2,
            density.shape[1],
        )

        # Target grid (axis values stay as provided by the user, in meters)
        xnew = np.linspace(xmin, xmax, nx_pic, endpoint=True)
        ynew = np.linspace(ymin, ymax, ny_pic, endpoint=True)
        xx, yy = np.meshgrid(xnew, ynew, indexing="ij")

        # Rotation of the map by theta (degrees, counterclockwise) about the
        # rotation center. To rotate the *data* by +theta on the fixed target
        # grid, we sample the source data at the inverse-rotated coordinates:
        #   (x_src, y_src) = R(-theta) . (x_tgt - c) + c
        if rot_center_x is None:
            rot_center_x = 0.5 * (xmin + xmax)
        if rot_center_y is None:
            rot_center_y = 0.5 * (ymin + ymax)

        theta_rad = np.deg2rad(theta)
        cos_t = np.cos(theta_rad)
        sin_t = np.sin(theta_rad)

        dx_c = xx - rot_center_x
        dy_c = yy - rot_center_y
        x_src = rot_center_x + cos_t * dx_c + sin_t * dy_c
        y_src = rot_center_y - sin_t * dx_c + cos_t * dy_c

        if theta != 0.0:
            print(
                "Rotating map by theta = "
                + str(theta)
                + " deg about center ("
                + str(rot_center_x)
                + ", "
                + str(rot_center_y)
                + ") m"
            )

        # Interpolate to new grid (z=0 slice), sampling at the rotated coordinates.
        # Points outside the FLASH domain are filled with 0.
        f_dens = interpolate.interpn(
            (x, y), density[:, :, 0].v, (x_src, y_src),
            method="linear", bounds_error=False, fill_value=0.0
        )
        f_tele = interpolate.interpn(
            (x, y), t_electron[:, :, 0].v, (x_src, y_src),
            method="linear", bounds_error=False, fill_value=0.0
        )
        f_tion = interpolate.interpn(
            (x, y), t_ion[:, :, 0].v, (x_src, y_src),
            method="linear", bounds_error=False, fill_value=0.0
        )
        f_trad = interpolate.interpn(
            (x, y), t_rad[:, :, 0].v, (x_src, y_src),
            method="linear", bounds_error=False, fill_value=0.0
        )
        f_ye = interpolate.interpn(
            (x, y), ye[:, :, 0].v, (x_src, y_src),
            method="linear", bounds_error=False, fill_value=0.0
        )

        # Optional density cutoff
        if dens_cutoff:
            # dens_limit is in g/cm^3 -> convert to kg/m^3
            mask = f_dens <= dens_limit * 1e3
            f_dens[mask] = 0.0
            f_tele[mask] = 0.0
            f_tion[mask] = 0.0
            f_trad[mask] = 0.0
            f_ye[mask] = 0.0

        # Optional CPL layer: zero density (and optionally other fields) near boundaries
        if add_cpl and cpl_cells > 0:
            # Normalize sides
            if isinstance(cpl_sides, str):
                sides = {s.strip().lower() for s in cpl_sides.split(',')}
            else:
                sides = set()
            if "all" in sides or not sides:
                sides = {"left", "right", "top", "bottom"}

            # Clamp depth to array size
            nx, ny = f_dens.shape
            d = max(0, min(cpl_cells, nx, ny))

            def apply_cpl(arr, zero_this_field):
                if not zero_this_field or d <= 0:
                    return arr
                if "left" in sides:
                    arr[:d, :] = 0.0
                if "right" in sides:
                    arr[-d:, :] = 0.0
                if "bottom" in sides:
                    arr[:, :d] = 0.0
                if "top" in sides:
                    arr[:, -d:] = 0.0
                return arr

            # Always zero density in the CPL region
            apply_cpl(f_dens, True)

            # Optionally zero all other fields in the same region
            if cpl_zero_all:
                apply_cpl(f_tele, True)
                apply_cpl(f_tion, True)
                apply_cpl(f_trad, True)
                apply_cpl(f_ye, True)

        # Grid spacing written to openPMD, in meters (xmin/xmax are already in m)
        grid_spacing_x = (xmax - xmin) / (nx_pic - 1)
        grid_spacing_y = (ymax - ymin) / (ny_pic - 1)

        print("fdens range: [" + str(np.min(f_dens[:, :])) + ", " + str(np.max(f_dens[:, :])) + "]")
        print("ftele range: [" + str(np.min(f_tele[:, :])) + ", " + str(np.max(f_tele[:, :])) + "]")
        print("ftion range: [" + str(np.min(f_tion[:, :])) + ", " + str(np.max(f_tion[:, :])) + "]")
        print("fye range: [" + str(np.min(f_ye[:, :])) + ", " + str(np.max(f_ye[:, :])) + "]")

        return (
            f_dens.astype("float32"),
            f_tele.astype("float32"),
            f_tion.astype("float32"),
            f_trad.astype("float32"),
            f_ye.astype("float32"),
            grid_spacing_x,
            grid_spacing_y,
        )

    def write2openpmd(
        self,
        density,
        t_electron,
        t_ion,
        t_rad,
        ye,
        xmin,
        ymin,
        grid_spacing_x,
        grid_spacing_y,
        output_filename,
        output_folder,
    ):
        """Write the converted data to an openPMD-compatible HDF5 file."""
        os.makedirs(output_folder, exist_ok=True)
        output_path = os.path.join(output_folder, output_filename + ".h5")

        series_out = io.Series(output_path, io.Access.create)
        it = series_out.iterations[0]

        def store_mesh_data(name, data, unit_dim):
            mesh = it.meshes[name]
            mesh_mrc = mesh[io.Mesh_Record_Component.SCALAR]
            mesh_mrc.reset_dataset(io.Dataset(data.dtype, data.shape))
            mesh.set_attribute("dataOrder", "C")
            mesh.set_attribute("axisLabels", ["x", "y"])
            mesh.set_attribute("unitDimension", unit_dim)
            mesh.set_attribute("gridGlobalOffset", [xmin, ymin])
            mesh.set_attribute("gridSpacing", [grid_spacing_x, grid_spacing_y])
            mesh.set_attribute("gridUnitSI", 1.0)
            mesh_mrc.store_chunk(data)

        # Diagnostics on the final data that will be written
        print("dens range: [" + str(np.min(density)) + ", " + str(np.max(density)) + "]")
        print("tele range: [" + str(np.min(t_electron)) + ", " + str(np.max(t_electron)) + "]")
        print("tion range: [" + str(np.min(t_ion)) + ", " + str(np.max(t_ion)) + "]")
        print("ye range: [" + str(np.min(ye)) + ", " + str(np.max(ye)) + "]")

        # Density (kg/m^3)
        store_mesh_data(
            "dens",
            np.asfortranarray(density),
            [-3.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        )
        # Temperatures (thermodynamic temperature dimension)
        store_mesh_data(
            "tele",
            np.asfortranarray(t_electron),
            [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
        )
        store_mesh_data(
            "tion",
            np.asfortranarray(t_ion),
            [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
        )
        store_mesh_data(
            "trad",
            np.asfortranarray(t_rad),
            [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
        )
        # Electron fraction (dimensionless)
        store_mesh_data(
            "ye",
            np.asfortranarray(ye),
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        )

        series_out.flush()
        del series_out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert FLASH output files to openPMD format, with optional in-plane rotation."
    )
    parser.add_argument(
        "-input_file",
        required=True,
        help="Path to the FLASH output file.",
    )
    parser.add_argument(
        "-nx_pic",
        type=int,
        default=4800,
        help="Resolution in x-direction of the OpenPMD output file.",
    )
    parser.add_argument(
        "-ny_pic",
        type=int,
        default=4800,
        help="Resolution in y-direction of the OpenPMD output file.",
    )

    parser.add_argument(
        "-xmin",
        type=float,
        default=-20e-6,
        help="Minimum x boundary of the OpenPMD output file (in meters).",
    )
    parser.add_argument(
        "-ymin",
        type=float,
        default=-20e-6,
        help="Minimum y boundary of the OpenPMD output file (in meters).",
    )
    parser.add_argument(
        "-xmax",
        type=float,
        default=20e-6,
        help="Maximum x boundary of the OpenPMD output file (in meters).",
    )
    parser.add_argument(
        "-ymax",
        type=float,
        default=20e-6,
        help="Maximum y boundary of the OpenPMD output file (in meters).",
    )

    parser.add_argument(
        "-dens_cutoff",
        action="store_true",
        help="Apply density cutoff.",
    )
    parser.add_argument(
        "-dens_limit",
        type=float,
        default=1e-5,
        help="Density limit for cutoff (in g/cm^3).",
    )

    parser.add_argument(
        "-output_file",
        required=True,
        help="Output filename for openPMD data (without extension).",
    )
    parser.add_argument(
        "-output_folder",
        type=str,
        default="output",
        help="Output folder for openPMD data.",
    )

    parser.add_argument(
        "-add_cpl",
        action="store_true",
        help="Enable CPL layer (zero-out density near selected boundaries).",
    )
    parser.add_argument(
        "-cpl_cells",
        type=int,
        default=0,
        help="CPL depth in cells (integer).",
    )
    parser.add_argument(
        "-cpl_sides",
        type=str,
        default="all",
        help="Comma-separated sides to apply CPL: left,right,top,bottom,all.",
    )
    parser.add_argument(
        "-cpl_zero_all",
        action="store_true",
        help="Also zero tele/tion/trad/ye within CPL region.",
    )

    # Rotation options: rotate the map by theta on the fixed target grid
    parser.add_argument(
        "-theta",
        type=float,
        default=0.0,
        help="Rotation angle of the map in degrees (counterclockwise, "
             "mathematically positive). The target grid axes stay unchanged; "
             "only the data is rotated.",
    )
    parser.add_argument(
        "-rot_center_x",
        type=float,
        default=None,
        help="x-coordinate of the rotation center in meters "
             "(default: center of the target grid, (xmin+xmax)/2).",
    )
    parser.add_argument(
        "-rot_center_y",
        type=float,
        default=None,
        help="y-coordinate of the rotation center in meters "
             "(default: center of the target grid, (ymin+ymax)/2).",
    )

    args = parser.parse_args()

    run_directory, filename = os.path.split(args.input_file)
    ts = Convert(run_directory, filename)

    density, t_electron, t_ion, t_rad, ye, grid_spacing_x, grid_spacing_y = ts.read4Flash(
        args.nx_pic,
        args.ny_pic,
        args.xmin,
        args.ymin,
        args.xmax,
        args.ymax,
        args.dens_cutoff,
        args.dens_limit,
        add_cpl=args.add_cpl,
        cpl_cells=args.cpl_cells,
        cpl_sides=args.cpl_sides,
        cpl_zero_all=args.cpl_zero_all,
        theta=args.theta,
        rot_center_x=args.rot_center_x,
        rot_center_y=args.rot_center_y,
    )

    ts.write2openpmd(
        density,
        t_electron,
        t_ion,
        t_rad,
        ye,
        args.xmin,
        args.ymin,
        grid_spacing_x,
        grid_spacing_y,
        args.output_file,
        args.output_folder,
    )

    print("Conversion completed successfully.")