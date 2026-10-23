# flash2openpmd

`flash2openpmd.py` goes the other way from [`openpmd2txt`](../openpmd2txt/README.md).
It converts a 2D FLASH checkpoint (or plot file) into an openPMD HDF5 file of
mesh records, so that a hydrodynamic profile computed with FLASH can be used
to initialise a PIConGPU (or other PIC) simulation.

## Requirements

- Python 3
- `numpy`, `scipy`
- [`yt`](https://yt-project.org/) to read the FLASH file
- [`openpmd-api`](https://github.com/openPMD/openPMD-api) with HDF5 support

```bash
pip install -r requirements.txt
# or: conda install -c conda-forge openpmd-api yt numpy scipy
```

## Input

A FLASH 2D HDF5 file that yt can load (`*_hdf5_chk_*` or `*_hdf5_plt_cnt_*`)
containing the fields `dens`, `tele`, `tion`, `trad` and `ye`. All five must be
present. For a plot file this means they must be in the `plot_var_*` list,
which is the case in `Test_case/flash.par`. FLASH coordinates are assumed to
be in cm.

## Usage

```bash
python flash2openpmd.py -input_file <FLASH file> -output_file <name> [options]
```

Example: resample a 40 µm × 40 µm region on 4800 × 4800 points, remove the
low-density background, clear a 20-cell layer at the left and right
boundaries and rotate the profile by 30°:

```bash
python flash2openpmd.py \
    -input_file /path/to/run/lasslab_hdf5_chk_0010 \
    -output_file profile_t10 \
    -nx_pic 4800 -ny_pic 4800 \
    -xmin -20e-6 -xmax 20e-6 -ymin -20e-6 -ymax 20e-6 \
    -dens_cutoff -dens_limit 1e-5 \
    -add_cpl -cpl_cells 20 -cpl_sides left,right \
    -theta 30
```

This writes `output/profile_t10.h5`. The folder `output/` is **relative to the
current working directory**, unless `-output_folder` is given.

### Options

| Option | Default | Description |
|---|---|---|
| `-input_file` | *(required)* | FLASH checkpoint / plot file |
| `-output_file` | *(required)* | Output file name, without `.h5` |
| `-output_folder` | `output` | Output directory (created if missing) |
| `-nx_pic`, `-ny_pic` | `4800`, `4800` | Number of points of the target grid in x / y |
| `-xmin`, `-xmax` | `-20e-6`, `20e-6` | x extent of the target grid (m, FLASH coordinates converted to m) |
| `-ymin`, `-ymax` | `-20e-6`, `20e-6` | y extent of the target grid (m) |
| `-dens_cutoff` | off | Set **all** fields to 0 where the density is ≤ `dens_limit` |
| `-dens_limit` | `1e-5` | Cutoff threshold, in **g/cm³** |
| `-add_cpl` | off | Set the density to 0 in a layer along the selected boundaries (for the PIC absorbing layer) |
| `-cpl_cells` | `0` | Depth of that layer, in target-grid cells |
| `-cpl_sides` | `all` | Comma-separated list of `left`, `right`, `bottom`, `top`, `all` (left/right = x, bottom/top = y) |
| `-cpl_zero_all` | off | Also set `tele`, `tion`, `trad` and `ye` to 0 in the layer |
| `-theta` | `0` | Rotate the data by this angle (degrees, counterclockwise). The grid axes stay fixed |
| `-rot_center_x`, `-rot_center_y` | grid centre | Rotation centre (m) |

## What the conversion does

1. **Load** the file with yt and build a uniform covering grid at AMR level 4
   (`level` argument of `read4Flash`, not exposed on the command line).
2. **Convert units**: density g/cm³ → kg/m³, FLASH coordinates cm → m.
   Temperatures stay in K, `ye` is dimensionless.
3. **Resample** the `z = 0` slice bilinearly onto the `nx_pic × ny_pic` grid
   covering `[xmin, xmax] × [ymin, ymax]` (end points included).
   With `-theta`, the target points are rotated back by `theta` around the rotation
   centre before sampling, so the data is rotated in a single interpolation step.
   Points outside the FLASH domain are set to 0.
4. **Optionally** apply the density cutoff, then the boundary layer.
5. **Write** iteration 0 of an openPMD series with five `float32` scalar meshes:

   | Mesh | FLASH field | Unit |
   |---|---|---|
   | `dens` | `dens` | kg/m³ |
   | `tele` | `tele` | K |
   | `tion` | `tion` | K |
   | `trad` | `trad` | K |
   | `ye` | `ye` | dimensionless |

   Each mesh has `gridGlobalOffset = [xmin, ymin]`,
   `gridSpacing = [(xmax - xmin) / (nx_pic - 1), (ymax - ymin) / (ny_pic - 1)]`
   (m), `gridUnitSI = 1`, `axisLabels = ["x", "y"]` and `dataOrder = "C"`.

> **Older files:** earlier versions of this script applied an extra cm → m
> factor to values already in meters, so their `gridSpacing` is 100× too
> small. Regenerate them, or multiply `gridSpacing` by 100, if your PIC
> reader uses it.

## Known limitations

Check these before relying on the output:

- **Axis order of the stored data.** The arrays are passed to openPMD-api as
  Fortran-ordered (`np.asfortranarray`), but openPMD-api writes the memory
  buffer as C order. As a result:
  - for a **square** grid (`nx_pic == ny_pic`) the stored array is the
    transpose of the computed one: the first index is y, the second x,
    although `axisLabels` says `["x", "y"]`;
  - for a **non-square** grid the data is scrambled. **Use
    `nx_pic == ny_pic`.**
- **AMR level: keep `level` at 4 (or 2).** yt builds the covering grid with
  `domain_dimensions × level**2` cells, while the cell size at that level is
  set by the real refinement factor `2**level`. The two agree only for levels
  2 and 4 (and 0). For any other level the covering grid does not match the
  FLASH domain: it is too large (e.g. level 3: 9× instead of 8×) and wraps
  around periodically, or too small, and the x/y coordinates assigned to the
  data are wrong. `level` is an argument of `read4Flash` and is not exposed on
  the command line, so this only matters if you change it in the code.
- **Half-cell shift.** The covering-grid values are cell-centred, but the
  script places them on `linspace(domain_left, domain_right, n)`, i.e. the
  first and last values sit on the domain edges instead of half a cell inside.
  The profile is therefore shifted and stretched by up to half a cell of the
  covering grid (`Δx_base / 2**(level+1)` at the domain edges, 0 in the
  centre). At level 4 this is 1/32 of a base-grid cell, usually negligible
  compared with the PIC grid, but noticeable for sharp gradients near the
  domain boundaries.
- Only the `z = 0` slice is used.
