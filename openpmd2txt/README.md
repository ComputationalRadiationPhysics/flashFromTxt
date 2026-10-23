# openpmd2txt

`openpmd2txt.py` converts a 2D PIConGPU openPMD output into the text tables
read by the `fromtxt` initialisation of the FLASH test case
(see [../Test_case/README.md](../Test_case/README.md)).

## Requirements

- Python 3
- `numpy`, `scipy`
- [`openpmd-api`](https://github.com/openPMD/openPMD-api) with the backend your
  data uses (ADIOS2 for `.bp`, HDF5 for `.h5`)

```bash
pip install -r requirements.txt
# or: conda install -c conda-forge openpmd-api numpy scipy
```

## Input

The openPMD series must contain these 2D mesh records (PIConGPU derived fields
for an electron species `e` and a hydrogen species `H`):

| Record | Meaning |
|---|---|
| `H_all_density` | ion number density |
| `e_all_density` | electron number density |
| `H_all_energyDensity` | ion kinetic energy density |
| `e_all_energyDensity` | electron kinetic energy density |

Each record is multiplied by its `unit_SI`. If the series holds several
iterations, **only the last one is used**.

`-input_file` accepts any path that `openpmd_api.Series` can open, e.g.
`simData_%T.bp` for a file-based series or `simData.bp` for a group-based one.

## Usage

```bash
python openpmd2txt.py -input_file <series> [options]
```

Example: a 40 µm × 40 µm region centred on the PIC point (10 µm, 10 µm),
resampled on 2400 × 800 points. This matches the shape and extent of the example
tables in `Test_case/input`.

```bash
python openpmd2txt.py -input_file input/simData_%T.bp \
       -output_name Collision_Laser_25fs_4P16_L_100nm_smooth \
       -xmin -10e-6 -xmax 30e-6 -ymin -10e-6 -ymax 30e-6 \
       -xcenter 10e-6 -ycenter 10e-6 \
       -nx_flash 2400 -ny_flash 800 -smooth -debug
```

> These box limits are only an illustration. Set them for your own PIC geometry.

Output is written to `./output/` **relative to the current working directory**:

```
output/<output_name>_xfile.txt   x coordinates [cm], nx_flash values, increasing
output/<output_name>_yfile.txt   y coordinates [cm], ny_flash values, increasing
output/<output_name>_dens.txt    mass density [g/cm^3]   nx_flash rows × ny_flash columns
output/<output_name>_tele.txt    electron temperature [K] nx_flash rows × ny_flash columns
output/<output_name>_tion.txt    ion temperature [K]      nx_flash rows × ny_flash columns
```

Row *i* of each 2D table corresponds to `xfile[i]` and column *j* to `yfile[j]`.
This is the layout the FLASH reader expects.

### Options

| Option | Default | Description |
|---|---|---|
| `-input_file` | *(required)* | openPMD series to read |
| `-output_name` | `Jet_test` | Base name of the output files |
| `-nx_flash`, `-ny_flash` | `800`, `800` | Number of points of the output table in x / y |
| `-xmin`, `-xmax` | `5e-6`, `15e-6` | x range (m, PIC coordinates) to resample |
| `-ymin`, `-ymax` | `5e-6`, `15e-6` | y range (m, PIC coordinates) to resample |
| `-xcenter`, `-ycenter` | `10e-6`, `10e-6` | PIC point (m) mapped to the FLASH origin (0, 0) |
| `-delta_x_pic` | `0.8e-6/96` | PIC cell size (m), assumed equal in x and y |
| `-dens_cutoff` | `10` | Cells with ion density below `dens_cutoff × n_c` are treated as vacuum |
| `-min_dens` | `1e-6` | Floor / vacuum density (g/cm³) |
| `-min_temp` | `400` | Floor / vacuum temperature (K) |
| `-density_conversion_factor` | `1.00784 × 1.66e-24` | Ion mass (g) used to convert cm⁻³ to g/cm³. **Change it for other ion species.** |
| `-smooth` | off | Apply a Gaussian filter to dens, tele and tion |
| `-smooth_level` | `2.0` | Gaussian σ, in output grid points |
| `-debug` | off | Print value ranges at each stage |

## What the conversion does

1. **Read** the four records of the last iteration and convert them to SI.
2. **Resample** onto a regular `nx_flash × ny_flash` grid covering
   `[xmin, xmax] × [ymin, ymax]` (nearest-neighbour interpolation). The PIC grid
   is rebuilt as `linspace(0, n·delta_x_pic, n)`. Points outside the PIC box get
   the minimum of the field.
3. **Shift and convert the axes**: `x_flash = (x − xcenter) × 100` (m → cm), and the same for y.
4. **Compute the temperatures** from the energy per particle:

   T_e = (2 / 3k_B) · ε_e / n_e,   T_i = (2 / 3k_B) · ε_i / n_i

   - Where `n_i < dens_cutoff · n_c` (n_c = 1.742×10²⁷ m⁻³, the critical density
     at 800 nm), density and temperatures are set to the floors `min_dens` /
     `min_temp`.
   - Temperatures are floored at `min_temp`.
   - T_i is capped at T_e.
5. **Optionally smooth** dens, tele and tion with a Gaussian filter
   (σ = `smooth_level`), keeping vacuum cells at the floor values.
6. **Convert the density** to g/cm³: `ρ = n_i[m⁻³] × 1e-6 × density_conversion_factor`,
   floored at `min_dens`.
7. **Write** the five text files.

## Feeding the result to FLASH

Copy (or link) the five files into the FLASH run directory, typically under
`input/`, and set in `flash.par`:

```
sim_initGeom          = "fromtxt"
sim_inputFilenameRoot = "input/<output_name>"
```

The FLASH domain (`xmin`…`ymax` in `flash.par`, in cm) does not have to match
the table extent exactly. Cells outside the table take the value at the
nearest edge (see the Test_case README).
