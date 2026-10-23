# Test_case: FLASH initialised from text tables

A FLASH 4 2D, three-temperature (3T) simulation derived from the `LaserSlab`
example. Besides `slab` and `circle2dpolar`, it adds a third initial geometry,
**`sim_initGeom = "fromtxt"`**. This mode reads density, electron temperature and
ion temperature from text tables (produced by
[`openpmd2txt`](../openpmd2txt/README.md)) and bilinearly interpolates them onto
every AMR cell.

## Files

| File | Role |
|---|---|
| `Config` | Setup requirements and runtime parameters. Adds `"fromtxt"` to `sim_initGeom` and adds `sim_inputFilenameRoot`. |
| `Simulation_data.F90` | Module holding the runtime parameters **and the tabulated initial state** (axes and 2D tables). |
| `Simulation_init.F90` | Reads the tables **once** on the master rank and broadcasts them with MPI. |
| `Simulation_initBlock.F90` | Fills each block by bilinear interpolation from the in-memory tables. |
| `Simulation_freeUserArrays.F90` | Frees the tables once the initial grid has been built. |
| `Simulation_adjustEvolution.F90`, `Simulation_customizeProlong.F90` | Unchanged from LaserSlab (empty hook / standard prolongation). |
| `Makefile` | Adds `Simulation_data.o` to the Simulation unit. |
| `flash.par` | Runtime parameters for the test (H target, 40 µm × 40 µm box, 200 fs). |
| `H_long.cn4` | IONMIX4 EOS / opacity table for hydrogen. |
| `input/*.txt.gz` | Example tables (gzipped). |
| `start_rosi`, `start_rosi_parallel` | Build-and-submit scripts for the HZDR *rosi* cluster (serial / MPI). |
| `submit_slurm_rosi.sh`, `submit_slurm_multinode.sh` | Matching Slurm job scripts. |

## Installation

1. Copy this folder into your FLASH 4 source tree as a simulation directory, e.g.

   ```bash
   cp -r Test_case $FLASH_DIR/source/Simulation/SimulationMain/FromTxt
   ```

   (The `start_rosi*` scripts assume this depth: they `cd ../../../..` to reach
   the FLASH root.)

2. Decompress the example tables, or put your own tables in `input/`:

   ```bash
   gunzip -k input/*.txt.gz
   ```

3. Run setup and compile from the FLASH root:

   ```bash
   ./setup FromTxt -auto -2d +cartesian -nxb=16 -nyb=16 -maxblocks=22000 \
           +hdf5typeio species=cham,targ +mtmmmt +uhd3t +mgd mgd_meshgroups=1
   cd object && make -j
   ```

   `species=cham,targ` is required. `-maxblocks` must hold the full refined
   tree on each rank (≈ 22 000 for one rank with `nblockx = nblocky = 16` and
   `lrefine_max = 4`; much less per rank in parallel).

4. Copy `flash4`, `flash.par`, `H_long.cn4` and `input/` into a run directory
   and start with `mpirun -np <N> ./flash4`. The `start_rosi*` scripts do steps
   3 and 4 automatically and submit to Slurm. **Edit `run_dir` for your system.**
   They source `~/FLASH_rome.sh` and need the rosi `Makefile.h`. See
   [../rosi/README.md](../rosi/README.md).

## Runtime parameters for `fromtxt`

```
sim_initGeom          = "fromtxt"
sim_inputFilenameRoot = "input/Collision_Laser_25fs_4P16_L_100nm_smooth"
```

FLASH reads `<root>_xfile.txt`, `<root>_yfile.txt`, `<root>_dens.txt`,
`<root>_tele.txt` and `<root>_tion.txt`, relative to the run directory.

The chamber values (`sim_rhoCham`, `sim_teleCham`, `sim_tionCham`) still
matter: they set the lower bound of the "outside" values (see below). Make
the FLASH domain (`xmin`, `xmax`, `ymin`, `ymax`, in cm) cover the region of interest.

## Table format

| File | Content |
|---|---|
| `_xfile.txt` | `nx` x coordinates in cm, one per line, **strictly increasing** |
| `_yfile.txt` | `ny` y coordinates in cm, one per line, **strictly increasing** |
| `_dens.txt`  | density (g/cm³): `nx` lines of `ny` values; line *i* ↔ `x(i)`, column *j* ↔ `y(j)` |
| `_tele.txt`  | electron temperature (K), same layout |
| `_tion.txt`  | ion temperature (K), same layout |

The axes need not be uniformly spaced. The example tables are 2400 × 800 and
cover −20 µm … +20 µm in both directions.

## The modifications in detail

### 1. Load once, broadcast: `Simulation_init.F90`

When `sim_initGeom == "fromtxt"` and the run is **not a restart**,
`Simulation_init` calls the internal routine `sim_loadTxtTables`:

- **Master rank only** reads the files:
  - `sim_readVector` reads the x and y axes. It first counts the lines, so the
    table size does not have to be set anywhere; it is taken from the axis files.
  - `sim_checkAxis` aborts if an axis is not strictly increasing, because the
    interpolation relies on a binary search.
  - `sim_readMatrix` reads each 2D table as `nx` records of `ny` values and
    aborts with the failing row number on a read error.
- The table sizes are broadcast with `MPI_Bcast`. The other ranks allocate
  their arrays, then the axes and the three tables are broadcast.
- One-time preprocessing on all ranks:
  - negative values are clamped to 0;
  - the outside values are set: `sim_densOutside = max(minval(dens), sim_rhoCham)`,
    and likewise for tele and tion.
- The master rank writes `loaded <nx> x <ny> tables from <root>` to stdout and to the log file.

`fromtxt` is implemented for 2D only (`NDIM /= 2` aborts).

On a **restart**, the initial state comes from the checkpoint, so the tables
are not read.

### 2. Shared storage: `Simulation_data.F90`

New module variables:

```fortran
integer :: sim_nxFile, sim_nyFile
real, allocatable :: sim_xFile(:), sim_yFile(:)                 ! axes (cm)
real, allocatable :: sim_densFile(:,:), sim_teleFile(:,:), sim_tionFile(:,:)  ! (nx, ny)
real :: sim_densOutside, sim_teleOutside, sim_tionOutside
logical, parameter :: sim_txtFloorOutside = .false.
character(len=MAX_STRING_LENGTH) :: sim_inputFilenameRoot
```

### 3. Bilinear interpolation: `Simulation_initBlock.F90`

For a `fromtxt` block:

1. For each **column** *i* of the block, `locate` finds by binary search the
   table index `ix(i)` and the weight `tx(i) ∈ [0, 1]` such that
   `x ≈ (1 − t)·x(ix) + t·x(ix+1)`. The same is done for each **row** *j*
   (`iy`, `ty`). This costs O((nxb + nyb) log N) per block instead of a search
   per cell.
2. Each cell is then a 4-point weighted sum:

   ```
   f(x, y) = (1−tx)(1−ty) f(ix, iy)   + tx(1−ty) f(ix+1, iy)
           + (1−tx) ty    f(ix, iy+1) + tx ty    f(ix+1, iy+1)
   ```

   applied to dens, tele and tion.
3. **Cells outside the table extent** are handled by `sim_txtFloorOutside`
   (a compile-time parameter in `Simulation_data.F90`):
   - `.false.` (default): the coordinate is clamped to the edge, so the cell
     takes the value at the nearest table edge, with no discontinuity;
   - `.true.`: the cell gets `sim_densOutside`, `sim_teleOutside`,
     `sim_tionOutside`.
4. The values are written as `DENS_VAR`, `TEMP_VAR`, `TELE_VAR`, `TION_VAR`.
   Since there is no radiation table, **T_rad = T_e** (passed through
   `RadTrans_mgdEFromT`).
5. Every cell is assigned to the **target species** (`targ`), so the target
   material parameters (`ms_targA`, `ms_targZ`, `eos_targ*`) describe the
   whole domain.

The `slab` and `circle2dpolar` geometries behave as in LaserSlab.

### 4. Freeing memory: `Simulation_freeUserArrays.F90`

The FLASH Driver calls this routine after the initial AMR grid is complete,
including all refinement passes, each of which calls `Simulation_initBlock` again.
It deallocates the axes and tables, so they use no memory during the
evolution.

### Why this design

In a naive implementation every block (and every refinement pass) re-opens and
re-reads the text files, and every rank does this in parallel. With tables of
~24 MB per field and thousands of blocks, initialisation then dominates the run
time and floods the file system. Reading once on one rank and broadcasting
makes the I/O cost independent of the block and rank counts.

## Notes

- `Config` lists the example tables as `DATAFILES` and uses
  `input/Collision_Laser_25fs_4P16_L_100nm_smooth` as the default `sim_inputFilenameRoot`.
  If you use your own tables, update both (or set `sim_inputFilenameRoot` in
  `flash.par`).
- Decompress the example tables (`gunzip -k input/*.txt.gz`) **before** running
  `./setup`, since the `DATAFILES` entries point to the `.txt` files.
