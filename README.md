# flashFromTxt

Initialise a 2D [FLASH](https://flash.rochester.edu/) radiation-hydrodynamics
simulation from the output of a PIConGPU particle-in-cell simulation.

The workflow has two parts:

```
 PIConGPU openPMD output            openpmd2txt.py                 FLASH (Test_case)
 (H / e density, energy density) ──► text tables (dens, tele, ──► reads the tables once, broadcasts
                                     tion, x and y axes)          them and interpolates them onto
                                                                   every AMR block at initialisation
```

The opposite direction, from FLASH back to PIC, is covered by `flash2openpmd`:

```
 FLASH checkpoint / plot file       flash2openpmd.py               PIConGPU (or other PIC code)
 (dens, tele, tion, trad, ye)   ──► openPMD HDF5 meshes on a   ──► initial plasma profile
                                    uniform, optionally rotated
                                    2D grid (SI units)
```

Together they allow a PIC → FLASH → PIC chain, e.g. a short-pulse
interaction in PIConGPU, its hydrodynamic expansion in FLASH, and a later
PIC simulation on the expanded profile.

| Folder | Content |
|---|---|
| [`openpmd2txt/`](openpmd2txt/) | Python script that reads an openPMD series, converts it to FLASH units and writes the text tables. See [openpmd2txt/README.md](openpmd2txt/README.md). |
| [`Test_case/`](Test_case/) | A FLASH 4 simulation directory (derived from `LaserSlab`) with a new `sim_initGeom = "fromtxt"` mode that reads and interpolates the tables. See [Test_case/README.md](Test_case/README.md). |
| [`flash2openpmd/`](flash2openpmd/) | Python script that converts a 2D FLASH checkpoint into an openPMD file (dens, tele, tion, trad, ye) for initialising a PIC simulation, with optional density cutoff, boundary layer and rotation. See [flash2openpmd/README.md](flash2openpmd/README.md), including its known limitations. |
| [`rosi/`](rosi/) | Shell profile (`FLASH_rome.sh`, goes in `~`), `Makefile.h` (goes in `sites/Prototypes/Linux/`) and an HDF5 + HYPRE install script for building and running FLASH 4.8 on the rosi cluster. See [rosi/README.md](rosi/README.md). |

## Quick start

1. **Make the tables** from a PIConGPU run:

   ```bash
   cd openpmd2txt
   pip install -r requirements.txt
   python openpmd2txt.py -input_file /path/to/simData_%T.bp -output_name MyRun \
          -nx_flash 2400 -ny_flash 800 -smooth -debug
   # -> output/MyRun_{xfile,yfile,dens,tele,tion}.txt
   ```

2. **Install the test case** in your FLASH 4 tree, copy the tables into
   `Test_case/input/`, and set `sim_inputFilenameRoot = "input/MyRun"` in `flash.par`.

3. **Build and run** FLASH as usual (see `Test_case/README.md` for the setup line).
   On rosi, first install the profile and `Makefile.h` from [`rosi/`](rosi/README.md).

4. **Back to PIC** (optional): convert a FLASH checkpoint to openPMD:

   ```bash
   cd flash2openpmd
   pip install -r requirements.txt
   python flash2openpmd.py -input_file /path/to/run/lasslab_hdf5_chk_0010 \
          -output_file profile_t10 -nx_pic 4800 -ny_pic 4800
   # -> output/profile_t10.h5
   ```

There is an example set of tables in `Test_case/input/*.txt.gz`. It is gzipped
to keep the repository small, so decompress it before running:

```bash
gunzip -k Test_case/input/*.txt.gz
```

## Licence note

The files in `Test_case/` are modified FLASH source files. FLASH is distributed
under its own user licence, so make sure anyone you share this repository
with is a licensed FLASH user.
