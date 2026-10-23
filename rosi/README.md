# Building and running FLASH 4.8 on rosi

Files for compiling and running FLASH 4.8 on the HZDR **rosi** cluster
(`cpu-rome` partition) with GCC 7.3.0 and OpenMPI 4.1.8:

| File | Where it goes | Role |
|---|---|---|
| `FLASH_rome.sh` | your home directory, as `~/FLASH_rome.sh` (**not** in the FLASH tree) | Shell profile: loads the modules and sets the HDF5 / HYPRE paths |
| `Makefile.h` | `<FLASH root>/sites/Prototypes/Linux/Makefile.h` | FLASH compiler and library definitions used by `./setup` |
| `install_hdf5_hypre_mpif90.sh` | run once from the repository (installs into `$HOME/lib/...`) | Builds the parallel HDF5 and HYPRE libraries FLASH needs |

## Prerequisites

FLASH needs **parallel HDF5** and **HYPRE** built with the same compiler and MPI.
They are not provided as modules, so both files expect your own builds in:

```
$HOME/lib/gcc7.3.0-openmpi4.1.8/hdf5-1.14.6-parallel-mpif90
$HOME/lib/gcc7.3.0-openmpi4.1.8/hypre-2.30.0-mpif90
```

Build both with the `mpif90` / `mpicc` of the `openmpi` module loaded by
`FLASH_rome.sh`. `install_hdf5_hypre_mpif90.sh` downloads, builds and installs
them into exactly these directories:

```bash
source ~/FLASH_rome.sh                 # same compiler and MPI as FLASH
bash rosi/install_hdf5_hypre_mpif90.sh
```

It also writes `$HOME/lib/gcc7.3.0-openmpi4.1.8/env-hdf5-hypre-mpif90.sh` and
compiles two small test programs. Set `MAKE_JOBS=8` (or similar) on the login
node, and `RUN_CHECKS=1` only on a compute node. The other settings
(`HDF5_VERSION`, `HYPRE_VERSION`, `STACK`, `BUILD_SHARED`) are listed at the
top of the script.

If you install them elsewhere, change the paths in **both**
files (`HDF5_PREFIX` / `HYPRE_PREFIX` in `FLASH_rome.sh`, `HDF5_PATH` /
`HYPRE_PATH` in `Makefile.h`).

## 1. Shell profile: `FLASH_rome.sh`

```bash
cp rosi/FLASH_rome.sh ~/FLASH_rome.sh
source ~/FLASH_rome.sh
```

Source it **before running `./setup`/`make`** and **in every job script**, so
that the code runs with the same libraries it was compiled against
(all the `Test_case/start_rosi*` and `submit_slurm_*.sh` scripts already do this).

It:

- purges the modules and loads `rome`, `gcc/7.3.0`, `openmpi`, `zlib`, `lapack`,
  `python` and `gnuplot`;
- adds the HDF5 and HYPRE builds to `PATH`, `LD_LIBRARY_PATH`,
  `CMAKE_PREFIX_PATH` and `PKG_CONFIG_PATH`;
- defines helpers for interactive sessions on `cpu-rome`:

  | Command | Allocation |
  |---|---|
  | `getNode [n]` | *n* full nodes (default 1) for 1 h |
  | `getDevice [n]` | *n* ≤ 2 tasks × 20 cores for 1 h |
  | `get3HoursDevice [n]`, `get24HoursDevice [n]` | same, for 3 h / 24 h |
  | `sq` | `squeue` for your own jobs (`$USER`) |

## 2. `Makefile.h`

```bash
cp rosi/Makefile.h $FLASH_DIR/sites/Prototypes/Linux/Makefile.h
```

On Linux, `./setup` picks this prototype automatically when no site-specific
directory matches the host name. You can also select it explicitly with
`-site=Prototypes/Linux`.

The file sets:

- `MPI_PATH` to the OpenMPI 4.1.8 / GCC 7.3.0 install in
  `/data/rosi/shared/pkg/...`. It must be the same OpenMPI as the `openmpi`
  module, and the one HDF5 and HYPRE were built with.
- `HDF5_PATH`, `HYPRE_PATH` to the builds above, and `LAPACK_PATH` to the
  rosi LAPACK 3.8.0 module.
- the `mpif90` / `mpicc` / `mpicxx` compiler wrappers, `-fdefault-real-8
  -fdefault-double-8`, and the `_OPT` (`-O2`), `_DEBUG` (`-O0`, bounds and
  FPE checks) and `_TEST` flag sets;
- `-DH5_USE_16_API` for HDF5 1.14 compatibility with FLASH's HDF5 I/O.

## Checking the build

`Test_case/submit_slurm_multinode.sh` prints the libraries the executable is
linked against:

```bash
ldd ./flash4 | grep -Ei 'libmpi|hdf5|HYPRE'
```

All three should come from the OpenMPI / HDF5 / HYPRE described above. If
they don't (e.g. a system HDF5 or another OpenMPI is picked up), check that
`~/FLASH_rome.sh` was sourced both for the build and in the job script.
