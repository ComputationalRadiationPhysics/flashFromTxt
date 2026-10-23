#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Local install script: parallel HDF5 + hypre built with MPI wrappers (mpicc/mpif90)
#
# Usage (on rosi, load the same toolchain as FLASH first):
#   source ~/FLASH_rome.sh
#   bash rosi/install_hdf5_hypre_mpif90.sh
#
# Optional overrides (env vars):
#   HDF5_VERSION=1.14.6
#   HYPRE_VERSION=2.30.0
#   STACK=$HOME/lib/gcc7.3.0-openmpi4.1.8
#   MAKE_JOBS=8
#   RUN_CHECKS=1          # run "make check" for HDF5 (often needs compute node)
#   BUILD_SHARED=1        # 1 => build shared+static; 0 => static only
###############################################################################

# ---- user-tunable defaults ----
HDF5_VERSION="${HDF5_VERSION:-1.14.6}"
HYPRE_VERSION="${HYPRE_VERSION:-2.30.0}"
STACK="${STACK:-$HOME/lib/gcc7.3.0-openmpi4.1.8}"
SRC="${SRC:-$STACK/src}"
MAKE_JOBS="${MAKE_JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}"
RUN_CHECKS="${RUN_CHECKS:-0}"
BUILD_SHARED="${BUILD_SHARED:-1}"

# If your cluster uses Environment Modules, you can uncomment these:
# module purge
# module load gcc/7.3.0 openmpi/4.1.8

# ---- derived install prefixes ----
HDF5_PREFIX="${HDF5_PREFIX:-$STACK/hdf5-${HDF5_VERSION}-parallel-mpif90}"
HYPRE_PREFIX="${HYPRE_PREFIX:-$STACK/hypre-${HYPRE_VERSION}-mpif90}"

# ---- helpers ----
die(){ echo "ERROR: $*" >&2; exit 1; }
need_cmd(){ command -v "$1" >/dev/null 2>&1 || die "Missing required command: $1"; }

fetch() {
  local url="$1"
  local out="$2"
  if command -v wget >/dev/null 2>&1; then
    wget -O "$out" "$url"
  elif command -v curl >/dev/null 2>&1; then
    curl -L -o "$out" "$url"
  else
    die "Need wget or curl to download sources."
  fi
}

banner(){ echo; echo "=== $* ==="; }

# ---- sanity checks: MPI wrappers must exist ----
need_cmd mpicc
need_cmd mpif90

# ---- create dirs ----
mkdir -p "$SRC"

###############################################################################
# 1) Build & install parallel HDF5 with Fortran interface (h5pfc uses mpif90)
###############################################################################
banner "Building HDF5 ${HDF5_VERSION} (parallel + Fortran) -> $HDF5_PREFIX"

cd "$SRC"
HDF5_TGZ="hdf5-${HDF5_VERSION}.tar.gz"
HDF5_URL="https://github.com/HDFGroup/hdf5/releases/download/hdf5_${HDF5_VERSION}/hdf5-${HDF5_VERSION}.tar.gz"

if [[ ! -f "$HDF5_TGZ" ]]; then
  fetch "$HDF5_URL" "$HDF5_TGZ"
fi

rm -rf "hdf5-${HDF5_VERSION}"
tar -xzf "$HDF5_TGZ"
cd "hdf5-${HDF5_VERSION}"

rm -rf build && mkdir build && cd build

HDF5_SHARED_FLAGS=()
if [[ "$BUILD_SHARED" == "1" ]]; then
  HDF5_SHARED_FLAGS+=(--enable-shared --enable-static)
else
  HDF5_SHARED_FLAGS+=(--disable-shared --enable-static)
fi

# Use mpif90 explicitly as requested
CC=mpicc CXX=mpicxx FC=mpif90 F77=mpif90 \
../configure \
  --prefix="$HDF5_PREFIX" \
  --enable-parallel \
  --enable-fortran \
  --enable-hl \
  --with-pic \
  "${HDF5_SHARED_FLAGS[@]}"

make -j "$MAKE_JOBS"

if [[ "$RUN_CHECKS" == "1" ]]; then
  # Note: may need to run on a compute node; leave off if your login node disallows mpirun.
  make check || die "HDF5 make check failed (try RUN_CHECKS=0 if your system blocks MPI runs on login nodes)."
fi

make install

# quick checks
[[ -x "$HDF5_PREFIX/bin/h5pfc" ]] || die "Expected $HDF5_PREFIX/bin/h5pfc not found (Fortran wrapper)."
echo "HDF5 h5pfc uses:"
"$HDF5_PREFIX/bin/h5pfc" -show | head -n 2 || true

###############################################################################
# 2) Build & install hypre with MPI + Fortran support using mpif90
###############################################################################
banner "Building hypre v${HYPRE_VERSION} (MPI + Fortran) -> $HYPRE_PREFIX"

cd "$SRC"
HYPRE_TGZ="hypre-${HYPRE_VERSION}.tar.gz"
HYPRE_URL="https://github.com/hypre-space/hypre/archive/refs/tags/v${HYPRE_VERSION}.tar.gz"

if [[ ! -f "$HYPRE_TGZ" ]]; then
  fetch "$HYPRE_URL" "$HYPRE_TGZ"
fi

rm -rf "hypre-${HYPRE_VERSION}"
tar -xzf "$HYPRE_TGZ"
cd "hypre-${HYPRE_VERSION}/src"

HYPRE_SHARED_FLAGS=()
if [[ "$BUILD_SHARED" == "1" ]]; then
  HYPRE_SHARED_FLAGS+=(--enable-shared)
fi

CC=mpicc CXX=mpicxx FC=mpif90 F77=mpif90 \
./configure \
  --prefix="$HYPRE_PREFIX" \
  --with-MPI \
  --enable-fortran \
  ${HYPRE_SHARED_FLAGS[@]+"${HYPRE_SHARED_FLAGS[@]}"}

make -j "$MAKE_JOBS"
make install

###############################################################################
# 3) Write an env script to activate this stack
###############################################################################
banner "Writing environment setup -> $STACK/env-hdf5-hypre-mpif90.sh"

cat > "$STACK/env-hdf5-hypre-mpif90.sh" <<EOF
# Source this file to use the locally installed HDF5 + hypre stack (mpif90-based)
#   source "$STACK/env-hdf5-hypre-mpif90.sh"

# If you normally use modules, load the same toolchain you built with:
# module purge
# module load gcc/7.3.0 openmpi/4.1.8

export HDF5_PREFIX="$HDF5_PREFIX"
export HYPRE_PREFIX="$HYPRE_PREFIX"

# Common names used by build systems:
export HDF5_DIR="\$HDF5_PREFIX"
export HDF5_PATH="\$HDF5_PREFIX"
export HYPRE_DIR="\$HYPRE_PREFIX"

export PATH="\$HDF5_PREFIX/bin:\$PATH"
export LD_LIBRARY_PATH="\$HDF5_PREFIX/lib:\$HYPRE_PREFIX/lib:\$LD_LIBRARY_PATH"
export CMAKE_PREFIX_PATH="\$HDF5_PREFIX:\$HYPRE_PREFIX:\$CMAKE_PREFIX_PATH"
export PKG_CONFIG_PATH="\$HDF5_PREFIX/lib/pkgconfig:\$HYPRE_PREFIX/lib/pkgconfig:\$PKG_CONFIG_PATH"
EOF

###############################################################################
# 4) Smoke tests (compile only)
###############################################################################
banner "Smoke tests (compile/link)"

# HDF5 Fortran test program (compile with h5pfc which invokes mpif90)
cat > "$STACK/test_hdf5_mpi.f90" <<'EOF'
program test_hdf5_mpi
  use mpi
  use hdf5
  implicit none
  integer :: ierr, rank, nprocs
  integer(HID_T) :: file_id

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

  call h5open_f(ierr)
  call h5fcreate_f("test_parallel.h5", H5F_ACC_TRUNC_F, file_id, ierr)
  call h5fclose_f(file_id, ierr)
  call h5close_f(ierr)

  if (rank == 0) print *, "OK: created test_parallel.h5 with ", nprocs, " ranks"
  call MPI_Finalize(ierr)
end program
EOF

"$HDF5_PREFIX/bin/h5pfc" "$STACK/test_hdf5_mpi.f90" -o "$STACK/test_hdf5_mpi.exe"

# hypre link test (no calls; just ensure mpif90 can link libHYPRE)
cat > "$STACK/test_hypre_link.f90" <<'EOF'
program test_hypre_link
  implicit none
  print *, "OK: hypre link test"
end program
EOF

mpif90 "$STACK/test_hypre_link.f90" -o "$STACK/test_hypre_link.exe" \
  -L"$HYPRE_PREFIX/lib" -lHYPRE -lm

banner "Done"
echo "HDF5 : $HDF5_PREFIX"
echo "hypre: $HYPRE_PREFIX"
echo "Environment: source $STACK/env-hdf5-hypre-mpif90.sh"

