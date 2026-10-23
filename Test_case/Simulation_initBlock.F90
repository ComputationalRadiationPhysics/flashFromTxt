!!****if* source/Simulation/SimulationMain/LaserSlab/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_initBlock(integer(IN) :: blockID)
!!
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.
!!
!!  For sim_initGeom == "fromtxt" the tables are NOT read here: they
!!  were loaded once in Simulation_init and live in Simulation_data.
!!  This routine only does a bilinear interpolation from memory:
!!    - the bracketing table index and weight are found once per
!!      i column and once per j row of the block (binary search),
!!    - then every cell is a 4-point weighted sum.
!!
!! ARGUMENTS
!!
!!  blockID -        the number of the block to initialize
!!
!!
!!
!!***

subroutine Simulation_initBlock(blockId)
  use Simulation_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
       Grid_getCellCoords, Grid_putPointData, &
       Grid_getBlkPtr, Grid_releaseBlkPtr

  use Driver_interface, ONLY: Driver_abortFlash
  use RadTrans_interface, ONLY: RadTrans_mgdEFromT
  implicit none

#include "constants.h"
#include "Flash.h"

  ! Inputs
  integer, intent(in) :: blockId

  ! Local variables
  integer :: i, j, k, n
  integer :: blkLimits(2, MDIM)
  integer :: blkLimitsGC(2, MDIM)
  integer :: axis(MDIM)
  real, allocatable :: xcent(:), ycent(:)
  real :: tradActual
  real :: rho, tele, trad, tion, bound
  integer :: species
  real, pointer, dimension(:,:,:,:) :: facexData, faceyData
#if NDIM > 0
  real, pointer, dimension(:,:,:,:) :: facezData
#endif

  ! Per-column / per-row interpolation stencils (fromtxt only)
  logical :: fromTxt
  integer, allocatable :: ix(:), iy(:)
  real,    allocatable :: tx(:), ty(:)
  logical, allocatable :: inX(:), inY(:)
  real :: w00, w10, w01, w11

#ifndef CHAM_SPEC
  integer :: CHAM_SPEC = 1, TARG_SPEC = 2
#endif

  call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC)

  ! Coordinate and grid initialization
  allocate(xcent(blkLimitsGC(HIGH, IAXIS)))
  call Grid_getCellCoords(IAXIS, blockId, CENTER, .true., xcent, blkLimitsGC(HIGH, IAXIS))
  allocate(ycent(blkLimitsGC(HIGH, JAXIS)))
  call Grid_getCellCoords(JAXIS, blockId, CENTER, .true., ycent, blkLimitsGC(HIGH, JAXIS))

#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_getBlkPtr(blockID, facexData, FACEX)
     call Grid_getBlkPtr(blockID, faceyData, FACEY)
     if (NDIM > 2) call Grid_getBlkPtr(blockID, facezData, FACEZ)
  endif
#endif

  fromTxt = (trim(sim_initGeom) == "fromtxt")

  if (fromTxt) then
     if (.not. allocated(sim_densFile)) call Driver_abortFlash( &
          "Simulation_initBlock: fromtxt tables not loaded (Simulation_init)")

     ! Locate each column / row once in the tables: O((nx+ny) log N)
     allocate(ix(size(xcent)), tx(size(xcent)), inX(size(xcent)))
     allocate(iy(size(ycent)), ty(size(ycent)), inY(size(ycent)))
     do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
        call locate(sim_xFile, xcent(i), ix(i), tx(i), inX(i))
     end do
     do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
        call locate(sim_yFile, ycent(j), iy(j), ty(j), inY(j))
     end do
  end if

  ! Loop over cells and set the initial state
  do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
     do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
        do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)

           axis(IAXIS) = i
           axis(JAXIS) = j
           axis(KAXIS) = k
           bound = 1.0

           if (fromTxt) then
              species = TARG_SPEC

              if (sim_txtFloorOutside .and. .not. (inX(i) .and. inY(j))) then
                 rho  = sim_densOutside
                 tele = sim_teleOutside
                 tion = sim_tionOutside
              else
                 w00 = (1.0 - tx(i)) * (1.0 - ty(j))
                 w10 =        tx(i)  * (1.0 - ty(j))
                 w01 = (1.0 - tx(i)) *        ty(j)
                 w11 =        tx(i)  *        ty(j)
                 rho  = bilin(sim_densFile, ix(i), iy(j), w00, w10, w01, w11)
                 tele = bilin(sim_teleFile, ix(i), iy(j), w00, w10, w01, w11)
                 tion = bilin(sim_tionFile, ix(i), iy(j), w00, w10, w01, w11)
              end if
              trad = tele          ! no separate _trad file (as before)
              bound = -1.0
           else
              species = CHAM_SPEC
              if (sim_initGeom == "slab") then
                 if (xcent(i) <= sim_targetRadius .and. &
                     xcent(i) >= -1.*sim_targetRadius) then
                    species = TARG_SPEC
                 end if
              elseif(sim_initGeom == "circle2dpolar" .and. NDIM == 2) then
                 if (sqrt(xcent(i)**2 + ycent(j)**2) <= sim_targetRadius) then
                    species = TARG_SPEC
                 end if
              else
                 call Driver_abortFlash("Simulation_initBlock: sim_initGeom unknown!")
              end if
              if (species == TARG_SPEC) then
                 rho = sim_rhoTarg
                 tele = sim_teleTarg
                 tion = sim_tionTarg
                 trad = sim_tradTarg
                 bound=-1.0
              else
                 rho = sim_rhoCham
                 tele = sim_teleCham
                 tion = sim_tionCham
                 trad = sim_tradCham
                 bound=-1.0
              end if
           end if

           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho)
           call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, tele)

#ifdef FLASH_3T
           call Grid_putPointData(blockId, CENTER, TION_VAR, EXTERIOR, axis, tion)
           call Grid_putPointData(blockId, CENTER, TELE_VAR, EXTERIOR, axis, tele)

           ! Set up radiation energy density:
           call RadTrans_mgdEFromT(blockId, axis, trad, tradActual)
           call Grid_putPointData(blockId, CENTER, TRAD_VAR, EXTERIOR, axis, tradActual)
#endif

           if (NSPECIES > 0) then
              ! Fill mass fractions in solution array if we have any SPECIES defined
              do n = SPECIES_BEGIN, SPECIES_END
                 if (n == species) then
                    call Grid_putPointData(blockID, CENTER, n, EXTERIOR, axis, 1.0e0-(NSPECIES-1)*sim_smallX)
                 else
                    call Grid_putPointData(blockID, CENTER, n, EXTERIOR, axis, sim_smallX)
                 end if
              end do
           end if

#if NFACE_VARS > 0
           if (sim_killdivb) then
              if (NDIM == 2) then
                 facexData(MAG_FACE_VAR, i, j, k) = 0.0
                 faceyData(MAG_FACE_VAR, i, j, k) = 0.0
                 if (NDIM > 2) facezData(MAG_FACE_VAR, i, j, k) = 0.0
              end if
           end if
#endif
        end do
     end do
  end do

#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_releaseBlkPtr(blockID, facexData, FACEX)
     call Grid_releaseBlkPtr(blockID, faceyData, FACEY)
     if (NDIM > 2) call Grid_releaseBlkPtr(blockID, facezData, FACEZ)
  endif
#endif

  deallocate(xcent)
  deallocate(ycent)
  if (fromTxt) deallocate(ix, iy, tx, ty, inX, inY)

  return

contains

  !---------------------------------------------------------------------
  ! Binary search in a strictly increasing axis xs(1:n).
  ! Returns idx in [1, n-1] and weight t in [0, 1] such that
  !   x ~ (1-t)*xs(idx) + t*xs(idx+1).
  ! Points outside [xs(1), xs(n)] are clamped to the nearest edge
  ! (t = 0 or 1) and flagged with inside = .false.
  !---------------------------------------------------------------------
  pure subroutine locate(xs, x, idx, t, inside)
    real,    intent(in)  :: xs(:), x
    integer, intent(out) :: idx
    real,    intent(out) :: t
    logical, intent(out) :: inside
    integer :: lo, hi, mid, nx

    nx = size(xs)
    if (x <= xs(1)) then
       idx = 1;      t = 0.0;  inside = (x == xs(1))
       return
    else if (x >= xs(nx)) then
       idx = nx - 1; t = 1.0;  inside = (x == xs(nx))
       return
    end if

    lo = 1
    hi = nx
    do while (hi - lo > 1)
       mid = (lo + hi) / 2
       if (xs(mid) <= x) then
          lo = mid
       else
          hi = mid
       end if
    end do
    idx    = lo
    t      = (x - xs(lo)) / (xs(lo+1) - xs(lo))
    inside = .true.
  end subroutine locate

  pure real function bilin(f, i0, j0, a00, a10, a01, a11)
    real,    intent(in) :: f(:,:)
    integer, intent(in) :: i0, j0
    real,    intent(in) :: a00, a10, a01, a11
    bilin = a00 * f(i0,   j0  ) + a10 * f(i0+1, j0  ) + &
            a01 * f(i0,   j0+1) + a11 * f(i0+1, j0+1)
  end function bilin

end subroutine Simulation_initBlock