!!****f* source/Simulation/Simulation_adjustEvolution
!!
!! NAME
!!  Simulation_adjustEvolution
!!
!!
!! SYNOPSIS
!!  Simulation_adjustEvolution( integer(IN) :: blkcnt,
!!                              integer(IN) :: blklst(blkcnt),
!!                              integer(IN) :: nstep,
!!                              real(IN) :: dt,
!!                              real(IN) :: stime )
!!
!! DESCRIPTION
!!  This routine is called every cycle. It can be used to adjust
!!  the simulation while it is running.
!!  
!! ARGUMENTS
!!  blkcnt - number of blocks
!!  blklist - block list
!!  nstep - current cycle number
!!  dt - current time step length
!!  stime - current simulation time
!!
!!***
subroutine Simulation_adjustEvolution(blkcnt, blklst, nstep, dt, stime)

#include "constants.h"
#include "Flash.h"

  use Simulation_data

  use Grid_interface, ONLY: Grid_getBlkIndexLimits
  use Grid_interface, ONLY: Grid_getBlkPtr
  use Grid_interface, ONLY: Grid_releaseBlkPtr
  use Grid_interface, ONLY: Grid_getBlkBC

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
       Grid_getCellCoords, Grid_putPointData,&
                           Grid_getBlkPtr,         &
                           Grid_releaseBlkPtr, Grid_getPointData

  implicit none


  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)

  integer :: i, j, k, n
  integer :: blkLimits(2, MDIM)
  integer :: blkLimitsGC(2, MDIM)
  integer :: axis(MDIM)
  integer :: lb
  real, allocatable :: xcent(:), ycent(:), zcent(:)
  real :: tele
  real, pointer, dimension(:,:,:,:) :: facexData,faceyData
#if NDIM > 0
  real, pointer, dimension(:,:,:,:) :: facezData
#endif
  real, pointer :: blkPtr(:,:,:,:)
  integer, intent(in) :: blkcnt
  integer, intent(in) :: blklst(blkcnt)
  integer, intent(in) :: nstep
  real, intent(in) :: dt
  real, intent(in) :: stime



  !------------------------------------------------------------------------------


  return

end subroutine Simulation_adjustEvolution




