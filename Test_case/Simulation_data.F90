!!****if* source/Simulation/SimulationMain/LaserSlab/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!  Use Simulation_data
!!
!! DESCRIPTION
!!
!!  Store the simulation data
!!
!!
!!***
module Simulation_data

  implicit none

#include "constants.h"

  !! *** Runtime Parameters *** !!
  real, save :: sim_targetRadius

  real,    save :: sim_rhoTarg
  real,    save :: sim_teleTarg
  real,    save :: sim_tionTarg
  real,    save :: sim_tradTarg
  real,    save :: sim_zminTarg
  integer, save :: sim_eosTarg

  real,    save :: sim_rhoCham
  real,    save :: sim_teleCham
  real,    save :: sim_tionCham
  real,    save :: sim_tradCham
  integer, save :: sim_eosCham

  logical, save :: sim_killdivb = .FALSE.
  real, save :: sim_smallX
  character(len=MAX_STRING_LENGTH), save :: sim_initGeom

  character(len=MAX_STRING_LENGTH), save :: sim_inputFilenameRoot

  !! *** Tabulated initial state (sim_initGeom = "fromtxt") *** !!
  !!
  !! Read ONCE in Simulation_init (by the master rank, then broadcast),
  !! used by every call of Simulation_initBlock, and released in
  !! Simulation_freeUserArrays once the initial grid is built.
  !!
  !! Layout (unchanged from the original files):
  !!   sim_xFile(nx), sim_yFile(ny)       : strictly increasing coordinates
  !!   sim_xxxFile(nx, ny)                : row i of the text file = x index i
  !!
  !! Negative values are already clamped to 0 when loaded.
  !! Trad uses the tele table (there is no separate _trad file),
  !! exactly as in the original implementation.
  integer, save :: sim_nxFile = 0, sim_nyFile = 0
  real, allocatable, save :: sim_xFile(:), sim_yFile(:)
  real, allocatable, save :: sim_densFile(:,:), sim_teleFile(:,:), sim_tionFile(:,:)

  ! Values used for cells lying outside the tabulated region when
  ! sim_txtFloorOutside is .true. (see Simulation_initBlock):
  ! max(min(table), chamber value), as in the original code.
  real, save :: sim_densOutside, sim_teleOutside, sim_tionOutside

  ! Behaviour for cells whose centre lies outside the x/y range of the files:
  !   .false. : take the value at the nearest edge of the table (no jump)
  !   .true.  : use sim_densOutside / sim_teleOutside / sim_tionOutside
  logical, parameter :: sim_txtFloorOutside = .false.

end module Simulation_data