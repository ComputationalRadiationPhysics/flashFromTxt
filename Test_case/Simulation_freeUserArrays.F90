!!****if* source/Simulation/SimulationMain/LaserSlab/Simulation_freeUserArrays
!!
!! NAME
!!
!!  Simulation_freeUserArrays
!!
!! SYNOPSIS
!!
!!  call Simulation_freeUserArrays()
!!
!! DESCRIPTION
!!
!!  Called by the Driver after all initialization steps (including the
!!  construction of the initial AMR grid) are done. Releases the
!!  tabulated initial state read in Simulation_init, which is no
!!  longer needed during the evolution.
!!
!!***

subroutine Simulation_freeUserArrays()
  use Simulation_data, ONLY : sim_xFile, sim_yFile, &
       sim_densFile, sim_teleFile, sim_tionFile
  implicit none

  if (allocated(sim_xFile))    deallocate(sim_xFile)
  if (allocated(sim_yFile))    deallocate(sim_yFile)
  if (allocated(sim_densFile)) deallocate(sim_densFile)
  if (allocated(sim_teleFile)) deallocate(sim_teleFile)
  if (allocated(sim_tionFile)) deallocate(sim_tionFile)

end subroutine Simulation_freeUserArrays
