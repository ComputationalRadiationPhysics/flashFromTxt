!!****if* source/Simulation/SimulationMain/LaserSlab/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!!
!! DESCRIPTION
!!
!!  Initializes all the parameters needed for a particular simulation.
!!
!!  When sim_initGeom == "fromtxt", the tabulated initial state is read
!!  here ONCE: the master rank reads the text files and broadcasts the
!!  tables to all other ranks. Simulation_initBlock then only
!!  interpolates from memory, and Simulation_freeUserArrays releases
!!  the tables once the initial AMR grid has been built.
!!
!!
!! ARGUMENTS
!!
!!
!!
!! PARAMETERS
!!
!!***

subroutine Simulation_init()
  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype, Driver_getComm

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  real :: xmin, xmax, ymin, ymax
  integer :: lrefine_max, nblockx, nblocky
  character(len=MAX_STRING_LENGTH) :: str
  logical :: restart

  call RuntimeParameters_get('sim_targetRadius', sim_targetRadius)

  call RuntimeParameters_get('sim_rhoTarg', sim_rhoTarg)
  call RuntimeParameters_get('sim_teleTarg', sim_teleTarg)
  call RuntimeParameters_get('sim_tionTarg', sim_tionTarg)
  call RuntimeParameters_get('sim_tradTarg', sim_tradTarg)

  call RuntimeParameters_get('sim_rhoCham', sim_rhoCham)
  call RuntimeParameters_get('sim_teleCham', sim_teleCham)
  call RuntimeParameters_get('sim_tionCham', sim_tionCham)
  call RuntimeParameters_get('sim_tradCham', sim_tradCham)

  call RuntimeParameters_get('smallX', sim_smallX)

  call RuntimeParameters_get('sim_initGeom', sim_initGeom)

  call RuntimeParameters_get('sim_inputFilenameRoot', sim_inputFilenameRoot)

#ifdef FLASH_USM_MHD
  call RuntimeParameters_get('killdivb', sim_killdivb)
#endif

  ! Load the tables once. On a restart Simulation_initBlock is not
  ! called (the state comes from the checkpoint), so skip the read.
  call RuntimeParameters_get('restart', restart)
  if (trim(sim_initGeom) == "fromtxt" .and. .not. restart) then
     call sim_loadTxtTables()
  end if

contains

  !---------------------------------------------------------------------
  ! Master reads all files, then everything is broadcast to all ranks.
  !---------------------------------------------------------------------
  subroutine sim_loadTxtTables()
    integer :: myPE, globalComm, ierr
    integer :: dims(2)
    character(len=MAX_STRING_LENGTH) :: root

    if (NDIM /= 2) call Driver_abortFlash( &
         "Simulation_init: sim_initGeom='fromtxt' is only implemented for 2D")

    call Driver_getMype(GLOBAL_COMM, myPE)
    call Driver_getComm(GLOBAL_COMM, globalComm)

    root = sim_inputFilenameRoot

    if (myPE == MASTER_PE) then
       call sim_readVector(trim(root)//'_xfile.txt', sim_xFile)
       call sim_readVector(trim(root)//'_yfile.txt', sim_yFile)
       sim_nxFile = size(sim_xFile)
       sim_nyFile = size(sim_yFile)
       call sim_checkAxis(sim_xFile, trim(root)//'_xfile.txt')
       call sim_checkAxis(sim_yFile, trim(root)//'_yfile.txt')

       call sim_readMatrix(trim(root)//'_dens.txt', sim_densFile, sim_nxFile, sim_nyFile)
       call sim_readMatrix(trim(root)//'_tele.txt', sim_teleFile, sim_nxFile, sim_nyFile)
       call sim_readMatrix(trim(root)//'_tion.txt', sim_tionFile, sim_nxFile, sim_nyFile)
       dims = (/ sim_nxFile, sim_nyFile /)
    end if

    ! --- broadcast sizes, allocate on the other ranks, broadcast data ---
    call MPI_Bcast(dims, 2, FLASH_INTEGER, MASTER_PE, globalComm, ierr)
    sim_nxFile = dims(1)
    sim_nyFile = dims(2)

    if (myPE /= MASTER_PE) then
       allocate(sim_xFile(sim_nxFile), sim_yFile(sim_nyFile))
       allocate(sim_densFile(sim_nxFile, sim_nyFile))
       allocate(sim_teleFile(sim_nxFile, sim_nyFile))
       allocate(sim_tionFile(sim_nxFile, sim_nyFile))
    end if

    call MPI_Bcast(sim_xFile,    sim_nxFile,            FLASH_REAL, MASTER_PE, globalComm, ierr)
    call MPI_Bcast(sim_yFile,    sim_nyFile,            FLASH_REAL, MASTER_PE, globalComm, ierr)
    call MPI_Bcast(sim_densFile, sim_nxFile*sim_nyFile, FLASH_REAL, MASTER_PE, globalComm, ierr)
    call MPI_Bcast(sim_teleFile, sim_nxFile*sim_nyFile, FLASH_REAL, MASTER_PE, globalComm, ierr)
    call MPI_Bcast(sim_tionFile, sim_nxFile*sim_nyFile, FLASH_REAL, MASTER_PE, globalComm, ierr)

    ! --- one-time preprocessing (previously redone for every block and field) ---
    sim_densFile = max(sim_densFile, 0.0)
    sim_teleFile = max(sim_teleFile, 0.0)
    sim_tionFile = max(sim_tionFile, 0.0)

    sim_densOutside = max(minval(sim_densFile), sim_rhoCham)
    sim_teleOutside = max(minval(sim_teleFile), sim_teleCham)
    sim_tionOutside = max(minval(sim_tionFile), sim_tionCham)

    if (myPE == MASTER_PE) then
       write(str,'(a,i0,a,i0,a)') "loaded ", sim_nxFile, " x ", sim_nyFile, &
            " tables from "//trim(root)
       call Logfile_stamp(trim(str), "[Simulation_init]")
       write(*,'(a)') " [Simulation_init] "//trim(str)
    end if
  end subroutine sim_loadTxtTables


  !---------------------------------------------------------------------
  ! One value per line.
  !---------------------------------------------------------------------
  subroutine sim_readVector(filename, vec)
    character(len=*), intent(in)    :: filename
    real, allocatable, intent(out)  :: vec(:)
    integer :: u, ios, n, i
    real    :: tmp

    open(newunit=u, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) call Driver_abortFlash("Simulation_init: cannot open "//filename)

    n = 0
    do
       read(u, *, iostat=ios) tmp
       if (ios /= 0) exit
       n = n + 1
    end do
    if (n < 2) call Driver_abortFlash("Simulation_init: need >= 2 points in "//filename)

    rewind(u)
    allocate(vec(n))
    do i = 1, n
       read(u, *, iostat=ios) vec(i)
       if (ios /= 0) call Driver_abortFlash("Simulation_init: read error in "//filename)
    end do
    close(u)
  end subroutine sim_readVector


  !---------------------------------------------------------------------
  ! nrows lines, each with ncols values (row i <-> x index i).
  !---------------------------------------------------------------------
  subroutine sim_readMatrix(filename, mat, nrows, ncols)
    character(len=*), intent(in)    :: filename
    integer, intent(in)             :: nrows, ncols
    real, allocatable, intent(out)  :: mat(:,:)
    integer :: u, ios, i, j
    character(len=32) :: istr

    open(newunit=u, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) call Driver_abortFlash("Simulation_init: cannot open "//filename)

    allocate(mat(nrows, ncols))
    do i = 1, nrows
       read(u, *, iostat=ios) (mat(i, j), j = 1, ncols)
       if (ios /= 0) then
          write(istr,'(i0)') i
          call Driver_abortFlash("Simulation_init: failed reading row "//trim(istr)// &
               " of "//filename)
       end if
    end do
    close(u)
  end subroutine sim_readMatrix


  !---------------------------------------------------------------------
  ! The interpolation uses a binary search, so axes must be increasing.
  !---------------------------------------------------------------------
  subroutine sim_checkAxis(ax, filename)
    real, intent(in)             :: ax(:)
    character(len=*), intent(in) :: filename
    integer :: i
    do i = 2, size(ax)
       if (.not. (ax(i) > ax(i-1))) call Driver_abortFlash( &
            "Simulation_init: coordinates must be strictly increasing in "//filename)
    end do
  end subroutine sim_checkAxis

end subroutine Simulation_init