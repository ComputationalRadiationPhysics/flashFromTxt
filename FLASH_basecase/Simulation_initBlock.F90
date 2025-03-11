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
!! ARGUMENTS
!!
!!  blockID -        the number of the block to initialize
!!  
!!
!!
!!***

module readTxtData
  use Simulation_data
  implicit none
  private
  public :: read_txt_grid, read_txt_mesh, interpolate_txt_to_grid, init_geom_from_txt, cubic_interpolate_clamped

contains

subroutine read_txt_grid(filename, grid_data)
  implicit none

  ! Inputs
  character(len=*), intent(in) :: filename

  ! Outputs
  real, allocatable, intent(out) :: grid_data(:)

  ! Local variables
  integer :: i, io_stat
  logical :: file_exists
  real :: temp
  integer :: count

  ! Check if the file exists
  inquire(file=filename, exist=file_exists)
  if (.not. file_exists) then
    print *, "Error: File not found:", trim(filename)
    stop
  end if

  ! Open the file and count lines
  open(unit=10, file=filename, status='old', action='read', iostat=io_stat)
  if (io_stat /= 0) then
    print *, "Error: Unable to open file:", trim(filename)
    stop
  end if

  count = 0
  do
    read(10, *, iostat=io_stat)
    if (io_stat /= 0) exit
    count = count + 1
  end do

  rewind(10)
  allocate(grid_data(count))

  do i = 1, count
    read(10, *, iostat=io_stat) temp
    if (io_stat /= 0) then
      print *, "Error: Issue reading data from file:", trim(filename)
      stop
    end if
    grid_data(i) = temp
  end do

  close(10)
end subroutine read_txt_grid




subroutine read_txt_mesh(filename, mesh_data, nrows, ncols)
  implicit none

  ! Inputs
  character(len=*), intent(in) :: filename
  integer, intent(in) :: nrows, ncols

  ! Outputs
  real, allocatable, intent(out) :: mesh_data(:,:)

  ! Local variables
  integer :: i, j, io_stat
  logical :: file_exists

  ! Check if the file exists
  inquire(file=filename, exist=file_exists)
  if (.not. file_exists) then
    print *, "Error: File not found:", trim(filename)
    stop
  end if

  !print *, "Debug: Opening file:", trim(filename)
  open(unit=10, file=filename, status='old', action='read', iostat=io_stat)
  if (io_stat /= 0) then
    print *, "Error: Unable to open file:", trim(filename)
    stop
  end if

  ! Allocate the array
  allocate(mesh_data(nrows, ncols))

  ! Read the data into the array
  !print *, "Debug: Reading mesh data with dimensions:", nrows, "x", ncols
  do i = 1, nrows
    read(10, *, iostat=io_stat) (mesh_data(i, j), j = 1, ncols)
    if (io_stat /= 0) then
      print *, "Error: Failed to read row", i, "from file:", trim(filename)
      stop
    end if
  end do

  close(10)
  !print *, "Debug: Successfully read mesh data."
end subroutine read_txt_mesh




subroutine cubic_interpolate_clamped(x_grid, y_grid, data, x_file, y_file, file_data, default_value)
  implicit none

  ! Inputs
  real, intent(in) :: x_grid(:), y_grid(:)
  real, intent(in) :: x_file(:), y_file(:)
  real, intent(in) :: file_data(:,:)
  real, intent(in) :: default_value

  ! Outputs
  real, intent(out) :: data(size(x_grid), size(y_grid))

  ! Local variables
  integer :: i, j, i_low, j_low
  real :: dx, dy, tx, ty
  real :: f00, f01, f10, f11, f20, f21, f02, f12, f22
  real :: interp_value
  real :: min_value, max_value

  ! Clamp file_data to avoid negative values and get min/max bounds
  real, allocatable :: clamped_file_data(:,:)
  allocate(clamped_file_data(size(file_data, 1), size(file_data, 2)))
  clamped_file_data = max(file_data, 0.0)  ! Replace negative values with zero
  min_value = minval(clamped_file_data)   ! Minimum value in the clamped data
  max_value = maxval(clamped_file_data)   ! Maximum value in the clamped data

  ! Initialize output data to default_value
  data = default_value

  ! Loop over each grid point
  do j = 1, size(y_grid)
     do i = 1, size(x_grid)

        ! Find the cell in the input grid
        i_low = 1
        do while (i_low < size(x_file) .and. x_file(i_low) < x_grid(i))
           i_low = i_low + 1
        end do
        i_low = max(2, min(i_low, size(x_file) - 2))  ! Ensure valid range for bicubic

        j_low = 1
        do while (j_low < size(y_file) .and. y_file(j_low) < y_grid(j))
           j_low = j_low + 1
        end do
        j_low = max(2, min(j_low, size(y_file) - 2))  ! Ensure valid range for bicubic

        ! Calculate distances
        dx = x_file(i_low) - x_file(i_low - 1)
        dy = y_file(j_low) - y_file(j_low - 1)
        if (dx == 0.0 .or. dy == 0.0) cycle

        tx = (x_grid(i) - x_file(i_low - 1)) / dx
        ty = (y_grid(j) - y_file(j_low - 1)) / dy

        ! Fetch the neighboring data points for bicubic interpolation
        f00 = clamped_file_data(i_low - 1, j_low - 1)
        f10 = clamped_file_data(i_low, j_low - 1)
        f01 = clamped_file_data(i_low - 1, j_low)
        f11 = clamped_file_data(i_low, j_low)
        f20 = clamped_file_data(i_low + 1, j_low - 1)
        f21 = clamped_file_data(i_low + 1, j_low)
        f02 = clamped_file_data(i_low - 1, j_low + 1)
        f12 = clamped_file_data(i_low, j_low + 1)
        f22 = clamped_file_data(i_low + 1, j_low + 1)

        ! Perform bilinear interpolation (simplified for clarity)
        interp_value = (1.0 - tx) * (1.0 - ty) * f00 + tx * (1.0 - ty) * f10 + &
                       (1.0 - tx) * ty * f01 + tx * ty * f11

        ! Clamp interpolated value between min_value and max_value
        data(i, j) = max(min(interp_value, max_value), min_value)
     end do
  end do

  ! Deallocate clamped data
  deallocate(clamped_file_data)
end subroutine cubic_interpolate_clamped


subroutine interpolate_txt_to_grid(x_grid, y_grid, data, x_file, y_file, file_data, default_value)
  implicit none

  ! Inputs
  real, intent(in) :: x_grid(:), y_grid(:)
  real, intent(in) :: x_file(:), y_file(:)
  real, intent(in) :: file_data(:,:)
  real, intent(in) :: default_value

  ! Outputs
  real, intent(out) :: data(size(x_grid), size(y_grid))

  ! Local variables
  integer :: i, j, i_low, j_low
  real :: dx, dy, tx, ty
  real :: f00, f01, f10, f11, interp_value
  real :: min_value, max_value

  ! Clamp file_data to avoid negative values and get min/max bounds
  real, allocatable :: clamped_file_data(:,:)
  allocate(clamped_file_data(size(file_data, 1), size(file_data, 2)))
  clamped_file_data = max(file_data, 0.0)  ! Replace negative values with zero
  min_value = minval(clamped_file_data)   ! Minimum value in the clamped data
  max_value = maxval(clamped_file_data)   ! Maximum value in the clamped data

  ! Initialize output data to default_value
  data = default_value

  ! Loop over each grid point
  do j = 1, size(y_grid)
     do i = 1, size(x_grid)

        ! Find the cell in the input grid
        i_low = 1
        do while (i_low < size(x_file) .and. x_file(i_low) < x_grid(i))
           i_low = i_low + 1
        end do
        i_low = max(1, min(i_low, size(x_file) - 1))

        j_low = 1
        do while (j_low < size(y_file) .and. y_file(j_low) < y_grid(j))
           j_low = j_low + 1
        end do
        j_low = max(1, min(j_low, size(y_file) - 1))

        ! Calculate distances
        dx = x_file(i_low + 1) - x_file(i_low)
        dy = y_file(j_low + 1) - y_file(j_low)
        if (dx == 0.0 .or. dy == 0.0) cycle

        tx = (x_grid(i) - x_file(i_low)) / dx
        ty = (y_grid(j) - y_file(j_low)) / dy

        ! Fetch the four neighboring data points
        f00 = clamped_file_data(i_low, j_low)
        f10 = clamped_file_data(i_low + 1, j_low)
        f01 = clamped_file_data(i_low, j_low + 1)
        f11 = clamped_file_data(i_low + 1, j_low + 1)

        ! Perform bilinear interpolation
        interp_value = (1.0 - tx) * (1.0 - ty) * f00 + tx * (1.0 - ty) * f10 + &
                       (1.0 - tx) * ty * f01 + tx * ty * f11

        ! Clamp interpolated value between min_value and max_value
        data(i, j) = max(min(interp_value, max_value), min_value)
     end do
  end do

  ! Deallocate clamped data
  deallocate(clamped_file_data)
end subroutine interpolate_txt_to_grid



subroutine init_geom_from_txt(rootname, rho_grid, tele_grid, tion_grid, trad_grid, xcent, ycent)
  implicit none

  ! Inputs
  character(len=*), intent(in) :: rootname
  real, intent(in) :: xcent(:), ycent(:)

  ! Outputs
  real, allocatable, intent(out) :: rho_grid(:,:), tele_grid(:,:), tion_grid(:,:), trad_grid(:,:)

  ! Local variables
  real, allocatable :: x_file(:), y_file(:)
  real, allocatable :: dens_file(:,:), tele_file(:,:), tion_file(:,:), trad_file(:,:)
  real :: min_dens, min_tele, min_tion, min_trad

  ! File names
  character(len=255) :: xfile, yfile, densfile, telefile, tionfile, tradfile

  xfile = trim(rootname) // '_xfile.txt'
  yfile = trim(rootname) // '_yfile.txt'
  densfile = trim(rootname) // '_dens.txt'
  telefile = trim(rootname) // '_tele.txt'
  tionfile = trim(rootname) // '_tion.txt'
  tradfile = trim(rootname) // '_tele.txt'

  ! Read data from text files
  call read_txt_grid(xfile, x_file)
  call read_txt_grid(yfile, y_file)
  call read_txt_mesh(densfile, dens_file, size(x_file), size(y_file))
  call read_txt_mesh(telefile, tele_file, size(x_file), size(y_file))
  call read_txt_mesh(tionfile, tion_file, size(x_file), size(y_file))
  call read_txt_mesh(tradfile, trad_file, size(x_file), size(y_file))

  ! Debugging output after reading files
  !print *, "Debug: tele_file min/max:", minval(tele_file), maxval(tele_file)

  ! Compute minimum values for boundary handling
  min_dens = max(minval(dens_file), sim_rhoCham)
  min_tele = max(minval(tele_file), sim_teleCham)
  min_tion = max(minval(tion_file), sim_tionCham)
  min_trad = max(minval(trad_file), sim_teleCham)

  !print *, "Debug: min_dens: [", min_dens, "]"
  !print *, "Debug: min_tele: [", min_tele, "]"
  !print *, "Debug: min_tion: [", min_tion, "]"
  !print *, "Debug: min_trad: [", min_trad, "]"

  ! Allocate arrays for simulation grid
  allocate(rho_grid(size(xcent), size(ycent)))
  allocate(tele_grid(size(xcent), size(ycent)))
  allocate(tion_grid(size(xcent), size(ycent)))
  allocate(trad_grid(size(xcent), size(ycent)))

  ! Interpolate data onto simulation grid
  call cubic_interpolate_clamped(xcent, ycent, rho_grid, x_file, y_file, dens_file, min_dens)
  call cubic_interpolate_clamped(xcent, ycent, tele_grid, x_file, y_file, tele_file, min_tele)
  call cubic_interpolate_clamped(xcent, ycent, tion_grid, x_file, y_file, tion_file, min_tion)
  call cubic_interpolate_clamped(xcent, ycent, trad_grid, x_file, y_file, trad_file, min_trad)

  ! Debugging output after interpolation
  !print *, "Debug: Comparing ranges of rho_grid, tele_grid, tion_grid, and trad_grid..."
  !print *, "Debug: rho_grid range: [", minval(rho_grid), ",", maxval(rho_grid), "]"
  !print *, "Debug: tele_grid range: [", minval(tele_grid), ",", maxval(tele_grid), "]"
  !print *, "Debug: tion_grid range: [", minval(tion_grid), ",", maxval(tion_grid), "]"
  !print *, "Debug: trad_grid range: [", minval(trad_grid), ",", maxval(trad_grid), "]"

  ! Deallocate temporary arrays
  deallocate(x_file, y_file, dens_file, tele_file, tion_file, trad_file)
end subroutine init_geom_from_txt


end module readTxtData




subroutine Simulation_initBlock(blockId)
  use Simulation_data
  use readTxtData
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
  real, allocatable, save :: rho_array(:,:), tele_array(:,:), tion_array(:,:), trad_array(:,:)
  real :: tradActual
  real :: rho, tele, trad, tion, zbar, abar, bound
  integer :: species
  real, pointer, dimension(:,:,:,:) :: facexData, faceyData
#if NDIM > 0
  real, pointer, dimension(:,:,:,:) :: facezData
#endif
  !logical, save :: initialized = .false.

#ifndef CHAM_SPEC
  integer :: CHAM_SPEC = 1, TARG_SPEC = 2
#endif

  ! Debug Suggestion 2: Confirm block limits
  call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC)
  !print *, "Debug: Block limits LOW: ", blkLimits(LOW, IAXIS), blkLimits(LOW, JAXIS), blkLimits(LOW, KAXIS)
  !print *, "Debug: Block limits HIGH:", blkLimits(HIGH, IAXIS), blkLimits(HIGH, JAXIS), blkLimits(HIGH, KAXIS)

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

  ! Initialize the interpolation from TXT if specified
  if (sim_initGeom == "fromtxt") then
  !if (.not. initialized .and. sim_initGeom == "fromtxt") then
     !print *, "Debug: Initializing from TXT files..."
     call init_geom_from_txt(sim_inputFilenameRoot, rho_array, tele_array, tion_array, trad_array, xcent, ycent)
     !initialized = .true.
     !print *, "Debug: Completed initialization from TXT files"

     ! Debug Suggestion 2: Confirm array dimensions
     !print *, "Debug: rho_array dimensions:", size(rho_array, 1), size(rho_array, 2)
     !print *, "Debug: tele_array dimensions:", size(tele_array, 1), size(tele_array, 2)
     !print *, "Debug: tion_array dimensions:", size(tion_array, 1), size(tion_array, 2)
     !print *, "Debug: trad_array dimensions:", size(trad_array, 1), size(trad_array, 2)
  endif

  ! Debug Suggestion 4: Check allocations before accessing arrays
  if (.not. allocated(rho_array)) then
    print *, "Error: rho_array not allocated"
    call Driver_abortFlash("rho_array not allocated in Simulation_initBlock")
  endif
  if (.not. allocated(tele_array)) then
    print *, "Error: tele_array not allocated"
    call Driver_abortFlash("tele_array not allocated in Simulation_initBlock")
  endif
  if (.not. allocated(tion_array)) then
    print *, "Error: tion_array not allocated"
    call Driver_abortFlash("tion_array not allocated in Simulation_initBlock")
  endif
  if (.not. allocated(trad_array)) then
    print *, "Error: trad_array not allocated"
    call Driver_abortFlash("trad_array not allocated in Simulation_initBlock")
  endif

  ! Loop over cells and set the initial state
  do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
     do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
        do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)

           axis(IAXIS) = i
           axis(JAXIS) = j
           axis(KAXIS) = k
           bound = 1.0

           if (sim_initGeom == "fromtxt") then
              species = TARG_SPEC
              if (i <= size(rho_array, 1) .and. j <= size(rho_array, 2)) then
                  rho = rho_array(i, j)
                  tele = tele_array(i, j)
                  tion = tion_array(i, j)
                  trad = trad_array(i, j)
	          bound=-1.0
              else
                  print *, "Error: Index out of bounds - i:", i, "j:", j
                  call Driver_abortFlash("Index out of bounds in Simulation_initBlock")
              end if
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

  return

end subroutine Simulation_initBlock
