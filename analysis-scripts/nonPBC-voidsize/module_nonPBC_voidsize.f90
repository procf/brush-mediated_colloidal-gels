! Fortran module for analyzing experimental data converted to a BD-style GSD file
! NOTE: requires compilation with compile-module
! NOTE: is run from a matching sim-analysis-BDexp.py Python script
! NOTE: ASSUMES 1 colloid type
! (Deepak Mangal and Rob Campbell)

! Calculates:
! * void size calculation for the final frame or a selection of frames, using two methods: 
!     Torquato's Pore Size Distribution and Gubbibn's Pore Size Distribution  
!     ! NOTE: This is the most complicated code, it has MULTIPLE subroutines
!     !       uses a linked-list implementation, and REQUIRES an additional module
!     !       generated with solvopt.f90 (which is also compiled in compile module)

  ! NOTE: f2py only understands limited number of KIND parameters
  ! use real(kind=8) not real64

  ! NOTE: fortran indexes from 1 by default, switch to base 0 to match C/Python
  ! for indexing from 0, 0:ncolloid-1 is of length ncolloid, 0:2 is length 3

  ! NOTE: must tell f2py about any variable dependencies or you will get Value Errors


!###########################################
! Calculate the void size for all frames using two methods: 
!   - Torquato’s Pore Size Distribution
!   - Gubbins’s Pore Size Distribution
!! NOTE: This is the most complicated code, it has MULTIPLE
!!       subroutines and uses a linked-list implementation
module void_size_calculation
implicit none

! Total number of colloids in simulation box and number of frames in trajectory
integer, public :: ncolloid

!Framewise particles positions and simulation box size
real(8), allocatable, save, public :: rxi(:)         ! framewise x-coordinates
real(8), allocatable, save, public :: ryi(:)         ! framewise y-cordinates
real(8), allocatable, save, public :: rzi(:)         ! framewise z-coordinates
real(8), allocatable, save, public :: box_length(:)  ! framewise box size
real(8), allocatable, save, public :: radii(:)         ! colloid radii, AKA the hard-sphere distance to avoid overlap between probe and colloid

real(8), public :: box_size(3)                         ! simulation box size in a frame
real(8), public :: inv_box_size(3)                     ! inverse of simulation box size in a frame
real(8), public :: rp(3)                               ! void-size calculation location point
real(8), allocatable, public :: rpos(:,:)              ! colloid positions in a given frame

! Variables for linked-list method
real(8), public :: dcell_init                          ! initial cell-size
integer, save, public :: tot_cell                      ! total number of cells
integer, save, public :: ncell(3)                      ! number of cells in each direction
real(8), save, public:: dcell(3)                       ! cell size in each direction
integer, allocatable, save, public :: head(:,:,:)
integer, allocatable, save, public :: list(:)
integer, allocatable, save, public :: colloid_id(:,:)

contains
!------------------------------------------------------------------------

!subroutine for initial random seed

SUBROUTINE init_random_seed()
implicit none
integer :: i, n, clock
integer, allocatable:: seed(:)
CALL RANDOM_SEED(size = n)
ALLOCATE(seed(n))
CALL SYSTEM_CLOCK(COUNT=clock)
seed = clock + 37 * (/ (i - 1, i = 1, n) /)
CALL RANDOM_SEED(PUT = seed)
DEALLOCATE(seed)
END SUBROUTINE init_random_seed
!--------------------------------------------------------------------------

subroutine void_size_calc(data_outpath, frame, ncolloids_py, nprobe, rhs_py, &
dcell_init_py, rxi_py, ryi_py, rzi_py, box_length_py)
implicit none

! Variables imported from Python file
character*256, intent(in) :: data_outpath
integer, intent(in) :: ncolloids_py, nprobe, frame
real(8), intent(in) :: dcell_init_py
real(8), intent(in) :: rhs_py(ncolloids_py)
real(8), intent(in) :: rxi_py(0:ncolloids_py-1)
real(8), intent(in) :: ryi_py(0:ncolloids_py-1)
real(8), intent(in) :: rzi_py(0:ncolloids_py-1)
real(8), intent(in) :: box_length_py(3)
!f2py depend(ncolloids_py) rxi_py,ryi_py,rzi_py
!f2py depend(ncolloids_py) rhs_py

! Internal Variables
integer :: j                         ! index for particle (j)
real(8) :: f_T, f_G                  ! Fraction of void volume (as calculated by Torquato's PSD and Gubbin's PSD)
real(8) :: voidsize_T, voidsize_G    ! void diameter as (as calculated by Torquato's PSD and Gubbin's PSD)
real(8) :: ri(3)                     ! reference location 
character(len=50) :: filename_void   ! output filename/path
logical :: file_exists               ! a logical check to add to existing files

integer,pointer :: atemp => null()
call init_random_seed()

! assign imported variables to global variables
dcell_init = dcell_init_py
ncolloid = ncolloids_py

allocate( rxi(0:ncolloid-1) )
allocate( ryi(0:ncolloid-1) )
allocate( rzi(0:ncolloid-1) )
allocate( box_length(3) )
allocate( radii(0:ncolloid-1) )

radii(:) = rhs_py(:)


do j=0,ncolloid-1
  rxi(j) = rxi_py(j)
  ryi(j) = ryi_py(j)
  rzi(j) = rzi_py(j)
enddo
do j=1,3
  box_length(j) = box_length_py(j)
enddo



! open the output file and write labels
! append if it exists, or create new
filename_void = trim(adjustl(data_outpath))
inquire(file=trim(filename_void), exist=file_exists)
if (file_exists) then
  open(unit=14, file=trim(filename_void), status='old', action='write', position='append')
else
  open(unit=14, file=trim(filename_void), status='new', action='write')
  write(14,fmt="(*(g0:','))") 'frame', 'probe_posx','probe_posy', 'probe_posz', &
'voidcenter_x','voidcenter_y','voidcenter_z','voiddiameter_T','voiddiameter_G'
endif
! calculate the void size

! read box_size in a given frame
box_size(:) = box_length(:)
inv_box_size(:) = 1.d0/box_size(:)

! read colloid positions for the given frame
allocate(rpos(3,0:ncolloid-1))
rpos(:,:)=0.d0
do j=0,ncolloid-1 
  rpos(1,j)=rxi(j)
  rpos(2,j)=ryi(j)
  rpos(3,j)=rzi(j)
enddo

! Linked list initialization, formation, and check
call init_list()
call link_list()
call check_list()

! calculate void size at every random prob point in the box 
do j = 1,nprobe
  f_T=0.d0
  ! select random point in simulation box and check overlap with colloids
  do while( f_T == 0.d0)
    call random_number(rp(1))
    call random_number(rp(2))
    call random_number(rp(3))
    ! rescale position from (0,L) to (-L/2,L/2) coordinates to match sim box
    rp(1)=0.5d0*box_size(1)*(2.d0*rp(1)-1.d0)
    rp(2)=0.5d0*box_size(2)*(2.d0*rp(2)-1.d0)
    rp(3)=0.5d0*box_size(3)*(2.d0*rp(3)-1.d0)
    ! use the current position as the center of a void
    ri(:) = rp(:)
    ! get the fraction of void volume with Torquato's method 
    call fun(ri,f_T)
  enddo

  ! get the voidsize with Torquato's method 
  voidsize_T = 2.d0*dsqrt(-f_T)

  ! get the fraction of void volume with Gubbin's method 
  f_G = 0.d0
  call solvopt(3, ri, f_G, fun, .false., atemp, .true., func, .false., atemp)
  ! get the voidsize with Gubbin's method 
  voidsize_G = 2.d0*dsqrt(-f_G)

  ! save the data for both voidsize calculations
  write(14,fmt="(*(g0:','))") frame, rp, ri, voidsize_T, voidsize_G
  !write(14,'(I8, 8f16.5)')  framechoice(i), rp, ri, voidsize_T, voidsize_G
enddo
deallocate(rpos)
deallocate(rxi)
deallocate(ryi)
deallocate(rzi)
deallocate(box_length)
deallocate(radii)
call finalize_list()

close(14)
end subroutine void_size_calc
!--------------------------------------------

!Linked-list initialization - allocation of head and list arrays

subroutine init_list()
implicit none

! calculate the number of cells
ncell(:) = floor(box_size(:) / dcell_init)
! check number of cells
if (any(ncell(:) < 3)) then
  print *, 'system is too small to use cell links'
  stop
endif
! calculate the total number of cells
tot_cell = ncell(1) * ncell(2) * ncell(3)
! calculate the cell dimensions
dcell(:) = box_size(:) / dble(ncell(:))
! initialize variables for head position, list of colloids, and colloid IDs
allocate( head(0 : ncell(1)-1, 0 : ncell(2)-1, 0 : ncell(3)-1) )
allocate( list(0:ncolloid-1) )
allocate( colloid_id(3,0:ncolloid-1) )
head(:,:,:) = 0
list(:) = 0
colloid_id(:,:) = 0
return
end subroutine init_list
!--------------------------------------------------------------------------

!Compute cell numbers of monomers

subroutine list_cell_i(ri,cell_i)
implicit none
real(8), intent(in) :: ri(3)         ! colloid position
integer, intent(out) :: cell_i(3)    ! colloid index, cell at this index

!if (any(dabs(ri(:)/box_size(:)) > 0.5d0 )) then
!  print *, 'colloid not in the main-box'
!  stop
!end if

cell_i(:) = 0
cell_i(:) = floor( (ri(:)/box_size(:)+0.5d0) * dble(ncell(:)) )
cell_i(:) = modulo( cell_i(:), ncell(:) )
return
end subroutine list_cell_i
!---------------------------------------------------------------------------

!Formation of head and list arrays

subroutine  link_list()
implicit none
real(8) :: ri(3)
integer :: i, cell_i(3)

do i=0,ncolloid-1
  ri(:) = rpos(:,i)
  cell_i(:) = 0
  call list_cell_i(ri,cell_i)
  list(i) = head(cell_i(1),cell_i(2),cell_i(3))
  head(cell_i(1),cell_i(2),cell_i(3)) = i
  colloid_id(:,i) = cell_i(:)
enddo
return
end subroutine link_list
!----------------------------------------------------------------------

!Check formation of head and list arrays

subroutine check_list()
implicit none
integer :: i, j, k      ! indices for x (i), y (j), and z (k) dimensions of each cell
integer :: c            ! the index of the current cell
integer :: cell_i(3)    ! the current cell
real(8) :: ri(3)        ! individual colloid position 
do i=0,ncolloid-1
  ! read in the colloid position
  ri(:) = rpos(:,i)
  ! check the current cell has the correct contents
  cell_i(:) = 0
  call list_cell_i(ri,cell_i)
  if(any(cell_i(:) .ne. colloid_id(:,i))) then
    print *, 'inconsistency1 found, list_cell_i does not returned assigned colloidIDs:', i, cell_i, colloid_id(:,i)
    stop
  endif
enddo
do i=0, ncell(1)-1
  do j=0, ncell(2)-1
    do k=0, ncell(3)-1
      cell_i(:) = (/i,j,k/)
      c = head(i,j,k)
      do while ( c .ne. 0)
        if(any( cell_i(:) .ne. colloid_id(:,c))) then
          print *, 'inconsistency2 found, c=head() does not return assigned colloidIDs', c, cell_i(:), colloid_id(:,c)
          stop
        endif
        c = list(c)
      enddo
    enddo
  enddo
enddo
return
end subroutine check_list
!-----------------------------------------------------------------------

!Deallocation of arrays

subroutine finalize_list()
implicit none
deallocate(list)
deallocate(head)
deallocate(colloid_id)
end subroutine finalize_list
!------------------------------------------------------------------------

!Subroutine to compute void size
subroutine fun(rin,f)
implicit none

! Inputs
real(8),intent(in) :: rin(3)      ! center of the void being calculated 

! Outputs
real(8), intent(out) :: f         ! the calculated fraction of void volume

! Internal Variables
integer :: c                      ! the index for the particle tag in the current cell
integer :: i, j, k                ! indices for x,y,z cell coordinates
integer :: new_i, new_j, new_k    ! updated location of the cell
integer :: cell_i(3)              ! the current cell
real(8) :: dist                   ! scalar r_ij distance between the void center and the nearest particles
real(8) :: rij(3)                 ! the center-center distance between the void center and the nearest particles
real(8) :: ri(3)                  ! location of the center of the void (AKA rc_in)
real(8) :: void_size              ! current void_size estimate (updated during the calculation)
real(8) :: max_size               ! maximum void size (the size of the simulation box) AKA initial void_size estimate

f=0.d0
! initialize the void_size as it's max possible value (sim box volume)
max_size = dmax1(box_size(1),box_size(2),box_size(3))
void_size = max_size

! read the current void location being checked
ri = rin
! adjust the position to include interactions across the periodic boundaries
!ri(:) = ri(:) - box_size(:) * dnint(ri(:)*inv_box_size(:))

! initialize the current cell at zero
cell_i(:) = 0

! find the cell that contains the current void being checked
call list_cell_i(ri,cell_i)

! check/update cell information with neighbors
do i = cell_i(1)-1, cell_i(1)+1
  do j = cell_i(2)-1, cell_i(2)+1
    do k = cell_i(3)-1, cell_i(3)+1
      new_i = modulo(i,ncell(1))
      new_j = modulo(j,ncell(2))
      new_k = modulo(k,ncell(3))
      c = head(new_i,new_j,new_k)

      ! for each colloid in the cell
      do while(c .ne. 0)
        ! calculate the distance between the center of the void and the colloid 
        rij = rpos(:,c) - ri
        ! update the position for interactions across the sim's periodic boundaries
        !rij(:) = rij(:) - box_size(:)* dnint(rij(:)*inv_box_size(:))
        ! convert to scalar and subtract the colloid radius
        dist = dsqrt(sum(rij*rij)) - radii(c)
        ! if the probe is inside a colloid, then dist<0 (and there is no void)
        if(dist .lt. 0.d0) then
          return
        ! else, if the distance is less than the current void size (initialized at max), update the size
        else if (dist .lt. void_size) then
          void_size=dist     
        end if
        ! update c and continue checking through the remaining particles in the cell
        c = list(c)
      enddo
    enddo
  enddo
enddo
! if the void size is still larger than the size of a cell, check against all colloids 
if(void_size .ge. dcell_init) then
  do i = 0,ncolloid-1
    rij = rpos(:,i) - ri
    !rij(:) = rij(:) - box_size(:) * dnint(rij(:)*inv_box_size(:))
    dist = dsqrt(sum(rij*rij)) - radii(i)
    if (dist .lt. void_size) then
      void_size = dist
    endif
  enddo
endif
! use void_size to update f
if(void_size .gt. 0) f = -(void_size) * (void_size)
return
end subroutine fun
!----------------------------------------------------------------------------

!Subroutine to compute constraint value
subroutine func(rin,f)
implicit none

! Inputs
real(8),intent(in) :: rin(3)  ! the center of the current void

! Outputs
real(8), intent(out) :: f     ! the calculated fraction of void volume

! Internal Values
real(8) :: f_local, dr2       ! the input fraction of void volume, the constraint value

! get the Torquato's void size
call fun(rin,f_local)

! check the constraint value (should be less than or equal to 0)
dr2 = sum((rp-rin)*(rp-rin)) + f_local
f = max(0.d0,dr2)
return

end subroutine func
!----------------------------------------------------------------------------

end module void_size_calculation
                 
