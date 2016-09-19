!----------------------------------------------------------------------------
!   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of SPECFEM3D_Cartesian version 3.0 and ASKI version 1.2.
!
!   SPECFEM3D_Cartesian version 3.0 and ASKI version 1.2 are free software: 
!   you can redistribute it and/or modify it under the terms of the GNU 
!   General Public License as published by the Free Software Foundation, 
!   either version 2 of the License, or (at your option) any later version.
!
!   SPECFEM3D_Cartesian version 3.0 and ASKI version 1.2 are distributed in 
!   the hope that they will be useful, but WITHOUT ANY WARRANTY; without 
!   even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
!   PURPOSE.  See the GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with SPECFEM3D_Cartesian version 3.0 and ASKI version 1.2.
!   If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
subroutine prepare_timerun_ASKI()

  use constants,only: IN_DATA_FILES
  use specfem_par,only: CUSTOM_REAL,SIZE_REAL,PRINT_SOURCE_TIME_FUNCTION,NPROC,myrank
  use specfem_par_ASKI

  implicit none

  integer :: iproc
  integer, dimension(1) :: i_array_one_value

  ! before doing anything else, nullify all involved pointers (just to be save for deallocation below)
  nullify(ASKI_Goertzel_U0_local_double,ASKI_Goertzel_U1_local_double,ASKI_Goertzel_U2_local_double,&
       ASKI_Goertzel_U0_local_single,ASKI_Goertzel_U1_local_single,ASKI_Goertzel_U2_local_single)

  call read_Par_file_ASKI()

  if (.not.COMPUTE_ASKI_OUTPUT) return

  if (CUSTOM_REAL /= SIZE_REAL) call exit_MPI_without_rank('only single-precision SPECFEM simulations '//&
       'supported for ASKI output (i.e. CUSTOM_MPI_TYPE must be MPI_REAL in precision.h)')

  if (ASKI_DECONVOLVE_STF.and.(.not.PRINT_SOURCE_TIME_FUNCTION)) call exit_MPI_without_rank('PRINT_SOURCE_TIME_FUNCTION '//&
       'must be set to .true. in Par_file in case of ASKI_DECONVOLVE_STF = .true.')

  ! depending on type of inversion grid, collect information on GLL points
  ! where ASKI output should be produced (like number of points, their indices, 
  ! coordinates, model values)

  ASKI_np_local = 0
  select case(ASKI_type_inversion_grid)
  case(2,3) ! ASKI internal, but SPECFEM independent inversion grids, (scartInversionGrid,ecartInversionGrid) so store at all inner GLL points
     call search_ASKI_wavefield_points_type_invgrid_2_3()
  case(4) ! use SPECFEM elements as inversion grid (specfem3dInversionGrid)
     call search_ASKI_wavefield_points_type_invgrid_4()
  case default
     call exit_MPI_without_rank('values for ASKI_type_inversion_grid other than 2,3,4 not supported yet')
  end select ! ASKI_type_inversion_grid

  ! in the routines search_ASKI_wavefield_points*, the following variables are defined:
  !   ASKI_np_local (number of ASKI wavefield points for this proc)
  !   ASKI_indx_local

  ! gather ASKI_np_local from everybody on proc 0
  call synchronize_all()
  if(myrank == 0) then
     allocate(ASKI_np_local_all(NPROC))
     ASKI_np_local_all(1) = ASKI_np_local ! this is me, rank 0
     do iproc = 1,NPROC-1
        ! receive ASKI_np_local from rank iproc
        call recv_i_t(i_array_one_value,1,iproc)
        ASKI_np_local_all(iproc+1) = i_array_one_value(1)
     end do ! iproc

     if(sum(ASKI_np_local_all) .le. 0) &
          call exit_MPI_without_rank('no ASKI wavefield points found, so no ASKI output can be computed')

  else ! (myrank == 0)
     ! send ASKI_np_local to rank 0
     i_array_one_value(1) = ASKI_np_local
     call send_i_t(i_array_one_value,1,0)

  end if ! (myrank == 0)

  ! define discrete fourier transform factors exp(...) once here, before time loop, plus checks, plus allcoation
  call prepare_ASKI_output()

  if(myrank == 0) call write_ASKI_log_start()

  ! write wavefield points and kernel reference model (and jacobian, in case type_invgrid = 4)
  ! and frequency info to file
  call write_ASKI_main_file()

  if(ASKI_MAIN_FILE_ONLY) then
     ! wait until the main parfile has been written
     call synchronize_all()
     ! abort this run
     call exit_MPI_without_rank("logical parameter 'ASKI_MAIN_FILE_ONLY' in "//trim(IN_DATA_FILES)//&
          "Par_file_ASKI requests to only "//&
          "write the main ASKI output file, hence aborting this run")
  end if

  call synchronize_all()

end subroutine prepare_timerun_ASKI
!
! ----------------------------------------------------------------------------------------------------------
!
subroutine read_Par_file_ASKI()

  use specfem_par,only: myrank
  use constants,only: IN_DATA_FILES,IMAIN
  use specfem_par_ASKI

  implicit none

  character(len=500), dimension(:), allocatable :: val_parfile
  character(len=100), dimension(:), allocatable :: key_parfile
  character(len=601) :: line
  character(len=500) :: val
  integer :: npar,ios,IOASKI,eqindx

  ! open Par_file_ASKI and find number of valid lines
  call get_file_unit_ASKI(IOASKI)
  open(unit=IOASKI,file=trim(IN_DATA_FILES)//'Par_file_ASKI',&
       form='formatted',status='old',action='read',iostat=ios)
  if(ios/=0) then
     close(IOASKI)
     COMPUTE_ASKI_OUTPUT = .false.
     if(myrank==0) call write_ASKI_log('LOG_ASKI_start.txt',"could not open file '"//trim(IN_DATA_FILES)//&
          "Par_file_ASKI', so no ASKI output is produced")
     return
  end if
  ! number of valid lines
  npar = 0
  do while(ios==0)
     read(IOASKI,"(a601)",iostat=ios) line
     if( len_trim(line) > 0 ) then
        line = adjustl(line)
        if(line(1:1) /= '#') then ! ignore comment lines
           eqindx = index(line,'=') ! only allow lines with at least one character in front of '='
           if(eqindx>1) npar = npar + 1
        end if
     end if
  end do
  close(IOASKI)

  if(npar == 0) call exit_MPI_without_rank("no valid lines in file '"//trim(IN_DATA_FILES)//"Par_file_ASKI'")
  allocate(key_parfile(npar),val_parfile(npar))

  ! now open again and store key,val pairs of valid lines
  call get_file_unit_ASKI(IOASKI)
  open(unit=IOASKI,file=trim(IN_DATA_FILES)//'Par_file_ASKI',&
       form='formatted',status='old',action='read',iostat=ios)
  npar = 0
  do while(ios==0)
     read(IOASKI,"(a601)",iostat=ios) line
     if( len_trim(line) > 0 ) then
        line = adjustl(line)
        if(line(1:1) /= '#') then ! ignore comment lines
           eqindx = index(line,'=') ! only allow lines with at least one character in front of '='
           if(eqindx>1) then
                npar = npar + 1
                key_parfile(npar) = line(1:eqindx-1)
                val_parfile(npar) = adjustl(line(eqindx+1:))
           end if
        end if
     end if
  end do
  close(IOASKI)

  ! now set values of variables in module specfem3D_par_ASKI according to content of Par_file_ASKI

  ! COMPUTE_ASKI_OUTPUT
  call get_value_Par_file_ASKI('COMPUTE_ASKI_OUTPUT',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) COMPUTE_ASKI_OUTPUT
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'COMPUTE_ASKI_OUTPUT' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  if(.not.COMPUTE_ASKI_OUTPUT) then
     if(myrank == 0) then
        call write_ASKI_log('LOG_ASKI_start.txt',"in '"//trim(IN_DATA_FILES)//&
             "Par_file_ASKI': COMPUTE_ASKI_OUTPUT is .false., so no ASKI output is produced")
        write(IMAIN,*) "in '"//trim(IN_DATA_FILES)//&
             "Par_file_ASKI': COMPUTE_ASKI_OUTPUT is .false., so no ASKI output is produced"
     end if
     deallocate(key_parfile,val_parfile)
     return
  end if

  ! ASKI_MAIN_FILE_ONLY
  call get_value_Par_file_ASKI('ASKI_MAIN_FILE_ONLY',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_MAIN_FILE_ONLY
  if(ios/=0) call exit_MPI_without_rank("invalid value for logical parameter 'ASKI_MAIN_FILE_ONLY' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! OVERWRITE_ASKI_OUTPUT
  call get_value_Par_file_ASKI('OVERWRITE_ASKI_OUTPUT',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) OVERWRITE_ASKI_OUTPUT
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'OVERWRITE_ASKI_OUTPUT' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! ASKI_outfile
  call get_value_Par_file_ASKI('ASKI_outfile',ASKI_outfile,key_parfile,val_parfile,npar)

  ! ASKI_output_ID
  call get_value_Par_file_ASKI('ASKI_output_ID',val,key_parfile,val_parfile,npar)
  ASKI_output_ID = val(1:length_ASKI_output_ID)

  ! ASKI_DECONVOLVE_STF
  call get_value_Par_file_ASKI('ASKI_DECONVOLVE_STF',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_DECONVOLVE_STF
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_DECONVOLVE_STF' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! ASKI_df
  call get_value_Par_file_ASKI('ASKI_df',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_df
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_df' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")
  if(ASKI_df<0.d0) call exit_MPI_without_rank("value for 'ASKI_df' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI' must be positive")

  ! ASKI_nf
  call get_value_Par_file_ASKI('ASKI_nf',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_nf
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_nf' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")
  if(ASKI_nf<1) call exit_MPI_without_rank("value for 'ASKI_nf' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI' must be positive")

  allocate(ASKI_jf(ASKI_nf))
  ! ASKI_jf
  call get_value_Par_file_ASKI('ASKI_jf',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_jf
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_jf' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! ASKI_DFT_method
  call get_value_Par_file_ASKI('ASKI_DFT_method',val,key_parfile,val_parfile,npar)
  ASKI_DFT_method = val(1:length_ASKI_DFT_method)
  select case(ASKI_DFT_method)
  case('EXPLICIT_SUMMATION','GOERTZEL_STANDARD')
     ! OK, do nothing
  case default
     call exit_MPI_without_rank("invalid value '"//trim(ASKI_DFT_method)//"' for parameter 'ASKI_DFT_method' in '"//&
       trim(IN_DATA_FILES)//"Par_file_ASKI': only values 'EXPLICIT_SUMMATION' and 'GOERTZEL_STANDARD' are supported")
  end select

  ! ASKI_DFT_double
  call get_value_Par_file_ASKI('ASKI_DFT_double',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_DFT_double
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_DFT_double' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! ASKI_DFT_apply_taper
  call get_value_Par_file_ASKI('ASKI_DFT_apply_taper',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_DFT_apply_taper
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_DFT_apply_taper' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! ASKI_DFT_taper_percentage
  call get_value_Par_file_ASKI('ASKI_DFT_taper_percentage',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_DFT_taper_percentage
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_DFT_taper_percentage' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")
  if(ASKI_DFT_taper_percentage<0.d0 .or. ASKI_DFT_taper_percentage>1.d0) &
       call exit_MPI_without_rank("value for 'ASKI_DFT_taper_percentage' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI' must be between 0.0 and 1.0")

  ! ASKI_type_inversion_grid
  call get_value_Par_file_ASKI('ASKI_type_inversion_grid',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_type_inversion_grid
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_type_inversion_grid' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! ASKI_wx
  call get_value_Par_file_ASKI('ASKI_wx',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_wx
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_wx' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! ASKI_wy
  call get_value_Par_file_ASKI('ASKI_wy',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_wy
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_wy' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! ASKI_wz
  call get_value_Par_file_ASKI('ASKI_wz',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_wz
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_wz' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! ASKI_rot_X
  call get_value_Par_file_ASKI('ASKI_rot_X',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_rot_X
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_rot_X' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! ASKI_rot_Y
  call get_value_Par_file_ASKI('ASKI_rot_Y',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_rot_Y
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_rot_Y' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! ASKI_rot_Z
  call get_value_Par_file_ASKI('ASKI_rot_Z',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_rot_Z
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_rot_Z' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! ASKI_cx
  call get_value_Par_file_ASKI('ASKI_cx',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_cx
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_cx' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! ASKI_cy
  call get_value_Par_file_ASKI('ASKI_cy',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_cy
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_cy' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! ASKI_cz
  call get_value_Par_file_ASKI('ASKI_cz',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_cz
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_cz' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")


!!IOASKI  COMPUTE_ASKI_OUTPUT = .true.
!!$  ASKI_outfile = &
!!$'/rscratch/minos27/Kernel/specfem3D/inversions/test/new_specfem3d/iteration_step_001/kernel_displacements/S001'
!!$  ASKI_df = 10.0
!!$  ASKI_nf = 23
!!$  ASKI_jf = (/11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33/)
!!$  ASKI_type_inversion_grid = 2
!!$  ASKI_cx = 0.0  
!!$  ASKI_cy = 0.0  
!!$  ASKI_cz = 155.0
!!$  ASKI_wx = 128.0
!!$  ASKI_wy = 128.0
!!$  ASKI_wz = 128.0


  if(allocated(key_parfile)) deallocate(key_parfile)
  if(allocated(val_parfile)) deallocate(val_parfile)
end subroutine read_Par_file_ASKI
!
! ----------------------------------------------------------------------------------------------------------
!
subroutine get_value_Par_file_ASKI(key,val,key_parfile,val_parfile,npar)
  use constants,only: IN_DATA_FILES
  character(len=*), intent(in) :: key
  integer, intent(in) :: npar
  character(len=*), dimension(npar), intent(in) :: key_parfile,val_parfile
  character(len=500), intent(out) :: val
  integer :: ipar
  logical :: found
  found = .false.
  do ipar = 1,size(key_parfile)
     if(key == key_parfile(ipar)) then
        val = val_parfile(ipar)
        found = .true.
        exit
     end if
  end do ! ipar
  if(.not.found) call exit_MPI_without_rank("definition of parameter '"//trim(key)//"' not found in '"&
          //trim(IN_DATA_FILES)//"Par_file_ASKI'")
end subroutine get_value_Par_file_ASKI
!
! ----------------------------------------------------------------------------------------------------------
!
subroutine search_ASKI_wavefield_points_type_invgrid_2_3()

  use specfem_par,only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NSPEC_AB,ibool,xstore,ystore,zstore
  use specfem_par_ASKI

  implicit none

  ! local variables
  integer, dimension(:,:), allocatable :: ASKI_indx_local_tmp
  real(kind=CUSTOM_REAL) :: xtmp,ytmp,ztmp
  integer :: ispec,i,j,k,iglob
  real, dimension(3,3) :: Mrot_tmp,Mrot
  real, dimension(3) :: xyz_rotated
  logical :: apply_rotation
  real, parameter :: deg2rad = 0.017453292519943295

  ! define inverse transformation matrix Mrot, which applies the inverse (i.e. clockwise) rotations defined by angles ASKI_rot_(XYZ)
  ! Mrot is used to back-transform a point in real coordinates to non-rotated block
  apply_rotation = .false.
  select case(ASKI_type_inversion_grid)
  case (2)
     if(ASKI_rot_Z /= 0.0) then
        apply_rotation = .true.
        Mrot(1,:) = (/  cos(ASKI_rot_Z*deg2rad), sin(ASKI_rot_Z*deg2rad), 0. /)
        Mrot(2,:) = (/ -sin(ASKI_rot_Z*deg2rad), cos(ASKI_rot_Z*deg2rad), 0. /)
        Mrot(3,:) = (/ 0., 0. , 1. /)
     end if
  case (3)
     if(ASKI_rot_X /= 0.0 .or. ASKI_rot_Y /= 0.0 .or. ASKI_rot_Z /= 0.0) then
        apply_rotation = .true.
        Mrot(1,:) = (/ 1., 0., 0. /)
        Mrot(2,:) = (/ 0., 1., 0. /)
        Mrot(3,:) = (/ 0., 0., 1. /)
        if(ASKI_rot_X /= 0.0) then
           Mrot_tmp(1,:) = (/ 1., 0. , 0. /)
           Mrot_tmp(2,:) = (/ 0.,  cos(ASKI_rot_X*deg2rad), sin(ASKI_rot_X*deg2rad) /)
           Mrot_tmp(3,:) = (/ 0., -sin(ASKI_rot_X*deg2rad), cos(ASKI_rot_X*deg2rad) /)
           Mrot = matmul(Mrot,Mrot_tmp)
        end if
        if(ASKI_rot_Y /= 0.0) then
           Mrot_tmp(1,:) = (/  cos(ASKI_rot_Y*deg2rad), 0., sin(ASKI_rot_Y*deg2rad) /)
           Mrot_tmp(2,:) = (/ 0., 1. , 0. /)
           Mrot_tmp(3,:) = (/ -sin(ASKI_rot_Y*deg2rad), 0., cos(ASKI_rot_Y*deg2rad) /)
           Mrot = matmul(Mrot,Mrot_tmp)
        end if
        if(ASKI_rot_Z /= 0.0) then
           Mrot_tmp(1,:) = (/  cos(ASKI_rot_Z*deg2rad), sin(ASKI_rot_Z*deg2rad), 0. /)
           Mrot_tmp(2,:) = (/ -sin(ASKI_rot_Z*deg2rad), cos(ASKI_rot_Z*deg2rad), 0. /)
           Mrot_tmp(3,:) = (/ 0., 0. , 1. /)
           Mrot = matmul(Mrot,Mrot_tmp)
        end if
     end if
  end select

  allocate(ASKI_indx_local_tmp((NGLLX-2)*(NGLLY-2)*(NGLLZ-2)*NSPEC_AB,4))
  ASKI_indx_local_tmp(:,:) = 0

  ASKI_np_local = 0

  ! loop only on points inside the element. That results in a sufficiently uniform scatter of gridpoints in case of NGLL = 5
  do ispec=1,NSPEC_AB
     do k=2,NGLLZ-1
        do j=2,NGLLY-1
           do i=2,NGLLX-1

              iglob = ibool(i,j,k,ispec)
              xtmp = xstore(iglob) - ASKI_cx ! shift x,y,z coordinates back to original center of ASKI volume x=y=z=0
              ytmp = ystore(iglob) - ASKI_cy
              ztmp = zstore(iglob) - ASKI_cz

              ! if there is any rotation of the ASKI volume, apply the inverse rotation here ...
              if(apply_rotation) then
                 xyz_rotated = matmul(Mrot,(/xtmp, ytmp, ztmp/))
                 xtmp = xyz_rotated(1)
                 ytmp = xyz_rotated(2)
                 ztmp = xyz_rotated(3)
              end if

              ! ... before checking if the shifted and back-rotated point lies in standard block with defined width
              if (         xtmp .ge. ( - ASKI_wx/2.0) .and. xtmp .le. (ASKI_wx/2.0) &
                   & .and. ytmp .ge. ( - ASKI_wy/2.0) .and. ytmp .le. (ASKI_wy/2.0) &
                   & .and. ztmp .ge. ( - ASKI_wz/2.0) .and. ztmp .le. (ASKI_wz/2.0) ) then

                 ! increment index of points found in kernel chunk
                 ASKI_np_local = ASKI_np_local + 1 

                 ! store index of element
                 ASKI_indx_local_tmp(ASKI_np_local,1) = ispec
                 ! store index of x - gridpoint in that element
                 ASKI_indx_local_tmp(ASKI_np_local,2) = i
                 ! store index of y - gridpoint in that element
                 ASKI_indx_local_tmp(ASKI_np_local,3) = j
                 ! store index of z - gridpoint in that element
                 ASKI_indx_local_tmp(ASKI_np_local,4) = k  

              end if ! current point lies in ASKI output volume

           end do ! i
        end do ! j
     enddo ! k
  enddo ! ispec

  if(ASKI_np_local .gt. 0) then
     allocate(ASKI_indx_local(ASKI_np_local,4))
     ASKI_indx_local(:,:) = ASKI_indx_local_tmp(1:ASKI_np_local,:)
  end if

  if(allocated(ASKI_indx_local_tmp)) deallocate(ASKI_indx_local_tmp)

end subroutine search_ASKI_wavefield_points_type_invgrid_2_3
!
! ----------------------------------------------------------------------------------------------------------
!
subroutine search_ASKI_wavefield_points_type_invgrid_4()

  use specfem_par,only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NSPEC_AB,ibool,xstore,ystore,zstore
  use specfem_par_ASKI

  implicit none

  ! local variables
  real(kind=CUSTOM_REAL) :: xtmp,ytmp,ztmp
  integer :: ispec,jspec,i,j,k,iglob,nspec_in_ASKI_volume,ip
  real, dimension(3,3) :: Mrot_tmp,Mrot
  real, dimension(3) :: xyz_rotated
  logical :: apply_rotation
  real, parameter :: deg2rad = 0.017453292519943295
  integer, dimension(:), allocatable :: ispec_in_ASKI_volume
  logical :: kji_loop_exit

  ! define inverse transformation matrix Mrot, which applies the inverse (i.e. clockwise) rotations defined by angles ASKI_rot_(XYZ)
  ! Mrot is used to back-transform a point in real coordinates to non-rotated block
  apply_rotation = .false.
  if(ASKI_rot_X /= 0.0 .or. ASKI_rot_Y /= 0.0 .or. ASKI_rot_Z /= 0.0) then
     apply_rotation = .true.
     Mrot(1,:) = (/ 1., 0., 0. /)
     Mrot(2,:) = (/ 0., 1., 0. /)
     Mrot(3,:) = (/ 0., 0., 1. /)
     if(ASKI_rot_X /= 0.0) then
        Mrot_tmp(1,:) = (/ 1., 0. , 0. /)
        Mrot_tmp(2,:) = (/ 0.,  cos(ASKI_rot_X*deg2rad), sin(ASKI_rot_X*deg2rad) /)
        Mrot_tmp(3,:) = (/ 0., -sin(ASKI_rot_X*deg2rad), cos(ASKI_rot_X*deg2rad) /)
        Mrot = matmul(Mrot,Mrot_tmp)
     end if
     if(ASKI_rot_Y /= 0.0) then
        Mrot_tmp(1,:) = (/  cos(ASKI_rot_Y*deg2rad), 0., sin(ASKI_rot_Y*deg2rad) /)
        Mrot_tmp(2,:) = (/ 0., 1. , 0. /)
        Mrot_tmp(3,:) = (/ -sin(ASKI_rot_Y*deg2rad), 0., cos(ASKI_rot_Y*deg2rad) /)
        Mrot = matmul(Mrot,Mrot_tmp)
     end if
     if(ASKI_rot_Z /= 0.0) then
        Mrot_tmp(1,:) = (/  cos(ASKI_rot_Z*deg2rad), sin(ASKI_rot_Z*deg2rad), 0. /)
        Mrot_tmp(2,:) = (/ -sin(ASKI_rot_Z*deg2rad), cos(ASKI_rot_Z*deg2rad), 0. /)
        Mrot_tmp(3,:) = (/ 0., 0. , 1. /)
        Mrot = matmul(Mrot,Mrot_tmp)
     end if
  end if

  allocate(ispec_in_ASKI_volume(NSPEC_AB))
  ! loop on all points inside an element. If any point is contained in the ASKI volume, remember the element index and
  ! later add THE WHOLE element (and all points contained in it) to the list of inversion grid cells / wavefield points
  nspec_in_ASKI_volume = 0
  do ispec=1,NSPEC_AB
     kji_loop_exit = .false.
     do k=1,NGLLZ
        do j=1,NGLLY
           do i=1,NGLLX

              iglob = ibool(i,j,k,ispec)
              xtmp = xstore(iglob) - ASKI_cx ! shift x,y,z coordinates back to original center of ASKI volume x=y=z=0
              ytmp = ystore(iglob) - ASKI_cy
              ztmp = zstore(iglob) - ASKI_cz

              ! if there is any rotation of the ASKI volume, apply the inverse rotation here ...
              if(apply_rotation) then
                 xyz_rotated = matmul(Mrot,(/xtmp, ytmp, ztmp/))
                 xtmp = xyz_rotated(1)
                 ytmp = xyz_rotated(2)
                 ztmp = xyz_rotated(3)
              end if

              ! ... before checking if the shifted and back-rotated point lies in standard block with defined width
              if (         xtmp .ge. ( - ASKI_wx/2.0) .and. xtmp .le. (ASKI_wx/2.0) &
                   & .and. ytmp .ge. ( - ASKI_wy/2.0) .and. ytmp .le. (ASKI_wy/2.0) &
                   & .and. ztmp .ge. ( - ASKI_wz/2.0) .and. ztmp .le. (ASKI_wz/2.0) ) then

                 nspec_in_ASKI_volume = nspec_in_ASKI_volume + 1
                 ispec_in_ASKI_volume(nspec_in_ASKI_volume) = ispec
                 kji_loop_exit = .true.

              end if ! current point lies in ASKI output volume

              if(kji_loop_exit) exit
           end do ! i
           if(kji_loop_exit) exit
        end do ! j
        if(kji_loop_exit) exit
     enddo ! k
  enddo ! ispec

  if(nspec_in_ASKI_volume > 0) then
     ! store wavefield point information
     ASKI_np_local = NGLLX*NGLLY*NGLLZ*nspec_in_ASKI_volume
     allocate(ASKI_indx_local(ASKI_np_local,4))

     ip = 0
     do jspec=1,nspec_in_ASKI_volume; ispec = ispec_in_ASKI_volume(jspec)
        do k=1,NGLLZ
           do j=1,NGLLY
              do i=1,NGLLX
                 ! increment point index
                 ip = ip + 1 

                 ! store index of element
                 ASKI_indx_local(ip,1) = ispec
                 ! store index of x - gridpoint in that element
                 ASKI_indx_local(ip,2) = i
                 ! store index of y - gridpoint in that element
                 ASKI_indx_local(ip,3) = j
                 ! store index of z - gridpoint in that element
                 ASKI_indx_local(ip,4) = k  
              end do ! i
           end do ! j
        enddo ! k
     enddo ! ispec

  else ! nspec_in_ASKI_volume > 0
     ASKI_np_local = 0
  end if ! nspec_in_ASKI_volume > 0

  if(allocated(ispec_in_ASKI_volume)) deallocate(ispec_in_ASKI_volume)
end subroutine search_ASKI_wavefield_points_type_invgrid_4
!
! ----------------------------------------------------------------------------------------------------------
!
subroutine prepare_ASKI_output()

  use specfem_par,only: DT,NSTEP,PI,TWO_PI,myrank,USE_RICKER_TIME_FUNCTION,USE_FORCE_POINT_SOURCE
  use specfem_par_ASKI

  implicit none

  double precision :: wtaper!,nefactors
  integer :: jt,jf,IOASKI,ios
  logical :: file_exists
  character(len=509) :: filename

  if(ASKI_DECONVOLVE_STF .and. USE_RICKER_TIME_FUNCTION) call exit_MPI_without_rank("ASKI_DECONVOLVE_STF = .true. "//&
       "is only supported for case USE_RICKER_TIME_FUNCTION = .false.")

  if(.not.ASKI_DECONVOLVE_STF) then
     ! always store displacement in this case (and not velocity)
     ASKI_store_veloc = .false.
  else
     ! if the stf should be deconvolved, we must distinguish between point source and moment tensor source
     ! (since in the SPECFEM3D_Cartesian release by June 2015 either a Gaussian or an error function are used, 
     !  dependent on the type of source mechanism)
     if(USE_FORCE_POINT_SOURCE) then
        ! a thin gaussian is used, so store displacement and deconvolve the source time function directly
        ASKI_store_veloc = .false.
     else
        ! a steep error function is used (integral of thin gaussian), so store velocity and deconvolve the
        ! differentiated source time function (i.e. gaussian). this is done for reasons of numerical stability, 
        ! deconvolving a spectrum which is nearly constant 1.0 )
        ASKI_store_veloc = .true.
     end if
  end if

  if(myrank == 0) then

     do jf=0,ASKI_nf
        if(jf==0) then
           filename = trim(ASKI_outfile)//".main"
        else
           write(filename,"(a,i6.6)") trim(ASKI_outfile)//".jf",ASKI_jf(jf)
        end if

        inquire(file=filename,exist=file_exists)

        ! check here, according to OVERWRITE_ASKI_OUTPUT: is there any conflict with existing files?
        ! if so, abort the current run
        if(file_exists .and. .not.OVERWRITE_ASKI_OUTPUT) &
             call exit_MPI(myrank,"output file '"//trim(filename)//"' "//&
             "already exists and must not be overwritten according to OVERWRITE_ASKI_OUTPUT")

        ! check if files can be opened to write
        ! open file to write with status='unknown' in order not to delete the files when opening
        call get_file_unit_ASKI(IOASKI)
        open(unit=IOASKI, file=filename, status='unknown', access="stream", &
             form="unformatted", action='write',iostat=ios)
        if(ios/=0) call exit_MPI(myrank,"will not be able to write ASKI output file '"//trim(filename)//"'")
        if(file_exists) then
           close(IOASKI,status='keep')
        else
           close(IOASKI,status='delete')
        end if

     end do ! jf        

  end if ! myrank == 0

  ! check if frequency discretization as defined in Par_file_ASKI is valid:
  if(ASKI_df < 0.d0) call exit_MPI_without_rank('ASKI_df in  must not be smaller than zero')
  if(ASKI_nf < 1) call exit_MPI_without_rank('ASKI_nf in Par_file_ASKI must be a strictly positive number')
  if(size(ASKI_jf) /= ASKI_nf) call exit_MPI_without_rank('size(ASKI_jf) must equal ASKI_nf in Par_file_ASKI')

  ! the variables ASKI_efactors and ASKI_Goertzel_Wr are only used by procs which compute any ASKI output 
  ! at local wavefield points 
  ! ADDITIONALLY they are used by rank 0 below in deconvolve_stf_from_ASKI_output even if it does 
  ! not have local wavefield points (therefore DO NOT incorporate tapering into efactors!)
  if (ASKI_np_local .gt. 0 .or. myrank == 0) then

     ! if you plan on doing an inverse fourier transform afterwards, ASKI_df should be chosen 
     ! in a way, that ASKI_df = 1/(NSTEP-1)*DT (1/length_of_timeseries), 
     ! resulting in N = NSTEP-1 in the formula if exp^(-i2pi(k)(n)/N)
     ! which matches the general rule of computing the discrete fourier transform.
     ! If N is not integer, the forward fourier transform works out fine, but it is in general
     ! problematic to do an inverse fourier transform afterwards, as the exp^(...) are no roots of 1 anymore
!!$     nefactors = 1./(ASKI_df*DT)

     select case(ASKI_DFT_method)
     case('EXPLICIT_SUMMATION')
        ! allocate and compute efactors
        allocate(ASKI_efactors(ASKI_nf,NSTEP))
        do jt = 1,NSTEP
           do jf = 1,ASKI_nf
!!$           ASKI_efactors(jf,jt) = cexp(-sngl(TWO_PI)*cmplx(0.,1.)*ASKI_jf(jf)*(jt-1)/sngl(nefactors))*DT
              ASKI_efactors(jf,jt) = exp(-(0.d0,1.d0)*TWO_PI*(ASKI_jf(jf)*ASKI_df)*(dble(jt-1)*DT))*DT
           end do ! jf
        end do ! jt
     case('GOERTZEL_STANDARD')
        allocate(ASKI_Goertzel_Wr(ASKI_nf),ASKI_Goertzel_Wi(ASKI_nf))
        do jf = 1,ASKI_nf
           ! Also compare ASKI developers manual, section 4.5:
           ! In order to account for the time-reversal in Goertzel's algorithm, choose
           ! x = +omega*dt for computing Wr = 2*cos(x)
           ASKI_Goertzel_Wr(jf) = 2.d0 * cos( TWO_PI*ASKI_jf(jf)*ASKI_df*DT )
           ASKI_Goertzel_Wi(jf) = sin( TWO_PI*ASKI_jf(jf)*ASKI_df*DT )
        end do
     end select

     if(ASKI_DFT_apply_taper) then
        if(ASKI_DFT_taper_percentage<0.d0 .or. ASKI_DFT_taper_percentage>1.d0) &
             call exit_MPI_without_rank('ASKI_DFT_taper_percentage must be between 0.0 and 1.0 in Par_file_ASKI')
        ! set hanning taper parameters, taper the last ASKI_DFT_taper_percentage of the time series
        wtaper = DT*dble(NSTEP-1)*ASKI_DFT_taper_percentage
        ASKI_DFT_ntaper_start = NSTEP - int(wtaper/DT) + 1 ! 2 <= ASKI_DFT_ntaper_start <= NSTEP+1

        if(ASKI_DFT_ntaper_start<=NSTEP) then
           allocate(ASKI_DFT_taper_values(ASKI_DFT_ntaper_start:NSTEP))
           do jt = ASKI_DFT_ntaper_start,NSTEP
              ASKI_DFT_taper_values(jt) = 0.5d0*(1.d0-cos(PI*DT*dble(NSTEP-jt)/wtaper))
           end do ! jt
        else ! ASKI_DFT_ntaper_start<=NSTEP
           ASKI_DFT_apply_taper = .false.
        end if ! ASKI_DFT_ntaper_start<=NSTEP
     end if ! ASKI_DFT_apply_taper
  end if ! ASKI_np_local .gt. 0 .or. myrank == 0

  ! the following allocations are only needed for procs which compute any ASKI output at local wavefield points 
  if (ASKI_np_local .gt. 0) then
     select case(ASKI_DFT_method)
     case('EXPLICIT_SUMMATION')
        ! allocate for spectra, first rank: 3 underived components, plus 6 strains = 9
        if(ASKI_DFT_double) then
           allocate(ASKI_spectra_local_double(9,ASKI_nf,ASKI_np_local))
           ASKI_spectra_local_double = (0.d0,0.d0)
        else
           allocate(ASKI_spectra_local_single(9,ASKI_nf,ASKI_np_local))
           ASKI_spectra_local_single = (0.0,0.0)
        end if
     case('GOERTZEL_STANDARD')
        ! allocate for Goertzel's algorithm, first rank: 3 underived components, plus 6 strains = 9
        if(ASKI_DFT_double) then
           allocate(ASKI_Goertzel_U0_local_double(9,ASKI_nf,ASKI_np_local),&
                ASKI_Goertzel_U1_local_double(9,ASKI_nf,ASKI_np_local),&
                ASKI_Goertzel_U2_local_double(9,ASKI_nf,ASKI_np_local))
           ASKI_Goertzel_U1_local_double = 0.d0
           ASKI_Goertzel_U2_local_double = 0.d0
        else
           allocate(ASKI_Goertzel_U0_local_single(9,ASKI_nf,ASKI_np_local),&
                ASKI_Goertzel_U1_local_single(9,ASKI_nf,ASKI_np_local),&
                ASKI_Goertzel_U2_local_single(9,ASKI_nf,ASKI_np_local))
           ASKI_Goertzel_U1_local_single = 0.0
           ASKI_Goertzel_U2_local_single = 0.0
        end if
     end select
  end if ! ASKI_np_local .gt. 0

end subroutine prepare_ASKI_output
!
!-------------------------------------------------------------------------------------------
!
subroutine write_ASKI_main_file()

  use specfem_par,only: CUSTOM_REAL,ibool,xstore,ystore,zstore,kappastore,mustore,rhostore,&
       NGLLX,NGLLY,NGLLZ,jacobian,myrank,NPROC,itag
  use specfem_par_ASKI

  implicit none

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: xyz_local,xyz
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: model_local,model
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: jacob_local,jacob
  integer :: ip,iglob,ispec,i,j,k
  integer :: iproc,specfem_version
  integer :: IOASKI
  logical :: file_exists
  character (len=7) :: open_status
  integer, dimension(:,:), allocatable :: neighbours
  integer :: nnb,ncell

  ! first, locally collect all coordinates of wavefield points and compute kernel reference model
  if(ASKI_np_local > 0) then
     allocate(xyz_local(ASKI_np_local,3),model_local(ASKI_np_local,3))

     do ip = 1,ASKI_np_local
        ! set ispec,i,j,k
        ispec = ASKI_indx_local(ip,1)
        i = ASKI_indx_local(ip,2)
        j = ASKI_indx_local(ip,3)
        k = ASKI_indx_local(ip,4)
        iglob = ibool(i,j,k,ispec)

        ! x, y, z
        xyz_local(ip,1) = xstore(iglob)
        xyz_local(ip,2) = ystore(iglob)
        xyz_local(ip,3) = zstore(iglob)

        ! rho
        model_local(ip,1) = rhostore(i,j,k,ispec)
        ! vp
        model_local(ip,2) = sqrt((kappastore(i,j,k,ispec) + 4.*mustore(i,j,k,ispec)/3.)/rhostore(i,j,k,ispec))
        ! vs
        model_local(ip,3) = sqrt(mustore(i,j,k,ispec)/rhostore(i,j,k,ispec))

     end do ! ip

     ! in case of specfem3dInversionGrid, additionally store jacobian at all points
     select case(ASKI_type_inversion_grid)
     case(4)
        allocate(jacob_local(ASKI_np_local))
        do ip = 1,ASKI_np_local
           ! set ispec,i,j,k
           ispec = ASKI_indx_local(ip,1)
           i = ASKI_indx_local(ip,2)
           j = ASKI_indx_local(ip,3)
           k = ASKI_indx_local(ip,4)

           ! jacobian
           jacob_local(ip) = jacobian(i,j,k,ispec)
        end do ! ip
     end select

  end if ! ASKI_np_local > 0

  ! gather wavefield points on proc 0
  call synchronize_all()
  if(myrank == 0) then
     allocate(xyz(sum(ASKI_np_local_all),3))
     ip = 0

     ! this is me, rank 0
     if(ASKI_np_local > 0) then
        xyz(1:ASKI_np_local,:) = xyz_local(:,:)
        ip = ip + ASKI_np_local
     end if
     do iproc = 1,NPROC-1
        if(ASKI_np_local_all(iproc+1) > 0) then
           call recvv_cr(xyz(ip+1:ip+ASKI_np_local_all(iproc+1),:), ASKI_np_local_all(iproc+1)*3, iproc, itag)
           ip = ip + ASKI_np_local_all(iproc+1)
        end if
     end do ! iproc

  else ! myrank == 0

     if(ASKI_np_local > 0) call sendv_cr(xyz_local, ASKI_np_local*3, 0, itag)

  end if ! myrank == 0

  select case(ASKI_type_inversion_grid)
  case(4)
     ! gather jacobian on proc 0 and compute neighbours
     call synchronize_all()
     if(myrank == 0) then
        allocate(jacob(sum(ASKI_np_local_all)))
        ip = 0

        ! this is me, rank 0
        if(ASKI_np_local > 0) then
           jacob(1:ASKI_np_local) = jacob_local(:)
           ip = ip + ASKI_np_local
        end if
        do iproc = 1,NPROC-1
           if(ASKI_np_local_all(iproc+1) > 0) then
              call recvv_cr(jacob(ip+1:ip+ASKI_np_local_all(iproc+1)), ASKI_np_local_all(iproc+1), iproc, itag)
              ip = ip + ASKI_np_local_all(iproc+1)
           end if
        end do ! iproc

        ! find neighbouring elements
        ncell = sum(ASKI_np_local_all)/(NGLLX*NGLLY*NGLLZ)
        allocate(neighbours(7,ncell)); neighbours = 0
        call find_ASKI_neighbours_type_invgrid_4(neighbours,ncell,xyz,sum(ASKI_np_local_all))

     else ! myrank == 0

        if(ASKI_np_local > 0) &
             call sendv_cr(jacob_local, ASKI_np_local, 0, itag)

     end if ! myrank == 0
  end select ! ASKI_type_inversion_grid

  ! gather model on proc 0
  call synchronize_all()
  if(myrank == 0) then
     allocate(model(sum(ASKI_np_local_all),3))
     ip = 0

     ! this is me, rank 0
     if(ASKI_np_local > 0) then
        model(1:ASKI_np_local,:) = model_local(:,:)
        ip = ip + ASKI_np_local
     end if
     do iproc = 1,NPROC-1
        if(ASKI_np_local_all(iproc+1) > 0) then
           call recvv_cr(model(ip+1:ip+ASKI_np_local_all(iproc+1),:), ASKI_np_local_all(iproc+1)*3, iproc, itag)
           ip = ip + ASKI_np_local_all(iproc+1)
        end if
     end do ! iproc

  else ! myrank == 0

     if(ASKI_np_local > 0) &
          call sendv_cr(model_local, ASKI_np_local*3, 0, itag)

  end if ! myrank == 0

  if(myrank == 0) then
     ! do not worry about whether to overwrite or not, if program comes here, it already has been checked 
     ! if the value of OVERWRITE_ASKI_OUTPUT is in conflict with existing files
     inquire(file=trim(ASKI_outfile)//".main",exist=file_exists)
     if(file_exists) then
        open_status = 'replace'
     else
        open_status = 'new'
     end if

     specfem_version = 2
     ! master proc writes out file containing general info and wavefield points and kernel reference model
     ! also write NPROC to file, as a safety feature in order to assure that wavefield points (written only in ASKI_main_file
     ! and kernel values  (written only in kernel files) are not confused
     call get_file_unit_ASKI(IOASKI)
     open(unit=IOASKI, file=trim(ASKI_outfile)//".main", status=open_status, access="stream", form="unformatted")
     write(IOASKI) specfem_version,length_ASKI_output_ID,ASKI_output_ID,NPROC,ASKI_type_inversion_grid,&
          sum(ASKI_np_local_all),ASKI_df,ASKI_nf
     write(IOASKI) ASKI_jf
     write(IOASKI) xyz(:,1)
     write(IOASKI) xyz(:,2)
     write(IOASKI) xyz(:,3)
     write(IOASKI) model(:,1)
     write(IOASKI) model(:,2)
     write(IOASKI) model(:,3)
     select case(ASKI_type_inversion_grid)
     case(4)
        write(IOASKI) NGLLX,NGLLY,NGLLZ
        write(IOASKI) jacob
        write(IOASKI) ncell,sum(neighbours(1,:))+ncell !sum(neighbours(1,:))+ncell is the total number of integers following now, before model values
        do ispec = 1,ncell
           nnb = neighbours(1,ispec)
           write(IOASKI) nnb
           if(nnb > 0) write(IOASKI) neighbours(2:1+nnb,ispec)
        end do
     end select
     close(IOASKI)
  end if ! myrank == 0

  ! deallocate local stuff
  if(allocated(model)) deallocate(model)
  if(allocated(model_local)) deallocate(model_local)
  if(allocated(xyz)) deallocate(xyz)
  if(allocated(xyz_local)) deallocate(xyz_local)
  if(allocated(jacob)) deallocate(jacob)
  if(allocated(jacob_local)) deallocate(jacob_local)
  if(allocated(neighbours)) deallocate(neighbours)

end subroutine write_ASKI_main_file
!
!-------------------------------------------------------------------------------------------
!
subroutine find_ASKI_neighbours_type_invgrid_4(neighbours,ncell,xyz,np)

  use specfem_par,only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

  implicit none

  integer :: ncell,np
  integer, dimension(7,ncell) :: neighbours
  real(kind=CUSTOM_REAL), dimension(np,3) :: xyz
  ! local
  integer :: icell,jcell,iface,nnb,cell_shift,np_cell,&
       nb1_shift,nb2_shift,nb3_shift,nb4_shift,nb5_shift,nb6_shift
  real(kind=CUSTOM_REAL), dimension(6*ncell) :: xnb,ynb,znb
  real(kind=CUSTOM_REAL), dimension(3) :: v
  real :: eps

  if(mod(NGLLX,2)/=1 .or. mod(NGLLY,2)/=1 .or. mod(NGLLZ,2)/=1) &
       call exit_MPI_without_rank("Neighbour search for specfem3dInversionGrid for "//&
       "ASKI only supported for odd NGLLX,NGLLY,NGLLZ")

  neighbours = 0

  ! loop on all cells and store their 6 face-center points in arrays xnb,ynb,znb

  nb1_shift = NGLLX*NGLLY*(NGLLZ-1)/2 + (NGLLX+1)/2 ! point-index of ymin-face relative to cell
  nb2_shift = NGLLX*NGLLY*(NGLLZ-1)/2 + NGLLX*(NGLLY+1)/2 ! xmax-face
  nb3_shift = NGLLX*NGLLY*(NGLLZ-1)/2 + NGLLX*(NGLLY-1) + (NGLLX+1)/2 ! ymax-face
  nb4_shift = NGLLX*NGLLY*(NGLLZ-1)/2 + NGLLX*(NGLLY-1)/2 + 1 ! xmin-face
  nb5_shift = NGLLX*(NGLLY-1)/2 + (NGLLX+1)/2 ! zmin-face
  nb6_shift = NGLLX*NGLLY*(NGLLZ-1) + NGLLX*(NGLLY-1)/2 + (NGLLX+1)/2 ! zmax-face
  np_cell = NGLLX*NGLLY*NGLLZ

  cell_shift = 0
  do icell = 1,ncell
     ! coordinates of first face-center
     xnb((icell-1)*6+1) = xyz(cell_shift+nb1_shift,1)
     ynb((icell-1)*6+1) = xyz(cell_shift+nb1_shift,2)
     znb((icell-1)*6+1) = xyz(cell_shift+nb1_shift,3)

     ! coordinates of second face-center
     xnb((icell-1)*6+2) = xyz(cell_shift+nb2_shift,1)
     ynb((icell-1)*6+2) = xyz(cell_shift+nb2_shift,2)
     znb((icell-1)*6+2) = xyz(cell_shift+nb2_shift,3)

     ! coordinates of third face-center
     xnb((icell-1)*6+3) = xyz(cell_shift+nb3_shift,1)
     ynb((icell-1)*6+3) = xyz(cell_shift+nb3_shift,2)
     znb((icell-1)*6+3) = xyz(cell_shift+nb3_shift,3)

     ! coordinates of fourth face-center
     xnb((icell-1)*6+4) = xyz(cell_shift+nb4_shift,1)
     ynb((icell-1)*6+4) = xyz(cell_shift+nb4_shift,2)
     znb((icell-1)*6+4) = xyz(cell_shift+nb4_shift,3)

     ! coordinates of fifth face-center
     xnb((icell-1)*6+5) = xyz(cell_shift+nb5_shift,1)
     ynb((icell-1)*6+5) = xyz(cell_shift+nb5_shift,2)
     znb((icell-1)*6+5) = xyz(cell_shift+nb5_shift,3)

     ! coordinates of sixth face-center
     xnb((icell-1)*6+6) = xyz(cell_shift+nb6_shift,1)
     ynb((icell-1)*6+6) = xyz(cell_shift+nb6_shift,2)
     znb((icell-1)*6+6) = xyz(cell_shift+nb6_shift,3)

     ! go to next NGLLX*NGLLY*NGLLZ points in array xyz, which form the next cell
     cell_shift = cell_shift + np_cell
  end do ! icell


  ! again loop on all cells and check for all faces if any other face-cell has the same coordinates (i.e. the faces are a match)

  do icell = 1,ncell
     ! define threashold for this cell to test distances to other face-centers as the minimum over all dimensions of
     ! face width per number of GLL points (multiplied scaled by 1e-3)
     v(1) = xnb((icell-1)*6+1)-xnb((icell-1)*6+3) ! vector v contains xmin-face-center minus xmax-face-center
     v(2) = ynb((icell-1)*6+1)-ynb((icell-1)*6+3)
     v(3) = znb((icell-1)*6+1)-znb((icell-1)*6+3)
     eps = sqrt( v(1)*v(1) + v(2)*v(2) + v(3)*v(3) )/real(NGLLX)
     v(1) = xnb((icell-1)*6+2)-xnb((icell-1)*6+4) ! vector v contains ymin-face-center minus ymax-face-center
     v(2) = ynb((icell-1)*6+2)-ynb((icell-1)*6+4)
     v(3) = znb((icell-1)*6+2)-znb((icell-1)*6+4)
     eps = min(eps, sqrt( v(1)*v(1) + v(2)*v(2) + v(3)*v(3) )/real(NGLLY) )
     v(1) = xnb((icell-1)*6+5)-xnb((icell-1)*6+6) ! vector v contains zmin-face-center minus zmax-face-center
     v(2) = ynb((icell-1)*6+5)-ynb((icell-1)*6+6)
     v(3) = znb((icell-1)*6+5)-znb((icell-1)*6+6)
     eps = min(eps, sqrt( v(1)*v(1) + v(2)*v(2) + v(3)*v(3) )/real(NGLLZ) )
     eps = eps * 0.001

     ! loop on all faces of cell icell
     do iface = 1,6
        ! v contains current face-center
        v(1) = xnb((icell-1)*6+iface)
        v(2) = ynb((icell-1)*6+iface)
        v(3) = znb((icell-1)*6+iface)

        ! loop on all other cells now and compare all face centers of those faces to current face center v
        do jcell = 1,ncell
           if(jcell==icell) cycle ! don't compare to own face-centers

           nnb = count( sqrt( (xnb((jcell-1)*6+1:jcell*6)-v(1))**2 + &
                (ynb((jcell-1)*6+1:jcell*6)-v(2))**2 + &
                (znb((jcell-1)*6+1:jcell*6)-v(3))**2) < eps)

           if(nnb>1) then
              write(*,*) "ASKI ERROR: for ",iface,"'th face of cell ",icell," there are ",nnb,&
                   " neighbouring faces in cell ",jcell
              call exit_MPI_without_rank("ASKI ERROR finding neighbours for specfem3dInversionGrid: "//&
                   "more than 1 face neighbours for one face")
           end if
           if(nnb==1) then
              ! increase the count of found neighbours
              neighbours(1,icell) = neighbours(1,icell) + 1
              ! store cell index of this found neighbour
              neighbours(neighbours(1,icell)+1,icell) = jcell
              exit ! loop over jcell, goto next face
           end if ! nnb==1

        end do ! jcell
     end do ! iface
  end do ! icell

end subroutine find_ASKI_neighbours_type_invgrid_4
!
!-------------------------------------------------------------------------------------------
!
subroutine write_ASKI_output()

  use constants,only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ
  use specfem_par,only: it,it_begin,it_end,NSTEP,DT,TWO_PI,ibool, &
       xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
       hprime_xx,hprime_yy,hprime_zz       
  use specfem_par_elastic,only: veloc,displ
  use specfem_par_ASKI

  implicit none

  integer :: ip,ispec,i,j,k,l,iglob,jf
  real(kind=CUSTOM_REAL) :: xix_ip,xiy_ip,xiz_ip, &
       etax_ip,etay_ip,etaz_ip, &
       gammax_ip,gammay_ip,gammaz_ip
  real(kind=CUSTOM_REAL) :: sumxxi,sumyxi,sumzxi, &
       sumxeta,sumyeta,sumzeta, &
       sumxgamma,sumygamma,sumzgamma
  real(kind=CUSTOM_REAL) :: ux,uy,uz,&
       uxdx,uxdy,uxdz, &
       uydx,uydy,uydz, &
       uzdx,uzdy,uzdz, &
       e11,e22,e33,e23,e13,e12
  double precision :: taper_value,cos_phi,sin_phi
  complex(kind=kind(1.d0)) :: efactor_tapered !variable is also used in case there is no tapering!
  double precision, dimension(:,:,:), pointer :: ASKI_Goertzel_U_tmp_double
  real, dimension(:,:,:), pointer :: ASKI_Goertzel_U_tmp_single

  if(.not.COMPUTE_ASKI_OUTPUT) return

  ! it_begin = 1 and it_end = NSTEP are assumed in this subroutine, as well as for deconvolving the source time 
  ! function below, so check this here (not sure wheter at all there would be SPECFEM standard functionality
  ! which uses it_begin /= 1 or it_end /= NSTEP ?!)
  if(it_begin /= 1 .or. it_end /= NSTEP) call exit_MPI_without_rank('for ASKI output it is assumed that shared '//&
       'parameters it_begin == 1 and it_end == NSTEP (at least one is violated): someone messed around with the '//&
       'time loop in subroutine iterate_time() or we are not running a forward simulation (?)')

  ! fourier transform to spectra only, if there are local wavefield points
  if(ASKI_np_local > 0) then

     ! reduce array access by reading from array ASKI_DFT_taper_values only once before the following loop
     if(ASKI_DFT_apply_taper .and. it >= ASKI_DFT_ntaper_start) taper_value = ASKI_DFT_taper_values(it)

     do ip=1,ASKI_np_local

        ! set ispec,i,j,k
        ispec = ASKI_indx_local(ip,1)
        i = ASKI_indx_local(ip,2)
        j = ASKI_indx_local(ip,3)
        k = ASKI_indx_local(ip,4)

        ! calculate derivatives

        ! set derivatives of xi,eta,gamma
        xix_ip = xix(i,j,k,ispec)
        xiy_ip = xiy(i,j,k,ispec)
        xiz_ip = xiz(i,j,k,ispec)
        etax_ip = etax(i,j,k,ispec)
        etay_ip = etay(i,j,k,ispec)
        etaz_ip = etaz(i,j,k,ispec)
        gammax_ip = gammax(i,j,k,ispec)
        gammay_ip = gammay(i,j,k,ispec)
        gammaz_ip = gammaz(i,j,k,ispec)

        ! calculate sum with h'_i(xi_l)=hprime_xx(i,l)
        sumxxi = 0.0
        sumyxi = 0.0
        sumzxi = 0.0

        do l=1,NGLLX

           iglob = ibool(l,j,k,ispec)
           if(ASKI_store_veloc) then
              sumxxi = sumxxi + veloc(1,iglob)*hprime_xx(i,l)
              sumyxi = sumyxi + veloc(2,iglob)*hprime_xx(i,l)
              sumzxi = sumzxi + veloc(3,iglob)*hprime_xx(i,l)
           else
              sumxxi = sumxxi + displ(1,iglob)*hprime_xx(i,l)
              sumyxi = sumyxi + displ(2,iglob)*hprime_xx(i,l)
              sumzxi = sumzxi + displ(3,iglob)*hprime_xx(i,l)
           end if

        end do ! l

        ! calculate sum with h'_j(eta_l)=hprime_yy(j,l)
        sumxeta = 0.0
        sumyeta = 0.0
        sumzeta = 0.0

        do l=1,NGLLY

           iglob = ibool(i,l,k,ispec)
           if(ASKI_store_veloc) then
              sumxeta = sumxeta + veloc(1,iglob)*hprime_yy(j,l)
              sumyeta = sumyeta + veloc(2,iglob)*hprime_yy(j,l)
              sumzeta = sumzeta + veloc(3,iglob)*hprime_yy(j,l)
           else
              sumxeta = sumxeta + displ(1,iglob)*hprime_yy(j,l)
              sumyeta = sumyeta + displ(2,iglob)*hprime_yy(j,l)
              sumzeta = sumzeta + displ(3,iglob)*hprime_yy(j,l)
           end if

        end do ! l

        ! calculate sum with h'_k(gamma_l)=hprime_zz(k,l)
        sumxgamma = 0.0
        sumygamma = 0.0
        sumzgamma = 0.0

        do l=1,NGLLZ

           iglob = ibool(i,j,l,ispec)
           if(ASKI_store_veloc) then
              sumxgamma = sumxgamma + veloc(1,iglob)*hprime_zz(k,l)
              sumygamma = sumygamma + veloc(2,iglob)*hprime_zz(k,l)
              sumzgamma = sumzgamma + veloc(3,iglob)*hprime_zz(k,l)
           else
              sumxgamma = sumxgamma + displ(1,iglob)*hprime_zz(k,l)
              sumygamma = sumygamma + displ(2,iglob)*hprime_zz(k,l)
              sumzgamma = sumzgamma + displ(3,iglob)*hprime_zz(k,l)
           end if

        end do ! l

        ! now calculate the derivative of veloc (displ) w.r.t. x, y and z with help of the sums calculated above
        ! also call it u if veloc is stored, since in the context of ASKI, this velocity field w.r.t Heaviside excitation is 
        ! interpreted as the DISPLACEMENT field w.r.t Delta-impulse excitation!

        ! derivative by x
        uxdx = xix_ip*sumxxi + etax_ip*sumxeta + gammax_ip*sumxgamma
        uydx = xix_ip*sumyxi + etax_ip*sumyeta + gammax_ip*sumygamma
        uzdx = xix_ip*sumzxi + etax_ip*sumzeta + gammax_ip*sumzgamma

        ! derivative by y
        uxdy = xiy_ip*sumxxi + etay_ip*sumxeta + gammay_ip*sumxgamma
        uydy = xiy_ip*sumyxi + etay_ip*sumyeta + gammay_ip*sumygamma
        uzdy = xiy_ip*sumzxi + etay_ip*sumzeta + gammay_ip*sumzgamma

        ! derivative by z
        uxdz = xiz_ip*sumxxi + etaz_ip*sumxeta + gammaz_ip*sumxgamma
        uydz = xiz_ip*sumyxi + etaz_ip*sumyeta + gammaz_ip*sumygamma
        uzdz = xiz_ip*sumzxi + etaz_ip*sumzeta + gammaz_ip*sumzgamma

        ! store underived velocity wavefield
        ! call it u, since in the context of ASKI, this velocity field w.r.t Heaviside excitation is 
        ! interpreted as the DISPLACEMENT field w.r.t Delta-impulse excitation!
        iglob = ibool(i,j,k,ispec)
        if(ASKI_store_veloc) then
           ux = veloc(1,iglob)
           uy = veloc(2,iglob)
           uz = veloc(3,iglob)
        else
           ux = displ(1,iglob)
           uy = displ(2,iglob)
           uz = displ(3,iglob)
        end if

        ! strains
        e11 = uxdx
        e22 = uydy
        e33 = uzdz
        e23 = 0.5*(uydz + uzdy)
        e13 = 0.5*(uxdz + uzdx)
        e12 = 0.5*(uxdy + uydx)

        ! conduct DFT
        if(ASKI_DFT_double) then

           select case(ASKI_DFT_method)

           case('EXPLICIT_SUMMATION')
              do jf = 1,ASKI_nf
                 ! reduce array access by reading from array efactors_tapered only once and taper, if necessary
                 if(ASKI_DFT_apply_taper .and. it >= ASKI_DFT_ntaper_start) then
                    efactor_tapered = ASKI_efactors(jf,it)*taper_value
                 else
                    efactor_tapered = ASKI_efactors(jf,it)
                 end if

                 ! underived wavefiled spectra
                 ASKI_spectra_local_double(1,jf,ip) = ASKI_spectra_local_double(1,jf,ip) + ux*efactor_tapered
                 ASKI_spectra_local_double(2,jf,ip) = ASKI_spectra_local_double(2,jf,ip) + uy*efactor_tapered
                 ASKI_spectra_local_double(3,jf,ip) = ASKI_spectra_local_double(3,jf,ip) + uz*efactor_tapered

                 ! strain spectra in order e11,e22,e33,e23,e13,e12
                 ASKI_spectra_local_double(4,jf,ip) = ASKI_spectra_local_double(4,jf,ip) + e11*efactor_tapered
                 ASKI_spectra_local_double(5,jf,ip) = ASKI_spectra_local_double(5,jf,ip) + e22*efactor_tapered
                 ASKI_spectra_local_double(6,jf,ip) = ASKI_spectra_local_double(6,jf,ip) + e33*efactor_tapered
                 ASKI_spectra_local_double(7,jf,ip) = ASKI_spectra_local_double(7,jf,ip) + e23*efactor_tapered
                 ASKI_spectra_local_double(8,jf,ip) = ASKI_spectra_local_double(8,jf,ip) + e13*efactor_tapered
                 ASKI_spectra_local_double(9,jf,ip) = ASKI_spectra_local_double(9,jf,ip) + e12*efactor_tapered
              end do ! jf

           case('GOERTZEL_STANDARD')

              ! In Goertzel's algorithm, the very last time step is treated differently compared with the others
              if(it < NSTEP) then

                 if(ASKI_DFT_apply_taper .and. it >= ASKI_DFT_ntaper_start) then

                    do jf = 1,ASKI_nf
                       ! underived wavefiled spectra
                       ASKI_Goertzel_U0_local_double(1,jf,ip) = &
                            ux*taper_value + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(1,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(1,jf,ip)
                       ASKI_Goertzel_U0_local_double(2,jf,ip) = &
                            uy*taper_value + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(2,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(2,jf,ip)
                       ASKI_Goertzel_U0_local_double(3,jf,ip) = &
                            uz*taper_value + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(3,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(3,jf,ip)

                       ! strain spectra in order e11,e22,e33,e23,e13,e12
                       ASKI_Goertzel_U0_local_double(4,jf,ip) = &
                            e11*taper_value + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(4,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(4,jf,ip)
                       ASKI_Goertzel_U0_local_double(5,jf,ip) = &
                            e22*taper_value + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(5,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(5,jf,ip)
                       ASKI_Goertzel_U0_local_double(6,jf,ip) = &
                            e33*taper_value + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(6,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(6,jf,ip)
                       ASKI_Goertzel_U0_local_double(7,jf,ip) = &
                            e23*taper_value + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(7,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(7,jf,ip)
                       ASKI_Goertzel_U0_local_double(8,jf,ip) = &
                            e13*taper_value + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(8,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(8,jf,ip)
                       ASKI_Goertzel_U0_local_double(9,jf,ip) = &
                            e12*taper_value + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(9,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(9,jf,ip)
                    end do ! jf

                 else ! taper

                    do jf = 1,ASKI_nf
                       ! underived wavefiled spectra
                       ASKI_Goertzel_U0_local_double(1,jf,ip) = &
                            ux + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(1,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(1,jf,ip)
                       ASKI_Goertzel_U0_local_double(2,jf,ip) = &
                            uy + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(2,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(2,jf,ip)
                       ASKI_Goertzel_U0_local_double(3,jf,ip) = &
                            uz + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(3,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(3,jf,ip)

                       ! strain spectra in order e11,e22,e33,e23,e13,e12
                       ASKI_Goertzel_U0_local_double(4,jf,ip) = &
                            e11 + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(4,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(4,jf,ip)
                       ASKI_Goertzel_U0_local_double(5,jf,ip) = &
                            e22 + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(5,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(5,jf,ip)
                       ASKI_Goertzel_U0_local_double(6,jf,ip) = &
                            e33 + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(6,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(6,jf,ip)
                       ASKI_Goertzel_U0_local_double(7,jf,ip) = &
                            e23 + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(7,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(7,jf,ip)
                       ASKI_Goertzel_U0_local_double(8,jf,ip) = &
                            e13 + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(8,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(8,jf,ip)
                       ASKI_Goertzel_U0_local_double(9,jf,ip) = &
                            e12 + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(9,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(9,jf,ip)
                    end do ! jf

                 end if ! taper

              else ! it < NSTEP

                 ! AT THE LAST TIME STEP OF GOERTZEL'S ALGORITHM:
                 ! - COMPUTE THE REAL PART OF THE OUTPUT SPECTRA AND STORE TO U0 (rename below as U1)
                 ! - COMPUTE THE IMAGINARY PART OF THE OUTPUT SPECTRA AND STORE TO U2

                 if(ASKI_DFT_apply_taper .and. it >= ASKI_DFT_ntaper_start) then

                    do jf = 1,ASKI_nf
                       ! REAL PART

                       ! underived wavefiled spectra
                       ASKI_Goertzel_U0_local_double(1,jf,ip) = DT * ( &
                            ux*taper_value + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(1,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(1,jf,ip) )
                       ASKI_Goertzel_U0_local_double(2,jf,ip) = DT * ( &
                            uy*taper_value + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(2,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(2,jf,ip) )
                       ASKI_Goertzel_U0_local_double(3,jf,ip) = DT * ( &
                            uz*taper_value + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(3,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(3,jf,ip) )

                       ! strain spectra in order e11,e22,e33,e23,e13,e12
                       ASKI_Goertzel_U0_local_double(4,jf,ip) = DT * ( &
                            e11*taper_value + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(4,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(4,jf,ip) )
                       ASKI_Goertzel_U0_local_double(5,jf,ip) = DT * ( &
                            e22*taper_value + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(5,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(5,jf,ip) )
                       ASKI_Goertzel_U0_local_double(6,jf,ip) = DT * ( &
                            e33*taper_value + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(6,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(6,jf,ip) )
                       ASKI_Goertzel_U0_local_double(7,jf,ip) = DT * ( &
                            e23*taper_value + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(7,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(7,jf,ip) )
                       ASKI_Goertzel_U0_local_double(8,jf,ip) = DT * ( &
                            e13*taper_value + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(8,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(8,jf,ip) )
                       ASKI_Goertzel_U0_local_double(9,jf,ip) = DT * ( &
                            e12*taper_value + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(9,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(9,jf,ip) )

                       ! IMAGINARY PART

                       ! underived wavefiled spectra
                       ASKI_Goertzel_U2_local_double(1,jf,ip) = DT * ( &
                            ux*taper_value + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_double(1,jf,ip) )
                       ASKI_Goertzel_U2_local_double(2,jf,ip) = DT * ( &
                            uy*taper_value + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_double(2,jf,ip) )
                       ASKI_Goertzel_U2_local_double(3,jf,ip) = DT * ( &
                            uz*taper_value + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_double(3,jf,ip) )

                       ! strain spectra in order e11,e22,e33,e23,e13,e12
                       ASKI_Goertzel_U2_local_double(4,jf,ip) = DT * ( &
                            e11*taper_value + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_double(4,jf,ip) )
                       ASKI_Goertzel_U2_local_double(5,jf,ip) = DT * ( &
                            e22*taper_value + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_double(5,jf,ip) )
                       ASKI_Goertzel_U2_local_double(6,jf,ip) = DT * ( &
                            e33*taper_value + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_double(6,jf,ip) )
                       ASKI_Goertzel_U2_local_double(7,jf,ip) = DT * ( &
                            e23*taper_value + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_double(7,jf,ip) )
                       ASKI_Goertzel_U2_local_double(8,jf,ip) = DT * ( &
                            e13*taper_value + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_double(8,jf,ip) )
                       ASKI_Goertzel_U2_local_double(9,jf,ip) = DT * ( &
                            e12*taper_value + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_double(9,jf,ip) )
                    end do ! jf

                 else ! taper

                    do jf = 1,ASKI_nf
                       ! REAL PART

                       ! underived wavefiled spectra
                       ASKI_Goertzel_U0_local_double(1,jf,ip) = DT * ( &
                            ux + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(1,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(1,jf,ip) )
                       ASKI_Goertzel_U0_local_double(2,jf,ip) = DT * ( &
                            uy + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(2,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(2,jf,ip) )
                       ASKI_Goertzel_U0_local_double(3,jf,ip) = DT * ( &
                            uz + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(3,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(3,jf,ip) )

                       ! strain spectra in order e11,e22,e33,e23,e13,e12
                       ASKI_Goertzel_U0_local_double(4,jf,ip) = DT * ( &
                            e11 + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(4,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(4,jf,ip) )
                       ASKI_Goertzel_U0_local_double(5,jf,ip) = DT * ( &
                            e22 + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(5,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(5,jf,ip) )
                       ASKI_Goertzel_U0_local_double(6,jf,ip) = DT * ( &
                            e33 + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(6,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(6,jf,ip) )
                       ASKI_Goertzel_U0_local_double(7,jf,ip) = DT * ( &
                            e23 + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(7,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(7,jf,ip) )
                       ASKI_Goertzel_U0_local_double(8,jf,ip) = DT * ( &
                            e13 + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(8,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(8,jf,ip) )
                       ASKI_Goertzel_U0_local_double(9,jf,ip) = DT * ( &
                            e12 + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(9,jf,ip) &
                            - ASKI_Goertzel_U2_local_double(9,jf,ip) )

                       ! IMAGINARY PART

                       ! underived wavefiled spectra
                       ASKI_Goertzel_U2_local_double(1,jf,ip) = DT * ( &
                            ux + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_double(1,jf,ip) )
                       ASKI_Goertzel_U2_local_double(2,jf,ip) = DT * ( &
                            uy + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_double(2,jf,ip) )
                       ASKI_Goertzel_U2_local_double(3,jf,ip) = DT * ( &
                            uz + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_double(3,jf,ip) )

                       ! strain spectra in order e11,e22,e33,e23,e13,e12
                       ASKI_Goertzel_U2_local_double(4,jf,ip) = DT * ( &
                            e11 + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_double(4,jf,ip) )
                       ASKI_Goertzel_U2_local_double(5,jf,ip) = DT * ( &
                            e22 + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_double(5,jf,ip) )
                       ASKI_Goertzel_U2_local_double(6,jf,ip) = DT * ( &
                            e33 + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_double(6,jf,ip) )
                       ASKI_Goertzel_U2_local_double(7,jf,ip) = DT * ( &
                            e23 + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_double(7,jf,ip) )
                       ASKI_Goertzel_U2_local_double(8,jf,ip) = DT * ( &
                            e13 + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_double(8,jf,ip) )
                       ASKI_Goertzel_U2_local_double(9,jf,ip) = DT * ( &
                            e12 + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_double(9,jf,ip) )
                    end do ! jf

                 end if ! taper

              end if ! it < NSTEP

           end select

        else ! ASKI_DFT_double

           select case(ASKI_DFT_method)

           case('EXPLICIT_SUMMATION')
              do jf = 1,ASKI_nf
                 ! reduce array access by reading from array efactors_tapered only once and taper, if necessary
                 if(ASKI_DFT_apply_taper .and. it >= ASKI_DFT_ntaper_start) then
                    efactor_tapered = ASKI_efactors(jf,it)*taper_value
                 else
                    efactor_tapered = ASKI_efactors(jf,it)
                 end if

                 ! underived wavefiled spectra
                 ASKI_spectra_local_single(1,jf,ip) = ASKI_spectra_local_single(1,jf,ip) + ux*efactor_tapered
                 ASKI_spectra_local_single(2,jf,ip) = ASKI_spectra_local_single(2,jf,ip) + uy*efactor_tapered
                 ASKI_spectra_local_single(3,jf,ip) = ASKI_spectra_local_single(3,jf,ip) + uz*efactor_tapered

                 ! strain spectra in order e11,e22,e33,e23,e13,e12
                 ASKI_spectra_local_single(4,jf,ip) = ASKI_spectra_local_single(4,jf,ip) + e11*efactor_tapered
                 ASKI_spectra_local_single(5,jf,ip) = ASKI_spectra_local_single(5,jf,ip) + e22*efactor_tapered
                 ASKI_spectra_local_single(6,jf,ip) = ASKI_spectra_local_single(6,jf,ip) + e33*efactor_tapered
                 ASKI_spectra_local_single(7,jf,ip) = ASKI_spectra_local_single(7,jf,ip) + e23*efactor_tapered
                 ASKI_spectra_local_single(8,jf,ip) = ASKI_spectra_local_single(8,jf,ip) + e13*efactor_tapered
                 ASKI_spectra_local_single(9,jf,ip) = ASKI_spectra_local_single(9,jf,ip) + e12*efactor_tapered
              end do ! jf

           case('GOERTZEL_STANDARD')

              ! In Goertzel's algorithm, the very last time step is treated differently compared with the others
              if(it < NSTEP) then

                 if(ASKI_DFT_apply_taper .and. it >= ASKI_DFT_ntaper_start) then

                    do jf = 1,ASKI_nf
                       ! underived wavefiled spectra
                       ASKI_Goertzel_U0_local_single(1,jf,ip) = &
                            ux*taper_value + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(1,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(1,jf,ip)
                       ASKI_Goertzel_U0_local_single(2,jf,ip) = &
                            uy*taper_value + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(2,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(2,jf,ip)
                       ASKI_Goertzel_U0_local_single(3,jf,ip) = &
                            uz*taper_value + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(3,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(3,jf,ip)

                       ! strain spectra in order e11,e22,e33,e23,e13,e12
                       ASKI_Goertzel_U0_local_single(4,jf,ip) = &
                            e11*taper_value + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(4,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(4,jf,ip)
                       ASKI_Goertzel_U0_local_single(5,jf,ip) = &
                            e22*taper_value + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(5,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(5,jf,ip)
                       ASKI_Goertzel_U0_local_single(6,jf,ip) = &
                            e33*taper_value + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(6,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(6,jf,ip)
                       ASKI_Goertzel_U0_local_single(7,jf,ip) = &
                            e23*taper_value + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(7,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(7,jf,ip)
                       ASKI_Goertzel_U0_local_single(8,jf,ip) = &
                            e13*taper_value + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(8,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(8,jf,ip)
                       ASKI_Goertzel_U0_local_single(9,jf,ip) = &
                            e12*taper_value + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(9,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(9,jf,ip)
                    end do ! jf

                 else ! taper

                    do jf = 1,ASKI_nf
                       ! underived wavefiled spectra
                       ASKI_Goertzel_U0_local_single(1,jf,ip) = &
                            ux + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(1,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(1,jf,ip)
                       ASKI_Goertzel_U0_local_single(2,jf,ip) = &
                            uy + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(2,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(2,jf,ip)
                       ASKI_Goertzel_U0_local_single(3,jf,ip) = &
                            uz + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(3,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(3,jf,ip)

                       ! strain spectra in order e11,e22,e33,e23,e13,e12
                       ASKI_Goertzel_U0_local_single(4,jf,ip) = &
                            e11 + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(4,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(4,jf,ip)
                       ASKI_Goertzel_U0_local_single(5,jf,ip) = &
                            e22 + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(5,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(5,jf,ip)
                       ASKI_Goertzel_U0_local_single(6,jf,ip) = &
                            e33 + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(6,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(6,jf,ip)
                       ASKI_Goertzel_U0_local_single(7,jf,ip) = &
                            e23 + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(7,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(7,jf,ip)
                       ASKI_Goertzel_U0_local_single(8,jf,ip) = &
                            e13 + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(8,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(8,jf,ip)
                       ASKI_Goertzel_U0_local_single(9,jf,ip) = &
                            e12 + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(9,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(9,jf,ip)
                    end do ! jf

                 end if ! taper

              else ! it < NSTEP

                 ! AT THE LAST TIME STEP OF GOERTZEL'S ALGORITHM:
                 ! - COMPUTE THE REAL PART OF THE OUTPUT SPECTRA AND STORE TO U0 (rename below as U1)
                 ! - COMPUTE THE IMAGINARY PART OF THE OUTPUT SPECTRA AND STORE TO U2

                 if(ASKI_DFT_apply_taper .and. it >= ASKI_DFT_ntaper_start) then

                    do jf = 1,ASKI_nf
                       ! REAL PART

                       ! underived wavefiled spectra
                       ASKI_Goertzel_U0_local_single(1,jf,ip) = DT * ( &
                            ux*taper_value + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(1,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(1,jf,ip) )
                       ASKI_Goertzel_U0_local_single(2,jf,ip) = DT * ( &
                            uy*taper_value + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(2,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(2,jf,ip) )
                       ASKI_Goertzel_U0_local_single(3,jf,ip) = DT * ( &
                            uz*taper_value + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(3,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(3,jf,ip) )

                       ! strain spectra in order e11,e22,e33,e23,e13,e12
                       ASKI_Goertzel_U0_local_single(4,jf,ip) = DT * ( &
                            e11*taper_value + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(4,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(4,jf,ip) )
                       ASKI_Goertzel_U0_local_single(5,jf,ip) = DT * ( &
                            e22*taper_value + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(5,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(5,jf,ip) )
                       ASKI_Goertzel_U0_local_single(6,jf,ip) = DT * ( &
                            e33*taper_value + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(6,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(6,jf,ip) )
                       ASKI_Goertzel_U0_local_single(7,jf,ip) = DT * ( &
                            e23*taper_value + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(7,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(7,jf,ip) )
                       ASKI_Goertzel_U0_local_single(8,jf,ip) = DT * ( &
                            e13*taper_value + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(8,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(8,jf,ip) )
                       ASKI_Goertzel_U0_local_single(9,jf,ip) = DT * ( &
                            e12*taper_value + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(9,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(9,jf,ip) )

                       ! IMAGINARY PART

                       ! underived wavefiled spectra
                       ASKI_Goertzel_U2_local_single(1,jf,ip) = DT * ( &
                            ux*taper_value + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_single(1,jf,ip) )
                       ASKI_Goertzel_U2_local_single(2,jf,ip) = DT * ( &
                            uy*taper_value + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_single(2,jf,ip) )
                       ASKI_Goertzel_U2_local_single(3,jf,ip) = DT * ( &
                            uz*taper_value + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_single(3,jf,ip) )

                       ! strain spectra in order e11,e22,e33,e23,e13,e12
                       ASKI_Goertzel_U2_local_single(4,jf,ip) = DT * ( &
                            e11*taper_value + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_single(4,jf,ip) )
                       ASKI_Goertzel_U2_local_single(5,jf,ip) = DT * ( &
                            e22*taper_value + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_single(5,jf,ip) )
                       ASKI_Goertzel_U2_local_single(6,jf,ip) = DT * ( &
                            e33*taper_value + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_single(6,jf,ip) )
                       ASKI_Goertzel_U2_local_single(7,jf,ip) = DT * ( &
                            e23*taper_value + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_single(7,jf,ip) )
                       ASKI_Goertzel_U2_local_single(8,jf,ip) = DT * ( &
                            e13*taper_value + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_single(8,jf,ip) )
                       ASKI_Goertzel_U2_local_single(9,jf,ip) = DT * ( &
                            e12*taper_value + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_single(9,jf,ip) )
                    end do ! jf

                 else ! taper

                    do jf = 1,ASKI_nf
                       ! REAL PART

                       ! underived wavefiled spectra
                       ASKI_Goertzel_U0_local_single(1,jf,ip) = DT * ( &
                            ux + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(1,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(1,jf,ip) )
                       ASKI_Goertzel_U0_local_single(2,jf,ip) = DT * ( &
                            uy + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(2,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(2,jf,ip) )
                       ASKI_Goertzel_U0_local_single(3,jf,ip) = DT * ( &
                            uz + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(3,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(3,jf,ip) )

                       ! strain spectra in order e11,e22,e33,e23,e13,e12
                       ASKI_Goertzel_U0_local_single(4,jf,ip) = DT * ( &
                            e11 + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(4,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(4,jf,ip) )
                       ASKI_Goertzel_U0_local_single(5,jf,ip) = DT * ( &
                            e22 + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(5,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(5,jf,ip) )
                       ASKI_Goertzel_U0_local_single(6,jf,ip) = DT * ( &
                            e33 + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(6,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(6,jf,ip) )
                       ASKI_Goertzel_U0_local_single(7,jf,ip) = DT * ( &
                            e23 + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(7,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(7,jf,ip) )
                       ASKI_Goertzel_U0_local_single(8,jf,ip) = DT * ( &
                            e13 + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(8,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(8,jf,ip) )
                       ASKI_Goertzel_U0_local_single(9,jf,ip) = DT * ( &
                            e12 + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(9,jf,ip) &
                            - ASKI_Goertzel_U2_local_single(9,jf,ip) )

                       ! IMAGINARY PART

                       ! underived wavefiled spectra
                       ASKI_Goertzel_U2_local_single(1,jf,ip) = DT * ( &
                            ux + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_single(1,jf,ip) )
                       ASKI_Goertzel_U2_local_single(2,jf,ip) = DT * ( &
                            uy + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_single(2,jf,ip) )
                       ASKI_Goertzel_U2_local_single(3,jf,ip) = DT * ( &
                            uz + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_single(3,jf,ip) )

                       ! strain spectra in order e11,e22,e33,e23,e13,e12
                       ASKI_Goertzel_U2_local_single(4,jf,ip) = DT * ( &
                            e11 + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_single(4,jf,ip) )
                       ASKI_Goertzel_U2_local_single(5,jf,ip) = DT * ( &
                            e22 + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_single(5,jf,ip) )
                       ASKI_Goertzel_U2_local_single(6,jf,ip) = DT * ( &
                            e33 + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_single(6,jf,ip) )
                       ASKI_Goertzel_U2_local_single(7,jf,ip) = DT * ( &
                            e23 + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_single(7,jf,ip) )
                       ASKI_Goertzel_U2_local_single(8,jf,ip) = DT * ( &
                            e13 + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_single(8,jf,ip) )
                       ASKI_Goertzel_U2_local_single(9,jf,ip) = DT * ( &
                            e12 + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_single(9,jf,ip) )
                    end do ! jf

                 end if ! taper

              end if ! it < NSTEP

           end select

        end if ! ASKI_DFT_double

     end do ! ip

     select case(ASKI_DFT_method)
     case('GOERTZEL_STANDARD')
        if(ASKI_DFT_double) then
           if(it < NSTEP) then
              ! rename U2 = U1 and U1 = U0 by re-assigning the pointers of all three arrays accordingly
              ASKI_Goertzel_U_tmp_double => ASKI_Goertzel_U2_local_double
              ASKI_Goertzel_U2_local_double => ASKI_Goertzel_U1_local_double
              ASKI_Goertzel_U1_local_double => ASKI_Goertzel_U0_local_double
              ! use the memory of U2 (pointed to by U_tmp) for setting the values U0 in next time step
              ASKI_Goertzel_U0_local_double => ASKI_Goertzel_U_tmp_double
           else ! it < NSTEP
              ! Correct for phase shift in time-reversed goertzel algorithm by multiplying the spectral values
              ! by exp(-i*omega*(NSTEP-1)*DT) = cos(-omega*(NSTEP-1)*DT) + i*sin(-omega*(NSTEP-1)*DT)
              ! Equivalently, modify real and imaginary part of the spectral value explicitely.
              ! Moreover: for efficiency reasons, above the real part was stored to U0, but below it is assumed 
              ! to be contained in U1. Hence, move real part values here to U1.
              do jf = 1,ASKI_nf
                 cos_phi = cos(-TWO_PI*ASKI_jf(jf)*ASKI_df*(NSTEP-1)*DT)
                 sin_phi = sin(-TWO_PI*ASKI_jf(jf)*ASKI_df*(NSTEP-1)*DT)
                 ! REAL PART
                 ASKI_Goertzel_U1_local_double(:,jf,:) = &
                      ASKI_Goertzel_U0_local_double(:,jf,:)*cos_phi - ASKI_Goertzel_U2_local_double(:,jf,:)*sin_phi
                 ! IMAGINARY PART
                 ASKI_Goertzel_U2_local_double(:,jf,:) = &
                      ASKI_Goertzel_U2_local_double(:,jf,:)*cos_phi + ASKI_Goertzel_U0_local_double(:,jf,:)*sin_phi
              end do ! jf
              ! U0 is not needed anymore, so deallocate (as early as possible)
              deallocate(ASKI_Goertzel_U0_local_double)
              nullify(ASKI_Goertzel_U0_local_double)
           end if ! it < NSTEP
        else ! ASKI_DFT_double
           if(it < NSTEP) then
              ! rename U2 = U1 and U1 = U0 by re-assigning the pointers of all three arrays accordingly
              ASKI_Goertzel_U_tmp_single => ASKI_Goertzel_U2_local_single
              ASKI_Goertzel_U2_local_single => ASKI_Goertzel_U1_local_single
              ASKI_Goertzel_U1_local_single => ASKI_Goertzel_U0_local_single
              ! the third array (U0) is not needed anymore, so deallocate (as early as possible)
              ASKI_Goertzel_U0_local_single => ASKI_Goertzel_U_tmp_single
           else ! it < NSTEP
              ! Correct for phase shift in time-reversed goertzel algorithm by multiplying the spectral values
              ! by exp(-i*omega*(NSTEP-1)*DT) = cos(-omega*(NSTEP-1)*DT) + i*sin(-omega*(NSTEP-1)*DT)
              ! Equivalently, modify real and imaginary part of the spectral value explicitely.
              ! Moreover: for efficiency reasons, above the real part was stored to U0, but below it is assumed 
              ! to be contained in U1. Hence, move real part values here to U1.
              do jf = 1,ASKI_nf
                 cos_phi = cos(-TWO_PI*ASKI_jf(jf)*ASKI_df*(NSTEP-1)*DT)
                 sin_phi = sin(-TWO_PI*ASKI_jf(jf)*ASKI_df*(NSTEP-1)*DT)
                 ! REAL PART
                 ASKI_Goertzel_U1_local_single(:,jf,:) = &
                      ASKI_Goertzel_U0_local_single(:,jf,:)*cos_phi - ASKI_Goertzel_U2_local_single(:,jf,:)*sin_phi
                 ! IMAGINARY PART
                 ASKI_Goertzel_U2_local_single(:,jf,:) = &
                      ASKI_Goertzel_U2_local_single(:,jf,:)*cos_phi + ASKI_Goertzel_U0_local_single(:,jf,:)*sin_phi
              end do ! jf
              ! U0 is not needed anymore, so deallocate (as early as possible)
              deallocate(ASKI_Goertzel_U0_local_single)
              nullify(ASKI_Goertzel_U0_local_single)
           end if ! it < NSTEP
        end if ! ASKI_DFT_double
     end select

  end if ! ASKI_np_local > 0

end subroutine write_ASKI_output
!
!-------------------------------------------------------------------------------------------
!
subroutine save_ASKI_output()
  use specfem_par_ASKI

  implicit none

  if(.not.COMPUTE_ASKI_OUTPUT) return

  if(ASKI_DECONVOLVE_STF) call deconvolve_stf_from_ASKI_output()

  call write_ASKI_output_files()

  ! deallocate everything, simulation is over
  if(allocated(ASKI_np_local_all)) deallocate(ASKI_np_local_all)
  if(allocated(ASKI_indx_local)) deallocate(ASKI_indx_local)
  if(allocated(ASKI_efactors)) deallocate(ASKI_efactors)
  if(allocated(ASKI_Goertzel_Wr)) deallocate(ASKI_Goertzel_Wr)
  if(allocated(ASKI_Goertzel_Wi)) deallocate(ASKI_Goertzel_Wi)
  if(allocated(ASKI_DFT_taper_values)) deallocate(ASKI_DFT_taper_values)
  if(allocated(ASKI_stf_spectrum_double)) deallocate(ASKI_stf_spectrum_double)
  if(allocated(ASKI_spectra_local_double)) deallocate(ASKI_spectra_local_double)
  if(allocated(ASKI_spectra_local_single)) deallocate(ASKI_spectra_local_single)
  if(associated(ASKI_Goertzel_U0_local_double)) deallocate(ASKI_Goertzel_U0_local_double)
  if(associated(ASKI_Goertzel_U1_local_double)) deallocate(ASKI_Goertzel_U1_local_double)
  if(associated(ASKI_Goertzel_U2_local_double)) deallocate(ASKI_Goertzel_U2_local_double)
  if(associated(ASKI_Goertzel_U0_local_single)) deallocate(ASKI_Goertzel_U0_local_single)
  if(associated(ASKI_Goertzel_U1_local_single)) deallocate(ASKI_Goertzel_U1_local_single)
  if(associated(ASKI_Goertzel_U2_local_single)) deallocate(ASKI_Goertzel_U2_local_single)
  if(allocated(ASKI_jf)) deallocate(ASKI_jf)
end subroutine save_ASKI_output
!
!-------------------------------------------------------------------------------------------
!
subroutine deconvolve_stf_from_ASKI_output()
  use constants,only: OUTPUT_FILES_BASE
  use specfem_par,only: myrank,NSTEP,DT,TWO_PI
  use specfem_par_ASKI
  implicit none
  double precision, external :: comp_source_time_function
  double precision, dimension(:), allocatable :: stf_deconv
  double precision, dimension(:), pointer :: U0,U1,U2,U2_tmp
  double precision :: dummy_time,stf_value_left,stf_value_tmp
  integer :: jt,jf,IOASKI,ios
  character(len=41) :: filename_stf_deconv,filename_stf_spectrum_dconv

  nullify(U0,U1,U2,U2_tmp)

  allocate(stf_deconv(NSTEP))

  ! make this routine independent of the knowledge about the chosen stf 
  ! by reading in output OUTPUT_FILES/plot_source_time_function.txt (assuming it is produced)
  ! THIS WAY, THE CODE SHOULD WORK EVEN IN CASE THE CONVENTION ABOUT GAUSSIAN/ERROR FUNCTION
  ! HAS CHANGED!
  if(myrank==0) then
     call get_file_unit_ASKI(IOASKI)
     open(unit=IOASKI,file=trim(OUTPUT_FILES_BASE)//'plot_source_time_function.txt',&
          form='formatted',status='old',action='read',iostat=ios)
     if(ios/=0) call exit_MPI(myrank,"could not open source time function file '"//trim(OUTPUT_FILES_BASE)//&
          'plot_source_time_function.txt'//"' to read")
     do jt = 1,NSTEP
        read(IOASKI,*) dummy_time,stf_deconv(jt)
     end do ! jt
     close(IOASKI)
  end if

  ! Since it is a bit more complicated to broadcast the stf only to those ranks which hold any ASKI output,
  ! let all ranks enter this routine to this point here and broadcast to all.
  call bcast_all_dp(stf_deconv,NSTEP)

  ! Those ranks which do not hold any ASKI output can leave this routine now, EXCEPT rank 0 which writes the logs below
  if(ASKI_np_local <= 0 .and. myrank/=0) then
     deallocate(stf_deconv)
     return
  end if

  ! if velocity was stored above, differentiate source time function by central finite differences
  if(ASKI_store_veloc) then
     stf_value_left = stf_deconv(1)
     do jt = 2,NSTEP-1
        stf_value_tmp = stf_deconv(jt) ! need to memorize the original value (will be the left one for the next step)
        stf_deconv(jt) = (stf_deconv(jt+1)-stf_value_left)/(2.d0*DT)
        stf_value_left = stf_value_tmp
     end do
     ! make the derivative continuous at the beginning and the end
     ! however, when using this properly, it should be zero anyway AND ADDITIONALLY it should be tapered to zero at the end
     stf_deconv(1) = stf_deconv(2)
     stf_deconv(NSTEP) = stf_deconv(NSTEP-1)
  end if ! ASKI_store_veloc

  ! now transform stf_deconv to spectrum at the chosen discrete frequencies

  allocate(ASKI_stf_spectrum_double(ASKI_nf))
  select case(ASKI_DFT_method)
  case('EXPLICIT_SUMMATION')
!!$     ASKI_stf_spectrum_double = (0.d0,0.d0)
!!$     do jt = 1,NSTEP
!!$        do jf = 1,ASKI_nf
!!$           ASKI_stf_spectrum_double(jf) = ASKI_stf_spectrum_double(jf) + stf_deconv(jt)*ASKI_efactors(jf,jt)
!!$        end do ! jf
!!$     end do ! jt
     ASKI_stf_spectrum_double = matmul(ASKI_efactors,stf_deconv)
  case('GOERTZEL_STANDARD')
     allocate(U0(ASKI_nf),U1(ASKI_nf),U2(ASKI_nf))
     U1 = 0.d0
     U2 = 0.d0
     do jt = 1,NSTEP-1
        do jf = 1,ASKI_nf
           U0(jf) = stf_deconv(jt) + ASKI_Goertzel_Wr(jf)*U1(jf) - U2(jf)
        end do ! jf
        ! rename U2 = U1 and U1 = U0 by re-assigning the pointers of all three arrays accordingly
        U2_tmp => U2
        U2 => U1
        U1 => U0
        U0 => U2_tmp ! use the memory of U2 (pointed to by U2_tmp) for setting the values U0 in next time step
     end do ! jt
     ! In Goertzel's algorithm, the very last time step is treated differently compared with the others
     do jf = 1,ASKI_nf
        ASKI_stf_spectrum_double(jf) = cmplx( &
             DT * ( stf_deconv(NSTEP)+0.5d0*ASKI_Goertzel_Wr(jf)*U1(jf)-U2(jf) ) , &
             DT * ASKI_Goertzel_Wi(jf) * U1(jf) , kind=kind(1.d0) )
        ! After treating the very last time step, correct for phase shift in time-reversed Goertzel algorithm
        ! by multiplying the spectral values by exp(-i*omega*(NSTEP-1)*DT)
        ASKI_stf_spectrum_double(jf) = ASKI_stf_spectrum_double(jf) * exp(-(0.d0,1.d0)*TWO_PI*ASKI_jf(jf)*ASKI_df*(NSTEP-1)*DT)
     end do
     nullify(U2_tmp)
     deallocate(U0,U1,U2)
  end select

  ! RANK 0 WRITES OUT LOGS CONTAINING THE (DIFFERENTIATED) STF AND THE SPECTRUM WHICH IS DECONVOLVED
  if(myrank==0) then
     if(ASKI_store_veloc) then
        filename_stf_deconv = 'LOG_ASKI_DECONVOLVE_stf_diff.dat'
        filename_stf_spectrum_dconv = 'LOG_ASKI_DECONVOLVE_stf_diff_spectrum.dat'
        write(*,*) "WILL DECONVOLVE DIFFERENTIATED SOURCE-TIME-FUNCTION FROM ASKI KERNEL SPECTRA (velocity field was stored). ",&
             "LOGS CONTAINING THIS TIME-SERIES AND ITS SPECTRUM ARE WRITTEN NOW TO OUTPUT_FILES/"
        call write_ASKI_log("LOG_ASKI_DECONVOLVE_stf.txt",&
             "In Par_file_ASKI: ASKI_DECONVOLVE_STF is .true.; so source-time-function is "//&
             "deconvolved from spectral ASKI kernel output. A quasi-Heaviside stf should have "//&
             "been used by SPECFEM for this moment tensor source. For reasons of numerical stability: "//&
             "spectra of particle velocity were stored as kernel output, which now will be "//&
             "deconvolved by a gaussian (i.e. the differentiated stf which has spectrum close "//&
             "to 1.0, so can be deconvolved in a numerically stable way). Final kernel output, "//&
             "hence, will be spectra of particle displacement w.r.t. an actual dirac stf! Files "//&
             "'LOG_ASKI_DECONVOLVE_stf_diff.dat', 'LOG_ASKI_DECONVOLVE_stf_diff_spectrum.dat' "//&
             "contain the deconvolved stf and its spectrum at the ASKI frequencies")
     else ! ASKI_store_veloc
        filename_stf_deconv = 'LOG_ASKI_DECONVOLVE_stf.dat'
        filename_stf_spectrum_dconv = 'LOG_ASKI_DECONVOLVE_stf_spectrum.dat'
        write(*,*) "WILL DECONVOLVE SOURCE-TIME-FUNCTION FROM ASKI KERNEL SPECTRA (displacement field was stored). ",&
             "LOGS CONTAINING THIS TIME-SERIES AND ITS SPECTRUM ARE WRITTEN NOW TO OUTPUT_FILES/"
        call write_ASKI_log("LOG_ASKI_DECONVOLVE_stf.txt",&
             "In Par_file_ASKI: ASKI_DECONVOLVE_STF is .true.; so source-time-function is "//&
             "deconvolved from spectral ASKI kernel output. A thin Gaussian stf should have "//&
             "been used by SPECFEM for this single force source. "//&
             "Spectra of particle displacement were stored as kernel output, which now will be "//&
             "deconvolved by the gaussian stf. Final kernel output, "//&
             "hence, will be spectra of particle displacement w.r.t. an actual dirac stf! Files "//&
             "'LOG_ASKI_DECONVOLVE_stf.dat', 'LOG_ASKI_DECONVOLVE_stf_spectrum.dat' "//&
             "contain the deconvolved stf and its spectrum at the ASKI frequencies")
     end if ! ASKI_store_veloc

     call get_file_unit_ASKI(IOASKI)
     open(unit=IOASKI,file=trim(OUTPUT_FILES_BASE)//trim(filename_stf_deconv),&
          form='formatted',status='unknown',action='write')
     do jt = 1,NSTEP
        write(IOASKI,*) real(dble(jt-1)*DT),real(stf_deconv(jt))
     end do ! jt
     close(IOASKI)
     deallocate(stf_deconv)
     call get_file_unit_ASKI(IOASKI)
     open(unit=IOASKI,file=trim(OUTPUT_FILES_BASE)//trim(filename_stf_spectrum_dconv),&
          form='formatted',status='unknown',action='write')
     do jf = 1,ASKI_nf
        write(IOASKI,*) ASKI_jf(jf)*ASKI_df, ASKI_stf_spectrum_double(jf), &
             atan2(aimag(ASKI_stf_spectrum_double(jf)),real(ASKI_stf_spectrum_double(jf))), abs(ASKI_stf_spectrum_double(jf))
     end do ! jf
     close(IOASKI)

     ! if myrank==0 is only here to write the deconvolve log output, but does not have any aski output, 
     ! then it can deallocate the deconvolve-spectrum already here
     if(ASKI_np_local <= 0) deallocate(ASKI_stf_spectrum_double)
  end if ! myrank == 0

end subroutine deconvolve_stf_from_ASKI_output
!
!-------------------------------------------------------------------------------------------
!
subroutine write_ASKI_output_files()

  use specfem_par,only: myrank,NPROC
  use specfem_par_ASKI

  implicit none

  integer :: jf,ip,iproc,IOASKI,specfem_version
  complex, dimension(:,:), allocatable :: spectrum_one_frequency,spectrum_one_frequency_local
  character(len=509) :: filename
  logical :: file_exists
  character (len=7) :: open_status
  
  call synchronize_all()

  if(myrank == 0) then

     allocate(spectrum_one_frequency(sum(ASKI_np_local_all),9))

     do jf = 1,ASKI_nf

        ip = 0

        ! this is me, rank 0
        if(ASKI_np_local > 0) then
           if(ASKI_DFT_double) then
              select case(ASKI_DFT_method)
              case('EXPLICIT_SUMMATION')
                 spectrum_one_frequency(1:ASKI_np_local,:) = transpose(ASKI_spectra_local_double(:,jf,:))
              case('GOERTZEL_STANDARD')
                 spectrum_one_frequency(1:ASKI_np_local,:) = cmplx( &
                      transpose(ASKI_Goertzel_U1_local_double(:,jf,:)) , & ! real part
                      transpose(ASKI_Goertzel_U2_local_double(:,jf,:)) )   ! imaginary part
              end select
           else
              select case(ASKI_DFT_method)
              case('EXPLICIT_SUMMATION')
                 spectrum_one_frequency(1:ASKI_np_local,:) = transpose(ASKI_spectra_local_single(:,jf,:))
              case('GOERTZEL_STANDARD')
                 spectrum_one_frequency(1:ASKI_np_local,:) = cmplx( &
                      transpose(ASKI_Goertzel_U1_local_single(:,jf,:)) , & ! real part
                      transpose(ASKI_Goertzel_U2_local_single(:,jf,:)) )   ! imaginary part
              end select
           end if
           if(ASKI_DECONVOLVE_STF) spectrum_one_frequency(1:ASKI_np_local,:) = &
                spectrum_one_frequency(1:ASKI_np_local,:) / ASKI_stf_spectrum_double(jf)
           ip = ip + ASKI_np_local
        end if
        do iproc = 1,NPROC-1
           if(ASKI_np_local_all(iproc+1) > 0) then
              call recv_c_t(spectrum_one_frequency(ip+1:ip+ASKI_np_local_all(iproc+1),:), ASKI_np_local_all(iproc+1)*9, iproc)
              ip = ip + ASKI_np_local_all(iproc+1)
           end if
        end do ! iproc

        ! write spectrum_one_frequency to file file basename.jf###### , where '######' contains the frequency index
        write(filename,"(a,i6.6)") trim(ASKI_outfile)//".jf",ASKI_jf(jf)
        ! do not worry about whether to overwrite or not, if program comes here, it already has been checked
        ! if the value of OVERWRITE_ASKI_OUTPUT is in conflict with existing files
        inquire(file=filename,exist=file_exists)
        if(file_exists) then
           open_status = 'replace'
        else
           open_status = 'new'
        end if

        specfem_version = 2
        ! also write NPROC to file, as a safety feature in order to assure that wavefield points (written only in ASKI_main_file
        ! and kernel values  (written only in kernel files) are not confused
        call get_file_unit_ASKI(IOASKI)
        open(unit=IOASKI, file=filename, status=open_status, access="stream", form="unformatted", action = 'write')
        write(IOASKI) specfem_version,length_ASKI_output_ID,ASKI_output_ID,NPROC,sum(ASKI_np_local_all),&
             ASKI_df,ASKI_jf(jf)
        write(IOASKI) spectrum_one_frequency
        close(IOASKI)

     end do ! jf

     deallocate(spectrum_one_frequency)

     ! write LOG_ASKI_finish.txt
     call write_ASKI_log('LOG_ASKI_finish.txt',"successfully created ASKI output, as specified in 'LOG_ASKI_start.txt'")

  else ! myrank == 0
     
     if(ASKI_np_local > 0) then
        allocate(spectrum_one_frequency_local(ASKI_np_local,9))
        do jf = 1,ASKI_nf
           if(ASKI_DFT_double) then
              select case(ASKI_DFT_method)
              case('EXPLICIT_SUMMATION')
                 spectrum_one_frequency_local = transpose(cmplx(ASKI_spectra_local_double(:,jf,:)))
              case('GOERTZEL_STANDARD')
                 spectrum_one_frequency_local = cmplx( &
                      transpose(ASKI_Goertzel_U1_local_double(:,jf,:)) , & ! real part
                      transpose(ASKI_Goertzel_U2_local_double(:,jf,:)))    ! imaginary part
              end select
           else
              select case(ASKI_DFT_method)
              case('EXPLICIT_SUMMATION')
                 spectrum_one_frequency_local = transpose(ASKI_spectra_local_single(:,jf,:))
              case('GOERTZEL_STANDARD')
                 spectrum_one_frequency_local = cmplx( &
                      transpose(ASKI_Goertzel_U1_local_single(:,jf,:)) , & ! real part
                      transpose(ASKI_Goertzel_U2_local_single(:,jf,:)))    ! imaginary part
              end select
           end if
           if(ASKI_DECONVOLVE_STF) spectrum_one_frequency_local = spectrum_one_frequency_local / ASKI_stf_spectrum_double(jf)
           call send_c_t(spectrum_one_frequency_local, ASKI_np_local*9, 0)
        end do
        deallocate(spectrum_one_frequency_local)
     end if

  end if ! myrank == 0

end subroutine write_ASKI_output_files
!
!-------------------------------------------------------------------------------------------
!
subroutine write_ASKI_log(filename,log_message)
  use constants,only: OUTPUT_FILES_BASE
  implicit none
  character(len=*) :: filename,log_message
  integer :: fu
  call get_file_unit_ASKI(fu)
  open(unit=fu,file=trim(OUTPUT_FILES_BASE)//trim(filename),&
       form='formatted',status='unknown',action='write')
  write(fu,*) trim(log_message)
  close(fu)
end subroutine write_ASKI_log
!
!-------------------------------------------------------------------------------------------
!
subroutine write_ASKI_log_start()
  use constants,only: OUTPUT_FILES_BASE
  use specfem_par,only: NPROC,NGLLX,NGLLY,NGLLZ
  use specfem_par_ASKI
  integer :: fu,i
  character(len=509) :: filename
  character(len=100) :: overwrite_message
  logical :: file_exists
  real :: size_main_file,size_jf_file
  call get_file_unit_ASKI(fu)
  open(unit=fu,file=trim(OUTPUT_FILES_BASE)//'LOG_ASKI_start.txt',&
       form='formatted',status='unknown',action='write')

  if(ASKI_MAIN_FILE_ONLY) then
     write(fu,*) "ONLY THE MAIN ASKI OUTPUT FILE WILL BE PRODUCED, as indicated by the logical parameter "//&
          "'ASKI_MAIN_FILE_ONLY' in DATA/Par_file_ASKI"
     write(fu,*) "HENCE, NO FREQUENCY KERNEL OUTPUT FILES WILL BE WRITTEN, EVEN IF INDICATED BELOW IN THIS LOGFILE!"
     write(fu,*) "For reasons of debugging and checking, this output was kept nevertheless."
     write(fu,*) ""
  end if

  write(fu,*) "Hello, this is SPECFEM3D_Cartesian for ASKI"
  write(fu,*) ""
  write(fu,*) "computing ASKI output now on ",NPROC," procs with following parameters:"
  write(fu,*) ""
  write(fu,*) "ASKI_type_inversion_grid = ",ASKI_type_inversion_grid
  select case(ASKI_type_inversion_grid)
  case(2,3,4)
     write(fu,*) "   ASKI_wx = ",ASKI_wx
     write(fu,*) "   ASKI_wy = ",ASKI_wy
     write(fu,*) "   ASKI_wz = ",ASKI_wz
  end select
  select case(ASKI_type_inversion_grid)
  case(3,4)
     write(fu,*) "   ASKI_rot_X = ",ASKI_rot_X
     write(fu,*) "   ASKI_rot_Y = ",ASKI_rot_Y
     write(fu,*) "   ASKI_rot_Z = ",ASKI_rot_Z
  case(2)
     write(fu,*) "   ASKI_rot_Z = ",ASKI_rot_Z
  end select
  select case(ASKI_type_inversion_grid)
  case(2,3,4)
     write(fu,*) "   ASKI_cx = ",ASKI_cx
     write(fu,*) "   ASKI_cy = ",ASKI_cy
     write(fu,*) "   ASKI_cz = ",ASKI_cz
  end select
  select case(ASKI_type_inversion_grid)
  case(4)
     write(fu,*) "   NGLLX = ",NGLLX
     write(fu,*) "   NGLLY = ",NGLLY
     write(fu,*) "   NGLLZ = ",NGLLZ
     write(fu,*) "   total number of inversion grid cells = ",sum(ASKI_np_local_all)/(NGLLX*NGLLY*NGLLZ)
  end select
  write(fu,*) ""
  write(fu,*) "local number of wavefield points (at which ASKI output is computed):"
  do i = 1,NPROC
     write(fu,*) "   proc ",i-1," : ",ASKI_np_local_all(i)
  end do ! iproc
  write(fu,*) "in total : ",sum(ASKI_np_local_all)
  write(fu,*) ""
  write(fu,*) "output spectra will be computed for ",ASKI_nf," frequencies f = df*jf defined by "
  write(fu,*) "   df = ",ASKI_df
  write(fu,*) "   jf = ",ASKI_jf
  write(fu,*) ""
  size_main_file = (length_ASKI_output_ID + 4.*(7+ASKI_nf+6*sum(ASKI_np_local_all)))/1048576.
  size_jf_file = (length_ASKI_output_ID + 4.*(6+2*9*sum(ASKI_np_local_all)))/1048576.
  write(fu,*) "in total ",ASKI_nf*size_jf_file+size_main_file,&
       " MiB output will be written to ",ASKI_nf+1," files with base_filename = "
  write(fu,*) "'"//trim(ASKI_outfile)//"' :"

  overwrite_message = ''
  if(OVERWRITE_ASKI_OUTPUT) then
     inquire(file=trim(ASKI_outfile)//".main",exist=file_exists)
     if(file_exists) then
        overwrite_message = 'exists and will be overwritten'
     else
        overwrite_message = 'does not exist and will be newly created'
     end if
  else
     ! if the file existed, this would have been detected above already, so if program comes here, the file does not exist
     overwrite_message = 'does not exist and will be newly created'
  end if
  write(fu,*) "   base_filename.main  (",size_main_file," MiB)  "//trim(overwrite_message)

  do i = 1,ASKI_nf
     write(filename,"(a,i6.6)") trim(ASKI_outfile)//".jf",ASKI_jf(i)
     overwrite_message = ''
     if(OVERWRITE_ASKI_OUTPUT) then
        inquire(file=trim(filename),exist=file_exists)
        if(file_exists) then
           overwrite_message = 'exists and will be overwritten'
        else
           overwrite_message = 'does not exist and will be newly created'
        end if
     else
        ! if the file existed, this would have been detected above already, so if program comes here, the file does not exist
        overwrite_message = 'does not exist and will be newly created'
     end if
     write(filename,"('   base_filename.jf',i6.6,'  ')") ASKI_jf(i)
     write(fu,*) trim(filename),"  (",size_jf_file," MiB)  "//trim(overwrite_message)
  end do ! i
  close(fu)
end subroutine write_ASKI_log_start
!
!-------------------------------------------------------------------------------------------
!
  subroutine get_file_unit_ASKI(unit_out)
   implicit none
   integer :: unit_out,fu
   logical :: is_open
   integer, parameter :: min_unit = 20
   integer, parameter :: max_unit = 99

   unit_out = -1
   do fu = min_unit, max_unit
      inquire(unit = fu, opened = is_open)
      if (.not. is_open) then
         unit_out = fu
         return
      end if
   end do
   call exit_MPI_without_rank('no file unit between 20 and 99 available')
 end subroutine get_file_unit_ASKI
