!----------------------------------------------------------------------------
!   Copyright 2013 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!   2012 Main authors: Dimitri Komatitsch and Jeroen Tromp
!
!   This file is part of SPECFEM3D_Cartesian version 2.1 and ASKI version 0.3.
!
!   SPECFEM3D_Cartesian version 2.1 and ASKI version 0.3 are free software: 
!   you can redistribute it and/or modify it under the terms of the GNU 
!   General Public License as published by the Free Software Foundation, 
!   either version 2 of the License, or (at your option) any later version.
!
!   SPECFEM3D_Cartesian version 2.1 and ASKI version 0.3 are distributed in 
!   the hope that they will be useful, but WITHOUT ANY WARRANTY; without 
!   even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
!   PURPOSE.  See the GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with SPECFEM3D_Cartesian version 2.1 and ASKI version 0.3.
!   If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------
!
! generic model file
!
! note: the idea is to super-impose velocity model values on the GLL points,
!          additional to the ones assigned on the CUBIT mesh
!
! most of the routines here are place-holders, please add/implement your own routines
!
!--------------------------------------------------------------------------------------------------



  module model_ASKI

!---
!
! ADD YOUR MODEL HERE
!
!---

  ! only here to illustrate an example
  !  type model_external_variables
  !    sequence
  !    double precision dvs(0:dummy_size)
  !  end type model_external_variables
  !  type (model_external_variables) MEXT_V
    

    integer :: model_ASKI_myrank


    type model_ASKI_cells
       ! inversion grid cells
       integer :: ncell
       real, dimension(:,:), pointer :: cc
       real, dimension(:), pointer :: r
       integer :: max_nnb
       integer, dimension(:,:), pointer :: nb ! of size (max_nnb+1,ncell)
    end type model_ASKI_cells
    type (model_ASKI_cells) :: mAc

    type model_ASKI_isoLame
       ! model values for parametrization isoLame
       real :: maxr_rho,maxr_lambda,maxr_mu
       integer :: nval_rho,nval_lambda,nval_mu
       integer, dimension(:), pointer :: idx_rho,idx_lambda,idx_mu
       real, dimension(:), pointer :: rho,lambda,mu
    end type model_ASKI_isoLame
    type (model_ASKI_isoLame) :: mAisoL

    type model_ASKI_isoVelocity
       ! model values for parametrization isoVelocity
       real :: maxr_rho,maxr_vp,maxr_vs
       integer :: nval_rho,nval_vp,nval_vs
       integer, dimension(:), pointer :: idx_rho,idx_vp,idx_vs
       real, dimension(:), pointer :: rho,vp,vs
    end type model_ASKI_isoVelocity
    type (model_ASKI_isoVelocity) :: mAisoV

    ! model_ASKI_pmtrz has values according to integer parameters 
    ! ipmtrz_isoLame,ipmtrz_isoVelocity,... below (selects which variable is used: mAisoV,mAisoL,...)
    integer :: model_ASKI_pmtrz

    ! type interpolation between of external model that is used
    ! type = 1 : case "shepard_standard" on first line in file "model_external_ASKI"
    !            interpolates model values by modified 3D Shepard interpolation with standard factor for influence radius
    ! type = 2 : case "shepard_factor_radius" on first line in file "model_external_ASKI"
    !            interpolates model values by modified 3D Shepard interpolation with given factor for influence radius
    integer :: model_ASKI_interpolation_type

    ! additional parameters dependent on model_ASKI_interpolation_type
    real :: model_ASKI_factor_shepard_radius


    ! other definitions
    integer, parameter :: ipmtrz_isoLame = 1
    integer, parameter :: ipmtrz_isoVelocity = 2

  contains

! must put this subroutine in the module, as dummy variables require an explicit interface
!---------------------------------------------------------------------------------
  subroutine shepard_interpolation_model_ASKI(x,y,z,factor_radius,nidx_cell,idx_cell,&
       n_C,C_prime,w)

    implicit none

    real :: factor_radius,x,y,z
    integer :: nidx_cell,n_C
    integer, dimension(nidx_cell) :: idx_cell
    integer, dimension(:), pointer :: C_prime
    real, dimension(:), pointer :: w

    real, dimension(nidx_cell) :: d
    logical, dimension(nidx_cell) :: larray
    integer :: iclose,i,j,n_C_tmp
    real :: r,r_prime,h,sum_s,smax_2
    real, dimension(:), pointer :: s,s_tmp
    real, dimension(:), allocatable :: t
    logical, dimension(:), allocatable :: larray2
    integer, dimension(:), pointer :: C_prime_tmp

!
! BEWARE: THIS ROUTINE MIGHT WELL BE OF NOT VERY GOOD PERFORMANCE!
! PLEASE DON'T HESITATE TO IMPROVE!
!

! nomenclature and method as in paper
!    Donald Shepard, "A two-dimensional interpolation function for irregularly-spaced data", 
!    Proceedings-1968 ACM National Conference

    n_C = 0
    nullify(C_prime,w)

    ! for given point P=(x,y,z) compute the distances d_i to all cell centers D_i
    d = sqrt( (mAc%cc(1,idx_cell)-x)**2 + (mAc%cc(2,idx_cell)-y)**2 + (mAc%cc(3,idx_cell)-z)**2 )

    ! select those cells, for which P is within their radius
    larray = d <= mAc%r(idx_cell)
    if(count(larray) == 0) return

    ! among all cells, for which P is within their radius, select the one with cell center closest to P
    iclose = minloc(d,1,larray)

    ! use the closest cell to P found above (has index iclose in incoming array idx_cell), 
    ! to define a total influence radius r= factor_radius*mAc%r(idx_cell(iclose)) about P 
    ! all cells with centers within radius r (n_C many) are chosen as collection
    ! C_prime, among which the interpolation will be done
    r = factor_radius*mAc%r(idx_cell(iclose))
    larray = d <= r
    n_C = count(larray)
    if(n_C == 0) return
    if(n_C == 1) then
       allocate(C_prime(1),w(1))
       C_prime = minloc(d,1,larray)
       w = 1.
       return
    end if

    ! if the minimum distance to a cell is exactly zero, we have to stop here, as we 
    ! cannot conduct the computations below and return the closest cell as only interpolation point
    if(minval(d,larray) == 0.) then
       n_C = 1
       allocate(C_prime(1),w(1))
       C_prime(1) = minloc(d,1,larray)
       w(1) = 1.
       return
    end if

    ! otherwise, pack preliminary collection C_prime_tmp
    n_C_tmp = n_C
    allocate(C_prime_tmp(n_C_tmp))
    C_prime_tmp = pack( (/ (i,i=1,nidx_cell) /) , larray)

    ! choose radius r_prime, beyond which the influence of a cell center will be absolutely zero, 
    ! as the shortest distance to P of a cell center which is not in collection C_prime_tmp, 
    ! or r if all cell centers are in the collection
    if(n_C_tmp == nidx_cell) then
       r_prime = r
    else
       r_prime = minval(d,.not.larray)
    end if

    ! compute preliminary values s_i
    allocate(s_tmp(n_C_tmp))
    do i = 1,n_C_tmp
       if(d(C_prime_tmp(i)) .le. r_prime/3.) then
          s_tmp(i) = 1./d(C_prime_tmp(i)) ! this should not be a problem, since it was checked above that all d(C_prime_tmp) > 0.
       else
          h = d(C_prime_tmp(i))/r_prime - 1.
          s_tmp(i) = 27.*h*h/(4.*r_prime)
       end if
    end do ! i

    ! for reasons of numerical stability, we need to assure that the values (s_i)^2 do not
    ! cause any significant roundoff when summing them up
    ! for this reason, we neglect all cell centers i in the collection for which (s_i)^2/(s_max)^2 < epsilon(1.0)
    ! (as those will behave as zero-weight points)
    smax_2 = maxval(s_tmp); smax_2 = smax_2*smax_2
    allocate(larray2(n_C_tmp))
    larray2 = s_tmp*s_tmp/smax_2 >= epsilon(1.0)
    n_C = count(larray2) ! note that n_C > 0 as s_tmp*s_tmp/smax_2 == 1 for maxval(s_tmp)
    if(n_C < n_C_tmp) then
       allocate(C_prime(n_C),s(n_C))
       C_prime = pack(C_prime_tmp,larray2)
       s = pack(s_tmp,larray2)
       deallocate(C_prime_tmp,s_tmp)
    else
       C_prime => C_prime_tmp
       s => s_tmp
       nullify(C_prime_tmp,s_tmp)
    end if
    deallocate(larray2)

    ! if there is only a single cell center left (very close to point P), the only weight is 1.
    if(n_C == 1) then
       allocate(w(1))
       w(1) = 1.
       deallocate(s)
       return
    end if

    ! otherwise define interpolation weights for final collection C_prime from s and 
    ! direction factors t (interpolation function f_3 in the paper)

    sum_s = sum(s)

    allocate(t(n_C))
    do i = 1,n_C
       t(i) = 0.
       do j = 1,n_C
          if(j==i) cycle
          h = (mAc%cc(1,C_prime(i)) - x)*(mAc%cc(1,C_prime(j)) - x) + &
              (mAc%cc(2,C_prime(i)) - y)*(mAc%cc(2,C_prime(j)) - y) + &
              (mAc%cc(3,C_prime(i)) - z)*(mAc%cc(3,C_prime(j)) - z)
          t(i) = t(i) + s(j)*(1.-h/(d(C_prime(i))*d(C_prime(j))))
       end do ! j
       t(i) = t(i)/sum_s
    end do ! i

    allocate(w(n_C))
    w = s*s*(1.+t)
    w = w/sum(w)

    deallocate(s,t)
  end subroutine shepard_interpolation_model_ASKI
!---------------------------------------------------------------------------------
! must put this subroutine in the module, as dummy variables require an explicit interface

  end module model_ASKI

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_external_broadcast(myrank)

! standard routine to setup model

  use model_ASKI

  implicit none

  include "constants.h"

  integer :: myrank

  ! local parameters
  character(len=800) :: error_message

  ! dummy to ignore compiler warnings
  model_ASKI_myrank = myrank

!---
!
! ADD YOUR MODEL HERE
!
!---

  ! the variables read are declared and stored in structure MEXT_V
  !if(myrank == 0) call read_external_model()

  ! broadcast the information read on the master to the nodes
  !call bcast_all_dp(MEXT_V%dvs, size(MEXT_V%dvs))

  if(myrank == 0) call read_external_model()

  call sync_all()

  call bcast_all_i(model_ASKI_interpolation_type,1)
  select case(model_ASKI_interpolation_type)
  case( 2 )
     call bcast_all_r(model_ASKI_factor_shepard_radius,1)
  end select

  call bcast_all_i(mAc%ncell,1)
  call bcast_all_i(mAc%max_nnb,1)
  if(myrank .ne. 0) then
     allocate(mAc%cc(3,mAc%ncell),mAc%r(mAc%ncell),&
          mAc%nb(mAc%max_nnb+1,mAc%ncell))
  end if
  call bcast_all_r(mAc%cc,size(mAc%cc))
  call bcast_all_r(mAc%r,size(mAc%r))
  call bcast_all_i(mAc%nb,size(mAc%nb))

  call bcast_all_i(model_ASKI_pmtrz,1)
  select case (model_ASKI_pmtrz)
     case ( ipmtrz_isoLame )

        call bcast_all_r(mAisoL%maxr_rho,1)
        call bcast_all_r(mAisoL%maxr_lambda,1)
        call bcast_all_r(mAisoL%maxr_mu,1)
        call bcast_all_i(mAisoL%nval_rho,1)
        call bcast_all_i(mAisoL%nval_lambda,1)
        call bcast_all_i(mAisoL%nval_mu,1)

        if(myrank .ne. 0) then
           if(mAisoL%nval_rho>0) allocate(mAisoL%idx_rho(mAisoL%nval_rho),&
                mAisoL%rho(mAisoL%nval_rho))
           if(mAisoL%nval_lambda>0) allocate(mAisoL%idx_lambda(mAisoL%nval_lambda),&
                mAisoL%lambda(mAisoL%nval_lambda))
           if(mAisoL%nval_mu>0) allocate(mAisoL%idx_mu(mAisoL%nval_mu),&
                mAisoL%mu(mAisoL%nval_mu))
        end if

        if(mAisoL%nval_rho>0) then
           call bcast_all_i(mAisoL%idx_rho,size(mAisoL%idx_rho))
           call bcast_all_r(mAisoL%rho,size(mAisoL%rho))
        end if
        if(mAisoL%nval_lambda>0) then
           call bcast_all_i(mAisoL%idx_lambda,size(mAisoL%idx_lambda))
           call bcast_all_r(mAisoL%lambda,size(mAisoL%lambda))
        end if
        if(mAisoL%nval_mu>0) then
           call bcast_all_i(mAisoL%idx_mu,size(mAisoL%idx_mu))
           call bcast_all_r(mAisoL%mu,size(mAisoL%mu))
        end if

     case ( ipmtrz_isoVelocity )

        call bcast_all_r(mAisoV%maxr_rho,1)
        call bcast_all_r(mAisoV%maxr_vp,1)
        call bcast_all_r(mAisoV%maxr_vs,1)
        call bcast_all_i(mAisoV%nval_rho,1)
        call bcast_all_i(mAisoV%nval_vp,1)
        call bcast_all_i(mAisoV%nval_vs,1)

        if(myrank .ne. 0) then
           if(mAisoV%nval_rho>0) allocate(mAisoV%idx_rho(mAisoV%nval_rho),&
                mAisoV%rho(mAisoV%nval_rho))
           if(mAisoV%nval_vp>0) allocate(mAisoV%idx_vp(mAisoV%nval_vp),&
                mAisoV%vp(mAisoV%nval_vp))
           if(mAisoV%nval_vs>0) allocate(mAisoV%idx_vs(mAisoV%nval_vs),&
                mAisoV%vs(mAisoV%nval_vs))
        end if

        if(mAisoV%nval_rho>0) then
           call bcast_all_i(mAisoV%idx_rho,size(mAisoV%idx_rho))
           call bcast_all_r(mAisoV%rho,size(mAisoV%rho))
        end if
        if(mAisoV%nval_vp>0) then
           call bcast_all_i(mAisoV%idx_vp,size(mAisoV%idx_vp))
           call bcast_all_r(mAisoV%vp,size(mAisoV%vp))
        end if
        if(mAisoV%nval_vs>0) then
           call bcast_all_i(mAisoV%idx_vs,size(mAisoV%idx_vs))
           call bcast_all_r(mAisoV%vs,size(mAisoV%vs))
        end if

     case default

        ! write error message to file and stop
        write(error_message,*) 'in model_external_broadcast: model_ASKI_pmtrz = '&
             ,model_ASKI_pmtrz,';  this parametrization index is not known: routines '//&
             'in model_external_values.f90 are inconsistent!'
        call stop_error_model_ASKI(error_message)

   end select

  end subroutine model_external_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_external_model()
!---
!
! ADD YOUR MODEL HERE
!
!---

  use model_ASKI

  implicit none
  
  include "constants.h"

  integer :: ier,IOASKI
  character(len=400) :: modelfile,line
  character(len=800) :: error_message
  character(len=50) :: char_type_external_model

  ier = -1

!
! READ FILE "model_external_ASKI" CONTAINING INTERPOLATION TYPE AND MODEL FILE NAME
!

  call get_file_unit_model_ASKI(IOASKI)
  open(unit=IOASKI,file=trim(IN_DATA_FILES_PATH)//'model_external_ASKI', &
       status='unknown',action='read',iostat=ier)
  if(ier .ne. 0) then
     close(IOASKI) ! ?? necessary? sensible?
     ! write error message to file and stop
     write(error_message,*) "in read_external_model: could not open file '"//trim(IN_DATA_FILES_PATH)//&
          'model_external_ASKI'//"' to read"
     call stop_error_model_ASKI(error_message)
  end if
  
  read(IOASKI,"(a400)",iostat=ier) line
  if(ier .ne. 0) then
     close(IOASKI)
     ! write error message to file and stop
     write(error_message,*) "in read_external_model: could not read first "//&
          "line of file '"//trim(IN_DATA_FILES_PATH)//'model_external_ASKI'//"'"
     call stop_error_model_ASKI(error_message)
  end if
  read(line,*,iostat=ier) char_type_external_model
  if(ier .ne. 0) then
     close(IOASKI)
     ! write error message to file and stop
     write(error_message,*) "in read_external_model: could not read interpolation type as first "//&
          "value on first line of file '"//trim(IN_DATA_FILES_PATH)//'model_external_ASKI'//"'"
     call stop_error_model_ASKI(error_message)
  end if
  select case(char_type_external_model)
     case('shepard_standard')
        model_ASKI_interpolation_type = 1
     case('shepard_factor_radius')
        model_ASKI_interpolation_type = 2
        read(line,*,iostat=ier) char_type_external_model,model_ASKI_factor_shepard_radius
        if(ier .ne. 0) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model: could not read factor for shepard "//&
                "interpolation radius (real number) as second value on first line of file '"//&
                trim(IN_DATA_FILES_PATH)//'model_external_ASKI'//"'"
           call stop_error_model_ASKI(error_message)
        end if
     case default
        close(IOASKI)
        ! write error message to file and stop
        write(error_message,*) "in read_external_model: interpolation type '"//trim(char_type_external_model)//&
             "' (first value on first line of file '"//trim(IN_DATA_FILES_PATH)//'model_external_ASKI'//&
             "') is not supported"
        call stop_error_model_ASKI(error_message)
  end select

  read(IOASKI,"(a400)",iostat=ier) line
  if(ier .ne. 0) then
     close(IOASKI)
     ! write error message to file and stop
     write(error_message,*) "in read_external_model: could not read second "//&
          "line of file '"//trim(IN_DATA_FILES_PATH)//'model_external_ASKI'//"'"
     call stop_error_model_ASKI(error_message)
  end if
  line = adjustl(trim(line))
  select case(line)
  case('kim_export')
     read(IOASKI,*,iostat=ier) modelfile
     if(ier .ne. 0) then
        close(IOASKI)
        ! write error message to file and stop
        write(error_message,*) "in read_external_model: could not read filename of ASKI kim_export file "//&
             "from third line of file '"//trim(IN_DATA_FILES_PATH)//'model_external_ASKI'//"'"
        call stop_error_model_ASKI(error_message)
     end if
     close(IOASKI)
     modelfile = adjustl(modelfile)
     call read_model_ASKI_kim_export(modelfile)
  case default
     close(IOASKI)
     ! write error message to file and stop
     write(error_message,*) "in read_external_model: model file type '"//trim(line)//&
          "' (value on second line of file '"//trim(IN_DATA_FILES_PATH)//'model_external_ASKI'//&
          "' is not supported, only 'kim_export' supported so far"
     call stop_error_model_ASKI(error_message)
  end select

  ! write logfile about this model to OUTPUT_FILES
  call get_file_unit_model_ASKI(IOASKI)
  open(unit=IOASKI,file=trim(OUTPUT_FILES_PATH)//'LOG_ASKI_model_external.txt',&
       form='formatted',status='unknown',action='write')
  write(IOASKI,*) "successfully read in ASKI external model"
  write(IOASKI,*) "ncell = ",mAc%ncell
  write(IOASKI,*) "maximal number of neighbours = ",mAc%max_nnb
  select case (model_ASKI_pmtrz)
     case ( ipmtrz_isoLame )
        write(IOASKI,*) "parametrization is isoLame"
        write(IOASKI,*) "nval_rho,nval_lambda,nval_mu = ",mAisoL%nval_rho,mAisoL%nval_lambda,mAisoL%nval_mu
        write(IOASKI,*) "maxr_rho,maxr_lambda,maxr_mu = ",mAisoL%maxr_rho,mAisoL%maxr_lambda,mAisoL%maxr_mu
     case ( ipmtrz_isoVelocity )
        write(IOASKI,*) "parametrization is isoVelocity"
        write(IOASKI,*) "nval_rho,nval_vp,nval_vs = ",mAisoV%nval_rho,mAisoV%nval_vp,mAisoV%nval_vs
        write(IOASKI,*) "maxr_rho,maxr_vp,maxr_vs = ",mAisoV%maxr_rho,mAisoV%maxr_vp,mAisoV%maxr_vs
  end select
  close(IOASKI)

  end subroutine read_external_model

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_model_ASKI_kim_export(modelfile)

  use model_ASKI

  implicit none

  include "constants.h"

  character(len=*) :: modelfile

  integer :: ier,IOASKI,nparam,iline
  character(len=400) :: line
  character(len=800) :: error_message

!
!NOW READ KIM EXPORT FILE
!

  call get_file_unit_model_ASKI(IOASKI)
  open(unit=IOASKI,file=trim(IN_DATA_FILES_PATH)//trim(modelfile),form='formatted', &
       status='old',action='read',iostat=ier)
  if(ier .ne. 0) then
     close(IOASKI) ! ?? necessary? sensible?
     ! write error message to file and stop
     write(error_message,*) "in read_model_ASKI_kim_export: could not open file '"//trim(IN_DATA_FILES_PATH)//&
          trim(modelfile)//"' to read (was given on second line of file '"//trim(IN_DATA_FILES_PATH)//&
          'model_external_ASKI'//"')"
     call stop_error_model_ASKI(error_message)
  end if

  iline = 0

  call check_pmtrz_ASKI_kim_export(IOASKI,iline,trim(IN_DATA_FILES_PATH)//trim(modelfile))

  call read_model_ASKI_kim_export_cells(IOASKI,iline,trim(IN_DATA_FILES_PATH)//trim(modelfile))

  select case (model_ASKI_pmtrz)
     case ( ipmtrz_isoLame )
        call read_external_model_ASKI_kim_export_isoLame(IOASKI,iline,trim(IN_DATA_FILES_PATH)//trim(modelfile))

     case ( ipmtrz_isoVelocity )
        call read_external_model_ASKI_kim_export_isoVelocity(IOASKI,iline,trim(IN_DATA_FILES_PATH)//trim(modelfile))

     case default
        close(IOASKI)
        ! write error message to file and stop
        write(error_message,*) 'in read_model_ASKI_kim_export: model_ASKI_pmtrz = '&
             ,model_ASKI_pmtrz,';  this parametrization index is not known: routines '//&
             'in model_external_values.f90 are inconsistent!'
        call stop_error_model_ASKI(error_message)
  end select

  close(IOASKI)

  end subroutine read_model_ASKI_kim_export

!
!-------------------------------------------------------------------------------------------------
!

  subroutine check_pmtrz_ASKI_kim_export(IOASKI,iline,filename)

  use model_ASKI

  implicit none

  integer :: IOASKI,iline
  character(len=*) :: filename

  integer :: ier,nparam
  character(len=11) :: pmtrz
  character(len=400) :: line
  character(len=800) :: error_message
  character(len=6), dimension(:), allocatable :: param

  ier = -1

  read(IOASKI,"(a400)",iostat=ier) line  ; iline = iline + 1
  if(ier .ne. 0) then
     close(IOASKI)
     ! write error message to file and stop
     write(error_message,*) "in check_pmtrz_ASKI_kim_export: could not read line ",iline,&
          "of file '"//trim(filename)//"'"
     call stop_error_model_ASKI(error_message)
  end if
  read(line,*,iostat=ier) pmtrz
  if(ier .ne. 0) then
     close(IOASKI)
     ! write error message to file and stop
     write(error_message,*) "in check_pmtrz_ASKI_kim_export: could not read model parametrization "//&
          "from line ",iline," of file '"//trim(filename)//"'"
     call stop_error_model_ASKI(error_message)
  end if

  select case(pmtrz)
  case('isoLame','isoVelocity') ! ok, do nothing
  case default
     close(IOASKI)
     ! write error message to file and stop
     write(error_message,*) "in check_pmtrz_ASKI_kim_export: model parametrization '"//&
          trim(pmtrz)//"' (on line ",iline," of file '"//trim(filename)//"')"//&
          " is not supported, only 'isoLame','isoVelocity' supported so far"
     call stop_error_model_ASKI(error_message)
  end select

  read(IOASKI,"(a400)",iostat=ier) line   ; iline = iline + 1
  if(ier .ne. 0) then
     close(IOASKI)
     ! write error message to file and stop
     write(error_message,*) "in check_pmtrz_ASKI_kim_export: could not read line ",iline,&
          "of file '"//trim(filename)//"'"
     call stop_error_model_ASKI(error_message)
  end if
  read(line,*,iostat=ier) nparam
  if(ier .ne. 0) then
     close(IOASKI)
     ! write error message to file and stop
     write(error_message,*) "in check_pmtrz_ASKI_kim_export: could not read number of parameters "//&
          "as first entry from line ",iline," of file '"//trim(filename)//"'"
     call stop_error_model_ASKI(error_message)
  end if
  if(nparam < 1) then
     close(IOASKI)
     ! write error message to file and stop
     write(error_message,*) "in check_pmtrz_ASKI_kim_export: number of parameters '",nparam,&
          "(first entry on line ",iline," of file '"//trim(filename)//&
          "') should be positive"
     call stop_error_model_ASKI(error_message)
  end if
  allocate(param(nparam))
  read(line,*,iostat=ier) nparam,param
  if(ier .ne. 0) then
     close(IOASKI)
     ! write error message to file and stop
     write(error_message,*) "in check_pmtrz_ASKI_kim_export: could not read ",nparam," parameter names "//&
          "as second to ",nparam+1,"-th entry from line ",iline," of file '"//trim(filename)//"'"
     call stop_error_model_ASKI(error_message)
  end if

  select case(pmtrz)
  case('isoLame')
     model_ASKI_pmtrz = ipmtrz_isoLame
     if(nparam/=3 .or. .not.(any(param == 'rho') .and. any(param == 'lambda') .and. any(param == 'mu') )) then
        close(IOASKI)
        ! write error message to file and stop
        write(error_message,*) "in check_pmtrz_ASKI_kim_export: model parametrization isoLame "//&
             "expects the 3 parameter names (on line ",iline," of file '"//trim(filename)//&
             "'), namely 'rho','lambda','mu' (in any order)"
        call stop_error_model_ASKI(error_message)
     end if
  case('isoVelocity')
     model_ASKI_pmtrz = ipmtrz_isoVelocity
     if(nparam/=3 .or. .not.(any(param == 'rho') .and. any(param == 'vp') .and. any(param == 'vs') )) then
        close(IOASKI)
        ! write error message to file and stop
        write(error_message,*) "in check_pmtrz_ASKI_kim_export: model parametrization isoVelocity "//&
             "expects the 3 parameter names (on line ",iline," of file '"//trim(filename)//&
             "'), namely 'rho','vp','vs' (in any order)"
        call stop_error_model_ASKI(error_message)
     end if
  case default
  end select

  deallocate(param)

  end subroutine check_pmtrz_ASKI_kim_export

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_model_ASKI_kim_export_cells(IOASKI,iline,filename)

  use model_ASKI

  implicit none

  integer :: IOASKI,iline
  character(len=*) :: filename

  integer :: ier,icell,nnb
  character(len=400) :: line
  character(len=800) :: error_message
  real :: c1,c2,c3,r
  integer, dimension(:,:), pointer :: tmp

  ier = -1

  read(IOASKI,*,iostat=ier) mAc%ncell   ; iline = iline + 1
  if(ier .ne. 0) then
     close(IOASKI)
     ! write error message to file and stop
     write(error_message,*) "in read_model_ASKI_kim_export_cells: could not read number of invgrid cells "//&
          "from line ",iline," of file '"//trim(filename)//"'"
     call stop_error_model_ASKI(error_message)
  end if
  if(mAc%ncell<1) then
     close(IOASKI)
     ! write error message to file and stop
     write(error_message,*) "in read_model_ASKI_kim_export_cells: number of invgrid cells ",mAc%ncell,&
          "(on line ",iline," of file '"//trim(filename)//"') must be positive"
     call stop_error_model_ASKI(error_message)
  end if

  mAc%max_nnb = 0
  allocate(mAc%cc(3,mAc%ncell),mAc%r(mAc%ncell),mAc%nb(mAc%max_nnb+1,mAc%ncell))
  
  do icell = 1,mAc%ncell

     read(IOASKI,"(a400)",iostat=ier) line   ; iline = iline+1
     if(ier .ne. 0) then
        close(IOASKI)
        ! write error message to file and stop
        write(error_message,*) "in read_model_ASKI_kim_export_cells: could not read "//&
             "line ",iline," of file '"//trim(filename)//"'"
        call stop_error_model_ASKI(error_message)
     end if
     read(line,*,iostat=ier) mAc%cc(:,icell),mAc%r(icell),nnb
     if(ier .ne. 0) then
        close(IOASKI)
        ! write error message to file and stop
        write(error_message,*) "in read_model_ASKI_kim_export_cells: could not read cell center, "//&
             "radius and number of neighbours from line ",iline," of file '"//trim(filename)//"'"
        call stop_error_model_ASKI(error_message)
     end if

     if(nnb>0) then
        mAc%nb(1,icell) = nnb
        if(nnb > mAc%max_nnb) then
           allocate(tmp(nnb+1,mAc%ncell))
           tmp(1:mAc%max_nnb+1,:) = mAc%nb
           deallocate(mAc%nb)
           mAc%nb => tmp
           nullify(tmp)
           mAc%max_nnb = nnb
        end if
        read(line,*,iostat=ier) c1,c2,c3,r,nnb,mAc%nb(2:nnb+1,icell)
        if(ier .ne. 0) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_model_ASKI_kim_export_cells: could not read ",nnb,&
                " neighbour indices after number of neighbours from line ",iline," of file '"//&
                trim(filename)//"'"
           call stop_error_model_ASKI(error_message)
        end if
     elseif(nnb==0) then
        mAc%nb(1,icell) = 0
     else
        write(error_message,*) "in read_model_ASKI_kim_export_cells: number of neighbours ",nnb,&
             "of cell ",icell," on line ",iline," of file '"//trim(filename)//&
             "' must not be negative: must be 0 if no neighbours"
        call stop_error_model_ASKI(error_message)
     end if

  end do ! icell

  end subroutine read_model_ASKI_kim_export_cells

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_external_model_ASKI_kim_export_isoLame(IOASKI,iline,filename)

  use model_ASKI

  implicit none

  integer :: IOASKI,iline
  character(len=*) :: filename

  integer :: ier,iparam
  character(len=400) :: line
  character(len=800) :: error_message
  character(len=6) :: one_param
  logical :: rho_found,lambda_found,mu_found

  ier = -1

  rho_found = .false.; lambda_found = .false.; mu_found = .false.

  ! read the three blocks containing model values for rho,lambda,mu
  do iparam = 1,3

     read(IOASKI,*,iostat=ier) one_param   ; iline = iline + 1
     if(ier .ne. 0) then
        close(IOASKI)
        ! write error message to file and stop
        write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: could not read parameter name"//&
             "from line ",iline," of file '"//trim(filename)//"'"
        call stop_error_model_ASKI(error_message)
     end if

     select case(one_param)
     case('rho')
        if(rho_found) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: parameter name 'rho' "//&
                "occurs for the second time on line ",iline," of file '"//trim(filename)//&
                "', only one block for 'rho' allowed"
           call stop_error_model_ASKI(error_message)
        end if
        rho_found = .true.
        read(IOASKI,*,iostat=ier) mAisoL%nval_rho   ; iline = iline + 1
        if(ier .ne. 0) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: could not read number of rho "//&
                "model values from line ",iline," of file '"//trim(filename)//"'"
           call stop_error_model_ASKI(error_message)
        end if
        if(mAisoL%nval_rho < 0) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: number of rho model values ",&
                mAisoL%nval_rho," on line ",iline," of file '"//trim(filename)//&
                "' must not be negative: must be 0 if no values"
           call stop_error_model_ASKI(error_message)
        elseif(mAisoL%nval_rho > 0) then
           allocate(mAisoL%idx_rho(mAisoL%nval_rho))
           read(IOASKI,*,iostat=ier) mAisoL%idx_rho   ; iline = iline + 1
           if(ier .ne. 0) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: could not read ",mAisoL%nval_rho,&
                   " cell indices for rho model values from line ",iline," of file '"//trim(filename)//"'"
              call stop_error_model_ASKI(error_message)
           end if
           if(any(mAisoL%idx_rho < 1 .or. mAisoL%idx_rho > mAc%ncell)) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: there are cell indices for "//&
                   "rho model values on line ",iline," of file '"//trim(filename)//&
                   "' which exceed the lower limit 1 or the upper limit ncell = ",mAc%ncell
              call stop_error_model_ASKI(error_message)
           end if
           allocate(mAisoL%rho(mAisoL%nval_rho))
           read(IOASKI,*,iostat=ier) mAisoL%rho   ; iline = iline + 1
           if(ier .ne. 0) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: could not read ",mAisoL%nval_rho,&
                   " rho model values from line ",iline," of file '"//trim(filename)//"'"
              call stop_error_model_ASKI(error_message)
           end if
        end if

     case('lambda')
        if(lambda_found) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: parameter name 'lambda' "//&
                "occurs for the second time on line ",iline," of file '"//trim(filename)//&
                "', only one block for 'lambda' allowed"
           call stop_error_model_ASKI(error_message)
        end if
        lambda_found = .true.
        read(IOASKI,*,iostat=ier) mAisoL%nval_lambda   ; iline = iline + 1
        if(ier .ne. 0) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: could not read number of lambda "//&
                "model values from line ",iline," of file '"//trim(filename)//"'"
           call stop_error_model_ASKI(error_message)
        end if
        if(mAisoL%nval_lambda < 0) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: number of lambda model values ",&
                mAisoL%nval_lambda," on line ",iline," of file '"//trim(filename)//&
                "' must not be negative: must be 0 if no values"
           call stop_error_model_ASKI(error_message)
        elseif(mAisoL%nval_lambda > 0) then
           allocate(mAisoL%idx_lambda(mAisoL%nval_lambda))
           read(IOASKI,*,iostat=ier) mAisoL%idx_lambda   ; iline = iline + 1
           if(ier .ne. 0) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: could not read ",mAisoL%nval_lambda,&
                   " cell indices for lambda model values from line ",iline," of file '"//trim(filename)//"'"
              call stop_error_model_ASKI(error_message)
           end if
           if(any(mAisoL%idx_lambda < 1 .or. mAisoL%idx_lambda > mAc%ncell)) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: there are cell indices for "//&
                   "lambda model values on line ",iline," of file '"//trim(filename)//&
                   "' which exceed the lower limit 1 or the upper limit ncell = ",mAc%ncell
              call stop_error_model_ASKI(error_message)
           end if
           allocate(mAisoL%lambda(mAisoL%nval_lambda))
           read(IOASKI,*,iostat=ier) mAisoL%lambda   ; iline = iline + 1
           if(ier .ne. 0) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: could not read ",mAisoL%nval_lambda,&
                   " lambda model values from line ",iline," of file '"//trim(filename)//"'"
              call stop_error_model_ASKI(error_message)
           end if
        end if

     case('mu')
        if(mu_found) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: parameter name 'mu' "//&
                "occurs for the second time on line ",iline," of file '"//trim(filename)//&
                "', only one block for 'mu' allowed"
           call stop_error_model_ASKI(error_message)
        end if
        mu_found = .true.
        read(IOASKI,*,iostat=ier) mAisoL%nval_mu   ; iline = iline + 1
        if(ier .ne. 0) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: could not read number of mu "//&
                "model values from line ",iline," of file '"//trim(filename)//"'"
           call stop_error_model_ASKI(error_message)
        end if
        if(mAisoL%nval_mu < 0) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: number of mu model values ",&
                mAisoL%nval_mu," on line ",iline," of file '"//trim(filename)//&
                "' must not be negative: must be 0 if no values"
           call stop_error_model_ASKI(error_message)
        elseif(mAisoL%nval_mu > 0) then
           allocate(mAisoL%idx_mu(mAisoL%nval_mu))
           read(IOASKI,*,iostat=ier) mAisoL%idx_mu   ; iline = iline + 1
           if(ier .ne. 0) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: could not read ",mAisoL%nval_mu,&
                   " cell indices for mu model values from line ",iline," of file '"//trim(filename)//"'"
              call stop_error_model_ASKI(error_message)
           end if
           if(any(mAisoL%idx_mu < 1 .or. mAisoL%idx_mu > mAc%ncell)) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: there are cell indices for "//&
                   "mu model values on line ",iline," of file '"//trim(filename)//&
                   "' which exceed the lower limit 1 or the upper limit ncell = ",mAc%ncell
              call stop_error_model_ASKI(error_message)
           end if
           allocate(mAisoL%mu(mAisoL%nval_mu))
           read(IOASKI,*,iostat=ier) mAisoL%mu   ; iline = iline + 1
           if(ier .ne. 0) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: could not read ",mAisoL%nval_mu,&
                   " mu model values from line ",iline," of file '"//trim(filename)//"'"
              call stop_error_model_ASKI(error_message)
           end if
        end if

     case default
        close(IOASKI)
        ! write error message to file and stop
        write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: parameter name '"//trim(one_param)//&
             "' on line ",iline," of file '"//trim(filename)//"' is not valid in parametrization isoVelocity"
        call stop_error_model_ASKI(error_message)

     end select

  end do ! iparam

  close(IOASKI)

  ! compute maxr_rho,maxr_lambda,maxr_mu
  if(mAisoL%nval_rho > 0) mAisoL%maxr_rho = maxval(mAc%r(mAisoL%idx_rho))
  if(mAisoL%nval_lambda > 0) mAisoL%maxr_lambda = maxval(mAc%r(mAisoL%idx_lambda))
  if(mAisoL%nval_mu > 0) mAisoL%maxr_mu = maxval(mAc%r(mAisoL%idx_mu))

  end subroutine read_external_model_ASKI_kim_export_isoLame

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_external_model_ASKI_kim_export_isoVelocity(IOASKI,iline,filename)

  use model_ASKI

  implicit none

  integer :: IOASKI,iline
  character(len=*) :: filename

  integer :: ier,iparam
  character(len=400) :: line
  character(len=800) :: error_message
  character(len=6) :: one_param
  logical :: rho_found,vp_found,vs_found

  ier = -1

  rho_found = .false.; vp_found = .false.; vs_found = .false.

  ! read the three blocks containing model values for rho,vp,vs
  do iparam = 1,3

     read(IOASKI,*,iostat=ier) one_param   ; iline = iline + 1
     if(ier .ne. 0) then
        close(IOASKI)
        ! write error message to file and stop
        write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: could not read parameter name"//&
             "from line ",iline," of file '"//trim(filename)//"'"
        call stop_error_model_ASKI(error_message)
     end if

     select case(one_param)
     case('rho')
        if(rho_found) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: parameter name 'rho' "//&
                "occurs for the second time on line ",iline," of file '"//trim(filename)//&
                "', only one block for 'rho' allowed"
           call stop_error_model_ASKI(error_message)
        end if
        rho_found = .true.
        read(IOASKI,*,iostat=ier) mAisoV%nval_rho   ; iline = iline + 1
        if(ier .ne. 0) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: could not read number of rho "//&
                "model values from line ",iline," of file '"//trim(filename)//"'"
           call stop_error_model_ASKI(error_message)
        end if
        if(mAisoV%nval_rho < 0) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: number of rho model values ",&
                mAisoV%nval_rho," on line ",iline," of file '"//trim(filename)//&
                "' must not be negative: must be 0 if no values"
           call stop_error_model_ASKI(error_message)
        elseif(mAisoV%nval_rho > 0) then
           allocate(mAisoV%idx_rho(mAisoV%nval_rho))
           read(IOASKI,*,iostat=ier) mAisoV%idx_rho   ; iline = iline + 1
           if(ier .ne. 0) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: could not read ",mAisoV%nval_rho,&
                   " cell indices for rho model values from line ",iline," of file '"//trim(filename)//"'"
              call stop_error_model_ASKI(error_message)
           end if
           if(any(mAisoV%idx_rho < 1 .or. mAisoV%idx_rho > mAc%ncell)) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: there are cell indices for "//&
                   "rho model values on line ",iline," of file '"//trim(filename)//&
                   "' which exceed the lower limit 1 or the upper limit ncell = ",mAc%ncell
              call stop_error_model_ASKI(error_message)
           end if
           allocate(mAisoV%rho(mAisoV%nval_rho))
           read(IOASKI,*,iostat=ier) mAisoV%rho   ; iline = iline + 1
           if(ier .ne. 0) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: could not read ",mAisoV%nval_rho,&
                   " rho model values from line ",iline," of file '"//trim(filename)//"'"
              call stop_error_model_ASKI(error_message)
           end if
        end if

     case('vp')
        if(vp_found) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: parameter name 'vp' "//&
                "occurs for the second time on line ",iline," of file '"//trim(filename)//&
                "', only one block for 'vp' allowed"
           call stop_error_model_ASKI(error_message)
        end if
        vp_found = .true.
        read(IOASKI,*,iostat=ier) mAisoV%nval_vp   ; iline = iline + 1
        if(ier .ne. 0) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: could not read number of vp "//&
                "model values from line ",iline," of file '"//trim(filename)//"'"
           call stop_error_model_ASKI(error_message)
        end if
        if(mAisoV%nval_vp < 0) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: number of vp model values ",&
                mAisoV%nval_vp," on line ",iline," of file '"//trim(filename)//&
                "' must not be negative: must be 0 if no values"
           call stop_error_model_ASKI(error_message)
        elseif(mAisoV%nval_vp > 0) then
           allocate(mAisoV%idx_vp(mAisoV%nval_vp))
           read(IOASKI,*,iostat=ier) mAisoV%idx_vp   ; iline = iline + 1
           if(ier .ne. 0) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: could not read ",mAisoV%nval_vp,&
                   " cell indices for vp model values from line ",iline," of file '"//trim(filename)//"'"
              call stop_error_model_ASKI(error_message)
           end if
           if(any(mAisoV%idx_vp < 1 .or. mAisoV%idx_vp > mAc%ncell)) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: there are cell indices for "//&
                   "vp model values on line ",iline," of file '"//trim(filename)//&
                   "' which exceed the lower limit 1 or the upper limit ncell = ",mAc%ncell
              call stop_error_model_ASKI(error_message)
           end if
           allocate(mAisoV%vp(mAisoV%nval_vp))
           read(IOASKI,*,iostat=ier) mAisoV%vp   ; iline = iline + 1
           if(ier .ne. 0) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: could not read ",mAisoV%nval_vp,&
                   " vp model values from line ",iline," of file '"//trim(filename)//"'"
              call stop_error_model_ASKI(error_message)
           end if
        end if

     case('vs')
        if(vs_found) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: parameter name 'vs' "//&
                "occurs for the second time on line ",iline," of file '"//trim(filename)//&
                "', only one block for 'vs' allowed"
           call stop_error_model_ASKI(error_message)
        end if
        vs_found = .true.
        read(IOASKI,*,iostat=ier) mAisoV%nval_vs   ; iline = iline + 1
        if(ier .ne. 0) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: could not read number of vs "//&
                "model values from line ",iline," of file '"//trim(filename)//"'"
           call stop_error_model_ASKI(error_message)
        end if
        if(mAisoV%nval_vs < 0) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: number of vs model values ",&
                mAisoV%nval_vs," on line ",iline," of file '"//trim(filename)//&
                "' must not be negative: must be 0 if no values"
           call stop_error_model_ASKI(error_message)
        elseif(mAisoV%nval_vs > 0) then
           allocate(mAisoV%idx_vs(mAisoV%nval_vs))
           read(IOASKI,*,iostat=ier) mAisoV%idx_vs   ; iline = iline + 1
           if(ier .ne. 0) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: could not read ",mAisoV%nval_vs,&
                   " cell indices for vs model values from line ",iline," of file '"//trim(filename)//"'"
              call stop_error_model_ASKI(error_message)
           end if
           if(any(mAisoV%idx_vs < 1 .or. mAisoV%idx_vs > mAc%ncell)) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: there are cell indices for "//&
                   "vs model values on line ",iline," of file '"//trim(filename)//&
                   "' which exceed the lower limit 1 or the upper limit ncell = ",mAc%ncell
              call stop_error_model_ASKI(error_message)
           end if
           allocate(mAisoV%vs(mAisoV%nval_vs))
           read(IOASKI,*,iostat=ier) mAisoV%vs   ; iline = iline + 1
           if(ier .ne. 0) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: could not read ",mAisoV%nval_vs,&
                   " vs model values from line ",iline," of file '"//trim(filename)//"'"
              call stop_error_model_ASKI(error_message)
           end if
        end if

     case default
        close(IOASKI)
        ! write error message to file and stop
        write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: parameter name '"//trim(one_param)//&
             "' on line ",iline," of file '"//trim(filename)//"' is not valid in parametrization isoVelocity"
        call stop_error_model_ASKI(error_message)

     end select

  end do ! iparam

  close(IOASKI)

  ! compute maxr_rho,maxr_vp,maxr_vs
  if(mAisoV%nval_rho > 0) mAisoV%maxr_rho = maxval(mAc%r(mAisoV%idx_rho))
  if(mAisoV%nval_vp > 0) mAisoV%maxr_vp = maxval(mAc%r(mAisoV%idx_vp))
  if(mAisoV%nval_vs > 0) mAisoV%maxr_vs = maxval(mAc%r(mAisoV%idx_vs))

  end subroutine read_external_model_ASKI_kim_export_isoVelocity

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_external_values(xmesh,ymesh,zmesh,rho,vp,vs,qkappa_atten,qmu_atten,iflag_aniso,idomain_id )

! given a GLL point, returns super-imposed velocity model values

!  use generate_databases_par,only: nspec => NSPEC_AB,ibool

  use create_regions_mesh_ext_par,only: CUSTOM_REAL

  use model_ASKI

  implicit none

  ! GLL point
  double precision, intent(in) :: xmesh,ymesh,zmesh

  ! density, Vp and Vs
  real(kind=CUSTOM_REAL) :: vp,vs,rho

  ! attenuation flag
  real(kind=CUSTOM_REAL) :: qkappa_atten,qmu_atten

  ! anisotropy flag
  integer :: iflag_aniso

  ! acoustic/elastic/.. domain flag ( 1 = acoustic / 2 = elastic / ... )
  integer :: idomain_id

  ! local parameters
  real :: x,y,z

!---
!
! ADD YOUR MODEL HERE
!
!---

  ! GLL point location converted to real
  x = xmesh
  y = ymesh
  z = zmesh

  select case (model_ASKI_pmtrz)
  case ( ipmtrz_isoLame )
     ! no attenuation, so simply leave values qkappa_atten,qmu_atten untouched
     call model_external_values_ASKI_isoLame(x,y,z,rho,vp,vs)
  case ( ipmtrz_isoVelocity )
     ! no attenuation, so simply leave values qkappa_atten,qmu_atten untouched
     call model_external_values_ASKI_isoVelocity(x,y,z,rho,vp,vs)
  end select

  end subroutine model_external_values

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_external_values_ASKI_isoLame(x,y,z,rho,vp,vs)

! given a GLL point, do neighbour interpolation between cell centers defined by values in mAc,mAisoL

  use model_ASKI

! dummy variables
  ! point coordinates
  real :: x,y,z
  ! density, Vp and Vs
  real :: vp,vs,rho

! local variables
  character(len=800) :: error_message
  real :: factor_shepard_radius,rho_interpolate,lambda_interpolate,mu_interpolate
  integer :: n_rho,n_lambda,n_mu,j
  integer, dimension(:), pointer :: i_rho,i_lambda,i_mu
  real, dimension(:), pointer :: w_rho,w_lambda,w_mu

  ! we cannot compute vp,vs from rho,lambda,mu if any of lambda,mu values are not present (if no rho values present, use incoming default rho)
  if(mAisoL%nval_lambda == 0 .or. mAisoL%nval_mu == 0) return

  nullify(i_rho,i_lambda,i_mu,w_rho,w_lambda,w_mu)

  select case(model_ASKI_interpolation_type)
  case(1); factor_shepard_radius = 2.
  case(2); factor_shepard_radius = model_ASKI_factor_shepard_radius 
  end select

  ! deal with lambda, already checked that mAisoL%nval_lambda /= 0
  call shepard_interpolation_model_ASKI(x,y,z,factor_shepard_radius,mAisoL%nval_lambda,mAisoL%idx_lambda,&
       n_lambda,i_lambda,w_lambda)
  if(n_lambda==0) goto 1
  lambda_interpolate = sum( w_lambda * mAisoL%lambda(i_lambda) )

  ! deal with mu, already checked that mAisoL%nval_mu /= 0
  call shepard_interpolation_model_ASKI(x,y,z,factor_shepard_radius,mAisoL%nval_mu,mAisoL%idx_mu,&
       n_mu,i_mu,w_mu)
  if(n_mu/=n_lambda) goto 1
  do j = 1,n_mu
     if(.not.any(mAisoL%idx_lambda(i_lambda) == mAisoL%idx_mu(i_mu(j)))) goto 1
  end do
  mu_interpolate = sum( w_mu * mAisoL%mu(i_mu) )

  ! deal with rho
  if(mAisoL%nval_rho == 0) then
     ! if there are no rho values in the model, use incoming default rho (background rho)
     rho_interpolate = rho
  else
     call shepard_interpolation_model_ASKI(x,y,z,factor_shepard_radius,mAisoL%nval_rho,mAisoL%idx_rho,&
          n_rho,i_rho,w_rho)
     if(n_rho==0) then
        ! if point x,y,z is outside the area where rho values are given in the model, use incoming default rho (background rho)
        rho_interpolate = rho
     else
        if(n_rho/=n_lambda) goto 1
        do j = 1,n_rho
           if(.not.any(mAisoL%idx_lambda(i_lambda) == mAisoL%idx_rho(i_rho(j)))) goto 1
        end do
        rho_interpolate = sum( w_rho * mAisoL%rho(i_rho) )
     end if
  end if

  ! if code comes here, for ALL parameters lambda,mu (maybe rho is the background rho)
  ! the SAME inversion grid cell centers were chose for the interpolation rule, 
  ! hence we can compute sensible seismic velocities from an interpolation of rho,lambda and mu

  vp = sqrt( (lambda_interpolate+2.*mu_interpolate) / rho_interpolate )
  vs = sqrt( (mu_interpolate) / rho_interpolate )
  rho = rho_interpolate

1 if(associated(i_rho)) deallocate(i_rho)
  if(associated(w_rho)) deallocate(w_rho)
  if(associated(i_lambda)) deallocate(i_lambda)
  if(associated(w_lambda)) deallocate(w_lambda)
  if(associated(i_mu)) deallocate(i_mu)
  if(associated(w_mu)) deallocate(w_mu)

  end subroutine model_external_values_ASKI_isoLame

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_external_values_ASKI_isoVelocity(x,y,z,rho,vp,vs)

! given a GLL point, do neighbour interpolation between cell centers defined by values in mAc,mAisoV

  use model_ASKI

! dummy variables
  ! point coordinates
  real :: x,y,z
  ! density, Vp and Vs
  real :: vp,vs,rho

! local variables
  character(len=800) :: error_message
  real :: factor_shepard_radius
  integer :: n
  integer, dimension(:), pointer :: interpolation_index
  real, dimension(:), pointer :: interpolation_weight

  nullify(interpolation_index,interpolation_weight)

  select case(model_ASKI_interpolation_type)
  case(1); factor_shepard_radius = 2.
  case(2); factor_shepard_radius = model_ASKI_factor_shepard_radius 
  end select

  ! deal with rho
  if(mAisoV%nval_rho > 0) then
     call shepard_interpolation_model_ASKI(x,y,z,factor_shepard_radius,mAisoV%nval_rho,mAisoV%idx_rho,&
          n,interpolation_index,interpolation_weight)
     if(n>0) then
        ! assume that in this case, interpolation_index,interpolation_weight are associated and have size n
        rho = sum( interpolation_weight * mAisoV%rho(interpolation_index) )
        deallocate(interpolation_index,interpolation_weight)
     end if
  end if

  ! deal with vp
  if(mAisoV%nval_vp > 0) then
     call shepard_interpolation_model_ASKI(x,y,z,factor_shepard_radius,mAisoV%nval_vp,mAisoV%idx_vp,&
          n,interpolation_index,interpolation_weight)
     if(n>0) then
        ! assume that in this case, interpolation_index,interpolation_weight are associated and have size n
        vp = sum( interpolation_weight * mAisoV%vp(interpolation_index) )
        deallocate(interpolation_index,interpolation_weight)
     end if
  end if

  ! deal with vs
  if(mAisoV%nval_vs > 0) then
     call shepard_interpolation_model_ASKI(x,y,z,factor_shepard_radius,mAisoV%nval_vs,mAisoV%idx_vs,&
          n,interpolation_index,interpolation_weight)
     if(n>0) then
        ! assume that in this case, interpolation_index,interpolation_weight are associated and have size n
        vs = sum( interpolation_weight * mAisoV%vs(interpolation_index) )
        deallocate(interpolation_index,interpolation_weight)
     end if
  end if

  end subroutine model_external_values_ASKI_isoVelocity

!
!---------------------------------------------------------------------------------
!

  subroutine stop_error_model_ASKI(error_message)
  use model_ASKI
  implicit none
  include "constants.h"

  character(len=*) :: error_message
  character(len=400) :: filename
  integer :: IOASKI

  write(filename,"(a,i6.6,a)") trim(OUTPUT_FILES_PATH)//'ERROR_model_external_ASKI_',model_ASKI_myrank,'.txt'

  call get_file_unit_model_ASKI(IOASKI)
  open(unit=IOASKI,file=filename,form='formatted',status='unknown',action='write')
  write(IOASKI,*) trim(error_message)
  close(IOASKI)
  call stop_all()

  end subroutine stop_error_model_ASKI
!
!---------------------------------------------------------------------------------
!
  subroutine get_file_unit_model_ASKI(unit_out)
   implicit none
   integer :: unit_out
   integer :: fu
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
 end subroutine get_file_unit_model_ASKI
