!----------------------------------------------------------------------------
!   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!   2012 Main authors: Dimitri Komatitsch and Jeroen Tromp
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

    use constants,only: MAX_STRING_LEN

    implicit none

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
    

    ! defined by flag "USE_ASKI_BACKGROUND_MODEL" in Par_file_ASKI, logical "use_ASKI_background_model" indicates
    ! whether a global 1D background model should be set.
    ! (this will replace the default background model everywhere, before possibly imposing the inverted model)
    ! If IMPOSE_ASKI_BACKGROUND_MODEL = .true. , file_ASKI_background_model then contains the value of
    ! FILE_ASKI_BACKGROUND_MODEL, i.e. the model file, relative to DATA/
    logical :: use_ASKI_background_model = .false.
    character(len=MAX_STRING_LEN) :: file_ASKI_background_model = ''
    
    ! Defined by flag "IMPOSE_ASKI_INVERTED_MODEL" in Par_file_ASKI, logical "impose_ASKI_inverted_model" indicates
    ! whether the ASKI external model should be imposed.
    ! If IMPOSE_ASKI_INVERTED_MODEL = .true. , file_ASKI_inverted_model then contains the value of
    ! FILE_ASKI_INVERTED_MODEL, i.e. the model file, relative to DATA/
    logical :: impose_ASKI_inverted_model = .false.
    character(len=MAX_STRING_LEN) :: file_ASKI_inverted_model = ''


    integer :: model_ASKI_myrank


    type model_ASKI_1Dbackground
       ! 1D spline-interpolated background model to replace the overall SPECFEM3D background model,
       ! before possibly imposing the ASKI inverted model 
       ! This mechanism was implemented in order to properly continue an inversion where the first
       ! iteration was done by a 1D method. Since the simulation domain is usually larger than the inversion
       ! domain, the inverted 3D model should be extendet to the rest of the inversion domain by the very 1D
       ! reference model that was used before by the 1D method.
       integer :: nlayers
       real :: zmax
       integer, dimension(:), pointer :: nnodes
       real, dimension(:,:), pointer :: depth,rho,vp,vs,Qmu,Qkappa  ! depth and parameter arrays
       real, dimension(:,:), pointer :: sprho,spvp,spvs,spQmu,spQkappa  ! spline parameters p = s''(xj)
    end type model_ASKI_1Dbackground
    type (model_ASKI_1Dbackground) :: mA1Db

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
    ! ipmtrz_isoLameSI,ipmtrz_isoVelocitySI,... below (selects which variable is used: mAisoV,mAisoL,...)
    integer :: model_ASKI_pmtrz

    ! type interpolation between of external model that is used
    ! type = 1 : case "shepard_standard"
    !            interpolates model values by modified 3D Shepard interpolation with standard factor for influence radius
    ! type = 2 : case "shepard_factor_radius"
    !            interpolates model values by modified 3D Shepard interpolation with given factor for influence radius
    integer :: model_ASKI_interpolation_type

    ! additional parameters dependent on model_ASKI_interpolation_type
    real :: model_ASKI_factor_shepard_radius


    ! other definitions
    integer, parameter :: ipmtrz_isoLameSI = 1
    integer, parameter :: ipmtrz_isoVelocitySI = 2

  contains

! must put this subroutine in the module, as dummy variables require an explicit interface
!---------------------------------------------------------------------------------
  subroutine shepard_interpolation_model_ASKI(x,y,z,factor_radius,nidx_cell,idx_cell,&
       n_C,C_prime,w,enlarge_influence,radius_factor_closest_cell)

    implicit none

    real :: factor_radius,x,y,z
    integer :: nidx_cell,n_C
    integer, dimension(nidx_cell) :: idx_cell
    integer, dimension(:), pointer :: C_prime
    real, dimension(:), pointer :: w
    real, optional :: enlarge_influence,radius_factor_closest_cell

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
    if(present(enlarge_influence)) then
       larray = d <= mAc%r(idx_cell)*enlarge_influence
    else
       larray = d <= mAc%r(idx_cell)
    end if
    if(count(larray) == 0) return

    ! among all cells, for which P is within their radius, select the one with cell center closest to P
    iclose = minloc(d,1,larray)

    ! if radius_factor_closest_cell is present, define it here
    if(present(radius_factor_closest_cell)) then
       radius_factor_closest_cell = d(iclose)/mAc%r(idx_cell(iclose))
    end if

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

  integer :: myrank

  ! local parameters
  character(len=800) :: error_message
  integer :: maxnnodes
  integer, dimension(1) :: i_array_one_value
  real, dimension(1) :: r_array_one_value
  logical, dimension(1) :: l_array_one_value

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

  if(myrank == 0) then
     call read_external_model()
     if(use_ASKI_background_model) maxnnodes = maxval(mA1Db%nnodes)
  end if

  call synchronize_all()

  if(myrank == 0) l_array_one_value(1) = use_ASKI_background_model
  call bcast_all_l(l_array_one_value,1)
  if(myrank /= 0) use_ASKI_background_model = l_array_one_value(1)

  if(myrank == 0) l_array_one_value(1) = impose_ASKI_inverted_model
  call bcast_all_l(l_array_one_value,1)
  if(myrank /= 0) impose_ASKI_inverted_model = l_array_one_value(1)

  if(use_ASKI_background_model) then

     call bcast_all_string_world(file_ASKI_background_model)

     if(myrank == 0) r_array_one_value(1) = mA1Db%zmax
     call bcast_all_r(r_array_one_value,1)
     if(myrank /= 0) mA1Db%zmax = r_array_one_value(1)

     if(myrank == 0) i_array_one_value(1) = mA1Db%nlayers
     call bcast_all_i(i_array_one_value,1)
     if(myrank /= 0) mA1Db%nlayers = i_array_one_value(1)

     if(myrank == 0) i_array_one_value(1) = maxnnodes
     call bcast_all_i(i_array_one_value,1)
     if(myrank /= 0) maxnnodes = i_array_one_value(1)

     ! allocate for model values if I'm not rank 0
     if(myrank .ne. 0) allocate(mA1Db%nnodes(mA1Db%nlayers),mA1Db%depth(mA1Db%nlayers,maxnnodes), &
          mA1Db%rho(mA1Db%nlayers,maxnnodes),mA1Db%vp(mA1Db%nlayers,maxnnodes),mA1Db%vs(mA1Db%nlayers,maxnnodes), &
          mA1Db%Qmu(mA1Db%nlayers,maxnnodes),mA1Db%Qkappa(mA1Db%nlayers,maxnnodes),&
          mA1Db%sprho(mA1Db%nlayers,maxnnodes),mA1Db%spvp(mA1Db%nlayers,maxnnodes),mA1Db%spvs(mA1Db%nlayers,maxnnodes),&
          mA1Db%spQmu(mA1Db%nlayers,maxnnodes),mA1Db%spQkappa(mA1Db%nlayers,maxnnodes))
 
     call bcast_all_i(mA1Db%nnodes,size(mA1Db%nnodes))
     call bcast_all_r(mA1Db%depth,size(mA1Db%depth))
     call bcast_all_r(mA1Db%rho,size(mA1Db%rho))
     call bcast_all_r(mA1Db%vp,size(mA1Db%vp))
     call bcast_all_r(mA1Db%vs,size(mA1Db%vs))
     call bcast_all_r(mA1Db%Qmu,size(mA1Db%Qmu))
     call bcast_all_r(mA1Db%Qkappa,size(mA1Db%Qkappa))
     call bcast_all_r(mA1Db%sprho,size(mA1Db%sprho))
     call bcast_all_r(mA1Db%spvp,size(mA1Db%spvp))
     call bcast_all_r(mA1Db%spvs,size(mA1Db%spvs))
     call bcast_all_r(mA1Db%spQmu,size(mA1Db%spQmu))
     call bcast_all_r(mA1Db%spQkappa,size(mA1Db%spQkappa))
  end if ! use_ASKI_background_model

  if(impose_ASKI_inverted_model) then

     call bcast_all_string_world(file_ASKI_inverted_model)

     if(myrank == 0) i_array_one_value(1) = model_ASKI_interpolation_type
     call bcast_all_i(i_array_one_value,1)
     if(myrank /= 0) model_ASKI_interpolation_type = i_array_one_value(1)

     select case(model_ASKI_interpolation_type)
     case( 2 )
        if(myrank == 0) r_array_one_value(1) = model_ASKI_factor_shepard_radius
        call bcast_all_r(r_array_one_value,1)
        if(myrank /= 0) model_ASKI_factor_shepard_radius = r_array_one_value(1)
     end select

     if(myrank == 0) i_array_one_value(1) = mAc%ncell
     call bcast_all_i(i_array_one_value,1)
     if(myrank /= 0) mAc%ncell = i_array_one_value(1)

     if(myrank == 0) i_array_one_value(1) = mAc%max_nnb
     call bcast_all_i(i_array_one_value,1)
     if(myrank /= 0) mAc%max_nnb = i_array_one_value(1)

     if(myrank .ne. 0) then
        allocate(mAc%cc(3,mAc%ncell),mAc%r(mAc%ncell),&
             mAc%nb(mAc%max_nnb+1,mAc%ncell))
     end if
     call bcast_all_r(mAc%cc,size(mAc%cc))
     call bcast_all_r(mAc%r,size(mAc%r))
     call bcast_all_i(mAc%nb,size(mAc%nb))

     if(myrank == 0) i_array_one_value(1) = model_ASKI_pmtrz
     call bcast_all_i(i_array_one_value,1)
     if(myrank /= 0) model_ASKI_pmtrz = i_array_one_value(1)

     select case (model_ASKI_pmtrz)
     case ( ipmtrz_isoLameSI )

        if(myrank == 0) r_array_one_value(1) = mAisoL%maxr_rho
        call bcast_all_r(r_array_one_value,1)
        if(myrank /= 0) mAisoL%maxr_rho = r_array_one_value(1)

        if(myrank == 0) r_array_one_value(1) = mAisoL%maxr_lambda
        call bcast_all_r(r_array_one_value,1)
        if(myrank /= 0) mAisoL%maxr_lambda = r_array_one_value(1)

        if(myrank == 0) r_array_one_value(1) = mAisoL%maxr_mu
        call bcast_all_r(r_array_one_value,1)
        if(myrank /= 0) mAisoL%maxr_mu = r_array_one_value(1)

        if(myrank == 0) i_array_one_value(1) = mAisoL%nval_rho
        call bcast_all_i(i_array_one_value,1)
        if(myrank /= 0) mAisoL%nval_rho = i_array_one_value(1)

        if(myrank == 0) i_array_one_value(1) = mAisoL%nval_lambda
        call bcast_all_i(i_array_one_value,1)
        if(myrank /= 0) mAisoL%nval_lambda = i_array_one_value(1)

        if(myrank == 0) i_array_one_value(1) = mAisoL%nval_mu
        call bcast_all_i(i_array_one_value,1)
        if(myrank /= 0) mAisoL%nval_mu = i_array_one_value(1)

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

     case ( ipmtrz_isoVelocitySI )

        if(myrank == 0) r_array_one_value(1) = mAisoV%maxr_rho
        call bcast_all_r(r_array_one_value,1)
        if(myrank /= 0) mAisoV%maxr_rho = r_array_one_value(1)

        if(myrank == 0) r_array_one_value(1) = mAisoV%maxr_vp
        call bcast_all_r(r_array_one_value,1)
        if(myrank /= 0) mAisoV%maxr_vp = r_array_one_value(1)

        if(myrank == 0) r_array_one_value(1) = mAisoV%maxr_vs
        call bcast_all_r(r_array_one_value,1)
        if(myrank /= 0) mAisoV%maxr_vs = r_array_one_value(1)

        if(myrank == 0) i_array_one_value(1) = mAisoV%nval_rho
        call bcast_all_i(i_array_one_value,1)
        if(myrank /= 0) mAisoV%nval_rho = i_array_one_value(1)

        if(myrank == 0) i_array_one_value(1) = mAisoV%nval_vp
        call bcast_all_i(i_array_one_value,1)
        if(myrank /= 0) mAisoV%nval_vp = i_array_one_value(1)

        if(myrank == 0) i_array_one_value(1) = mAisoV%nval_vs
        call bcast_all_i(i_array_one_value,1)
        if(myrank /= 0) mAisoV%nval_vs = i_array_one_value(1)

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

  end if ! impose_ASKI_inverted_model

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
  use constants,only: IMAIN,OUTPUT_FILES_BASE,IN_DATA_FILES
  
  implicit none
  

  integer :: ier,IOASKI
  character(len=800) :: error_message

  ier = -1

  ! read the values of variables use_ASKI_background_model, impose_ASKI_inverted_model
  call read_Par_file_ASKI()

  if(use_ASKI_background_model) call read_ASKI_external_background_model()
  if(impose_ASKI_inverted_model) call read_model_ASKI_kim_export()

  ! write logfile about this model to OUTPUT_FILES
  call get_file_unit_model_ASKI(IOASKI)
  open(unit=IOASKI,file=trim(OUTPUT_FILES_BASE)//'LOG_ASKI_model_external.txt',&
       form='formatted',status='unknown',action='write')
  if(.not.use_ASKI_background_model) then
     write(IOASKI,*) "according to flag 'USE_ASKI_BACKGROUND_MODEL' in file '"//trim(IN_DATA_FILES)//&
          "Par_file_ASKI' (or the fact that the file is missing), this simulation will NOT use an external "//&
          "ASKI background model"
     write(IMAIN,*) "according to flag 'USE_ASKI_BACKGROUND_MODEL' in file '"//trim(IN_DATA_FILES)//&
          "Par_file_ASKI' (or the fact that the file is missing), this simulation will NOT use an external "//&
          "ASKI background model"
  else
     write(IOASKI,*) "according to flag 'USE_ASKI_BACKGROUND_MODEL' in file '"//trim(IN_DATA_FILES)//&
          "Par_file_ASKI', this simulation will use an external ASKI background model given by file '"//&
          trim(IN_DATA_FILES)//trim(file_ASKI_background_model)//"'"
     write(IMAIN,*) "according to flag 'USE_ASKI_BACKGROUND_MODEL' in file '"//trim(IN_DATA_FILES)//&
          "Par_file_ASKI', this simulation will use an external ASKI background model given by file '"//&
          trim(IN_DATA_FILES)//trim(file_ASKI_background_model)//"'"
  end if
  if(.not.impose_ASKI_inverted_model) then
     write(IOASKI,*) "according to flag 'IMPOSE_ASKI_INVERTED_MODEL' in file '"//trim(IN_DATA_FILES)//&
          "Par_file_ASKI' (or the fact that the file is missing), this simulation will NOT impose an external "//&
          "ASKI inverted model"
     write(IMAIN,*) "according to flag 'IMPOSE_ASKI_INVERTED_MODEL' in file '"//trim(IN_DATA_FILES)//&
          "Par_file_ASKI' (or the fact that the file is missing), this simulation will NOT impose an external "//&
          "ASKI inverted model"
  else
     write(IOASKI,*) "according to flag 'IMPOSE_ASKI_INVERTED_MODEL' in file '"//trim(IN_DATA_FILES)//&
          "Par_file_ASKI', this simulation will impose an external ASKI inverted model given by file '"//&
          trim(IN_DATA_FILES)//trim(file_ASKI_inverted_model)//"'"
     write(IMAIN,*) "according to flag 'IMPOSE_ASKI_INVERTED_MODEL' in file '"//trim(IN_DATA_FILES)//&
          "Par_file_ASKI', this simulation will impose an external ASKI inverted model given by file '"//&
          trim(IN_DATA_FILES)//trim(file_ASKI_inverted_model)//"'"
!
     write(IMAIN,*) "successfully read model file '",trim(IN_DATA_FILES)//trim(file_ASKI_inverted_model),&
          "' containing the ASKI inverted model"
     select case (model_ASKI_interpolation_type)
     case(1)
        write(IMAIN,*) "will interpolate by standard shepard method"
     case(2)
        write(IMAIN,*) "will interpolate by shepard method applying an additional radius factor of ",&
             model_ASKI_factor_shepard_radius
     case default
        write(error_message,*) "in read_ASKI_external_inverted_model: model_ASKI_interpolation_type = ",&
             model_ASKI_interpolation_type,"; must be either 1 or 2, should have been checked before. ",&
             "THIS ERROR SHOULD NOT OCCURR, HENCE THIS MODULE IS INCONSISTENT!"
        call stop_error_model_ASKI(error_message)
     end select
     write(IMAIN,*) "ncell = ",mAc%ncell
     write(IMAIN,*) "maximal number of neighbours = ",mAc%max_nnb
     select case (model_ASKI_pmtrz)
     case ( ipmtrz_isoLameSI )
        write(IMAIN,*) "parametrization is isoLameSI"
        write(IMAIN,*) "nval_rho,nval_lambda,nval_mu = ",mAisoL%nval_rho,mAisoL%nval_lambda,mAisoL%nval_mu
        write(IMAIN,*) "maxr_rho,maxr_lambda,maxr_mu = ",mAisoL%maxr_rho,mAisoL%maxr_lambda,mAisoL%maxr_mu
     case ( ipmtrz_isoVelocitySI )
        write(IMAIN,*) "parametrization is isoVelocitySI"
        write(IMAIN,*) "nval_rho,nval_vp,nval_vs = ",mAisoV%nval_rho,mAisoV%nval_vp,mAisoV%nval_vs
        write(IMAIN,*) "maxr_rho,maxr_vp,maxr_vs = ",mAisoV%maxr_rho,mAisoV%maxr_vp,mAisoV%maxr_vs
     end select
  end if
  close(IOASKI)

  end subroutine read_external_model

!
!-------------------------------------------------------------------------------------------------
!
  subroutine read_ASKI_external_background_model()
  use model_ASKI
  use constants,only: IMAIN,IN_DATA_FILES
  implicit none
  integer :: IOASKI,ier
  character(len=800) :: error_message
  integer :: maxnnodes,ilayer,inode,i
  real, dimension(:), allocatable :: d,u,wrho,wvp,wvs,wQmu,wQkappa
  character(len=MAX_STRING_LEN) :: modelfile
  character(len=200) :: line
  logical :: file_exists
!
  modelfile = file_ASKI_background_model
!
  inquire(file=trim(IN_DATA_FILES)//trim(modelfile),exist=file_exists)
  if(.not.file_exists) then
     write(error_message,*) "in read_ASKI_external_background_model: file '"//trim(IN_DATA_FILES)//&
          trim(modelfile)//"' does not exist"
     call stop_error_model_ASKI(error_message)
  end if
!
  ier = -1
  call get_file_unit_model_ASKI(IOASKI)
  open(unit=IOASKI,file=trim(IN_DATA_FILES)//trim(modelfile), &
       status='unknown',action='read',iostat=ier)
  if(ier .ne. 0) then
     close(IOASKI) ! ?? necessary? sensible?
     ! write error message to file and stop
     write(error_message,*) "in read_ASKI_external_background_model: could not open file '"//trim(IN_DATA_FILES)//&
          trim(modelfile)//"' to read"
     call stop_error_model_ASKI(error_message)
  end if


  write(IMAIN,*) "opened model file '",trim(IN_DATA_FILES)//trim(modelfile),"' for layered spline gradients successfully"

  ! ignore first line!
  read(IOASKI,*) line
  write(IMAIN,*) 'ignoring first line (comment line)'

  read(IOASKI,*) mA1Db%zmax
  write(IMAIN,*) 'zmax = ',mA1Db%zmax

  read(IOASKI,*) mA1Db%nlayers
  write(IMAIN,*) 'nlayers = ',mA1Db%nlayers

  allocate(mA1Db%nnodes(mA1Db%nlayers))
  read(IOASKI,*) mA1Db%nnodes
  write(IMAIN,*) 'number of nodes for each layer = ',mA1Db%nnodes

  if(minval(mA1Db%nnodes) .le. 1) &
       call stop_error_model_ASKI('in read_ASKI_external_background_model: number of nodes in a layer must be at least 2')

  write(IMAIN,*) 'the model values are:'
  write(IMAIN,*) 'depth [m]   density [Kg/m^3]   vp [m/s]   vs [m/s]   Qmu   Qkappa'
  write(IMAIN,*) '********************************************************************'

  maxnnodes = maxval(mA1Db%nnodes)
  allocate(mA1Db%depth(mA1Db%nlayers,maxnnodes),mA1Db%rho(mA1Db%nlayers,maxnnodes), &
           mA1Db%vp(mA1Db%nlayers,maxnnodes),mA1Db%vs(mA1Db%nlayers,maxnnodes),mA1Db%Qmu(mA1Db%nlayers,maxnnodes),&
           mA1Db%Qkappa(mA1Db%nlayers,maxnnodes),&
           mA1Db%sprho(mA1Db%nlayers,maxnnodes),mA1Db%spvp(mA1Db%nlayers,maxnnodes),mA1Db%spvs(mA1Db%nlayers,maxnnodes),&
           mA1Db%spQmu(mA1Db%nlayers,maxnnodes),mA1Db%spQkappa(mA1Db%nlayers,maxnnodes))
  mA1Db%depth(:,:) = 0.
  mA1Db%rho(:,:) = 0.
  mA1Db%vp(:,:) = 0.
  mA1Db%vs(:,:) = 0.
  mA1Db%Qmu(:,:) = 0.
  mA1Db%Qkappa(:,:) = 0.
  mA1Db%sprho(:,:) = 0.
  mA1Db%spvp(:,:) = 0.
  mA1Db%spvs(:,:) = 0.
  mA1Db%spQmu(:,:) = 0.
  mA1Db%spQkappa(:,:) = 0.
  
  do ilayer = 1,mA1Db%nlayers
     do inode = 1,mA1Db%nnodes(ilayer)
        read(IOASKI,*) mA1Db%depth(ilayer,inode),mA1Db%rho(ilayer,inode),mA1Db%vp(ilayer,inode),mA1Db%vs(ilayer,inode),&
             mA1Db%Qmu(ilayer,inode),mA1Db%Qkappa(ilayer,inode)
           write(IMAIN,*) mA1Db%depth(ilayer,inode),mA1Db%rho(ilayer,inode),mA1Db%vp(ilayer,inode),mA1Db%vs(ilayer,inode),&
                mA1Db%Qmu(ilayer,inode),mA1Db%Qkappa(ilayer,inode)
     end do ! inode
  enddo ! ilayer

  close(IOASKI)

!!$  write(IMAIN,*) "Hello, this is read_external_model. initiating mA1Db now."
!!$
!!$  !find out how many layers
!!$  mA1Db%nlayers = 3
!!$
!!$  write(IMAIN,*) 'mA1Db%nlayers == ',mA1Db%nlayers
!!$
!!$  !find out how many nodes in each layer
!!$  allocate(mA1Db%nnodes(mA1Db%nlayers))
!!$  mA1Db%nnodes(1) = 6
!!$  mA1Db%nnodes(2) = 24
!!$  mA1Db%nnodes(3) = 1 ! 1 means halfspace
!!$
!!$  write(IMAIN,*) 'mA1Db%nnodes == ',mA1Db%nnodes
!!$
!!$  !find out maximum number of nodes in a layer (in this case 24)
!!$  maxnnodes = maxval(mA1Db%nnodes)
!!$
!!$  write(IMAIN,*) 'maxnnodes == ',maxnnodes
!!$
!!$  allocate(mA1Db%depth(mA1Db%nlayers,maxnnodes),mA1Db%rho(mA1Db%nlayers,maxnnodes), &
!!$           mA1Db%vp(mA1Db%nlayers,maxnnodes),mA1Db%vs(mA1Db%nlayers,maxnnodes),mA1Db%sprho(mA1Db%nlayers,maxnnodes), &
!!$           mA1Db%spvp(mA1Db%nlayers,maxnnodes),mA1Db%spvs(mA1Db%nlayers,maxnnodes))
!!$
!!$  mA1Db%depth(:,:) = 0.
!!$  mA1Db%rho(:,:) = 0.
!!$  mA1Db%vp(:,:) = 0.
!!$  mA1Db%vs(:,:) = 0.
!!$  mA1Db%sprho(:,:) = 0.
!!$  mA1Db%spvp(:,:) = 0.
!!$  mA1Db%spvs(:,:) = 0.
!!$
!!$  write(IMAIN,*) 'size(mA1Db%depth) == ',size(mA1Db%depth)
!!$  write(IMAIN,*) 'size(mA1Db%rho) == ',size(mA1Db%rho)
!!$  write(IMAIN,*) 'size(mA1Db%vp) == ',size(mA1Db%vp)
!!$  write(IMAIN,*) 'size(mA1Db%vs) == ',size(mA1Db%vs)
!!$
!!$  ! set parameters (read from files actually)
!!$  mA1Db%depth(1,1) =   0.0000 
!!$  mA1Db%depth(1,2) =   0.5329 
!!$  mA1Db%depth(1,3) =   1.0657 
!!$  mA1Db%depth(1,4) =   1.5986 
!!$  mA1Db%depth(1,5) =   2.1314 
!!$  mA1Db%depth(1,6) =   2.6643 
!!$  mA1Db%depth(2,1) =   2.6643 
!!$  mA1Db%depth(2,2) =   3.2607 
!!$  mA1Db%depth(2,3) =   3.8571 
!!$  mA1Db%depth(2,4) =   4.4535 
!!$  mA1Db%depth(2,5) =   5.0499 
!!$  mA1Db%depth(2,6) =   5.6463 
!!$  mA1Db%depth(2,7) =   6.2427 
!!$  mA1Db%depth(2,8) =   6.8391 
!!$  mA1Db%depth(2,9) =   7.4355 
!!$  mA1Db%depth(2,10) =  8.0319  
!!$  mA1Db%depth(2,11) =  8.6283  
!!$  mA1Db%depth(2,12) =  9.2247  
!!$  mA1Db%depth(2,13) =  9.8210  
!!$  mA1Db%depth(2,14) = 10.4174  
!!$  mA1Db%depth(2,15) = 11.0138  
!!$  mA1Db%depth(2,16) = 11.6102  
!!$  mA1Db%depth(2,17) = 12.2066  
!!$  mA1Db%depth(2,18) = 12.8030  
!!$  mA1Db%depth(2,19) = 13.3994  
!!$  mA1Db%depth(2,20) = 13.9958  
!!$  mA1Db%depth(2,21) = 14.5922  
!!$  mA1Db%depth(2,22) = 15.1886  
!!$  mA1Db%depth(2,23) = 15.7850  
!!$  mA1Db%depth(2,24) = 16.3814  
!!$  mA1Db%depth(3,1) =  16.3814  
!!$
!!$  mA1Db%rho(1,1) =  1600.
!!$  mA1Db%rho(1,2) =  1600.
!!$  mA1Db%rho(1,3) =  1600.
!!$  mA1Db%rho(1,4) =  1600.
!!$  mA1Db%rho(1,5) =  1600.
!!$  mA1Db%rho(1,6) =  1600.
!!$  mA1Db%rho(2,1) =  1600.
!!$  mA1Db%rho(2,2) =  1600.
!!$  mA1Db%rho(2,3) =  1600.
!!$  mA1Db%rho(2,4) =  1600.
!!$  mA1Db%rho(2,5) =  1600.
!!$  mA1Db%rho(2,6) =  1600.
!!$  mA1Db%rho(2,7) =  1600.
!!$  mA1Db%rho(2,8) =  1600.
!!$  mA1Db%rho(2,9) =  1600.
!!$  mA1Db%rho(2,10) = 1600.
!!$  mA1Db%rho(2,11) = 1600.
!!$  mA1Db%rho(2,12) = 1600
!!$  mA1Db%rho(2,13) = 1600. 
!!$  mA1Db%rho(2,14) = 1600. 
!!$  mA1Db%rho(2,15) = 1600. 
!!$  mA1Db%rho(2,16) = 1600. 
!!$  mA1Db%rho(2,17) = 1600. 
!!$  mA1Db%rho(2,18) = 1600. 
!!$  mA1Db%rho(2,19) = 1600. 
!!$  mA1Db%rho(2,20) = 1600.
!!$  mA1Db%rho(2,21) = 1600.
!!$  mA1Db%rho(2,22) = 1600.
!!$  mA1Db%rho(2,23) = 1600.
!!$  mA1Db%rho(2,24) = 1600.
!!$  mA1Db%rho(3,1) =  2300.
!!$                  
!!$  mA1Db%vp(1,1) =  131.2
!!$  mA1Db%vp(1,2) =  245.3
!!$  mA1Db%vp(1,3) =  337.8
!!$  mA1Db%vp(1,4) =  408.5
!!$  mA1Db%vp(1,5) =  457.5
!!$  mA1Db%vp(1,6) =  484.8
!!$  mA1Db%vp(2,1) =  483.0
!!$  mA1Db%vp(2,2) =  493.3
!!$  mA1Db%vp(2,3) =  503.6
!!$  mA1Db%vp(2,4) =  513.9
!!$  mA1Db%vp(2,5) =  524.2
!!$  mA1Db%vp(2,6) =  534.5
!!$  mA1Db%vp(2,7) =  544.8
!!$  mA1Db%vp(2,8) =  555.1
!!$  mA1Db%vp(2,9) =  565.4
!!$  mA1Db%vp(2,10) = 575.7 
!!$  mA1Db%vp(2,11) = 586.1 
!!$  mA1Db%vp(2,12) = 596.4 
!!$  mA1Db%vp(2,13) = 606.7 
!!$  mA1Db%vp(2,14) = 617.0 
!!$  mA1Db%vp(2,15) = 627.3 
!!$  mA1Db%vp(2,16) = 637.6 
!!$  mA1Db%vp(2,17) = 647.9 
!!$  mA1Db%vp(2,18) = 658.2 
!!$  mA1Db%vp(2,19) = 668.5 
!!$  mA1Db%vp(2,20) = 678.8 
!!$  mA1Db%vp(2,21) = 689.1 
!!$  mA1Db%vp(2,22) = 699.4 
!!$  mA1Db%vp(2,23) = 709.7 
!!$  mA1Db%vp(2,24) = 720.0 
!!$  mA1Db%vp(3,1) =  3702.4
!!$
!!$  mA1Db%vs(1,1) =  71.0
!!$  mA1Db%vs(1,2) =  156.3
!!$  mA1Db%vs(1,3) =  217.9
!!$  mA1Db%vs(1,4) =  255.8
!!$  mA1Db%vs(1,5) =  270.1
!!$  mA1Db%vs(1,6) =  260.8
!!$  mA1Db%vs(2,1) =  258.8
!!$  mA1Db%vs(2,2) =  267.6
!!$  mA1Db%vs(2,3) =  276.3
!!$  mA1Db%vs(2,4) =  285.0
!!$  mA1Db%vs(2,5) =  293.8
!!$  mA1Db%vs(2,6) =  302.5
!!$  mA1Db%vs(2,7) =  311.2
!!$  mA1Db%vs(2,8) =  320.0
!!$  mA1Db%vs(2,9) =  328.7
!!$  mA1Db%vs(2,10) = 337.4 
!!$  mA1Db%vs(2,11) = 346.2 
!!$  mA1Db%vs(2,12) = 354.9 
!!$  mA1Db%vs(2,13) = 363.6 
!!$  mA1Db%vs(2,14) = 372.4 
!!$  mA1Db%vs(2,15) = 381.1 
!!$  mA1Db%vs(2,16) = 389.8 
!!$  mA1Db%vs(2,17) = 398.6 
!!$  mA1Db%vs(2,18) = 407.3 
!!$  mA1Db%vs(2,19) = 416.1 
!!$  mA1Db%vs(2,20) = 424.8 
!!$  mA1Db%vs(2,21) = 433.5 
!!$  mA1Db%vs(2,22) = 442.3 
!!$  mA1Db%vs(2,23) = 451.0 
!!$  mA1Db%vs(2,24) = 459.7 
!!$  mA1Db%vs(3,1) =  2150.8

  ! now compute the splines for each layer (and each parameter) in form of values for sprho,spvp,spvs = p = s''(xj) 
  ! (procedure, as well as notation of p,s,xj as written in "Algorithms" by Robert Sedgewick, ADDISON-WESLEY 2002, Chapter 38)
  allocate(d(maxnnodes),u(maxnnodes),wrho(maxnnodes),wvp(maxnnodes),wvs(maxnnodes),wQmu(maxnnodes),wQkappa(maxnnodes))  ! local variables

  ! for each layer calculate the second derivative of the respective spline at all nodes
  do ilayer = 1,mA1Db%nlayers
  if(mA1Db%nnodes(ilayer) .ge. 3) then 
     ! in case of mA1Db%nnodes(ilayer) == 2, the spline interpolation (in that case linear interpolation) works, as sprho,spvp,spvs = 0. initially

     ! initiate temporary variables and calculate their values
     d(:) = 0.
     u(:) = 0.
     wrho(:) = 0.
     wvp(:) = 0.
     wvs(:) = 0.
     wQmu(:) = 0.
     wQkappa(:) = 0.
     do i = 2,mA1Db%nnodes(ilayer) - 1
        d(i) = 2.*(mA1Db%depth(ilayer,i+1)-mA1Db%depth(ilayer,i-1))
     end do

     do i = 1,mA1Db%nnodes(ilayer) - 1
        u(i) = mA1Db%depth(ilayer,i+1)-mA1Db%depth(ilayer,i)
     end do

     do i = 2,mA1Db%nnodes(ilayer) - 1
        wrho(i) = 6.*((mA1Db%rho(ilayer,i+1)-mA1Db%rho(ilayer,i))/u(i) - &
                      (mA1Db%rho(ilayer,i)-mA1Db%rho(ilayer,i-1))/u(i-1))
        wvp(i) = 6.*((mA1Db%vp(ilayer,i+1)-mA1Db%vp(ilayer,i))/u(i) - &
                      (mA1Db%vp(ilayer,i)-mA1Db%vp(ilayer,i-1))/u(i-1))
        wvs(i) = 6.*((mA1Db%vs(ilayer,i+1)-mA1Db%vs(ilayer,i))/u(i) - &
                      (mA1Db%vs(ilayer,i)-mA1Db%vs(ilayer,i-1))/u(i-1))
        wQmu(i) = 6.*((mA1Db%Qmu(ilayer,i+1)-mA1Db%Qmu(ilayer,i))/u(i) - &
                      (mA1Db%Qmu(ilayer,i)-mA1Db%Qmu(ilayer,i-1))/u(i-1))
        wQkappa(i) = 6.*((mA1Db%Qkappa(ilayer,i+1)-mA1Db%Qkappa(ilayer,i))/u(i) - &
                      (mA1Db%Qkappa(ilayer,i)-mA1Db%Qkappa(ilayer,i-1))/u(i-1))
     end do

     ! now calculate the second derivatives of the spline, assuming them being zero at the extremal nodes (natural boundary conditions)
     mA1Db%sprho(ilayer,1) = 0.; mA1Db%sprho(ilayer,mA1Db%nnodes(ilayer)) = 0.
     mA1Db%spvp(ilayer,1) = 0.; mA1Db%spvp(ilayer,mA1Db%nnodes(ilayer)) = 0.
     mA1Db%spvs(ilayer,1) = 0.; mA1Db%spvs(ilayer,mA1Db%nnodes(ilayer)) = 0.
     mA1Db%spQmu(ilayer,1) = 0.; mA1Db%spQmu(ilayer,mA1Db%nnodes(ilayer)) = 0.
     mA1Db%spQkappa(ilayer,1) = 0.; mA1Db%spQkappa(ilayer,mA1Db%nnodes(ilayer)) = 0.

     ! then calculate the others by solving a tridiagonal system of equations
     if(mA1Db%nnodes(ilayer) > 3) then
        do i = 2,mA1Db%nnodes(ilayer) - 2
           wrho(i+1) = wrho(i+1) - wrho(i)*u(i)/d(i)
           wvp(i+1) = wvp(i+1) - wvp(i)*u(i)/d(i)
           wvs(i+1) = wvs(i+1) - wvs(i)*u(i)/d(i)
           wQmu(i+1) = wQmu(i+1) - wQmu(i)*u(i)/d(i)
           wQkappa(i+1) = wQkappa(i+1) - wQkappa(i)*u(i)/d(i)
           d(i+1) = d(i+1) - (u(i)**2)/d(i)
        end do
     endif

     do i = mA1Db%nnodes(ilayer)-1,2,-1
        mA1Db%sprho(ilayer,i) = (wrho(i) - u(i)*mA1Db%sprho(ilayer,i+1))/d(i)
        mA1Db%spvp(ilayer,i) = (wvp(i) - u(i)*mA1Db%spvp(ilayer,i+1))/d(i)
        mA1Db%spvs(ilayer,i) = (wvs(i) - u(i)*mA1Db%spvs(ilayer,i+1))/d(i)
        mA1Db%spQmu(ilayer,i) = (wQmu(i) - u(i)*mA1Db%spQmu(ilayer,i+1))/d(i)
        mA1Db%spQkappa(ilayer,i) = (wQkappa(i) - u(i)*mA1Db%spQkappa(ilayer,i+1))/d(i)
     end do

   end if
   end do ! ilayer

   deallocate(d,u,wrho,wvp,wvs,wQmu,wQkappa)
   ! done calculating splines

  end subroutine read_ASKI_external_background_model

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_model_ASKI_kim_export()

  use model_ASKI
  use constants,only: IN_DATA_FILES
  
  implicit none

  character(len=MAX_STRING_LEN) :: modelfile

  integer :: ier,IOASKI,iline
  character(len=800) :: error_message

!
!NOW READ KIM EXPORT FILE
!
  modelfile = file_ASKI_inverted_model

  call get_file_unit_model_ASKI(IOASKI)
  open(unit=IOASKI,file=trim(IN_DATA_FILES)//trim(modelfile),form='formatted', &
       status='old',action='read',iostat=ier)
  if(ier .ne. 0) then
     close(IOASKI) ! ?? necessary? sensible?
     ! write error message to file and stop
     write(error_message,*) "in read_model_ASKI_kim_export: could not open file '"//trim(IN_DATA_FILES)//&
          trim(modelfile)//"' to read"
     call stop_error_model_ASKI(error_message)
  end if

  iline = 0

  call check_pmtrz_ASKI_kim_export(IOASKI,iline,trim(IN_DATA_FILES)//trim(modelfile))

  call read_model_ASKI_kim_export_cells(IOASKI,iline,trim(IN_DATA_FILES)//trim(modelfile))

  select case (model_ASKI_pmtrz)
     case ( ipmtrz_isoLameSI )
        call read_external_model_ASKI_kim_export_isoLame(IOASKI,iline,trim(IN_DATA_FILES)//trim(modelfile))

     case ( ipmtrz_isoVelocitySI )
        call read_external_model_ASKI_kim_export_isoVelocity(IOASKI,iline,trim(IN_DATA_FILES)//trim(modelfile))

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
  character(len=15) :: pmtrz
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
  case('isoLameSI','isoVelocitySI') ! ok, do nothing
  case default
     close(IOASKI)
     ! write error message to file and stop
     write(error_message,*) "in check_pmtrz_ASKI_kim_export: model parametrization '"//&
          trim(pmtrz)//"' (on line ",iline," of file '"//trim(filename)//"')"//&
          " is not supported, only 'isoLameSI','isoVelocitySI' supported so far"
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
  case('isoLameSI')
     model_ASKI_pmtrz = ipmtrz_isoLameSI
     if(nparam/=3 .or. .not.(any(param == 'rho') .and. any(param == 'lambda') .and. any(param == 'mu') )) then
        close(IOASKI)
        ! write error message to file and stop
        write(error_message,*) "in check_pmtrz_ASKI_kim_export: model parametrization isoLameSI "//&
             "expects the 3 parameter names (on line ",iline," of file '"//trim(filename)//&
             "'), namely 'rho','lambda','mu' (in any order)"
        call stop_error_model_ASKI(error_message)
     end if
  case('isoVelocitySI')
     model_ASKI_pmtrz = ipmtrz_isoVelocitySI
     if(nparam/=3 .or. .not.(any(param == 'rho') .and. any(param == 'vp') .and. any(param == 'vs') )) then
        close(IOASKI)
        ! write error message to file and stop
        write(error_message,*) "in check_pmtrz_ASKI_kim_export: model parametrization isoVelocitySI "//&
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

  ! acoustic/elastic/.. domain flag ( 1 = acoustic / 2 = elastic / 3 = poroelastic )
  integer :: idomain_id

  ! local parameters
  real :: x,y,z

!---
!
! ADD YOUR MODEL HERE
!
!---

  if(.not.(use_ASKI_background_model .or. impose_ASKI_inverted_model)) return

  ! GLL point location converted to real
  x = xmesh
  y = ymesh
  z = zmesh

  if(use_ASKI_background_model) then
     ! get the background velocities vp_real,vs_real, and rho_real, Qmu_real, Qkappa_real
     call model_1Dbackground_values_ASKI(z,rho,vp,vs,qmu_atten,qkappa_atten)
  end if

  if(impose_ASKI_inverted_model) then
     ! no attenuation treated so far for inverted models, so simply leave values Qkappa,Qmu untouched here

     !! properly smooth out to the outside (if you are outside the reach of a cell center, check if
     !! there is a cell center in reach for twice the shepard radius, if so use all in reach and smooth out 
     !! linearly to the background (where for distance = twice the radius, the background model has full 
     !! contribution, and for distance = 1 shepard radius, the background model has no contribution)
     select case (model_ASKI_pmtrz)
     case ( ipmtrz_isoLameSI )
        ! no attenuation, so simply leave values qkappa_atten,qmu_atten untouched
        call model_external_values_ASKI_isoLame(x,y,z,rho,vp,vs)
     case ( ipmtrz_isoVelocitySI )
        ! no attenuation, so simply leave values qkappa_atten,qmu_atten untouched
        call model_external_values_ASKI_isoVelocity(x,y,z,rho,vp,vs)
     end select
  end if ! impose_ASKI_inverted_model

  end subroutine model_external_values

!
!--------------------------------------------------------------------------------------------
!

  subroutine model_1Dbackground_values_ASKI(z,rho,vp,vs,qmu_atten,qkappa_atten)

  use model_ASKI

  implicit none

! given a GLL point, do spline interpolation at depth zmax - z

  ! density, Vp and Vs
  real :: vp,vs,rho

  ! attenuation
  real :: qmu_atten,qkappa_atten

  ! imaterial_id: (something like REGION_CODE), associated material flag, could be used as a dummy variable
  ! to indicate on which side of a discontinuity this point is
!  integer :: imaterial_id

  real :: z

  logical :: values_defined

  character(len=400) :: errmsg

  integer :: inode,ilayer
!  real :: layerthickness
  real :: depth,t

  ! FOR TEST/DEBUGGING OUTPUT ONLY:
  !real :: depth2,vp2,vs2,rho2,t2 
  !character(len=250) :: filename
  !integer :: iz


  depth = mA1Db%zmax - z  ! [m]

  ! CHECK IF THE REQUESTED DEPTH IS ABOVE THE VERY BOTTOM OF THIS MODEL DEFINITION AND BELOW THE UPPERMOST NODE
  ! IF NOT, RETURN WITHOUT ADDING BACKGROUND MODEL, LEAVING THE DEFAULT MODEL VALUES UNTOUCHED
  if(depth > mA1Db%depth(mA1Db%nlayers,mA1Db%nnodes(mA1Db%nlayers)) .or. &
       depth < mA1Db%depth(1,1) ) then
     return
  end if

!!$  layerthickness = mA1Db%depth(imaterial_id,mA1Db%nnodes(imaterial_id)) - mA1Db%depth(imaterial_id,1)
!!$
!!$  if( (depth < mA1Db%depth(imaterial_id,1) - 0.001*layerthickness) .or. &
!!$      (depth > mA1Db%depth(imaterial_id,mA1Db%nnodes(imaterial_id)) + 0.001*layerthickness) ) then
!!$     write(*,*) 'in model_external_values_layered_spline_gradients: at this depth, imaterial_id is wrong!   x,y,z=',x,y,z, &
!!$                ';  depth=',depth,';  imaterial_id=',imaterial_id,';  layerthickness=',layerthickness
!!$     call MPI_ABORT(MPI_COMM_WORLD, 30, ier)
!!$     stop 'in model_external_values_layered_spline_gradients: depth and material_id do not fit the layered gradients model!'
!!$  endif

  ! FIRST FIND OUT THE LAYER IN WHICH THE CURRENT POINT IS, I.E. WITHIN WHICH SPLINE INTERPOLATION SHOULD BE CONDUCTED
  do ilayer = 1,mA1Db%nlayers
     if(depth <= mA1Db%depth(ilayer,mA1Db%nnodes(ilayer))) exit
  end do
  ! after this loop, ilayer should be the index of the layer which contains the current point (each layer contains its bottom depth, the first layer contains the surface of the earth)

  values_defined = .false.
  do inode = 2,mA1Db%nnodes(ilayer)
     if(depth <= mA1Db%depth(ilayer,inode)) then
        ! interpolate values at current depth
        t = (depth - mA1Db%depth(ilayer,inode-1)) / &
            (mA1Db%depth(ilayer,inode) - mA1Db%depth(ilayer,inode-1))
        rho = t*mA1Db%rho(ilayer,inode) + (1.-t)*mA1Db%rho(ilayer,inode-1) + &
              (mA1Db%depth(ilayer,inode)-mA1Db%depth(ilayer,inode-1))**2 * &
              ((t**3-t)*mA1Db%sprho(ilayer,inode) + ((1-t)**3-(1-t))*mA1Db%sprho(ilayer,inode-1))/6.
        vp = t*mA1Db%vp(ilayer,inode) + (1.-t)*mA1Db%vp(ilayer,inode-1) + &
             (mA1Db%depth(ilayer,inode)-mA1Db%depth(ilayer,inode-1))**2 * &
             ((t**3-t)*mA1Db%spvp(ilayer,inode) + ((1-t)**3-(1-t))*mA1Db%spvp(ilayer,inode-1))/6.
        vs = t*mA1Db%vs(ilayer,inode) + (1.-t)*mA1Db%vs(ilayer,inode-1) + &
             (mA1Db%depth(ilayer,inode)-mA1Db%depth(ilayer,inode-1))**2 * &
             ((t**3-t)*mA1Db%spvs(ilayer,inode) + ((1-t)**3-(1-t))*mA1Db%spvs(ilayer,inode-1))/6.
        qmu_atten = t*mA1Db%Qmu(ilayer,inode) + (1.-t)*mA1Db%Qmu(ilayer,inode-1) + &
             (mA1Db%depth(ilayer,inode)-mA1Db%depth(ilayer,inode-1))**2 * &
             ((t**3-t)*mA1Db%spQmu(ilayer,inode) + ((1-t)**3-(1-t))*mA1Db%spQmu(ilayer,inode-1))/6.
        qkappa_atten = t*mA1Db%Qkappa(ilayer,inode) + (1.-t)*mA1Db%Qkappa(ilayer,inode-1) + &
             (mA1Db%depth(ilayer,inode)-mA1Db%depth(ilayer,inode-1))**2 * &
             ((t**3-t)*mA1Db%spQkappa(ilayer,inode) + ((1-t)**3-(1-t))*mA1Db%spQkappa(ilayer,inode-1))/6.
        ! position between the nodes is found, so leave the do loop
        values_defined = .true.
        exit
     end if
  end do ! inode

  ! after this loop, there should have been found an interpolation
  ! if not, raise an error
  if(.not.values_defined) then
     write(errmsg,*) "in routine model_1Dbackground_values_ASKI: at depth ",depth,&
          " there was no model value defined, although it should have been. This routine is erroneous"
     call stop_error_model_ASKI(errmsg)
  end if

!! FS FS ##################################  TEST INTERPOLATION OF MODEL  #####################################
!!$if(testlayers(1) == 0) then
!!$   write(filename,"('/data/Kernel/specfem3D/OUTPUT_FILES/test_model_spline_proc',i3.3,'.dat')") mA1Db%myrank
!!$   open(unit=30,file=trim(filename),status='unknown',action='write')
!!$   do iz=0,2800
!!$      depth2 = real(iz)*0.01
!!$      ! find out layer, in which testpoint lies
!!$      do ilayer = 1,mA1Db%nlayers
!!$         layerthickness = mA1Db%depth(ilayer,mA1Db%nnodes(ilayer)) - mA1Db%depth(ilayer,1)
!!$         if(depth2 < mA1Db%depth(ilayer,mA1Db%nnodes(ilayer)) + 0.001*layerthickness) exit
!!$      end do ! ilayer
!!$      do inode = 2,mA1Db%nnodes(ilayer)
!!$         if(depth2 < mA1Db%depth(ilayer,inode) + 0.001*layerthickness) then
!!$            ! interpolate values at current depth
!!$            t2 = (depth2 - mA1Db%depth(ilayer,inode-1)) / &
!!$                 (mA1Db%depth(ilayer,inode) - mA1Db%depth(ilayer,inode-1))
!!$            rho2 = t2*mA1Db%rho(ilayer,inode) + (1.-t2)*mA1Db%rho(ilayer,inode-1) + &
!!$                 (mA1Db%depth(ilayer,inode)-mA1Db%depth(ilayer,inode-1))**2 * &
!!$                 ((t2**3-t2)*mA1Db%sprho(ilayer,inode) + ((1-t2)**3-(1-t2))*mA1Db%sprho(ilayer,inode-1))/6.
!!$            vp2 = t2*mA1Db%vp(ilayer,inode) + (1.-t2)*mA1Db%vp(ilayer,inode-1) + &
!!$                 (mA1Db%depth(ilayer,inode)-mA1Db%depth(ilayer,inode-1))**2 * &
!!$                 ((t2**3-t2)*mA1Db%spvp(ilayer,inode) + ((1-t2)**3-(1-t2))*mA1Db%spvp(ilayer,inode-1))/6.
!!$            vs2 = t2*mA1Db%vs(ilayer,inode) + (1.-t2)*mA1Db%vs(ilayer,inode-1) + &
!!$                 (mA1Db%depth(ilayer,inode)-mA1Db%depth(ilayer,inode-1))**2 * &
!!$                 ((t2**3-t2)*mA1Db%spvs(ilayer,inode) + ((1-t2)**3-(1-t2))*mA1Db%spvs(ilayer,inode-1))/6.
!!$
!!$            write(30,*) depth2,rho2,vp2,vs2
!!$            ! position between the nodes is found, so leave the do loop
!!$            exit
!!$         end if
!!$      end do ! inode
!!$   end do ! iz
!!$   close(30)
!!$   testlayers(1) = 1
!!$endif
!! FS FS ##################################  TEST INTERPOLATION OF MODEL  #####################################

!!$if(imaterial_id == 1) then
!!$   if(testlayers(1) == 0) then
!!$      if((depth < maxval(mA1Db%depth(1,:)) + 0.0001) .and. (depth > minval(mA1Db%depth(1,:)) - 0.0001)) then
!!$         write(*,*) 'this is model_external_values, myrank == ',myrank,'. imaterial_id == 1. depth is OK' 
!!$      else
!!$         write(*,*) 'this is model_external_values, myrank == ',myrank,'. imaterial_id == 1. WRONG! depth ==',depth
!!$      endif
!!$      testlayers(1) = 1
!!$   endif
!!$elseif(imaterial_id == 2) then
!!$   if(testlayers(2) == 0) then
!!$      if((depth < maxval(mA1Db%depth(2,:)) + 0.0001) .and. (depth > minval(mA1Db%depth(2,:)) - 0.0001)) then
!!$         write(*,*) 'this is model_external_values, myrank == ',myrank,'. imaterial_id == 2. depth is OK' 
!!$      else
!!$         write(*,*) 'this is model_external_values, myrank == ',myrank,'. imaterial_id == 2. WRONG! depth ==',depth
!!$      endif
!!$      testlayers(2) = 1
!!$   endif
!!$elseif(imaterial_id == 3) then
!!$   if(testlayers(3) == 0) then
!!$      if(depth > minval(mA1Db%depth(3,:)) - 0.0001) then
!!$         write(*,*) 'this is model_external_values, myrank == ',myrank,'. imaterial_id == 3. depth is OK' 
!!$      else
!!$         write(*,*) 'this is model_external_values, myrank == ',myrank,'. imaterial_id == 3. WRONG! depth ==',depth
!!$      endif
!!$      testlayers(3) = 1
!!$   endif
!!$endif

  end subroutine model_1Dbackground_values_ASKI
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
  real :: rho_bg,lambda_bg,mu_bg
  real :: rho_interpolated,lambda_interpolated,mu_interpolated
  real :: factor_shepard_radius,enlarge_influence,f
  integer :: n
  integer, dimension(:), pointer :: interpolation_index
  real, dimension(:), pointer :: interpolation_weight

  ! compute lambda,mu,(rho) from the incoming background velocity model
  ! whenever values are missing below, use the background values
  mu_bg = vs*vs*rho
  lambda_bg = vp*vp*rho-2.d0*mu_bg
  rho_bg = rho

  nullify(interpolation_index,interpolation_weight)

  select case(model_ASKI_interpolation_type)
  case(1); factor_shepard_radius = 2.
  case(2); factor_shepard_radius = model_ASKI_factor_shepard_radius 
  end select

  ! deal with rho
  ! first assign the background value as the interpolated one (in case that no if-clause below is entered)
  rho_interpolated = rho_bg
  if(mAisoL%nval_rho > 0) then
     call shepard_interpolation_model_ASKI(x,y,z,factor_shepard_radius,mAisoL%nval_rho,mAisoL%idx_rho,&
          n,interpolation_index,interpolation_weight)
     if(n>0) then
        ! assume that in this case, interpolation_index,interpolation_weight are associated and have size n
        rho_interpolated = sum( interpolation_weight * mAisoL%rho(interpolation_index) )
        deallocate(interpolation_index,interpolation_weight)
     else
        ! check if there is a cell center in reach for twice the radius. If so, linearly interpolate
        ! with the background model: 
        !   get the fractional distance f*shepard_radius of the current point to the closest cell ( 1 <= f <= 2)
        !   do a linear interpolation with the background model, using f-1  (0 <= f-1 <= 1)
        !   assign the return model values as:  (1 - (f-1) )*shepard_interpolation + (f-1)*background
        !   i.e. for distance = twice the radius (f=2), the background model has full 
        !   contribution, and for distance = 1 shepard radius (f=1), the background model has no contribution at all

        enlarge_influence = 2.0
        call shepard_interpolation_model_ASKI(x,y,z,factor_shepard_radius,mAisoL%nval_rho,mAisoL%idx_rho,&
             n,interpolation_index,interpolation_weight,&
             enlarge_influence=enlarge_influence,radius_factor_closest_cell=f)
        if(n>0) then
           if(f<2.0) rho_interpolated &
                = (1 - (f-1) )*sum( interpolation_weight * mAisoL%rho(interpolation_index) ) &
                + (f-1)*rho_bg
           deallocate(interpolation_index,interpolation_weight)
        end if
     end if
  end if

  ! deal with lambda
  ! first assign the background value as the interpolated one (in case that no if-clause below is entered)
  lambda_interpolated = lambda_bg
  if(mAisoL%nval_lambda > 0) then
     call shepard_interpolation_model_ASKI(x,y,z,factor_shepard_radius,mAisoL%nval_lambda,mAisoL%idx_lambda,&
          n,interpolation_index,interpolation_weight)
     if(n>0) then
        ! assume that in this case, interpolation_index,interpolation_weight are associated and have size n
        lambda_interpolated = sum( interpolation_weight * mAisoL%lambda(interpolation_index) )
        deallocate(interpolation_index,interpolation_weight)
     else
        ! check if there is a cell center in reach for twice the radius. If so, linearly interpolate
        ! with the background model: 
        !   get the fractional distance f*shepard_radius of the current point to the closest cell ( 1 <= f <= 2)
        !   do a linear interpolation with the background model, using f-1  (0 <= f-1 <= 1)
        !   assign the return model values as:  (1 - (f-1) )*shepard_interpolation + (f-1)*background
        !   i.e. for distance = twice the radius (f=2), the background model has full 
        !   contribution, and for distance = 1 shepard radius (f=1), the background model has no contribution at all

        enlarge_influence = 2.0
        call shepard_interpolation_model_ASKI(x,y,z,factor_shepard_radius,mAisoL%nval_lambda,mAisoL%idx_lambda,&
             n,interpolation_index,interpolation_weight,&
             enlarge_influence=enlarge_influence,radius_factor_closest_cell=f)
        if(n>0) then
           if(f<2.0) lambda_interpolated &
                = (1 - (f-1) )*sum( interpolation_weight * mAisoL%lambda(interpolation_index) ) &
                + (f-1)*lambda_bg
           deallocate(interpolation_index,interpolation_weight)
        end if
     end if
  end if

  ! deal with mu
  ! first assign the background value as the interpolated one (in case that no if-clause below is entered)
  mu_interpolated = mu_bg
  if(mAisoL%nval_mu > 0) then
     call shepard_interpolation_model_ASKI(x,y,z,factor_shepard_radius,mAisoL%nval_mu,mAisoL%idx_mu,&
          n,interpolation_index,interpolation_weight)
     if(n>0) then
        ! assume that in this case, interpolation_index,interpolation_weight are associated and have size n
        mu_interpolated = sum( interpolation_weight * mAisoL%mu(interpolation_index) )
        deallocate(interpolation_index,interpolation_weight)
     else
        ! check if there is a cell center in reach for twice the radius. If so, linearly interpolate
        ! with the background model: 
        !   get the fractional distance f*shepard_radius of the current point to the closest cell ( 1 <= f <= 2)
        !   do a linear interpolation with the background model, using f-1  (0 <= f-1 <= 1)
        !   assign the return model values as:  (1 - (f-1) )*shepard_interpolation + (f-1)*background
        !   i.e. for distance = twice the radius (f=2), the background model has full 
        !   contribution, and for distance = 1 shepard radius (f=1), the background model has no contribution at all

        enlarge_influence = 2.0
        call shepard_interpolation_model_ASKI(x,y,z,factor_shepard_radius,mAisoL%nval_mu,mAisoL%idx_mu,&
             n,interpolation_index,interpolation_weight,&
             enlarge_influence=enlarge_influence,radius_factor_closest_cell=f)
        if(n>0) then
           if(f<2.0) mu_interpolated = (1 - (f-1) )*sum( interpolation_weight * mAisoL%mu(interpolation_index) ) &
                + (f-1)*mu_bg
           deallocate(interpolation_index,interpolation_weight)
        end if
     end if
  end if

  ! compute vp,vs,(rho) from the above interpolated isoLame values
  vp = sqrt( (lambda_interpolated+2.*mu_interpolated) / rho_interpolated )
  vs = sqrt( (mu_interpolated) / rho_interpolated )
  rho = rho_interpolated

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
  real :: factor_shepard_radius,enlarge_influence,f
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
     else
        ! check if there is a cell center in reach for twice the radius. If so, linearly interpolate
        ! with the background model: 
        !   get the fractional distance f*shepard_radius of the current point to the closest cell ( 1 <= f <= 2)
        !   do a linear interpolation with the background model, using f-1  (0 <= f-1 <= 1)
        !   assign the return model values as:  (1 - (f-1) )*shepard_interpolation + (f-1)*background
        !   i.e. for distance = twice the radius (f=2), the background model has full 
        !   contribution, and for distance = 1 shepard radius (f=1), the background model has no contribution at all

        enlarge_influence = 2.0
        call shepard_interpolation_model_ASKI(x,y,z,factor_shepard_radius,mAisoV%nval_rho,mAisoV%idx_rho,&
             n,interpolation_index,interpolation_weight,&
             enlarge_influence=enlarge_influence,radius_factor_closest_cell=f)
        if(n>0) then
           if(f<2.0) rho = (1 - (f-1) )*sum( dble(interpolation_weight) * mAisoV%rho(interpolation_index) ) &
                + (f-1)*rho
           deallocate(interpolation_index,interpolation_weight)
        else
           ! do nothing, leave incoming model values as they are (incoming values contain background model)
        end if
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
     else
        ! check if there is a cell center in reach for twice the radius. If so, linearly interpolate
        ! with the background model: 
        !   get the fractional distance f*shepard_radius of the current point to the closest cell ( 1 <= f <= 2)
        !   do a linear interpolation with the background model, using f-1  (0 <= f-1 <= 1)
        !   assign the return model values as:  (1 - (f-1) )*shepard_interpolation + (f-1)*background
        !   i.e. for distance = twice the radius (f=2), the background model has full 
        !   contribution, and for distance = 1 shepard radius (f=1), the background model has no contribution at all

        enlarge_influence = 2.0
        call shepard_interpolation_model_ASKI(x,y,z,factor_shepard_radius,mAisoV%nval_vp,mAisoV%idx_vp,&
             n,interpolation_index,interpolation_weight,&
             enlarge_influence=enlarge_influence,radius_factor_closest_cell=f)
        if(n>0) then
           if(f<2.0) vp = (1 - (f-1) )*sum( dble(interpolation_weight) * mAisoV%vp(interpolation_index) ) &
                + (f-1)*vp
           deallocate(interpolation_index,interpolation_weight)
        else
           ! do nothing, leave incoming model values as they are (incoming values contain background model)
        end if
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
     else
        ! check if there is a cell center in reach for twice the radius. If so, linearly interpolate
        ! with the background model: 
        !   get the fractional distance f*shepard_radius of the current point to the closest cell ( 1 <= f <= 2)
        !   do a linear interpolation with the background model, using f-1  (0 <= f-1 <= 1)
        !   assign the return model values as:  (1 - (f-1) )*shepard_interpolation + (f-1)*background
        !   i.e. for distance = twice the radius (f=2), the background model has full 
        !   contribution, and for distance = 1 shepard radius (f=1), the background model has no contribution at all

        enlarge_influence = 2.0
        call shepard_interpolation_model_ASKI(x,y,z,factor_shepard_radius,mAisoV%nval_vs,mAisoV%idx_vs,&
             n,interpolation_index,interpolation_weight,&
             enlarge_influence=enlarge_influence,radius_factor_closest_cell=f)
        if(n>0) then
           if(f<2.0) vs = (1 - (f-1) )*sum( dble(interpolation_weight) * mAisoV%vs(interpolation_index) ) &
                + (f-1)*vs
           deallocate(interpolation_index,interpolation_weight)
        else
           ! do nothing, leave incoming model values as they are (incoming values contain background model)
        end if
     end if
  end if

  end subroutine model_external_values_ASKI_isoVelocity
!
!---------------------------------------------------------------------------------
!
subroutine read_Par_file_ASKI()

  use model_ASKI
  use constants,only: IN_DATA_FILES
  
  implicit none

  character(len=500), dimension(:), allocatable :: val_parfile
  character(len=100), dimension(:), allocatable :: key_parfile
  character(len=601) :: line
  character(len=500) :: val,errstr
  integer :: npar,ios,IOASKI,eqindx

  ! open Par_file_ASKI and find number of valid lines
  call get_file_unit_model_ASKI(IOASKI)
  open(unit=IOASKI,file=trim(IN_DATA_FILES)//'Par_file_ASKI',&
       form='formatted',status='old',action='read',iostat=ios)
  if(ios/=0) then
     close(IOASKI)
     ! if there is no file Par_file_ASKI:
     call stop_error_model_ASKI("ERROR: could not open file '"//trim(IN_DATA_FILES)//&
          "Par_file_ASKI' for external model specifications, but MODEL = external is requested in DATA/Par_file")
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

  if(npar == 0) call stop_error_model_ASKI("no valid lines in file '"//trim(IN_DATA_FILES)//"Par_file_ASKI'")
  allocate(key_parfile(npar),val_parfile(npar))

  ! now open again and store key,val pairs of valid lines
  call get_file_unit_model_ASKI(IOASKI)
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

  ! now set values of variables use_ASKI_background_model, impose_ASKI_inverted_model

  ! USE_ASKI_BACKGROUND_MODEL
  call get_value_Par_file_ASKI('USE_ASKI_BACKGROUND_MODEL',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) use_ASKI_background_model
  if(ios/=0) call stop_error_model_ASKI("invalid logical value for parameter 'USE_ASKI_BACKGROUND_MODEL' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  if(use_ASKI_background_model) then
     ! FILE_ASKI_BACKGROUND_MODEL
     call get_value_Par_file_ASKI('FILE_ASKI_BACKGROUND_MODEL',val,key_parfile,val_parfile,npar)
     read(val,*,iostat=ios) file_ASKI_background_model
     if(ios/=0) call stop_error_model_ASKI("invalid character value for parameter 'FILE_ASKI_BACKGROUND_MODEL' in '"&
          //trim(IN_DATA_FILES)//"Par_file_ASKI'")
  end if
     
  ! IMPOSE_ASKI_INVERTED_MODEL
  call get_value_Par_file_ASKI('IMPOSE_ASKI_INVERTED_MODEL',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) impose_ASKI_inverted_model
  if(ios/=0) call stop_error_model_ASKI("invalid logical value for parameter 'IMPOSE_ASKI_INVERTED_MODEL' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  if(impose_ASKI_inverted_model) then
     ! FILE_ASKI_INVERTED_MODEL
     call get_value_Par_file_ASKI('FILE_ASKI_INVERTED_MODEL',val,key_parfile,val_parfile,npar)
     read(val,*,iostat=ios) file_ASKI_inverted_model
     if(ios/=0) call stop_error_model_ASKI("invalid character value for parameter 'FILE_ASKI_INVERTED_MODEL' in '"&
          //trim(IN_DATA_FILES)//"Par_file_ASKI'")

     call get_value_Par_file_ASKI('ASKI_INVERTED_MODEL_INTERPOLATION_TYPE',val,key_parfile,val_parfile,npar)
     select case (val)
     case ('shepard_standard')
        model_aski_interpolation_type = 1
     case('shepard_factor_radius')
        model_ASKI_interpolation_type = 2
     case default
        call stop_error_model_ASKI("value '"//trim(val)//"' of parameter 'ASKI_INVERTED_MODEL_INTERPOLATION_TYPE' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI' not supported: must be one of 'shepard_standard', 'shepard_factor_radius'")
     end select
  end if 
  if(model_ASKI_interpolation_type == 2) then
     call get_value_Par_file_ASKI('ASKI_INVERTED_MODEL_FACTOR_SHEPARD_RADIUS',val,key_parfile,val_parfile,npar)
     read(val,*,iostat=ios) model_ASKI_factor_shepard_radius
     if(ios/=0) call stop_error_model_ASKI("invalid real value for parameter 'ASKI_INVERTED_MODEL_FACTOR_SHEPARD_RADIUS'"//&
          " in '"//trim(IN_DATA_FILES)//"Par_file_ASKI'")
     if(model_ASKI_factor_shepard_radius <= 0.0) then
        write(errstr,*) "value ",model_ASKI_factor_shepard_radius," of parameter 'ASKI_INVERTED_MODEL_FACTOR_SHEPARD_RADIUS'"//&
             " in '"//trim(IN_DATA_FILES)//"Par_file_ASKI' is invalid: must be strictly positve"
        call stop_error_model_ASKI(errstr)
     end if
  end if

  if(allocated(key_parfile)) deallocate(key_parfile)
  if(allocated(val_parfile)) deallocate(val_parfile)
end subroutine read_Par_file_ASKI

!
! ----------------------------------------------------------------------------------------------------------
!
subroutine get_value_Par_file_ASKI(key,val,key_parfile,val_parfile,npar)
  use constants,only: IN_DATA_FILES
  implicit none
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
!---------------------------------------------------------------------------------
!

  subroutine stop_error_model_ASKI(error_message)
  use model_ASKI
  use constants,only: OUTPUT_FILES_BASE
  implicit none

  character(len=*) :: error_message
  character(len=400) :: filename
  integer :: IOASKI

  write(filename,"(a,i6.6,a)") trim(OUTPUT_FILES_BASE)//'ERROR_model_external_ASKI_',model_ASKI_myrank,'.txt'

  call get_file_unit_model_ASKI(IOASKI)
  open(unit=IOASKI,file=filename,form='formatted',status='unknown',action='write')
  write(IOASKI,*) trim(error_message)
  close(IOASKI)
  call abort_mpi()

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


!----------------------------------------------------------------
!! !! ================= VM VM CUSTOM SUBROUTINE FOR DSM COUPLING
!----------------------------------------------------------------

  subroutine FindLayer(x,y,z)
!! FS FS THIS IS AN EMPTY DUMMY ROUTINE, ONLY FOR BEING ABLE TO COMPILE get_model.f90
    double precision :: x,y,z
  end subroutine FindLayer

