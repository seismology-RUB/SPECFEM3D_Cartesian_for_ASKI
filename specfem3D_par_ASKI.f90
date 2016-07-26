!=====================================================================

!----------------------------------------------------------------------------
!   Copyright 2015 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!
!   module specfem_par_ASKI is part of ASKI version 1.0.
!
!   ASKI version 1.0 is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 2 of the License, or
!   (at your option) any later version.
!
!   ASKI version 1.0 is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with ASKI version 1.0.  If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------

module specfem_par_ASKI

! specific parameter module for ASKI extension package

  implicit none

! parameters for generating ASKI output
  integer :: ASKI_np_local
  integer, dimension(:), allocatable :: ASKI_np_local_all
  integer, dimension(:,:), allocatable :: ASKI_indx_local
  complex(kind=kind(1.d0)), dimension(:,:), allocatable :: ASKI_efactors_tapered !variable is also used in case there is no tapering!
  complex(kind=kind(1.d0)), dimension(:,:,:), allocatable :: ASKI_spectra_local_double
  complex, dimension(:,:,:), allocatable :: ASKI_spectra_local_single
  logical :: ASKI_store_veloc

! content of Par_file_ASKI is following now:

!------------------------------------------------------------------------------------
! ASKI OUTPUT
!
! producing kernel green tensor / kernel displacement files
!------------------------------------------------------------------------------------

! parameter COMPUTE_ASKI_OUTPUT controls if any ASKI output files (i.e. kernel green tensor / 
! kernel displacement files) are produced
  logical :: COMPUTE_ASKI_OUTPUT

! Decide, whether to ONLY produce the main ASKI output file (at the beginning of this simulation, filename 'ASKI_outfile.main') 
! and then abort the simulation.
! This can be useful if you want to first run initBasics to check the ASKI output volume and the background model etc.
 logical :: ASKI_MAIN_FILE_ONLY

! choose to overwrite existing ASKI output files, or to abort if files already exist
  logical :: OVERWRITE_ASKI_OUTPUT

! absolute output filename including absolute path (will be used to open the output file as is),
! i.e. in case of a green tensor simulation, outfile should contain the extension "_[green_tensor_component]", like "_X","_Y","_Z"
  character(len=500) :: ASKI_outfile

! id of ASKI output (e.g. eventID, station name plus component, etc.)
  integer, parameter :: length_ASKI_output_ID = 13
  character(len=length_ASKI_output_ID) :: ASKI_output_ID

! flag which deconvolution of a source time function should be applied
  logical :: ASKI_DECONVOLVE_STF

!------------------------------------------------------------------------------------
! FREQUENCY DISCRETIZATION
!
!   the terms df,jf have the following meaning:
!   the spectra are saved for all frequencies f = (jf)*df
!------------------------------------------------------------------------------------

! predefined df, that is used to evaluate spectrum (in case we want to do an inverse FT, we need to choose with care 
! df=1/length_of_time_series and suitably high frequency indices (dependent on frequency content), as
! we could lose periodicity (if in exp^(-i2pi(k)(n)/N) "N" is no integer, these are no roots of 1 anymore))
! save the spectra for frequencies f = (ASKI_jf)*ASKI_df (ASKI_nf many)
  double precision :: ASKI_df
  integer :: ASKI_nf
  integer, dimension(:), allocatable :: ASKI_jf

! choose precision of Discrete Fourier Transform, if there is enough memory available, it is highly recommended
! to use ASKI_DFT_double = .true. in which casedouble precision complex coefficients exp^(-i*2pi*f*t) are used
! and double complex spectra are hold in memory(single precision is written to file, though, but less roundoffs
! during transformation
! otherwise choose ASKI_DFT_double = .false. inwhich case single precision will be used
  logical :: ASKI_DFT_double

! decide whether the (oversampled, noisy, ...) time series should be tapered bya hanning taper(on tail)
! before (i.e. while) applying the discrete fourier transform
! if ASKI_DFT_apply_taper = .true. the value ofASKI_DFT_taper_percentage (between 0.0 and 1.0)definesthe amount of
! total time for which the hanning taper will be applied at thetail ofthe time series
  logical :: ASKI_DFT_apply_taper
  double precision :: ASKI_DFT_taper_percentage

!------------------------------------------------------------------------------------
! INVERSION GRID
!
! ASKI supports several types of inversion grids for FORWARD_METHOD = SPECFEM3D:
!
!   ASKI_type_inversion_grid = 1 ('ccsInversionGrid') NOT TO BE USED WITH SPECFEM3D Cartesian!
!      ASKI internal, but SPECFEM independent spherical inverison grid
!
!   ASKI_type_inversion_grid = 2 ('scartInversionGrid')
!      ASKI internal, but SPECFEM independent cartesian inversion grid:
!      the values for ASKI output are stored at all inner GLL points of spectral elements which lie
!      inside the block volume defined below by parameters ASKI_(cw)(xyz)
!      the coordinates of those points are located inside the inversion grid cells, and ASKI computes
!      integration weights for them
!
!   ASKI_type_inversion_grid = 3 ('ecartInversionGrid')
!      external inversion grid provided, e.g. by CUBIT, which may contain tetrahedrons, as well as hexahedrons
!      as in case of ASKI_type_inversion_grid = 2, ASKI output is stored at all inner GLL points of elements
!      which are inside the volume defined by ASKI_(cw)(xyz)
!      ASKI locates the wavefield points inside the inversion grid and computes weights
!
!   ASKI_type_inversion_grid = 4 ('specfem3dInversionGrid')
!      use SPECFEM elements as inversion grid:
!      wavefield points are ALL GLL points of an element for elements which are (at least partly) inside the 
!      volume defined by ASKI_(cw)(xyz)
!      additionally store the jacobians for all wavefield points
!      assume ncell = ntot_wp/(NGLLX*NGLLY*NGLLY) as the number of inversion grid cells, and the order of 
!      wavefield points accordingly (do k=1,NGLLZ; do j=1,NGLLY; do i=1,NGLLX; ip=ip+1 ....)
!
! Dependent on ASKI_type_inversion_grid, (a selection of) the following parameters may be used to define a volume 
! within which wavefield points are searched for:
!
! First, ASKI_wx,ASKI_wy,ASKI_wz define the total width of a block which is centered in x=y=z=0
!    e.g. the total block extension in x-direction covers all points with 
!    x >= - 0.5*ASKI_wx and
!    x <=  0.5*ASKI_wx
! Then, ASKI_rot_X,ASKI_rot_Y,ASKI_rot_Z define rotation angles in degrees by which the block is rotated (anti-clockwise)
! about the Z, Y and X coordinate axis,
! before ASKI_cx,ASKI_cy,ASKI_cz define a vector by which the rotated block is shifted (new center of block)
!
! BE AWARE:
! - the parameters for rotation angles ASKI_rot_(XYZ) MUST ALWAYS be assinged to values! set to 0. if not used
! - scartInversionGrid only supports ASKI_rot_Z and uses a different definintion of the z-coverage
! - ecartInversionGrid and specfem3dInversionGrid use ALL rotation angles ASKI_rot_(XYZ)
!------------------------------------------------------------------------------------

! type of inversion grid, in order to produce correct output
  integer :: ASKI_type_inversion_grid

! width of ASKI output volume
  real :: ASKI_wx
  real :: ASKI_wy
  real :: ASKI_wz
! rotation angles in degrees of block about the respective axes
  real :: ASKI_rot_X
  real :: ASKI_rot_Y
  real :: ASKI_rot_Z
! center of ASKI output volume
  real :: ASKI_cx
  real :: ASKI_cy
  real :: ASKI_cz

end module specfem_par_ASKI

