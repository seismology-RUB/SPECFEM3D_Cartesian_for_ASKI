!=====================================================================

!----------------------------------------------------------------------------
!   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!
!   module specfem_par_ASKI is part of ASKI version 1.2.
!
!   ASKI version 1.2 is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 2 of the License, or
!   (at your option) any later version.
!
!   ASKI version 1.2 is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with ASKI version 1.2.  If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------

module specfem_par_ASKI

! specific parameter module for ASKI extension package

  implicit none

! parameters for generating ASKI output
  integer :: ASKI_np_local
  integer, dimension(:), allocatable :: ASKI_np_local_all
  integer, dimension(:,:), allocatable :: ASKI_indx_local
  complex(kind=kind(1.d0)), dimension(:,:), allocatable :: ASKI_efactors
  double precision, dimension(:), allocatable :: ASKI_Goertzel_Wr,ASKI_Goertzel_Wi
  complex(kind=kind(1.d0)), dimension(:,:,:), allocatable :: ASKI_spectra_local_double
  complex, dimension(:,:,:), allocatable :: ASKI_spectra_local_single
  double precision, dimension(:,:,:), pointer :: ASKI_Goertzel_U0_local_double,ASKI_Goertzel_U1_local_double,&
       ASKI_Goertzel_U2_local_double
  real, dimension(:,:,:), pointer :: ASKI_Goertzel_U0_local_single,ASKI_Goertzel_U1_local_single,&
       ASKI_Goertzel_U2_local_single
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

! absolute output filename including absolute path (will be used to open the output file as is),
! i.e. in case of a green tensor simulation, outfile should contain the extension "_[green_tensor_component]"
! ASKI_outfile is handled as a basename, to which ".main" and for each frequency ".jf######" 
! (e.g. ".jf000013" for frequency index 13) is appended
  character(len=500) :: ASKI_outfile

! choose to overwrite existing ASKI output files, or to abort if files already exist
  logical :: OVERWRITE_ASKI_OUTPUT

! id of ASKI output (e.g. eventID, station name plus component, etc.)
  integer, parameter :: length_ASKI_output_ID = 13
  character(len=length_ASKI_output_ID) :: ASKI_output_ID

! choose to deconvolve the derivative of the source time function from the wavefield spectra before writing them to files.
! select .true. for any Green function computations!
! even if a heaviside source time function is used, the velocity field is not exactly a Green function ( = displacement wavefield
! w.r.t. an impulse source time function), since a steep error function is used by SPECFEM to resemble a quasi-heaviside function.
! this steep error function, furthermore, is dependent on timestep DT! 
! hence, only by deconvolution of (the derivative of) this quasi-heaviside source time function, the real Green function, 
! which is independent of the time step can be computed.
  logical :: ASKI_DECONVOLVE_STF
  complex(kind=kind(1.d0)), dimension(:), allocatable :: ASKI_stf_spectrum_double

!------------------------------------------------------------------------------------
! FREQUENCY DISCRETIZATION AND FOURIER TRANSFORM
!
!   the double precision df[Hz] and integer values jf have the following meaning:
!   the spectra are saved for all frequencies f = (jf)*df [Hz]
!------------------------------------------------------------------------------------

! Predefined frequency step df, that is used to evaluate spectrum (in case we want to do an inverse FT, we need to choose with care 
! df=1/length_of_time_series and suitably high frequency indices (dependent on frequency content), as
! we could lose periodicity (if in exp^(-i2pi(k)(n)/N) "N" is no integer, these are no roots of 1 anymore)).
! Save the spectra for frequencies f = (ASKI_jf)*ASKI_df (ASKI_nf many)
  double precision :: ASKI_df
  integer :: ASKI_nf
  integer, dimension(:), allocatable :: ASKI_jf

! Choose the method of Discrete Fourier Transform:
! EXPLICIT_SUMMATION : on-the-fly summation of complex values s(t)*exp^(-i*2pi*f*t) (slightly more memory efficient than GOERTZEL_STANDARD)
! GOERTZEL_STANDARD: (recommended) using Goertzel's algorithm (only half the number of multiplications compared with EXPLICIT_SUMMATION -> much more efficient)
  integer, parameter :: length_ASKI_DFT_method = 18
  character(len=length_ASKI_DFT_method) :: ASKI_DFT_method

! Choose precision of Discrete Fourier Transform. If there is enough memory available, it is highly recommended
! to use ASKI_DFT_double = .true. in which case double precision spectral values are hold in memory (single  
! precision is written to file, though, but less roundoffs during summation).
! Otherwise choose ASKI_DFT_double = .false. in which case single precision spectral values will be used in memory.
! The wavefield point-independent transformation coefficients exp^(-i*2pi*f*t) (or values Wr in case of Goertzel)
! are always in double precision!
  logical :: ASKI_DFT_double

! Decide whether the (oversampled, noisy, ...) time series should be tapered by a hanning taper (on tail)
! before (i.e. on-the-fly while) applying the discrete fourier transform.
! If ASKI_DFT_apply_taper = .true. the value of ASKI_DFT_taper_percentage (between 0.0 and 1.0) defines the amount of
! total time for which the hanning taper will be applied at the tail of the time series.
  logical :: ASKI_DFT_apply_taper
  double precision :: ASKI_DFT_taper_percentage
  integer :: ASKI_DFT_ntaper_start
  double precision, dimension(:), allocatable :: ASKI_DFT_taper_values

!------------------------------------------------------------------------------------
! INVERSION GRID
!
! ASKI supports several types of inversion grids for FORWARD_METHOD = SPECFEM3D:
!
!   ASKI_type_inversion_grid = 1, (TYPE_INVERSION_GRID = schunkInversionGrid) NOT TO BE USED WITH SPECFEM3D Cartesian!
!      ASKI internal, but SPECFEM independent simple spherical inverison grid
!
!   ASKI_type_inversion_grid = 2, (TYPE_INVERSION_GRID = scartInversionGrid)
!      ASKI internal, but SPECFEM independent cartesian inversion grid:
!      the values for ASKI output are stored at all inner GLL points of spectral elements which lie
!      inside the block volume defined below by parameters ASKI_(cw)(xyz)
!      ASKI loactes the coordinates of those points inside the inversion grid cells and computes
!      integration weights for them
!
!   ASKI_type_inversion_grid = 3, (TYPE_INVERSION_GRID = ecartInversionGrid)
!      external inversion grid provided e.g. by CUBIT, which may contain tetrahedra, as well as hexahedra
!      as in case of ASKI_type_inversion_grid = 2, ASKI output is stored at all inner GLL points of elements
!      which are inside the volume defined by ASKI_(cw)(xyz)
!      ASKI locates the wavefield points inside the inversion grid and computes weights
!
!   ASKI_type_inversion_grid = 4, (TYPE_INVERSION_GRID = specfem3dInversionGrid)
!      use SPECFEM elements as inversion grid:
!      wavefield points are ALL GLL points of an element for elements which are (at least partly) inside the 
!      volume defined by ASKI_(cw)(xyz)
!      additionally store the jacobians for all wavefield points
!      assume ncell = ntot_wp/(NGLLX*NGLLY*NGLLY) as the number of inversion grid cells, and the order of 
!      wavefield points accordingly (do k=1,NGLLZ; do j=1,NGLLY; do i=1,NGLLX; ip=ip+1 ....)
!
!   ASKI_type_inversion_grid = 5, (TYPE_INVERSION_GRID = chunksInversionGrid) NOT TO BE USED WITH SPECFEM3D Cartesian!
!      ASKI internal, but SPECFEM independent more elaborate spherical inverison grid
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

