!----------------------------------------------------------------------------
!   Copyright 2015 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of ASKI version 1.0.
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
program transformSpecfem3dCartesianMeasuredData
  use specfem3dForASKI_mod
  use inversionBasics
  use seismicEvent
  use seismicEventList
  use seismicStation
  use seismicNetwork
  use discreteFourierTransform
  use componentTransformation
  use asciiDataIO
  use argumentParser
  use string
  use fileUnitHandler
  use errorMessage
  use mathConstants

  implicit none

  ! command line
  type (argument_parser) :: ap
  character(len=max_length_string) :: parfile,str
  character(len=max_length_string), dimension(:), pointer :: str_vec

  ! basics
  type (file_unit_handler) :: fuh
  type (error_message) :: errmsg
  character(len=39) :: myname = 'transformSpecfem3dCartesianMeasuredData'
  type (inversion_basics) :: invbasics

  ! component transformation
  integer :: ncomp,jcomp
  character(len=character_length_component), dimension(:), allocatable :: comp
  logical :: orientation_is_NEZ
  double precision, dimension(:,:), pointer :: trans_coef
  real, dimension(:,:,:), allocatable :: trans_coef_all_transpose

  ! specfem seismograms
  character(len=2) :: band_instrument_code
  character(len=100) :: seisfile_extension
  logical :: force_seisfile_extension
  integer :: NSTEP
  real :: DT
  real, dimension(:,:), pointer :: traces
  real, dimension(:), pointer :: stf

  ! fourier transformation
  type (discrete_fourier_transform) :: DFT
  complex, dimension(:,:), allocatable :: spectra,spectrum_rotated
  complex, dimension(:), allocatable :: stf_spectrum
  real :: df
  integer :: nfreq
  integer, dimension(:), pointer :: jf

  double precision :: unit_factor,inverse_unit_factor

  ! filtering
  logical :: apply_filter,apply_event_filter,apply_station_filter
  character(len=400) :: path_event_filter,path_station_filter
  complex, dimension(:), pointer :: event_filter,station_comp_filter,filter_values

  ! other stuff
  integer :: istat,lu,ifreq
  logical :: one_event_only,seismograms_are_bin,deconvolve_stf,differentiate_stf,diff_time_series,scale_time_series,&
       use_imag_freq_gemini,print_usage_and_stop
  real :: ts_scale_factor
  type (seismic_event) :: event
  character(len=character_length_evid) :: evid,evid_one_event_only
  type (seismic_station) :: station
  character(len=400) :: path_specfem_seismograms,path_measured_data,file_measured_data

!------------------------------------------------------------------------
!  arguments
!
  call init(ap,myname,"Transform standard SPECFEM3D Cartesian 3.0 output to ASKI 1.0 spectral data in measured-data format")
  call addPosarg(ap,"main_parfile","sval","Main parameter file of inversion")
  call addOption(ap,"-bicode",.true.,"(mandatory) bandcode and instrument code: the first two characters before "//&
       "the component in seismogram filename, e.g. 'LH' if your filenames look like 'network.staname.LH*.semd'",&
       "sval","")
  call addOption(ap,"-ori",.true.,"(mandatory) orientation : either 'NEZ' or 'XYZ', indicating the component "//&
       "orientations following band_instrument_code","sval","")
  call addOption(ap,"-dt",.true.,"(mandatory) time step of the seismograms (as in SPECFEM3D Par_file)","rval","0.0")
  call addOption(ap,"-nstep",.true.,"(mandatory) number of samples NSTEP as in SPECFEM3D Par_file","ival","0")
  call addOption(ap,"-ocomp",.true.,"(mandatory) receiver components for which measured data is produced; valid "//&
       "components: '"//trim(all_valid_components)//"')","svec","")
  call addOption(ap,"-ext",.true.,"(optional) force a specific file extension: WILL IGNORE STANDARD FUNCTIONALITY!"//&
       "should be ANYTHING following the orientation character (including ALL dots etc.), e.g. '.semv' if your "//&
       "filenames look like 'network.staname.FX*.semv'. ","sval","")
  call addOption(ap,"-filter",.false.,"(optional) if set, event and/or station (component) "//&
       "filters will be applied to the spectra according to the respective flags in main parfile")
  call addOption(ap,'-evid',.true.,"(optional) indicates a single event for which measured data is produced, "//&
       "otherwise measured data is produced for all events (as defined in ASKI FILE_EVENT_LIST)",'sval','')
  call addOption(ap,"-gemini",.false.,"(optional) use to produce GEMINI consistent data: COMPLEX valued "//&
       "frequencies  f = jf*df + i*sigma  with real part jf*df and constant (!) imaginary part sigma = -5*df/2pi "//&
       "are used for DFT (filter values are assumed to be at those frequencies)")
  call addOption(ap,"-dconv",.false.,"(optional) if set, the source time function will be deconvolved from "//&
       "SPECFEM seismograms; consistend with 'ASKI_DECONVOLVE_STF = .true.' in Par_file_ASKI")
  call addOption(ap,"-diffts",.false.,"(optional) if set, the time series will be differentiated (by simple "//&
       "first order central differences) before further processing")
  call addOption(ap,"-scale",.true.,"(optional) factor (different from 0) by which the time series are scaled "//&
       "before further processing","rval","1.0")
  call addOption(ap,"-bin",.false.,"(optional) indicates whether SPECFEM trace files are binary files or not. "//&
       "For ascii output simply do not set option -bin")
!
  call parse(ap)
  if (.level.(.errmsg.ap) == 2) goto 3
!
  str = ap.sval."main_parfile"
  parfile = str
  if (.level.(.errmsg.ap) == 2) goto 3
!
  print_usage_and_stop = .false.
!
  ! -bicode
  if(.not.(ap.optset."-bicode")) then
     write(*,*) "ERROR: please indicate -bicode"
     print_usage_and_stop = .true.
  else
     str = ap.sval."-bicode"
     band_instrument_code = str
  end if
!
  ! -ori
  if(.not.(ap.optset."-ori")) then
     write(*,*) "ERROR: please indicate -ori"
     print_usage_and_stop = .true.
  else
     select case(ap.sval."-ori")
     case('XYZ'); orientation_is_NEZ = .false.
     case('NEZ'); orientation_is_NEZ = .true.
     case default
        write(*,*) "ERROR: invalid -ori value '"//trim(ap.sval."-ori")//"'"
        print_usage_and_stop = .true.
     end select
  end if
!
  ! -ext
  force_seisfile_extension = ap.optset."-ext"
  if(force_seisfile_extension) then
     str = ap.sval."-ext"
     seisfile_extension = str
  end if
!
  ! -bin
  seismograms_are_bin = (ap.optset."-bin")
!
  ! -evid
  one_event_only = (ap.optset."-evid")
  if(one_event_only) then
     str = ap.sval."-evid"
     evid_one_event_only = str
  end if
!
  ! -dt
  if(.not.(ap.optset."-dt")) then
     write(*,*) "ERROR: please indicate -dt"
     print_usage_and_stop = .true.
  else
     DT = ap.rval."-dt"
  end if  
!
  ! -nstep
  if(.not.(ap.optset."-nstep")) then
     write(*,*) "ERROR: please indicate -nstep"
     print_usage_and_stop = .true.
  else
     NSTEP = ap.ival."-nstep"
  end if  
!
  ! -filter
  apply_filter = (ap.optset."-filter")
!
  ! -ocomp
  if(.not.(ap.optset."-ocomp")) then
     print *, "ERROR: please indicate -ocomp"
     print_usage_and_stop = .true.
  else
     str_vec => ap.svec.'-ocomp'
     if (.level.(.errmsg.ap) == 2) goto 3
     if(.not.associated(str_vec)) then
        write(*,*) "ERROR: for some reason, there is no list of station components returned by argument parser, "//&
             "even though there was no error parsing argument -ocomp. This is strange..."
        write(*,*) ""
        goto 3
     end if
     ncomp = size(str_vec)
     allocate(comp(ncomp))
     do jcomp=1,ncomp
        comp(jcomp) = str_vec(jcomp)
     end do
     deallocate(str_vec)
     if(.not.allValidComponents(comp,i_invalid=jcomp)) then
        write(*,*) "ERROR: ",jcomp,"'th output component '"//trim(comp(jcomp))//"' not valid. Valid components are '"//&
             all_valid_components//"'"
        goto 3
     end if
  end if
!
  ! -dconv
  deconvolve_stf = (ap.optset."-dconv")
!
  ! -diffts
  diff_time_series = (ap.optset."-diffts")
!
  ! -scale
  scale_time_series = (ap.optset."-scale")
  if(scale_time_series) then
     ts_scale_factor = ap.rval."-scale"
     if (.level.(.errmsg.ap) == 2) goto 3
     if(ts_scale_factor == 0.0) then
        write(*,*) "ERROR: scaling factor for time series ('-scale' option) must not be zero!"
        print_usage_and_stop = .true.
     end if
  end if ! scale_time_series
!
  ! -gemini
  use_imag_freq_gemini = (ap.optset."-gemini")
!
  print_usage_and_stop = print_usage_and_stop .or. .level.(.errmsg.ap) == 2
!
  if(print_usage_and_stop) goto 3
!
  call document(ap)
  write(*,*) ""
!
  ! creat file unit handler  
  call createFileUnitHandler(fuh,100)
!
!------------------------------------------------------------------------
!  setup basics
!
  ! setup inversion basics
  call new(errmsg,myname)
  call init(invbasics,trim(parfile),get(fuh),errmsg)
  call undo(fuh)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
!
!------------------------------------------------------------------------
!  preliminary processing
!
  ! in case of one_event_only, check if evid_one_event_only is valid
  if(one_event_only) then
     errmsg = searchEventidSeismicEventList(.evlist.invbasics,evid_one_event_only)
     if(.level.errmsg /=0) call print(errmsg)
     if(.level.errmsg == 2) then
        write(*,*) "ERROR: eventID '"//trim(evid_one_event_only)//"' given by option -evid is not contained in event list"
        goto 3
     end if
  end if
!
  ! get all transformation matrices here for transformations from specfem seismograms orientation
  ! XYZ (or NEZ, dependent on value of orientation_is_NEZ) to the requested output components contained in array comp
  allocate(trans_coef_all_transpose(3,ncomp,.nstat.(.statlist.invbasics)))
  istat = 0
  do while (nextStationSeismicNetwork(.statlist.invbasics,station))
     istat = istat + 1
     ! transpose trans_coef here by shwitching XYZ(NEZ) with comp (i.e. coef_in = comp and coef_out = CX,CY,CZ(N,E,UP)),
     ! as we need the transpose of the transformation matrix in matmul operation below
     if(orientation_is_NEZ) then
        trans_coef => transform(.comptrans.invbasics,comp,(/'N ','E ','UP'/),.staname.station)
     else
        trans_coef => transform(.comptrans.invbasics,comp,(/'CX','CY','CZ'/),.staname.station)
     end if
     if(.not.associated(trans_coef)) then
        write(*,*) "ERROR: no transformation coefficients for ",istat,"'th station '"//trim(.staname.station)//"'"
        goto 1
     end if
     trans_coef_all_transpose(:,:,istat) = trans_coef
     deallocate(trans_coef)
  end do ! while next station
!
  if(apply_filter) then
     apply_event_filter = (.inpar.invbasics).lval.'APPLY_EVENT_FILTER'
     apply_station_filter = (.inpar.invbasics).lval.'APPLY_STATION_FILTER'
     apply_filter = apply_event_filter .or. apply_station_filter ! set apply_filter dependent on parameters in main parfile
  end if
!
  unit_factor = .ufmdata.invbasics
  inverse_unit_factor = 1.d0/unit_factor
!
  ! compute fourier transformation factors efactors here
  ! taper parameters like in SPECFEM3D-for-ASKI codes
  df = rval(.inpar.invbasics,'MEASURED_DATA_FREQUENCY_STEP')
  nfreq = ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')
  jf => ivecp(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',nfreq)
  if(.not.associated(jf)) then
     write(*,*) "ERROR: could not read ",nfreq," frequency indices from vector 'MEASURED_DATA_INDEX_OF_FREQ' in main parfile '"&
          //trim(parfile)//"'"
     goto 1
  end if
  call new(errmsg,myname)
  if(use_imag_freq_gemini) then
     call initiateForwardDFT(DFT,DT,0,NSTEP-1,cmplx(jf*df)+mc_ci*cmplx(-5.*df/mc_two_pi),errmsg,hanning_taper=0.05)
  else
     call initiateForwardDFT(DFT,DT,0,NSTEP-1,jf*df,errmsg,hanning_taper=0.05)
  end if
  if(.level.errmsg /=0) call print(errmsg)
  if(.level.errmsg==2) goto 1
  call dealloc(errmsg)
!
  allocate(spectra(nfreq,3*.nstat.(.statlist.invbasics)),spectrum_rotated(nfreq,ncomp))
  if(deconvolve_stf) allocate(stf_spectrum(nfreq))
!
!------------------------------------------------------------------------
!  write some info about this run now
!
  if(one_event_only) then
     write(*,*) "creating ASKI measured data from SPECFEM3D seismograms for one event and ",&
          .nstat.(.statlist.invbasics)," stations, "
  else
     write(*,*) "creating ASKI measured data from SPECFEM3D seismograms for ",.nev.(.evlist.invbasics)," events and ",&
          .nstat.(.statlist.invbasics)," stations, "
  end if
  write(*,*) "as of main parameter file '"//trim(parfile)//"'"
  write(*,*) ""
  write(*,*) "input SPECFEM3D seismograms: "
  write(*,*) "   NSTEP =  ",NSTEP
  write(*,*) "   DT =  ",DT
  write(*,*) "   band and instrument code = ",band_instrument_code
  if(orientation_is_NEZ) then
     write(*,*) "   seismogram orientation = NEZ"
  else
     write(*,*) "   seismogram orientation = XYZ"
  end if
  write(*,*) ""
  write(*,*) "output spectra: "
  write(*,*) "   number of frequencies = ",nfreq
  write(*,*) "   frequency step df = ",df
  write(*,*) "   frequency indices jf = ",jf
  if(use_imag_freq_gemini) then
     write(*,*) "   assuming gemini-like complex frequencies   f = jf*df + i*sigma   with constant imaginary part ",&
          "sigma = -5*df/2pi = ",-5.*df/mc_two_pi
  else
     write(*,*) "   assuming real-valued frequencies f = jf*df"
  end if
  if(scale_time_series) then
     write(*,*) "   time series will be scaled by factor ",ts_scale_factor," before further processing"
  else
     write(*,*) "   time series will NOT be scaled by a factor"
  end if
  if(diff_time_series) then
     write(*,*) "   time series will be differentiated by simple 1st order finite differences before further processing"
  else
     write(*,*) "   time series will NOT be differentiated"
  end if
  if(deconvolve_stf) then
     if(.not.force_seisfile_extension) then
        write(*,*) "   source time function will be deconvolved from synthetics (producing displacement w.r.t. dirac)"
     else
        write(*,*) "   source time function will be deconvolved from synthetics, stf will be differentiated in "//&
             "case of moment tensor sources. make sure you properly chose the forced file extension '",&
             trim(seisfile_extension),"'"
     end if
  else
     ! if no stf is deconvolved, simply read in displacement
     if(.not.force_seisfile_extension) then
        seisfile_extension = '.semd'
        write(*,*) "   NO deconvolution of source time function will be applied. Reading in displacement, ",&
             "using seismogram file extension ",seisfile_extension
     else
        write(*,*) "   NO deconvolution of source time function will be applied. Reading in ",&
             "seismogram files having forced extension ",seisfile_extension
     end if
  end if
  path_event_filter = (.inpar.invbasics).sval.'PATH_EVENT_FILTER'
  path_station_filter = (.inpar.invbasics).sval.'PATH_STATION_FILTER'
  if(apply_event_filter) then
     write(*,*) "   according to flag -filter and parameter 'APPLY_EVENT_FILTER' in main parfile, event filters "//&
          "will be applied, as given by filter files in path '"//trim(path_event_filter)//"'"
  else
     write(*,*) "   according to flag -filter or parameter 'APPLY_EVENT_FILTER' in main parfile, NO event "//&
          "filters will be applied"
  end if ! apply_event_filter
  if(apply_station_filter) then
     write(*,*) "   according to flag -filter and parameter 'APPLY_STATION_FILTER' in main parfile, station "//&
          "filters will be applied, as given by filter files in path '"//trim(path_station_filter)//"'"
  else
     write(*,*) "   according to flag -filter or parameter 'APPLY_STATION_FILTER' in main parfile, NO station "//&
          "filters will be applied"
  end if ! apply_station_filter
  write(*,*) "   the output spectra are computed according to the unit factor ",unit_factor,&
       ", as of parameter 'UNIT_FACTOR_MEASURED_DATA' in main parfile"
  path_measured_data = (.inpar.invbasics).sval.'PATH_MEASURED_DATA'
  write(*,*) "   measured data files will be written to path '",trim(path_measured_data),"'"
  write(*,*) ""
!
!
!------------------------------------------------------------------------
!  loop on all events (even if one_event_only, for simplicity of coding) 
!
  lu = get(fuh)
!
  do while (nextEventSeismicEventList(.evlist.invbasics,event))
!
     if(one_event_only) then
        evid = evid_one_event_only
     else
        evid = .evid.event
     end if
!
     ! in case of event filtering, read in event specific filter values
     if(apply_event_filter) then
        if(associated(event_filter)) deallocate(event_filter)
        errmsg = readAsciiData(trim(path_event_filter)//"filter_"//trim(evid),lu,event_filter,ndata=nfreq)
        if(.level.errmsg /=0) call print(errmsg)
        if(.level.errmsg == 2) then
           write(*,*) "ERROR: could not read source filter from ascii file '"//trim(path_event_filter)//&
                "filter_"//trim(evid)
           goto 1
        end if
        if(.not.apply_station_filter) then
           ! Already at this point, prepare the variable filter_values which will be finally used to filter the spectra.
           ! Below, there will be no code executed related to station filtering (definition of variable filter_values)
           if(associated(filter_values)) deallocate(filter_values)
           filter_values => event_filter
           nullify(event_filter)
        else
           ! Otherwise, do nothing here, but keep the values event_filter to be used below to compose the final filter
        end if ! .not.apply_station_filter
     end if ! apply_event_filter
!
     path_specfem_seismograms = trim(path_measured_data)//'data_'//trim(evid)//'_OUTPUT_FILES/'
!
     write(*,*) "read all traces of event '",trim(evid),"' from path '",trim(path_specfem_seismograms),"'"
!
     ! in case that the source time function should be deconvolved, it is assumed here that it was NOT a Ricker 
     ! wavelet (but also works in that case!). In the compatible SPECFEM3D_Cartesian code, different "impulsive" 
     ! source time functions are used, dependent on the source mechanism: thin Gaussian for single force source; 
     ! steep error function for moment tensor source. check here, what kind of source type this event has.
     ! if it is a single force:
     ! - read in displacement seismograms and deconvolve the plain source time function
     ! if it is a moment tensor source:
     ! - read in velocity seismograms and deconvolve the differentiated source time function (for reasons of 
     !   numerical stability at very small frequencies)
     if(deconvolve_stf) then
        select case(.styp.event)
        case(0)
           ! SINGLE FORCE SOURCE
           if(.not.force_seisfile_extension) then
              seisfile_extension = '.semd'
              write(*,*) "   this event has a single force source mechanism, so reading in displacement using ",&
                   "seismogram file extension '",trim(seisfile_extension),&
                   "' and will NOT differentiate the source time function before deconvolving it"
           else
              write(*,*) "   this event has a single force source mechanism, so reading in ",&
                   "seismogram files having forced extension '",trim(seisfile_extension),&
                   "' and will NOT differentiate the source time function before deconvolving it"
           end if
           differentiate_stf = .false.
        case(1)
           ! MOMENT TENSOR SOURCE
           if(.not.force_seisfile_extension) then
              seisfile_extension = '.semv'
              write(*,*) "   this event has a moment tensor source mechanism, so reading in velocity using ",&
                   "seismogram file extension '",trim(seisfile_extension),&
                   "' and WILL differentiate the source time function before deconvolving it"
           else
              write(*,*) "   this event has a moment tensor source mechanism, so reading in ",&
                   "seismogram files having forced extension '",trim(seisfile_extension),&
                   "' and WILL differentiate the source time function before deconvolving it"
           end if
           differentiate_stf = .true.
        case default
           write(*,*) "ERROR: type of source of this event is: ",.styp.event,&
                ". It must be either single force (0) or moment tensor (1). Nothing else is supported here."
           goto 1
        end select
     end if ! deconvolve_stf
!
     ! read in all traces
     if(associated(traces)) deallocate(traces)
     call readTraces(traces,NSTEP,.statlist.invbasics,path_specfem_seismograms,band_instrument_code,&
          orientation_is_NEZ,seisfile_extension,seismograms_are_bin,lu)
     if(.not.associated(traces)) then
        write(*,*) "no spectra produced for this event"
        goto 2
     end if
!
     ! scale time series, if requested
     if(scale_time_series) traces = traces*ts_scale_factor
!
     ! differentiate by central differences, if requested
     if(diff_time_series) call diffTraces(traces,DT)
!
     if(deconvolve_stf) then
        if(associated(stf)) deallocate(stf)
        if(differentiate_stf) then
           write(*,*) "read and differentiate source time function of event '",trim(evid),"' from file '",&
                trim(path_specfem_seismograms),"plot_source_time_function.txt'"
           call getDiffStf(stf,NSTEP,DT,trim(path_specfem_seismograms)//'plot_source_time_function.txt',lu)
        else
           write(*,*) "read source time function of event '",trim(evid),"' from file '",&
                trim(path_specfem_seismograms),"plot_source_time_function.txt'"
           call getStf(stf,NSTEP,trim(path_specfem_seismograms)//'plot_source_time_function.txt',lu)
        end if
        if(.not.associated(stf)) then
           write(*,*) "no spectra produced for this event"
           goto 2
        end if
     end if ! deconvolve_stf
!
     write(*,*) "compute spectra from traces"
     ! fourier transform to frequency domain of all traces at once
     call new(errmsg,myname)
     call transformForwardDFT(DFT,traces,spectra,errmsg)
     if(.level.errmsg /=0) call print(errmsg)
     if(.level.errmsg==2) goto 1
     call dealloc(errmsg)
     if(deconvolve_stf) then
        write(*,*) "  compute spectrum of source time function:"
        ! fourier transform to frequency domain of all traces at once
        call new(errmsg,myname)
        call transformForwardDFT(DFT,stf,stf_spectrum,errmsg)
        if(.level.errmsg /=0) call print(errmsg)
        if(.level.errmsg==2) goto 1
        call dealloc(errmsg)
        do ifreq = 1,nfreq
           write(*,*) jf(ifreq)*df," Hz : ",stf_spectrum(ifreq)
        end do
        ! deconvolve spectra
        !   THE ONLY DIFFERENCE TO THE CODE IN specfem3d_for_ASKI.f90 IS:
        !   HERE, spectra AND stf_spectrum ARE ALWAYS SINGLE PRECISION!
        !   IN specfem3d_for_ASKI.f90 , stf_spectrum IS ALWAYS DOULBLE PRECISION, AND BY OPTIONAL FLAG spectra MAY BE DOUBLE PRECISION, TOO
        write(*,*) "  deconvolve stf from seismograms"
        do ifreq = 1,nfreq
           spectra(ifreq,:) = spectra(ifreq,:) / stf_spectrum(ifreq)
        end do ! ifreq
     end if ! deconvolve_stf
!
     ! write spectra to files
     write(*,*) "write measured data files to path '",trim(path_measured_data),"'"
     istat = 0
     do while (nextStationSeismicNetwork(.statlist.invbasics,station))
        istat = istat + 1
!
        ! transform all 3 traces of this station (assume order X,Y,Z or N,E,Z coming from routine readTraces,
        ! dependent on orientation_is_NEZ) to the requested components in comp and apply the conversion factor to the requested unit.
        ! Note that the standard SPECFEM3D seismograms are in SI units (displacements in meters). Thus, multiply by inverse_unit_factor here.
        spectrum_rotated(:,:) = inverse_unit_factor * &
             matmul(spectra(:,(istat-1)*3+1:(istat-1)*3+3) , trans_coef_all_transpose(:,:,istat))
!
        do jcomp = 1,ncomp
!
           ! in case of station filtering, read in station_comp_filter for this path and component
           if(apply_station_filter) then
              if(associated(station_comp_filter)) deallocate(station_comp_filter)
              errmsg = readAsciiData(trim(path_station_filter)//"filter_"//trim(.staname.station)//"_"//trim(comp(jcomp)),&
                   lu,station_comp_filter,ndata=nfreq)
              if(.level.errmsg /=0) call print(errmsg)
              if(.level.errmsg == 2) then
                 write(*,*) "ERROR: could not read station component filter from ascii file '"//trim(path_station_filter)//&
                      "filter_"//trim(.staname.station)//"_"//trim(comp(jcomp))
                 goto 1
              end if
              if(associated(filter_values)) deallocate(filter_values)
              filter_values => station_comp_filter
              nullify(station_comp_filter)
              ! in case that also event filters are to be applied, we need to compose here final filter values from event AND station filters
              if(apply_event_filter) filter_values = filter_values * event_filter
           end if ! apply_station_filter
!
           ! define filename of output file
           file_measured_data = "data_"//trim(evid)//"_"//trim(.staname.station)//"_"//trim(comp(jcomp))
!
           write(*,*) "writing measured data file '",trim(file_measured_data),"'"
           if(apply_filter) then
              errmsg = writeAsciiData(trim(path_measured_data)//file_measured_data,lu,&
                   filter_values*spectrum_rotated(:,jcomp))
           else
              errmsg = writeAsciiData(trim(path_measured_data)//file_measured_data,lu,spectrum_rotated(:,jcomp))
           end if
           if(.level.errmsg/=0) call print(errmsg)
           if(.level.errmsg==2) goto 1
           call dealloc(errmsg)
!
        end do ! jcomp
     end do ! while next station
!
2    write(*,*) ""
     if(one_event_only) exit
  end do ! while next event
!
!------------------------------------------------------------------------
!  clean up
!
  write(*,*) "good bye"
!
1 call dealloc(invbasics)
  call add(fuh,lu); call dealloc(fuh)
  call dealloc(ap)
  call dealloc(DFT)
  call dealloc(errmsg)
  if(allocated(comp)) deallocate(comp)
  if(associated(trans_coef)) deallocate(trans_coef)
  if(allocated(trans_coef_all_transpose)) deallocate(trans_coef_all_transpose)
  if(allocated(spectra)) deallocate(spectra)
  if(allocated(spectrum_rotated)) deallocate(spectrum_rotated)
  if(allocated(stf_spectrum)) deallocate(stf_spectrum)
  if(associated(jf)) deallocate(jf)
  if(associated(traces)) deallocate(traces)
  if(associated(event_filter)) deallocate(event_filter)
  if(associated(station_comp_filter)) deallocate(station_comp_filter)
  if(associated(filter_values)) deallocate(filter_values)
  if(associated(str_vec)) deallocate(str_vec)
!
  stop
!
3 if(.level.(.errmsg.ap)>=1) call print(.errmsg.ap)
  call usage(ap)
  goto 1
end program transformSpecfem3dCartesianMeasuredData
!
!------------------------------------------------------------------------
!
! subroutine printhelp
!   use componentTransformation, only: all_valid_components
!   print '(50(1h-))'
!   print *,'Usage:'
!   print *,''
!   print *,'  transformSpecfem3dCartesianMeasuredData [-h] -bicode band_instrument_code -ori orientation -ext file_extension'
!   print *,'     -dt time_step -nstep number_of_time_samples -ocomp "ncomp comp" [-filter] [-evid eventID] [-gemini]'
!   print *,'    [-dconvd] [-uf unit_factor] [-diffts] [-scale ts_scale_factor] [-bin] parfile'
!   print *,''
!   print *,'Arguments:'
!   print *,''
!   print *,"  parfile: main parameter file of inversion"
!   print *,''
!   print *,'Mandatory options:'
!   print *,''
!   print *,'  -bicode band_instrument_code : band_instrument_code must be two characters, band code and instrument code'
!   print *,'                                 i.e. the first two characters before the component in seismogram filename'
!   print *,"                                 e.g. 'LH' if your filenames look like 'staname.network.LH*.semd.ascii'"
!   print *,''
!   print *,"  -ori orientation : either 'NEZ' or 'XYZ', indicating the component orientations following band_instrument_code"
!   print *,''
!   print *,'  -ext file_extension : file_extension should be ANYTHING following the orientation character (including ALL dots etc.)'
!   print *,"                                 e.g. if your filenames look like 'staname.network.FX*.semv', file_extension = '.semv'"
!   print *,''
!   print *,"  -dt time_step : time_step is the real number defining the time step of the seismograms (as in SPECFEM3D Par_file)"
!   print *,''
!   print *,"  -nstep number_of_time_samples : number_of_time_samples is the number of samples NSTEP as in SPECFEM3D Par_file"
!   print *,''
!   print *,"  -ocomp : defines the components for which measured data output is produced"
!   print *,"          ncomp: number of components,  comp: ncomp components (valid components: '"//&
!        trim(all_valid_components)//"')"
!   print *,''
!   print *,'Optional options:'
!   print *,''
!   print *,'  -filter : if set, the respective event filters and station (component) filters will be applied to the spectra'
!   print *,''
!   print *,'  -evid eventID : if set, eventID indicates the single event for which measured data is produced. otherwise,'
!   print *,'                  measured data is produced for all events (as defined in ASKI FILE_EVENT_LIST)'
!   print *,''
!   print *,'  -gemini : if set, COMPLEX valued frequencies  f = jf*df + i*sigma  with usual real part jf*df and'
!   print *,'            constant (!) imaginary part sigma = -5*df/2pi are used for the discrete Fourier transform of the SPECFEM3D'
!   print *,'            output seismograms and the filters are also assumed to be given at those frequencies, etc. Use this flag'
!   print *,'            when producing synthetically compputed data for inversion with GEMINI, which computes spectral synthetics'
!   print *,'            and kernels at those kind of complex frequencies.'
!   print *,''
!   print *,'  -dconvd : if set, the derivative of the heaviside source time function will be deconvolved from the seismograms.'
!   print *,"            It is assumed that it was written to file 'plot_source_time_function.txt', i.e. flag"
!   print *,'            PRINT_SOURCE_TIME_FUNCTION was set to .true. in SPECFEM3D Par_file. '
!   print *,"            -dconvd is consistend with 'ASKI_DECONVOLVE_STF = .true.' in Par_file_ASKI. In case ASKI_DECONVOLVE_STF"
!   print *,"            was .true. , flag -dconvd should be set, too!! In this case, velocity seismograms (i.e. extensions '.semv')"
!   print *,"            must be used!"
!   print *,"            The content of file 'plot_source_time_function.txt' is NORMALIZED to maximum value 1 (relevant for use "//&
!        "with SPECFEM3D GLOBE)"
!   print *,''
!   print *,"  -uf unit_factor : if set, displacement spectra are produced in the unit according to value unit_factor > 0."
!  print *,"                    Here, the very same value should be set as UNIT_FACTOR_MEASURED_DATA in the ASKI main parameter file!"
!   print *,"                    The unit factor in ASKI is defined as follows: Multiplication by value unit_factor has the effect to"
!   print *,"                    transform the resulting output spectra to SI units [ms]. In particular, output spectra in [ms] "
!   print *,"                    correspond to unit_factor = 1.0."
!   print *,"                    BY DEFAULT, unit_factor = 1.0 (no need to set -uf if you want to produce spectra in [ms])."
!   print *,"                    For example: inverting time-domain displacement data given in the unit of nano meters, unit_factor"
!   print *,"                    should be set to the value of 1.0e-9 here (and in the ASKI main parfile), since nano meters times"
!   print *,"                    1.0e-9 gives the SI unit meters."
!   print *,""
!   print *,'  -diffts : If set, the time series will be differentiated (by simple first order central differences) before '//&
!        'further processing.'
!   print *,'            This option is sensible to set when processing displacement seismograms which were calculated w.r.t. a '//&
!        'Heaviside stf'
!   print *,'            (e.g. SPECFEM3D GLOBE). Using additionally option -dconvd, this then results in displacement spectra '//&
!        'w.r.t. a dirac impulse stf'
!   print *,''
!   print *,'  -scale ts_scale_factor : if set, the time series are scaled with factor ts_scale_factor before further processing'
!   print *,''
!   print *,'  -bin : indicates whether SPECFEM trace files are binary files or not. for ascii output simply do not set '//&
!        'option -bin'
!   print *,''
!   print *,'  -h : print this help message'
!   print '(50(1h-))'
!   return
! end subroutine printhelp
