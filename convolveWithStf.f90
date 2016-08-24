!----------------------------------------------------------------------------
!   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of ASKI version 1.2.
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
program convolveWithStf

  use specfem3dForASKI_mod
  use inversionBasics
  use seismicEvent
  use seismicEventList
  use seismicStation
  use seismicNetwork
  use componentTransformation
  use asciiDataIO
  use fileUnitHandler
  use argumentParser
  use string
  use errorMessage
  use inputParameter

  implicit none

  ! command line
  type (argument_parser) :: ap
  character(len=max_length_string) :: parfile,str

  ! basics
  type (file_unit_handler) :: fuh
  type (error_message) :: errmsg
  character(len=15) :: myname = 'convolveWithStf'
  type (inversion_basics) :: invbasics

  type (input_parameter) :: iterpar
  character (len=80), dimension(1) :: iter_inpar_keys
  ! keywords for input parameter
  data iter_inpar_keys/'PATH_KERNEL_DISPLACEMENTS'/

  ! component transformation
  logical :: orientation_is_NEZ
  double precision, dimension(:,:), pointer :: trans_coef
  real, dimension(:,:,:), allocatable :: trans_coef_all

  ! old source time function
  logical :: deconvolve_stf
  real, dimension(:), pointer :: old_stf
  real, dimension(:,:), pointer :: old_stf_time
  double complex, dimension(:), allocatable :: old_stf_spectrum

  ! new source time function
  character(len=max_length_string) :: stf_file
  double precision :: T_stf
  integer :: NSTEP_stf,NSTEP_stf_old
  double precision :: DT_stf
  real, dimension(:,:), pointer :: stf_trace_time
  real, dimension(:), allocatable :: stf_trace
  double complex, dimension(:), allocatable :: stf_spectrum

  ! specfem seismograms
  double precision :: T_seis
  character(len=2) :: band_instrument_code
  character(len=100) :: seisfile_extension
  integer :: NSTEP,ntrace
  double precision :: DT
  real :: t0
  real, dimension(:,:), pointer :: traces
  double complex, dimension(:,:), allocatable :: seis_spectra

  ! FFT
  integer :: n_FFT,nlog_FFT

  ! other stuff
  double precision :: DT_min
  integer :: NSTEP_max
  integer :: ios,istat,itrace,lu,n
  logical :: one_event_only,seismograms_are_bin,resample_seismograms,print_usage_and_stop,force_seisfile_extension,&
       differentiate_stf,stf_was_resized
  type (seismic_event) :: event
  character(len=character_length_evid) :: evid,evid_one_event_only
  type (seismic_station) :: station
  character(len=400) :: path_kernel_displacements,path_specfem_seismograms,path_output_convolved,file_output_trace

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  PROGRAM STARTS HERE
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  nullify(trans_coef,old_stf,old_stf_time,stf_trace_time,traces)

!------------------------------------------------------------------------
!  arguments
!
  call init(ap,myname,"Convolve standard SPECFEM3D Cartesian 3.0 output with a source-time function")
  call addPosarg(ap,"stf_file","sval","(text-) file containing trace of new source time function")
  call addPosarg(ap,"main_parfile","sval","Main parameter file of inversion")
  call addOption(ap,"-bicode",.true.,"(mandatory) bandcode and instrument code: the first two characters before "//&
       "the component in seismogram filename, e.g. 'LH' if your filenames look like 'network.staname.LH*.semd'",&
       "sval","")
  call addOption(ap,"-ori",.true.,"(mandatory) orientation : either 'NEZ' or 'XYZ', indicating the component "//&
       "orientations following band_instrument_code","sval","")
  call addOption(ap,"-dt",.true.,"(mandatory) time step of the seismograms (as in SPECFEM3D Par_file)","rval","0.0")
  call addOption(ap,"-nstep",.true.,"(mandatory) number of samples NSTEP as in SPECFEM3D Par_file","ival","0")
  call addOption(ap,"-ext",.true.,"(optional) NOT NEEDED FOR STANDARD FUNCTIONALITY! force a specific file "//&
       "extension: should be ANYTHING following the orientation character (including ALL dots etc.), e.g. '.semv' "//&
       "if your filenames look like 'network.staname.FX*.semv'. ","sval","")
  call addOption(ap,'-evid',.true.,"(optional) indicates a single event for which synthetic data is produced, "//&
       "otherwise synthetic data is produced for all events (as defined in ASKI FILE_EVENT_LIST)",'sval','')
  call addOption(ap,"-bin",.false.,"(optional) indicates whether SPECFEM trace files are binary files or not. "//&
       "For ascii output simply do not set option -bin")
  call addOption(ap,"-dconv",.false.,"(optional) if set, the source time function will be deconvolved from "//&
       "SPECFEM seismograms; consistend with 'ASKI_DECONVOLVE_STF = .true.' in Par_file_ASKI")
  call addOption(ap,"-opath",.true.,"(optional) output path (ending on '/') where convolved output traces will be "//&
       "written, relative to the respective specfem seismograms path (default: 'convolved/')","sval","convolved/")
!
  call parse(ap)
  if (.level.(.errmsg.ap) == 2) goto 3
!
  str = ap.sval."main_parfile"
  parfile = str
  if (.level.(.errmsg.ap) == 2) goto 3
!
  str = ap.sval."stf_file"
  stf_file = str
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
  ! -dt
  if(.not.(ap.optset."-dt")) then
     write(*,*) "ERROR: please indicate -dt"
     print_usage_and_stop = .true.
  else
     DT = dble(ap.rval."-dt")
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
  ! -ext
  force_seisfile_extension = ap.optset."-ext"
  if(force_seisfile_extension) then
     str = ap.sval."-ext"
     seisfile_extension = str
  end if
!
  ! -evid
  one_event_only = (ap.optset."-evid")
  if(one_event_only) then
     str = ap.sval."-evid"
     evid_one_event_only = str
  end if
!
  ! -bin
  seismograms_are_bin = (ap.optset."-bin")
!
  ! -dconv
  deconvolve_stf = (ap.optset."-dconv")
!
  ! -opath
  str = ap.sval."-opath"
  path_output_convolved = str
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
  ! create parfile of iterbasics by hand (since this is the only information needed from iterbasics,
  ! hence this is more robust if you want to do this operation before setting up any iter objects)
  call createKeywordsInputParameter(iterpar,iter_inpar_keys)
  call new(errmsg,myname)
  call readSubroutineInputParameter(iterpar,get(fuh),trim(.iterpath.invbasics)//&
       trim((.inpar.invbasics).sval.'PARFILE_ITERATION_STEP'),errmsg)
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
  ! when creating this program, one_event_only must always be .true.
  ! keep this construct anyway for future modification: could do that also for all events (or any subset)
  if(one_event_only) then
     errmsg = searchEventidSeismicEventList(.evlist.invbasics,evid_one_event_only)
     if(.level.errmsg /=0) call print(errmsg)
     if(.level.errmsg == 2) then
        write(*,*) "ERROR: eventID '"//trim(evid_one_event_only)//"' given by option -evid is not contained in event list"
        goto 3
     end if
  end if
!
  ! if orientation_is_NEZ, get all transformation matrices here
  if(orientation_is_NEZ) then
     allocate(trans_coef_all(3,3,.nstat.(.statlist.invbasics)))
     istat = 0
     do while (nextStationSeismicNetwork(.statlist.invbasics,station))
        istat = istat + 1
        ! transpose trans_coef here by shwitching XYZ with NEZ (i.e. coef_in = CX,CY,CZ and coef_out = N,E,UP),
        ! as we need the transpose in matmul operation when actually transforming later on
        trans_coef => transform(.comptrans.invbasics,(/'CX','CY','CZ'/),(/'N ','E ','UP'/),.staname.station)
        if(.not.associated(trans_coef)) then
           write(*,*) "ERROR: no transformation coefficients for ",istat,"'th station '"//trim(.staname.station)//"'"
           goto 1
        end if
        trans_coef_all(:,:,istat) = trans_coef
        deallocate(trans_coef)
     end do ! while next station
  end if
!
  T_seis = (NSTEP-1)*DT
!
!------------------------------------------------------------------------
!  write some info about this run now
!
  if(one_event_only) then
     write(*,*) "convolving SPECFEM3D seismograms of event '",trim(evid_one_event_only),"' for ",&
          .nstat.(.statlist.invbasics)," stations, "
  else
     write(*,*) "convolving SPECFEM3D seismograms of all ",.nev.(.evlist.invbasics)," events and ",&
          .nstat.(.statlist.invbasics)," stations, "
  end if
  write(*,*) "as of main parameter file '"//trim(parfile)//"'"
  write(*,*) ""
  write(*,*) "input SPECFEM3D seismograms: "
  write(*,*) "   NSTEP =  ",NSTEP
  write(*,*) "   DT =  ",real(DT),"      --> seismogram length T =  ",T_seis
  write(*,*) "   band and instrument code = ",band_instrument_code
  if(orientation_is_NEZ) then
     write(*,*) "   seismogram orientation = NEZ"
  else
     write(*,*) "   seismogram orientation = XYZ"
  end if
  write(*,*) ""
  write(*,*) "new source time function will be read from file '",trim(stf_file),"'"
  write(*,*) "the convolution will be done in frequency domain after applying an FFT to the seismograms and the new stf;"
  if(deconvolve_stf) then
     if(.not.force_seisfile_extension) then
        write(*,*) "   source time function will be deconvolved from synthetics first (producing displacement ",&
             "w.r.t. dirac), before convolving with new stf"
     else
        write(*,*) "   source time function will be deconvolved from synthetics first, before convolving with ",&
             "new stf (stf will be differentiated in case of moment tensor sources, make sure you properly ",&
             "accounted for that when chosing the forced file extension '",trim(seisfile_extension),"')"
     end if
  else
     ! if no stf is deconvolved, simply read in displacement
     if(.not.force_seisfile_extension) then
        seisfile_extension = '.semd'
        write(*,*) "   NO deconvolution of source time function will be applied (before convolving with ",&
             "new stf). Reading in displacement, using seismogram file extension ",seisfile_extension
     else
        write(*,*) "   NO deconvolution of source time function will be applied (before convolving with ",&
             "new stf). Reading in seismogram files having forced extension ",seisfile_extension
     end if
  end if
  write(*,*) ""
  write(*,*) "after applying an inverse FFT, the new traces are written (possibly with a finer sampling) to path '",&
       trim(path_output_convolved),"' (relative to PATH_KERNEL_DISPLACEMENTS/kernel_displ_eventID_OUTPUT_FILES/ )"
  write(*,*) ""
!
  path_kernel_displacements = trim(.iterpath.invbasics)//trim(iterpar.sval.'PATH_KERNEL_DISPLACEMENTS')
!
!------------------------------------------------------------------------
!  process the new source time function (spectrum)
!
  write(*,*) "reading, and processing of new source time function"
!
  errmsg = readAsciiData(stf_file,get(fuh),stf_trace_time,2)
  call undo(fuh)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
  if(.not.associated(stf_trace_time)) then
     write(*,*) "ERROR: new source time function was not read from file, error should have been raised already"
     goto 1
  end if
  NSTEP_stf = size(stf_trace_time,1)
  T_stf = dble(stf_trace_time(NSTEP_stf,1)) - dble(stf_trace_time(1,1))
  DT_stf = T_stf / dble(NSTEP_stf-1)
  write(*,*) "   number of samples read from stf file = ",size(stf_trace_time,1)
  write(*,*) "   time step DT of stf (inferred from stf file) = ",real(DT_stf)
!
  ! check if new source time function has the same total duration T as the seismograms. (and the original 
  ! source time function) If not, cut the new stf (or append zeros) in order to get (approximately) the 
  ! same df = 1/T (only then the (de)convolution in the frequency domain can be done sensibly)
  if( abs(T_stf-T_seis) > DT_stf ) then
     stf_was_resized = .true.
     NSTEP_stf_old = NSTEP_stf
     NSTEP_stf = floor((T_seis / DT_stf) + 1.d0)
     ! properly upround the number, if remainding decimal fraction is larger than 0.5
     if( (T_seis / DT_stf) + 1.d0 - dble(NSTEP_stf) > 0.5d0) NSTEP_stf = NSTEP_stf + 1
     ! reallocate stf_trace_time
     stf_trace_time => reallocate(stf_trace_time,NSTEP_stf,2)
     if(NSTEP_stf > NSTEP_stf_old) stf_trace_time(NSTEP_stf_old+1:NSTEP_stf,:) = 0.0
     write(*,*) "   since the total duration of the new stf T = ",T_stf,&
          " differs significantly from the duration of the seismograms T = ",T_seis,", "
     write(*,*) "   a corrected number of samples have been chosen for the new stf, namely: ",NSTEP_stf,&
          " (either cut short or appended zeros)"
  else
     stf_was_resized = .false.
     write(*,*) "   since the total duration of the new stf T = ",T_stf,&
          " does not differ significantly from the duration of the seismograms T = ",T_seis,", "
     write(*,*) "   the length of the stf is not modified"
  end if
  write(*,*) ""
!
  ! check if the new stf and the SPECFEM3D seismograms have the same time sampling (in single precision)
  ! if not, resample the coarser sampled one to the finer sampling
  if(real(DT_stf) > real(DT)) then
     resample_seismograms = .false.
     DT_min = DT
     NSTEP_max = NSTEP
     write(*,*) "resampling the new source time function"
     write(*,*) "   in order to get the same frequency sampling after applying FFTs,"
     write(*,*) "   both stf and SPECFEM3D seismograms should not only have the same time length,"
     write(*,*) "   but the same time sampling as well (due to the need of extending the time"
     write(*,*) "   series to sample lengths of the form 2^n)."
     write(*,*) "   for stability reasons, the source time function (which is coarser sampled than"
     write(*,*) "   the SPECFEM3D seismograms) will be resampled at the finer seismogram sampling,"
     write(*,*) "   DT = ",real(DT_min)," now. The resampled stf will additionally be written to file"
     write(*,*) "   '",trim(path_specfem_seismograms)//trim(path_output_convolved)//"new_stf_resampled.txt'"
     call resampleByLinearInterpolation(stf_trace_time,DT_stf,DT_min,NSTEP_max) ! the first column (containing the time) is also resampled = interpolated
     errmsg = writeAsciiData(trim(path_specfem_seismograms)//trim(path_output_convolved)//&
          'new_stf_resampled.txt',get(fuh),stf_trace_time)
     call undo(fuh)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 1
     call dealloc(errmsg)
     if(.not.associated(stf_trace_time)) then
        write(*,*) "ERROR: new source time function (resampled) could not be written to file, ",&
             "error should have been raised already"
        goto 1
     end if
  else ! real(DT_stf) > real(DT)
     ! even if DT_stf == DT, define values for DT_min,NSTEP_max
     DT_min = DT_stf
     NSTEP_max = NSTEP_stf
     ! if strictly DT_stf < DT, the seismograms must be resampled below, otherwise nothing must be resampled
     if(real(DT_stf) < real(DT)) then
        resample_seismograms = .true.
        write(*,*) "resampling the SPECFEM3D seismograms"
        write(*,*) "   in order to get the same frequency sampling after applying FFTs,"
        write(*,*) "   both stf and SPECFEM3D seismograms should not only have the same time length,"
        write(*,*) "   but the same time sampling as well (due to the need of extending the time"
        write(*,*) "   series to sample lengths of the form 2^n)."
        write(*,*) "   for stability reasons, the SPECFEM3D seismograms and possibly the original stf (which have"
        write(*,*) "   a coarser sampling than the new source time function) will be resampled at the finer"
        write(*,*) "   stf sampling, DT = ",real(DT_min)," below. This will result in the output traces"
        write(*,*) "   having the finer time sampling."
     else
        resample_seismograms = .false.
     end if
  end if ! real(DT_stf) > real(DT)
  write(*,*) ""
!
  ! the fast fourier transforms must be done with a number of samples which is of form 2^n
  ! calculate here the next number of that form larger than NSTEP_max
  nlog_FFT = ceiling(alog(float(NSTEP_max))/alog(2.))
  n_FFT = 2**nlog_FFT
  write(*,*) "will conduct (inverse) Fast Fourier Transforms with an actual number of samples = ",n_FFT
  write(*,*) "which is the next number of form 2^n larger than the maximum number of samples ",NSTEP_max
  write(*,*) ""
!
  if(deconvolve_stf) allocate(old_stf_spectrum(n_FFT))
  ntrace = 3*.nstat.(.statlist.invbasics)
  allocate(seis_spectra(n_FFT,ntrace))
!
  allocate(stf_trace(NSTEP_max))
  stf_trace = stf_trace_time(:,2)
  write(*,*) "applying Fast Fourier Transform to new source time function"
  allocate(stf_spectrum(n_FFT))
  call applyDPFFTOneTrace(stf_trace,NSTEP_max,stf_spectrum,n_FFT,nlog_FFT)
  write(*,*) ""
!
!------------------------------------------------------------------------
!  loop on all events (even if one_event_only, for simplicity of coding) 
!
  do while (nextEventSeismicEventList(.evlist.invbasics,event))
!
     if(one_event_only) then
        evid = evid_one_event_only
     else
        evid = .evid.event
     end if
!
     path_specfem_seismograms = trim(path_kernel_displacements)//'kernel_displ_'//trim(evid)//'_OUTPUT_FILES/'
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
                   "' and will NOT differentiate the original source time function before deconvolving it"
           else
              write(*,*) "   this event has a single force source mechanism, so reading in ",&
                   "seismogram files having forced extension '",trim(seisfile_extension),&
                   "' and will NOT differentiate the original source time function before deconvolving it"
           end if
           differentiate_stf = .false.
        case(1)
           ! MOMENT TENSOR SOURCE
           if(.not.force_seisfile_extension) then
              seisfile_extension = '.semv'
              write(*,*) "   this event has a moment tensor source mechanism, so reading in velocity using ",&
                   "seismogram file extension '",trim(seisfile_extension),&
                   "' and WILL differentiate the original source time function before deconvolving it"
           else
              write(*,*) "   this event has a moment tensor source mechanism, so reading in ",&
                   "seismogram files having forced extension '",trim(seisfile_extension),&
                   "' and WILL differentiate the original source time function before deconvolving it"
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
          orientation_is_NEZ,seisfile_extension,seismograms_are_bin,get(fuh),t0=t0)
     call undo(fuh)
     if(.not.associated(traces)) then
        write(*,*) "no output produced for this event"
        goto 2
     end if
!
     if(deconvolve_stf) then
        if(associated(old_stf)) deallocate(old_stf)
        if(differentiate_stf) then
           write(*,*) "read and differentiate source time function of event '",trim(evid),"' from file '",&
                trim(path_specfem_seismograms),"plot_source_time_function.txt'"
           call getDiffStf(old_stf,NSTEP,real(DT),trim(path_specfem_seismograms)//'plot_source_time_function.txt',lu)
        else
           write(*,*) "read source time function of event '",trim(evid),"' from file '",&
                trim(path_specfem_seismograms),"plot_source_time_function.txt'"
           call getStf(old_stf,NSTEP,trim(path_specfem_seismograms)//'plot_source_time_function.txt',lu)
        end if
        if(.not.associated(old_stf)) then
           write(*,*) "no spectra produced for this event"
           goto 2
        end if
        ! also produce here BY ASSUMPTION the time column (start at time t0, such as the traces)
        if(associated(old_stf_time)) deallocate(old_stf_time)
        allocate(old_stf_time(NSTEP,2))
        old_stf_time(:,2) = old_stf
        old_stf_time(:,1) = (/ (sngl(t0 + (n-1)*DT), n=1,NSTEP) /)
     end if ! deconvolve_stf
!
     if(resample_seismograms) then
        ! resample the SPECFEM3D seismograms, as well as the originial sorce time function (if -dconv was set)
        if(deconvolve_stf) then
           call resampleByLinearInterpolation(old_stf_time,DT,DT_min,NSTEP_max) ! the first column (containing the time) is also resampled = interpolated
           errmsg = writeAsciiData(trim(path_specfem_seismograms)//trim(path_output_convolved)//&
                'orig_stf_resampled.txt',get(fuh),old_stf_time)
           call undo(fuh)
           if (.level.errmsg /= 0) call print(errmsg)
           !call print(errmsg)
           if (.level.errmsg == 2) goto 1
           call dealloc(errmsg)
           if(.not.associated(old_stf_time)) then
              write(*,*) "ERROR: original source time function (resampled) could not be written to file, ",&
                   "error should have been raised already"
              goto 1
           end if
        end if ! deconvolve_stf
!
        call resampleByLinearInterpolation(traces,DT,DT_min,NSTEP_max)
        if(.not.associated(traces)) then
           write(*,*) "ERROR: SPECFEM3D seismogram traces vanished in memory after resampling; ",&
                   "error should have been raised already"
           goto 2
        end if
     end if
!
     if(deconvolve_stf) then
        write(*,*) "applying Fast Fourier Transform to original source time function"
        call applyDPFFTOneTrace(old_stf,NSTEP_max,old_stf_spectrum,n_FFT,nlog_FFT)
     end if
     write(*,*) "applying Fast Fourier Transform to all ",ntrace," seismogram traces of this event"
     call applyDPFFTTraces(traces,NSTEP_max,ntrace,seis_spectra,n_FFT,nlog_FFT)
!
     ! now do the (de)convolution
     write(*,*) "looping on all traces applying the (de)convolution(s)"
     do itrace = 1,ntrace
        if(deconvolve_stf) then
           seis_spectra(:,itrace) = seis_spectra(:,itrace) * stf_spectrum / old_stf_spectrum
        else
           seis_spectra(:,itrace) = seis_spectra(:,itrace) * stf_spectrum 
        end if
     end do ! itrace
!
     write(*,*) "applying Inverse Fast Fourier Transform to all ",ntrace," seismogram traces of this event"
     call applyDPIFFTTraces(seis_spectra,n_FFT,nlog_FFT,ntrace,traces,NSTEP_max)
!
     write(*,*) "writing the convolved seismogram traces to output path '",trim(path_specfem_seismograms),&
          trim(path_output_convolved),"'"
     call writeTraces(traces,real(DT),NSTEP,ntrace,.statlist.invbasics,&
          trim(path_specfem_seismograms)//trim(path_output_convolved),&
          band_instrument_code,orientation_is_NEZ,seisfile_extension,seismograms_are_bin,get(fuh),t0=t0)
     call undo(fuh)

2    write(*,*) ""
     if(one_event_only) exit
  end do ! while next event
!
  write(*,*) "good bye"
!
!------------------------------------------------------------------------
!  clean up
!
1 call dealloc(invbasics)
  call dealloc(iterpar)
  call add(fuh,lu); call dealloc(fuh)
  call dealloc(ap)
  call dealloc(errmsg)
  if(associated(trans_coef)) deallocate(trans_coef)
  if(allocated(trans_coef_all)) deallocate(trans_coef_all)
  if(associated(old_stf)) deallocate(old_stf)
  if(allocated(old_stf_spectrum)) deallocate(old_stf_spectrum)
  if(associated(stf_trace_time)) deallocate(stf_trace_time)
  if(allocated(stf_trace)) deallocate(stf_trace)
  if(allocated(stf_spectrum)) deallocate(stf_spectrum)
  if(allocated(seis_spectra)) deallocate(seis_spectra)
  if(associated(traces)) deallocate(traces)
  if(allocated(seis_spectra)) deallocate(seis_spectra)
!
  stop
!
3 if(.level.(.errmsg.ap)>=1) call print(.errmsg.ap)
  call usage(ap)
  goto 1
end program convolveWithStf
!
!------------------------------------------------------------------------
!
! subroutine printhelp
!   print '(50(1h-))'
!   print *,'Usage:'
!   print *,''
!   print *,'  convolveWithStf [-h] -bicode band_instrument_code -ori orientation -ext file_extension -dt time_step'
!   print *,'     -nstep number_of_time_samples -evid eventID [-bin] [-dconvd] [-opath output_path] stf_file parfile'
!   print *,''
!   print *,'Arguments:'
!   print *,''
!   print *,"  stf_file: file name of (text-) file containing trace of new source time function"
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
!   print *,"  -dt time_step : time_step is the real number defining the time step of the seismograms (as in SPECFEM3D Par_file)"
!   print *,''
!   print *,"  -nstep number_of_time_samples : number_of_time_samples is the number of samples NSTEP as in SPECFEM3D Par_file"
!   print *,''
!   print *,'  -evid eventID : if set, eventID indicates the single event for which synthetic data is produced. otherwise,'
!   print *,'                  synthetic data is produced for all events (as defined in ASKI FILE_EVENT_LIST)'
!   print *,''
!   print *,'Optional options:'
!   print *,''
!   print *,'  -bin : indicates whether SPECFEM trace files are binary files or not. for ascii output simply do not set '//&
!        'option -bin'
!   print *,''
!   print *,'  -dconvd : if set, the derivative of the heaviside source time function will be deconvolved from the synthetics.'
!   print *,"            It is assumed that it was written to file 'plot_source_time_function.txt', i.e. flag"
!   print *,'            PRINT_SOURCE_TIME_FUNCTION was set to .true. in SPECFEM3D Par_file.'
!   print *,"            (-dconvd is consistend with 'ASKI_DECONVOLVE_STF = .true.' in Par_file_ASKI). In case ASKI_DECONVOLVE_STF"
!   print *,"            was .true. , flag -dconvd should be set, too!! In this case, velocity seismograms (i.e. extensions '.semv')"
!   print *,"            must be used!"
!   print *,''
!   print *,"'  -opath output_path : output_path must end on '/' (!) if set, the convolved output traces will be written to path"
!   print *,"                        output_path, which is assumed relative to the respective specfem seismograms path (usually"
!   print *,"                        subdirectory). If not set, the default is 'convolved/'. It is assumed that this path exists!!"
!   print *,''
!   print *,'  -h : print this help message'
!   print '(50(1h-))'
!   return
! end subroutine printhelp
