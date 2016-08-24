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
program transformSpecfem3dCartesianSyntheticData
  use specfem3dForASKI_mod
  use inversionBasics
  use iterationStepBasics
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
  character(len=40) :: myname = 'transformSpecfem3dCartesianSyntheticData'
  type (inversion_basics) :: invbasics
  type (iteration_step_basics) :: iterbasics

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
  integer :: ifreq,nfreq
  integer, dimension(:), pointer :: jf

  double precision :: unit_factor,inverse_unit_factor

  ! other stuff
  integer :: istat,lu
  logical :: one_event_only,seismograms_are_bin,deconvolve_stf,differentiate_stf,diff_time_series,print_usage_and_stop
  type (seismic_event) :: event
  character(len=character_length_evid) :: evid,evid_one_event_only
  type (seismic_station) :: station
  character(len=400) :: path_specfem_seismograms,path_synthetic_data,file_synthetic_data

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  PROGRAM STARTS HERE
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  nullify(str_vec,trans_coef,traces,stf,jf)

!------------------------------------------------------------------------
!  arguments
!
  call init(ap,myname,"Transform standard SPECFEM3D Cartesian 3.0 output to ASKI 1.0-1.2 spectral data in synthetic-data format")
  call addPosarg(ap,"main_parfile","sval","Main parameter file of inversion")
  call addOption(ap,"-bicode",.true.,"(mandatory) bandcode and instrument code: the first two characters before "//&
       "the component in seismogram filename, e.g. 'LH' if your filenames look like 'network.staname.LH*.semd'",&
       "sval","")
  call addOption(ap,"-ori",.true.,"(mandatory) orientation : either 'NEZ' or 'XYZ', indicating the component "//&
       "orientations following band_instrument_code","sval","")
  call addOption(ap,"-dt",.true.,"(mandatory) time step of the seismograms (as in SPECFEM3D Par_file)","rval","0.0")
  call addOption(ap,"-nstep",.true.,"(mandatory) number of samples NSTEP as in SPECFEM3D Par_file","ival","0")
  call addOption(ap,"-ocomp",.true.,"(mandatory) receiver components for which synthetic data is produced; valid "//&
       "components: '"//trim(all_valid_components)//"')","svec","")
  call addOption(ap,'-evid',.true.,"(optional) indicates a single event for which synthetic data is produced, "//&
       "otherwise synthetic data is produced for all events (as defined in ASKI FILE_EVENT_LIST)",'sval','')
  call addOption(ap,"-dconv",.false.,"(optional) if set, the source time function will be deconvolved from "//&
       "SPECFEM seismograms; consistend with 'ASKI_DECONVOLVE_STF = .true.' in Par_file_ASKI")
  call addOption(ap,"-bin",.false.,"(optional) indicates whether SPECFEM trace files are binary files or not. "//&
       "For ascii output simply do not set option -bin")
  call addOption(ap,"-ext",.true.,"(optional) NOT NEEDED FOR STANDARD FUNCTIONALITY! force a specific file "//&
       "extension: should be ANYTHING following the orientation character (including ALL dots etc.), e.g. '.semv' "//&
       "if your filenames look like 'network.staname.FX*.semv'. ","sval","")
  call addOption(ap,"-diffts",.false.,"(optional) NOT NEEDED FOR STANDARD FUNCTIONALITY! if set, the time "//&
       "series will be differentiated (by simple first order central differences) before further processing")
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
  ! -evid
  one_event_only = (ap.optset."-evid")
  if(one_event_only) then
     str = ap.sval."-evid"
     evid_one_event_only = str
  end if
!
  ! -dconv
  deconvolve_stf = (ap.optset."-dconv")
!
  ! -bin
  seismograms_are_bin = (ap.optset."-bin")
!
  ! -ext
  force_seisfile_extension = ap.optset."-ext"
  if(force_seisfile_extension) then
     str = ap.sval."-ext"
     seisfile_extension = str
  end if
!
  ! -diffts
  diff_time_series = (ap.optset."-diffts")
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
  ! setup iteration step basics
  call new(errmsg,myname)
  call init(iterbasics,invbasics,fuh,errmsg)
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
     ! transpose trans_coef here by shwitching XYZ with NEZ (i.e. coef_in = CX,CY,CZ and coef_out = N,E,UP),
     ! as we need the transpose in matmul operation when actually transforming later on
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
  unit_factor = .ufmdata.invbasics
  inverse_unit_factor = 1.d0/unit_factor
!
  ! compute fourier transformation factors efactors here
  ! taper parameters like in SPECFEM3D-for-ASKI codes
  df = rval(.inpar.invbasics,'MEASURED_DATA_FREQUENCY_STEP')
  nfreq = ival(.inpar.iterbasics,'ITERATION_STEP_NUMBER_OF_FREQ')
  jf => ivecp(.inpar.iterbasics,'ITERATION_STEP_INDEX_OF_FREQ',nfreq)
  if(.not.associated(jf)) then
     write(*,*) "ERROR: could not read ",nfreq," frequency indices from vector 'ITERATION_STEP_INDEX_OF_FREQ' "//&
          "in iteration step parfile"
     goto 1
  end if
  call new(errmsg,myname)
  call initiateForwardDFT(DFT,DT,0,NSTEP-1,jf*df,errmsg,hanning_taper=0.05)
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
     write(*,*) "creating ASKI synthetic data from SPECFEM3D seismograms for one event and ",&
          .nstat.(.statlist.invbasics)," stations, "
  else
     write(*,*) "creating ASKI synthetic data from SPECFEM3D seismograms for ",.nev.(.evlist.invbasics)," events and ",&
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
  write(*,*) "   frequency step = ",df
  write(*,*) "   frequency indices = ",jf
  write(*,*) "   output spectra will be produced for the ",ncomp," receiver components: ",comp//", "
  if(diff_time_series) then
     write(*,*) "   time series will be differentiated by simple 1st order finite differences before further "//&
          "processing. THIS IS NOT NEEDED FOR STANDARD FUNCTIONALITY! make sure you know what you are doing!"
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
  write(*,*) "   the output spectra are computed according to the unit factor ",unit_factor,&
       ", as of parameter 'UNIT_FACTOR_MEASURED_DATA' in main parfile"
  path_synthetic_data = trim(.iterpath.invbasics)//trim((.inpar.iterbasics).sval.'PATH_SYNTHETIC_DATA')
  write(*,*) "   synthetic data files will be written to path '",trim(path_synthetic_data),"'"
  write(*,*) ""
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
     path_specfem_seismograms = trim(.iterpath.invbasics)//trim((.inpar.iterbasics).sval.'PATH_KERNEL_DISPLACEMENTS')//&
          'kernel_displ_'//trim(evid)//'_OUTPUT_FILES/'
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
     write(*,*) "write synthetic data files to path '",trim(path_synthetic_data),"'"
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
           ! define filename of output file
           file_synthetic_data = "synthetics_"//trim(evid)//"_"//trim(.staname.station)//"_"//trim(comp(jcomp))

           write(*,*) "writing synthetic data file '",trim(file_synthetic_data),"'"
           errmsg = writeAsciiData(trim(path_synthetic_data)//file_synthetic_data,lu,spectrum_rotated(:,jcomp))
           if(.level.errmsg/=0) call print(errmsg)
           if(.level.errmsg==2) goto 1
           call dealloc(errmsg)
        end do ! jcomp
!
     end do ! while next station
!
2    write(*,*) ""
     if(one_event_only) exit
  end do ! while next event
!
  write(*,*) "good bye"
!
!------------------------------------------------------------------------
!  clean up
!
1 call dealloc(invbasics); call dealloc(iterbasics)
  call add(fuh,lu); call dealloc(fuh)
  call dealloc(ap)
  call dealloc(DFT)
  call dealloc(errmsg)
  if(allocated(comp)) deallocate(comp)
  if(associated(trans_coef)) deallocate(trans_coef)
  if(allocated(trans_coef_all_transpose)) deallocate(trans_coef_all_transpose)
  if(associated(traces)) deallocate(traces)
  if(associated(stf)) deallocate(stf)
  if(allocated(spectra)) deallocate(spectra)
  if(allocated(spectrum_rotated)) deallocate(spectrum_rotated)
  if(allocated(stf_spectrum)) deallocate(stf_spectrum)
  if(associated(jf)) deallocate(jf)
  if(allocated(stf_spectrum)) deallocate(stf_spectrum)
!
  stop
!
3 if(.level.(.errmsg.ap)>=1) call print(.errmsg.ap)
  call usage(ap)
  goto 1
end program transformSpecfem3dCartesianSyntheticData
!
!------------------------------------------------------------------------
!
! subroutine printhelp
!   use componentTransformation, only: all_valid_components
!   print '(50(1h-))'
!   print *,'Usage:'
!   print *,''
!   print *,'  transformSpecfem3dCartesianSyntheticData [-h] -bicode band_instrument_code -ori orientation -dt time_step'
!   print *,'     -nstep number_of_time_samples -ocomp "ncomp comp" [-evid eventID] [-bin] [-dconv] [-uf unit_factor] parfile'
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
!   print *,"  -dt time_step : time_step is the real number defining the time step of the seismograms (as in SPECFEM3D Par_file)"
!   print *,''
!   print *,"  -nstep number_of_time_samples : number_of_time_samples is the number of samples NSTEP as in SPECFEM3D Par_file"
!   print *,''
!   print *,"  -ocomp : defines the receiver components for which synthetic data output is produced."
!   print *,"           ncomp: number of components,  comp: ncomp components (valid components: '"//&
!        trim(all_valid_components)//"')"
!   print *,''
!   print *,'Optional options:'
!   print *,''
!   print *,'  -bin : indicates whether SPECFEM trace files are binary files or not. for ascii output simply do not set '//&
!        'option -bin'
!   print *,''
!   print *,'  -evid eventID : if set, eventID indicates the single event for which synthetic data is produced. otherwise,'
!   print *,'                  synthetic data is produced for all events (as defined in ASKI FILE_EVENT_LIST)'
!   print *,''
!   print *,'  -dconv : if set, the source time function will be deconvolved from the synthetics.'
!   print *,"            It is assumed that it was written to file 'plot_source_time_function.txt', i.e. flag"
!   print *,'            PRINT_SOURCE_TIME_FUNCTION was set to .true. in SPECFEM3D Par_file.'
!   print *,"            (-dconv is consistend with 'ASKI_DECONVOLVE_STF = .true.' in Par_file_ASKI). In case"
!   print *,"            ASKI_DECONVOLVE_STF was .true. , flag -dconv should be set, too!"
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
!   print *,'  -h : print this help message'
!   print '(50(1h-))'
!   return
! end subroutine printhelp
