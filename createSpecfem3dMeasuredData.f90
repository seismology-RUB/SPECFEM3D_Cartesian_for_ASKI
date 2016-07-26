!----------------------------------------------------------------------------
!   Copyright 2013 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of ASKI version 0.3.
!
!   ASKI version 0.3 is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 2 of the License, or
!   (at your option) any later version.
!
!   ASKI version 0.3 is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with ASKI version 0.3.  If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
program createSpecfem3dMeasuredData
  use inversionBasics
  use iterationStepBasics
  use seismicEvent
  use seismicEventList
  use seismicStation
  use seismicNetwork
  use discreteFourierTransform
  use componentTransformation
  use asciiDataIO
  use commandLine
  use fileUnitHandler
  use errorMessage
  use mathConstants

  implicit none

  interface
     subroutine readTraces(traces,NSTEP,statlist,path_specfem_seismograms,band_instrument_code,&
          orientation_is_NEZ,seisfile_extension,seismograms_are_bin,lu)
       use seismicNetwork
       implicit none
       real, dimension(:,:), pointer :: traces
       integer :: NSTEP,lu
       type (seismic_network) :: statlist
       character(len=*) :: path_specfem_seismograms,band_instrument_code,seisfile_extension
       logical :: orientation_is_NEZ,seismograms_are_bin
     end subroutine readTraces
  end interface

  ! command line
  type (cmdLine) :: cl
  character(len=300) :: parfile,string

  ! basics
  type (file_unit_handler) :: fuh
  type (error_message) :: errmsg
  character(len=27) :: myname = 'createSpecfem3dMeasuredData'
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
  integer :: NSTEP
  real :: DT
  real, dimension(:,:), pointer :: traces

  ! fourier transformation
  type (discrete_fourier_transform) :: DFT
  complex, dimension(:,:), allocatable :: spectra,spectrum_rotated
  real :: df
  integer :: nfreq
  integer, dimension(:), pointer :: jf

  ! filtering
  logical :: apply_filter
  character(len=400) :: path_event_filter,path_station_filter
  complex, dimension(:), pointer :: event_filter,station_comp_filter

  ! other stuff
  integer :: istat,ios,lu
  logical :: one_event_only,seismograms_are_bin
  type (seismic_event) :: event
  character(len=character_length_evid) :: evid,evid_one_event_only
  type (seismic_station) :: station
  character(len=400) :: path_specfem_seismograms,path_measured_data,file_measured_data

  external printhelp

!------------------------------------------------------------------------
!  preliminary processing
!
  ! process command line
  call new(cl,10,1,'h bicode ori ext bin evid dt nstep filter ocomp','0 1 1 1 0 1 1 1 0 1',printhelp)
  parfile = clManarg(cl,1)
!
  if(.not.clOptset(cl,2)) then
     write(*,*) "ERROR: please indicate -bicode"
     call printhelp
     goto 1
  else
     band_instrument_code = clOptarg(cl,2)
  end if
!
  if(.not.clOptset(cl,3)) then
     write(*,*) "ERROR: please indicate -ori"
     call printhelp
     goto 1
  else
     select case(clOptarg(cl,3))
     case('XYZ'); orientation_is_NEZ = .false.
     case('NEZ'); orientation_is_NEZ = .true.
     case default
        write(*,*) "ERROR: invalid -ori value '"//trim(clOptarg(cl,3))//"'"
        call printhelp
        goto 1
     end select
  end if
!
  if(.not.clOptset(cl,4)) then
     write(*,*) "ERROR: please indicate -ext"
     call printhelp
     goto 1
  else
     seisfile_extension = clOptarg(cl,4)
  end if
!
  seismograms_are_bin = clOptset(cl,5)
!
  one_event_only = clOptset(cl,6)
  if(one_event_only) then
     evid_one_event_only = clOptarg(cl,6)
  end if
!
  if(.not.clOptset(cl,7)) then
     write(*,*) "ERROR: please indicate -dt"
     call printhelp
     goto 1
  else
     string = clOptarg(cl,7)
     read(string,*,iostat=ios) DT
     if(ios/=0) then
        write(*,*) "ERROR: could not read real value for DT from -dt input '"//trim(string)//"'"
        call printhelp
        goto 1
     end if
  end if  
!
  if(.not.clOptset(cl,8)) then
     write(*,*) "ERROR: please indicate -nstep"
     call printhelp
     goto 1
  else
     string = clOptarg(cl,8)
     read(string,*,iostat=ios) NSTEP
     if(ios/=0) then
        write(*,*) "ERROR: could not read integer value for NSTEP from -nstep input '"//trim(string)//"'"
        call printhelp
        goto 1
     end if
  end if  
!
  apply_filter = clOptset(cl,9)
!
  ! handle -ocomp
  if(clOptset(cl,10)) then
     string = clOptarg(cl,10)
     read(string,*,iostat=ios) ncomp
     if(ios/=0) then
        write(*,*) "ERROR: could not read integer number of components as first word of '-ocomp' input string '"&
             //trim(string)//"'"
        call printhelp
        stop
     end if
     if(ncomp.le.0) then
        write(*,*) "ERROR: number of components ",ncomp," (first word of '-ocomp' input string '"&
             //trim(string)//"' must be positive"
        call printhelp
        stop
     end if
     allocate(comp(ncomp))
     read(string,*,iostat=ios) ncomp,comp
     if(ios/=0) then
        write(*,*) "ERROR: could not read ",ncomp," components from '-ocomp' input string '"//trim(string)//"'"
        call printhelp
        stop
     end if
     do jcomp=1,ncomp
        if(.not.validComponent(comp(jcomp))) then
           write(*,*) "ERROR: ",jcomp,"'th component '"//trim(comp(jcomp))//"' not valid. Valid components are '"//&
                all_valid_components//"'"
           call printhelp
           stop
        end if ! .not.validComponent(comp(jcomp))
     end do ! comp
  else ! clOptset(cl,10)
     print *, "ERROR: please indicate '-ocomp'"
     call printhelp
     stop
  end if ! clOptset(cl,10)
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
        call printhelp
        goto 1
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
  call initiateForwardDFT(DFT,DT,0,NSTEP-1,jf*df,errmsg,hanning_taper=0.05)
  if(.level.errmsg /=0) call print(errmsg)
  if(.level.errmsg==2) goto 1
  call dealloc(errmsg)
!
  allocate(spectra(nfreq,3*.nstat.(.statlist.invbasics)),spectrum_rotated(nfreq,ncomp))
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
  write(*,*) "   seismogram file extension = ",seisfile_extension
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
  path_event_filter = (.inpar.invbasics).sval.'PATH_EVENT_FILTER'
  path_station_filter = (.inpar.invbasics).sval.'PATH_STATION_FILTER'
  if(apply_filter) then
     write(*,*) "   event filters will be applied, as defined by filter files in path '"//trim(path_event_filter)//"'"
     write(*,*) "   station filters will be applied, as defined by filter files in path '"//trim(path_station_filter)//"'"
  else
     write(*,*) "   no filtering, the specfem time series will simply be transformed to frequency domain"
  end if
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
     ! in case of filtering, read in event specific filter values
     if(apply_filter) then
        if(associated(event_filter)) deallocate(event_filter)
        errmsg = readAsciiData(trim(path_event_filter)//"filter_"//trim(evid),lu,event_filter,ndata=nfreq)
        if(.level.errmsg /=0) call print(errmsg)
        if(.level.errmsg == 2) then
           write(*,*) "ERROR: could not read source filter from ascii file '"//trim(path_event_filter)//&
                "filter_"//trim(evid)
           goto 1
        end if
     end if
!
     path_specfem_seismograms = trim(path_measured_data)//'data_'//trim(evid)//'_OUTPUT_FILES/'
!
     write(*,*) "read all traces of event '",trim(evid),"' from path '",trim(path_specfem_seismograms),"'"
     ! read in all traces
     if(associated(traces)) deallocate(traces)
     call readTraces(traces,NSTEP,.statlist.invbasics,path_specfem_seismograms,band_instrument_code,&
          orientation_is_NEZ,seisfile_extension,seismograms_are_bin,lu)
     if(.not.associated(traces)) then
        write(*,*) "no spectra produced for this event"
        goto 2
     end if
!
     write(*,*) "compute spectra from traces"
     ! fourier transform to frequency domain of all traces at once
     call new(errmsg,myname)
     call transformForwardDFT(DFT,traces,spectra,errmsg)
     if(.level.errmsg /=0) call print(errmsg)
     if(.level.errmsg==2) goto 1
     call dealloc(errmsg)
!
     ! write spectra to files
     write(*,*) "write measured data files"
     istat = 0
     do while (nextStationSeismicNetwork(.statlist.invbasics,station))
        istat = istat + 1
!
        ! transform all 3 traces of this station (assume order X,Y,Z or N,E,Z coming from routine readTraces,
        ! dependent on orientation_is_NEZ) to the requested components in comp
        spectrum_rotated(:,:) = matmul(spectra(:,(istat-1)*3+1:(istat-1)*3+3) , trans_coef_all_transpose(:,:,istat))
!
        do jcomp = 1,ncomp
!
           ! in case of filtering, read in station_comp_filter for this path and component
           if(apply_filter) then
              if(associated(station_comp_filter)) deallocate(station_comp_filter)
              errmsg = readAsciiData(trim(path_station_filter)//"filter_"//trim(.staname.station)//"_"//trim(comp(jcomp)),&
                   lu,station_comp_filter,ndata=nfreq)
              if(.level.errmsg /=0) call print(errmsg)
              if(.level.errmsg == 2) then
                 write(*,*) "ERROR: could not read station component filter from ascii file '"//trim(path_station_filter)//&
                      "filter_"//trim(.staname.station)//"_"//trim(comp(jcomp))
                 goto 1
              end if
           end if
!
           ! define filename of output file
           file_measured_data = "data_"//trim(evid)//"_"//trim(.staname.station)//"_"//trim(comp(jcomp))
!
           write(*,*) "writing measured data file '",trim(file_measured_data),"'"
           if(apply_filter) then
              errmsg = writeAsciiData(trim(path_measured_data)//file_measured_data,lu,&
                   event_filter*station_comp_filter*spectrum_rotated(:,jcomp))
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
1 call dealloc(invbasics); call dealloc(iterbasics)
  call add(fuh,lu); call dealloc(fuh)
  call dealloc(cl)
  call dealloc(DFT)
  call dealloc(errmsg)
  if(associated(trans_coef)) deallocate(trans_coef)
  if(allocated(trans_coef_all_transpose)) deallocate(trans_coef_all_transpose)
  if(allocated(spectra)) deallocate(spectra)
  if(allocated(spectrum_rotated)) deallocate(spectrum_rotated)
  if(associated(jf)) deallocate(jf)
  if(associated(traces)) deallocate(traces)
  if(associated(event_filter)) deallocate(event_filter)
  if(associated(station_comp_filter)) deallocate(station_comp_filter)
!
  write(*,*) "good bye"
end program createSpecfem3dMeasuredData
!
!-----------------------------------------------------------------------------------------------------------------
!
subroutine readTraces(traces,NSTEP,statlist,path_specfem_seismograms,band_instrument_code,&
     orientation_is_NEZ,seisfile_extension,seismograms_are_bin,lu)
  use seismicNetwork
  implicit None
  real, dimension(:,:), pointer :: traces
  integer :: NSTEP,lu
  type (seismic_network) :: statlist
  character(len=*) :: path_specfem_seismograms,band_instrument_code,seisfile_extension
  logical :: orientation_is_NEZ,seismograms_are_bin
  ! local
  integer :: itrace,icomp,ios,isamp
  type(seismic_station) :: station
  character(len=500) :: specfem_seismogram_file
  character(len=1), dimension(3) :: orientation
  real :: rdummy
!
  allocate(traces(NSTEP,3*.nstat.statlist))
  if(orientation_is_NEZ) then
     orientation = (/'N','E','Z'/)
  else
     orientation = (/'X','Y','Z'/)
  end if
!
  itrace = 0
  do while(nextStationSeismicNetwork(statlist,station))
!
     do icomp = 1,3
!
        specfem_seismogram_file = &
             trim(.staname.station)//"."//&
             trim(.netcode.station)//"."//&
             trim(band_instrument_code)//orientation(icomp)//&
             trim(seisfile_extension)
!
        if(seismograms_are_bin) then
!
           write(*,*) "reading binary seismogram file '",trim(specfem_seismogram_file),"'"
           open(unit=lu,file=trim(path_specfem_seismograms)//trim(specfem_seismogram_file),&
                status='old',form='unformatted',access='direct',&
                recl=4*NSTEP,iostat=ios)
           if(ios/=0) then
              print *,"ERROR: could not open file"
              close(lu)
              deallocate(traces); nullify(traces)
              return
           endif ! ios/=0
           ! read seismogram from file
           itrace = itrace + 1
           read(lu,rec=1,iostat=ios) traces(:,itrace)
           if(ios/=0) then
              print *,"ERROR: could not read trace from file"
              close(lu)
              deallocate(traces); nullify(traces)
              return
           endif ! ios/=0
           close(lu)
!
        else ! seismograms_are_bin
!
           write(*,*) "reading ascii seismogram file '",trim(specfem_seismogram_file),"'"
           open(unit=lu,file=trim(path_specfem_seismograms)//trim(specfem_seismogram_file),&
                status='old',form='formatted',action='read',iostat=ios)
           if(ios/=0) then
              print *,"ERROR: could not open file"
              close(lu)
              deallocate(traces); nullify(traces)
              return
           endif ! ios/=0
           ! read seismogram from file
           itrace = itrace + 1
           do isamp=1,NSTEP
              read(lu,*,iostat=ios) rdummy,traces(isamp,itrace)
              if(ios/=0) then
                 print *,"ERROR: could not read sample ",isamp," of trace from file"
                 close(lu)
                 deallocate(traces); nullify(traces)
                 return
              endif ! ios/=0
           end do ! isamp
           close(lu)
!
        end if ! seismograms_are_bin
!
     end do ! icomp
!
  end do ! while next station
!  
end subroutine readTraces
!
!------------------------------------------------------------------------
!
subroutine printhelp
  use componentTransformation
  print '(50(1h-))'
  print *,'Usage:'
  print *,''
  print *,'  createSpecfem3dMeasuredData [-h] -bicode band_instrument_code -ori orientation -ext file_extension'
  print *,'     -dt time_step -nstep number_of_time_samples -ocomp "ncomp comp" [-filter] [-evid eventID] [-bin] parfile'
  print *,''
  print *,'Arguments:'
  print *,''
  print *,"  parfile: main parameter file of inversion"
  print *,''
  print *,'Mandatory options:'
  print *,''
  print *,'  -bicode band_instrument_code : band_instrument_code must be two characters, band code and instrument code'
  print *,'                                 i.e. the first two characters before the component in seismogram filename'
  print *,"                                 e.g. 'LH' if your filenames look like 'staname.network.LH*.semd.ascii'"
  print *,''
  print *,"  -ori orientation : either 'NEZ' or 'XYZ', indicating the component orientations following band_instrument_code"
  print *,''
  print *,'  -ext file_extension : file_extension should be ANYTHING following the orientation character (including ALL dots etc.)'
  print *,"                                 e.g. if your filenames look like 'staname.network.FX*.semv', file_extension = '.semv'"
  print *,"  -dt time_step : time_step is the real number defining the time step of the seismograms (as in SPECFEM3D Par_file)"
  print *,''
  print *,"  -nstep number_of_time_samples : number_of_time_samples is the number of samples NSTEP as in SPECFEM3D Par_file"
  print *,''
  print *,"  -ocomp : defines the components for which measured data output is produced"
  print *,"          ncomp: number of components,  comp: ncomp components (valid components: '"//&
       trim(all_valid_components)//"')"
  print *,''
  print *,'Optional options:'
  print *,''
  print *,'  -filter : if set, the respective event filters and station (component) filters will be applied to the spectra'
  print *,''
  print *,'  -evid eventID : if set, eventID indicates the single event for which measured data is produced. otherwise,'
  print *,'                  measured data is produced for all events (as defined in ASKI FILE_EVENT_LIST)'
  print *,''
  print *,'  -bin : indicates whether SPECFEM trace files are binary files or not. for ascii output simply do not set '//&
       'option -bin'
  print *,''
  print *,'  -h : print this help message'
  print '(50(1h-))'
  return
end subroutine printhelp
