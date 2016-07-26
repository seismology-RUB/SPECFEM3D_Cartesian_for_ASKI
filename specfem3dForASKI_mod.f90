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
module specfem3dForASKI_mod
!
  use realloc
!
  implicit none
!
  interface diffTraces
     module procedure diffOneTrace
     module procedure diffManyTraces
  end interface diffTraces
!
contains
!
!-----------------------------------------------------------------------------------------------------------------
!
subroutine diffOneTrace(trace,DT)
  real, dimension(:) :: trace
  real :: DT
  ! local
  double precision :: one_over_DT,one_over_2DT
  real ::  tmp_sample,rtmp
  integer :: isamp,NSTEP
!
  NSTEP = size(trace)
  if(NSTEP <= 1) then
     write(*,*) "ERROR in diffOneTrace: there is only ",NSTEP," incoming samples of time series, cannot differentiate"
     stop
  end if
!
  one_over_DT = 1.d0/dble(DT)
  one_over_2DT = 0.5d0/dble(DT)
!
  ! treat first sample separately (right finite difference)
  tmp_sample = trace(1) ! remember current sample for differentiation of the next sample
  trace(1) = real(  (dble(trace(2))-dble(trace(1)))*one_over_DT  )
!
  if(NSTEP <= 2) goto 1 ! if there are only two samples, immediately jump to last sample after loop
!
  ! differentiate time series trace by simple central differences
  do isamp = 2,NSTEP-1
     rtmp = real(  (dble(trace(isamp+1))-dble(tmp_sample))*one_over_2DT )
     tmp_sample = trace(isamp)
     trace(isamp) = rtmp
  end do ! isamp
!
  ! treat last sample separately (left finite difference)
1 trace(NSTEP) = real(  (dble(trace(NSTEP))-dble(tmp_sample))*one_over_DT  )
end subroutine diffOneTrace
!
!-----------------------------------------------------------------------------------------------------------------
!
subroutine diffManyTraces(traces,DT)
  real, dimension(:,:) :: traces
  real :: DT
  ! local
  integer :: NSTEP,itrace,ntrace
!
  NSTEP = size(traces,1)
  ntrace = size(traces,2)
!
  if(NSTEP <= 1) then
     write(*,*) "ERROR in diffManyTraces: there is only ",NSTEP," incoming samples of time series, cannot differentiate"
     stop
  end if
!
  if(ntrace < 1) then
     write(*,*) "ERROR in diffManyTraces: there is no incoming time series to differentiate"
     stop
  end if
!
  do itrace = 1,ntrace
     call diffOneTrace(traces(:,itrace),DT)
  end do ! itrace
end subroutine diffManyTraces
!
!-----------------------------------------------------------------------------------------------------------------
!
subroutine getStf(stf,NSTEP,filename,lu,normalize)
  implicit none
  real, dimension(:), pointer :: stf
  integer :: NSTEP,lu
  logical, optional :: normalize
  character(len=*) :: filename
  ! local
  integer :: ios,isamp
  real :: max_stf,rtmp
!
  allocate(stf(NSTEP))
!
  ! read source time function from file
!
  open(unit=lu,file=filename,status='old',form='formatted',action='read',iostat=ios)
  if(ios/=0) then
     write(*,*) "ERROR in getStf(): could not open file '"//trim(filename)//"'"
     close(lu)
     deallocate(stf); nullify(stf)
     return
  endif ! ios/=0
!
  do isamp=1,NSTEP
     read(lu,*,iostat=ios) rtmp,stf(isamp)
     if(ios/=0) then
        write(*,*) "ERROR in getStf(): could not read sample ",isamp," of trace from file '"//trim(filename)//"'"
        close(lu)
        deallocate(stf); nullify(stf)
        return
     endif ! ios/=0
  end do ! isamp
  close(lu)
!
  if(present(normalize)) then
     if(normalize) then
        ! normalize source time function 
        ! (needed in case of SPECFEM3D GLOBE, where the error function is additionally multiplied by the seismic moment, 
        ! in order to plot a real moment rate function (I guess). However this factor is NOT used in the internal source term!)
        ! In turn, this means that ANY other (self-implemented) sourrce time functions MUST be used in a normalized fashion with SPECFEM3D!
        max_stf = maxval(stf)
        if(max_stf == 0.0) then
           write(*,*) "ERROR in getStf(): maximum value of source time function in file '"//trim(filename)//&
                "' is 0.0 , cannot normalize"
           deallocate(stf); nullify(stf)
           return     
        else
           write(*,*) "in getStf(): maximum value of source time function in file '"//&
                trim(filename)//"' is ",max_stf," , normalize stf by division through this value"
        end if
        stf = stf/max_stf
     end if ! normalize
  end if ! present(normalize)
end subroutine getStf
!
!-----------------------------------------------------------------------------------------------------------------
!
subroutine getDiffStf(stf_diff,NSTEP,DT,filename,lu,normalize)
  implicit none
  real, dimension(:), pointer :: stf_diff
  integer :: NSTEP,lu
  real :: DT
  character(len=*) :: filename
  logical, optional :: normalize
  ! local
!
  call getStf(stf_diff,NSTEP,filename,lu,normalize)
  if(.not.associated(stf_diff)) return
!
  call diffOneTrace(stf_diff,DT)
end subroutine getDiffStf
!
!-----------------------------------------------------------------------------------------------------------------
!
subroutine readTraces(traces,NSTEP,statlist,path_specfem_seismograms,band_instrument_code,&
     orientation_is_NEZ,seisfile_extension,seismograms_are_bin,lu,t0)
  use seismicNetwork
  use seismicStation
  implicit None
  real, dimension(:,:), pointer :: traces
  real, optional :: t0
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
             trim(.netcode.station)//"."//&
             trim(.staname.station)//"."//&
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
              write(*,*) "ERROR: could not open file"
              close(lu)
              deallocate(traces); nullify(traces)
              return
           endif ! ios/=0
           ! read seismogram from file
           itrace = itrace + 1
           read(lu,rec=1,iostat=ios) traces(:,itrace)
           if(ios/=0) then
              write(*,*) "ERROR: could not read trace from file"
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
              write(*,*) "ERROR: could not open file"
              close(lu)
              deallocate(traces); nullify(traces)
              return
           endif ! ios/=0
           ! read seismogram from file
           itrace = itrace + 1
           do isamp=1,NSTEP
              read(lu,*,iostat=ios) rdummy,traces(isamp,itrace)
              if(ios/=0) then
                 write(*,*) "ERROR: could not read sample ",isamp," of trace from file"
                 close(lu)
                 deallocate(traces); nullify(traces)
                 return
              endif ! ios/=0
              if(itrace == 1 .and. isamp == 1) then
                 if(present(t0)) t0 = rdummy
              end if
           end do ! isamp
           close(lu)
!
        end if ! seismograms_are_bin
!
     end do ! icomp
!
  end do ! while next station
!
  if(seismograms_are_bin) t0 = 0.0
!  
end subroutine readTraces
!
!-----------------------------------------------------------------------------------------------------------------
!
subroutine writeTraces(traces,DT,NSTEP,ntrace,statlist,output_path,band_instrument_code,&
     orientation_is_NEZ,seisfile_extension,seismograms_are_bin,lu,t0)
  use seismicNetwork
  use seismicStation
  implicit None
  real, dimension(NSTEP,ntrace) :: traces
  real :: DT
  real, optional :: t0
  integer :: NSTEP,ntrace,lu
  type (seismic_network) :: statlist
  character(len=*) :: output_path,band_instrument_code,seisfile_extension
  logical :: orientation_is_NEZ,seismograms_are_bin
  ! local
  integer :: itrace,icomp,ios,isamp
  type(seismic_station) :: station
  character(len=500) :: specfem_seismogram_file
  character(len=1), dimension(3) :: orientation
  real :: t0_used
!
  if(ntrace /= 3*.nstat.statlist) then
     write(*,*) "ERROR in writeTraces: incoming number of traces = ",ntrace,&
          " does not equal the number of receiver components contained in the station list = ",3*.nstat.statlist
     return
  end if
!
  if(present(t0)) then
     t0_used = t0
  else
     t0_used = 0.0
  end if
!
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
             trim(.netcode.station)//"."//&
             trim(.staname.station)//"."//&
             trim(band_instrument_code)//orientation(icomp)//&
             trim(seisfile_extension)
!
        if(seismograms_are_bin) then
!
           write(*,*) "writing binary seismogram file '",trim(specfem_seismogram_file),"'"
           open(unit=lu,file=trim(output_path)//trim(specfem_seismogram_file),&
                status='new',form='unformatted',access='direct',&
                recl=4*NSTEP,iostat=ios)
           if(ios/=0) then
              write(*,*) "ERROR: could not open file"
              close(lu)
              return
           endif ! ios/=0
           ! write seismogram to file
           itrace = itrace + 1
           write(lu,rec=1,iostat=ios) traces(:,itrace)
           if(ios/=0) then
              write(*,*) "ERROR: could not write trace to file"
              close(lu)
              return
           endif ! ios/=0
           close(lu)
!
        else ! seismograms_are_bin
!
           write(*,*) "writing ascii seismogram file '",trim(specfem_seismogram_file),"'"
           open(unit=lu,file=trim(output_path)//trim(specfem_seismogram_file),&
                status='new',form='formatted',action='write',iostat=ios)
           if(ios/=0) then
              write(*,*) "ERROR: could not open file"
              close(lu)
              return
           endif ! ios/=0
           ! read seismogram from file
           itrace = itrace + 1
           do isamp=1,NSTEP
              write(lu,*,iostat=ios) real(dble(t0_used) + dble(isamp-1)*dble(DT)),' ',traces(isamp,itrace)
              if(ios/=0) then
                 write(*,*) "ERROR: could not write sample ",isamp," of trace to file"
                 close(lu)
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
end subroutine writeTraces
!
!-----------------------------------------------------------------------------------------------------------------
!
subroutine resampleByLinearInterpolation(traces,old_dt,new_dt,new_NSTEP)
  implicit none
!
  real, dimension(:,:), pointer :: traces
  double precision :: old_dt,new_dt
  integer :: new_NSTEP
  ! local
  integer :: ntrace,itrace,isamp,NSTEP
  integer, dimension(:), allocatable :: isamp_old_left
  double precision, dimension(:), allocatable :: w_interpl
  double precision :: h,new_dt_over_old_dt
!
!-----------------------------------------------------------------------
!  subroutine starts here
!
  NSTEP = size(traces,1)
  ntrace = size(traces,2)
  ! ASSUME HERE, THAT NSTEP > 0 and ntrace > 0 and associated(traces)

  new_dt_over_old_dt = new_dt / old_dt

  ! find the intervals of the old samples, in which the new samples fall
  ! calculate the interpolation weights
  allocate(isamp_old_left(new_NSTEP),w_interpl(new_NSTEP))
  do isamp = 1,new_NSTEP
     ! current time in new sampling
     ! t = (isamp-1)*new_dt

     ! left sample is floor( t / old_dt ) + 1   ;  be aware that array indices start at 1. so + 1 here
     ! equivalent but more performant calculation: 
     h = (isamp-1)*new_dt_over_old_dt
     isamp_old_left(isamp) = floor(h) + 1

     ! isamp_old_right is assumed to be isamp_old_left + 1 below, so don't need to define in additional array

     ! w_interpl(isamp) shall be the weight assigned to the left sample, which is equal to 
     ! 1.d0 - (t - (isamp_old_left(isamp)-1)*old_dt) / old_dt   , equivalently:
      w_interpl(isamp) = -h + dble(isamp_old_left(isamp))
   end do ! isamp

   ! finally calculate the new interpolated values

   ! if there are more new samples than old ones, we need to reallocate BEFORE assigning new values
   if(new_NSTEP > NSTEP) traces => reallocate(traces,new_NSTEP,ntrace)

   do itrace = 1,ntrace
      traces(1:new_NSTEP,itrace) = &
           traces(isamp_old_left,itrace)*w_interpl &
           + traces(isamp_old_left + 1,itrace) * (1.d0 - w_interpl)
   end do ! itrace

   ! if there are less new samples than old ones, we need to reallocate AFTER assigning new values, in order
   ! not to throw away information too soon
   if(new_NSTEP < NSTEP) traces => reallocate(traces,new_NSTEP,ntrace)

   deallocate(isamp_old_left,w_interpl)

 end subroutine resampleByLinearInterpolation
!  
end module specfem3dForASKI_mod
