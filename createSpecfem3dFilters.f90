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
!
!
!#############################################################################################
!#############################################################################################
!## 
!##   A T T E N T I O N
!## 
!##  This main program is not yet fully developed, in the sense that the required information
!##  is implemented in a HARD-CODED way in the first lines of the program (df, nf, ifreq,
!##  dt, nt1, nt2, name of source-time-function file = 'ricker.txt'). The program should be 
!##  improved in such a way, that these values are read from the command line, or infered from 
!##  parameter files etc. (do not hesitate to do so!).
!##
!##  OTHERWISE, THE PROGRAM DOES WHAT IT IS SUPPOSED TO DO! IT WORKS FINE 
!##  (as far as I applied it).
!##
!##  IT MIGHT BE SENSIBLE, IN THE LONG RUN, TO PUT THIS PROGRAM INTO THE ASKI MAIN PACKAGE
!##
!##  Florian Schumacher, Ruhr-Universitaet Bochum, November 2015
!##
!##############################################################################################
!##############################################################################################
!
!
program createSpecfem3dFilters
!  use inversionBasics
!  use iterationStepBasics
  use discreteFourierTransform
  use asciiDataIO
!  use argumentParser
!  use fileUnitHandler
  use errorMessage

  implicit None

  type (error_message) :: errmsg
  
  type (discrete_fourier_transform) :: DFT

  integer, dimension(:), pointer :: ifreq
  real :: df,dt

  real, dimension(:,:), pointer :: ricker_stf,errf_stf
  complex, dimension(:), allocatable :: ricker_spec,errf_spec

  integer :: it,nt1,nt2,nf


! HARDCODED VALUES FOR ASKI_inversion_cross_borehole

  df = 2.0


  nf = 21
  allocate(ifreq(nf),ricker_spec(nf),errf_spec(nf))
  ifreq = (/ 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38 /)

  dt = 2.114e-4
  nt1 = 0; nt2 = 1184-1


  errmsg = readAsciiData("ricker.txt",11,ricker_stf,2)
  if(.level.errmsg/=0) call print(errmsg)
  if(.level.errmsg==2) stop
  call dealloc(errmsg)
  write(*,*) "size of ricker: ",size(ricker_stf,1),size(ricker_stf,2)


  ! errmsg = readAsciiData("errf.txt",11,errf_stf,2)
  ! if(.level.errmsg/=0) call print(errmsg)
  ! if(.level.errmsg==2) stop
  ! call dealloc(errmsg)
  ! write(*,*) "size of errf: ",size(errf_stf,1),size(errf_stf,2)

  ! if(size(errf_stf,1) /= size(ricker_stf,1)) then
  !    write(*,*) "ERROR, sizes of wavelets do not agree, not supported for now"
  !    stop
  ! end if

  ! ! derive errf by simple finite differences
  ! do it = 1,size(errf_stf,1)-1
  !    errf_stf(it,2) = (errf_stf(it+1,2) - errf_stf(it,2)) / dt
  ! end do ! it
  ! errf_stf(size(errf_stf,1),2) = errf_stf(size(errf_stf,1)-1,2)
  

  call new(errmsg,'createSpecfem3dFilters')
  call initiateForwardDFT(DFT,dt,nt1,nt2,ifreq*df,errmsg,hanning_taper=0.05)
  if(.level.errmsg/=0) call print(errmsg)
  if(.level.errmsg==2) stop

  call transformForwardDFT(DFT,ricker_stf(:,2),ricker_spec,errmsg)
  if(.level.errmsg/=0) call print(errmsg)
  if(.level.errmsg==2) stop

  ! call transformForwardDFT(DFT,errf_stf(:,2),errf_spec,errmsg)
  ! if(.level.errmsg/=0) call print(errmsg)
  ! if(.level.errmsg==2) stop

  call dealloc(errmsg)

  errmsg = writeAsciiData('ricker_spec.txt',11,ricker_spec)
  if(.level.errmsg/=0) call print(errmsg)
  if(.level.errmsg==2) stop
  call dealloc(errmsg)

  ! errmsg = writeAsciiData('errf_spec.txt',11,errf_spec)
  ! if(.level.errmsg/=0) call print(errmsg)
  ! if(.level.errmsg==2) stop
  ! call dealloc(errmsg)

  ! errmsg = writeAsciiData('ricker_divided_by_derivative_of_errf_spec.txt',11,ricker_spec/errf_spec)
  ! if(.level.errmsg/=0) call print(errmsg)
  ! if(.level.errmsg==2) stop
  ! call dealloc(errmsg)

  call dealloc(DFT)
  if(allocated(ricker_spec)) deallocate(ricker_spec)
  ! if(allocated(errf_spec)) deallocate(errf_spec)
  if(associated(ricker_stf)) deallocate(ricker_stf)
  ! if(associated(errf_stf)) deallocate(errf_stf)
  if(associated(ifreq)) deallocate(ifreq)

end program createSpecfem3dFilters
