!----------------------------------------------------------------------------
!   Copyright 2016 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of ASKI version 1.2.
!
!   ASKI version 1.2 are free software: you can
!   redistribute it and/or modify it under the terms of the GNU
!   General Public License as published by the Free Software
!   Foundation, either version 2 of the License, or (at your option) 
!   any later version.
!
!   ASKI version 1.2 are distributed in the hope that they
!   will be useful, but WITHOUT ANY WARRANTY; without even the implied
!   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with ASKI version 1.2.
!   If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
!> \brief Fast Fourier transform routines
!!
!! \details Fast Fourier transform of the  2**n complex values x + iy
!!  to frequency domain or back to time domain. Different types of
!!  normalization supported.
!!
!! \author Wolfgang Friederich
!
 module fourierTransform
!
   use mathConstants
!
   implicit none
!
   interface fastFourierTransform
      module procedure realFastFourierTransform
      module procedure doubleFastFourierTransform
   end interface fastFourierTransform
!
 contains
!----------------------------------------------------------------------
!> \brief routine for single precision FFT
!! \param x real part of incoming / outgoing values
!! \param y imaginary part of incoming / outgoing values
!! \param n x and y are expected to hold 2**n values each
!! \param is controls the direction of transformation, as well as normalization:  
!!        negativ is -> forward transform (time -> frequency domain)
!!        positive is -> inverse transform (frequency -> time domain)
!!        abs(is) == 1 -> without normalization
!!        abs(is) == 2 -> normalization by 1/2**(n+1)
!!        abs(is) == 3 -> normalization by 1/sqrt(2**(n+1))
!!        abs(is) == 3 -> no normalization, no control expression (?!)
!
   subroutine realFastFourierTransform(x,y,n,is)
     real, dimension(:) :: x,y
     integer :: is,n,m
     integer, dimension(21) :: zh
     integer :: l,ng,nar,lar,larh,jr,ja,nr,jb,j,js,k,ny,nny
     real :: gn,alpha,piz,beta,excos,exsin,zx,zy
!
     piz=2.*mc_pi
!
!  table of powers of two (zh = "zwei hoch" (German for "two to the power of"))
!
     zh(1)=1
     do l=1,n
        zh(l+1)=2*zh(l)
     enddo
     ng=zh(n+1)
     gn=1./float(ng)
!
!  core routine, triple loop over step / index / partial time/spectral-series
!
     do m=1,n
        nar=zh(m)
        lar=ng/nar
        larh=lar/2
        alpha =  piz/float(isign(lar,is))
        do jr=1,larh
           beta=alpha*float(jr-1)
           excos = cos(beta)
           exsin = sin(beta)
           ja=jr-lar
           do nr=1,nar
              ja=ja+lar
              jb=ja+larh
              zx = x(ja)-x(jb)
              zy = y(ja)-y(jb)
              x(ja) = x(ja)+x(jb)
              y(ja) = y(ja)+y(jb)
              x(jb) = zx*excos-zy*exsin
              y(jb) = zx*exsin+zy*excos
           enddo
        enddo
     enddo
!
!     normalization
!
     if (iabs(is) == 3) gn=sqrt(gn)
     if (iabs(is) == 2 .or. iabs(is) == 3) then
        y = y*gn
        x = x*gn
     endif
!
!     re-ordering by "bitreversed" indizes
!
     do j=1,ng
        js=j-1
        k=1
        nny=n+1
        do ny=1,n
           nny=nny-1
           if (js.lt.zh(nny)) cycle
           js=js-zh(nny)
           k=k+zh(ny)
        enddo
        if (j-k < 0) then
           zx = x(j)
           zy = y(j)
           x(j) = x(k)
           y(j) = y(k)
           x(k) = zx
           y(k) = zy
        else
           cycle
        endif
     enddo
   end subroutine realFastFourierTransform
!----------------------------------------------------------------------
!> \brief routine for double precision FFT
!! \param x real part of incoming / outgoing values
!! \param y imaginary part of incoming / outgoing values
!! \param n x and y are expected to hold 2**n values each
!! \param is controls the direction of transformation, as well as normalization:  
!!        negativ is -> forward transform (time -> frequency domain)
!!        positive is -> inverse transform (frequency -> time domain)
!!        abs(is) == 1 -> without normalization
!!        abs(is) == 2 -> normalization by 1/2**(n+1)
!!        abs(is) == 3 -> normalization by 1/sqrt(2**(n+1))
!!        abs(is) == 3 -> no normalization, no control expression (?!)
!
   subroutine doubleFastFourierTransform(x,y,n,is)
     double precision, dimension(:) :: x,y
     integer :: is,n,m
     integer, dimension(21) :: zh
     integer :: l,ng,nar,lar,larh,jr,ja,nr,jb,j,js,k,ny,nny
     double precision :: gn,alpha,piz,beta,excos,exsin,zx,zy
!
     piz=2.*mc_pid
!
!  table of powers of two (zh = "zwei hoch" (German for "two to the power of"))
!
     zh(1)=1
     do l=1,n
        zh(l+1)=2*zh(l)
     enddo
     ng=zh(n+1)
     gn=1./dble(ng)
!
!  core routine, triple loop over step / index / partial time/spectral-series
!
     do m=1,n
        nar=zh(m)
        lar=ng/nar
        larh=lar/2
        alpha =  piz/dble(isign(lar,is))
        do jr=1,larh
           beta=alpha*dble(jr-1)
           excos = dcos(beta)
           exsin = dsin(beta)
           ja=jr-lar
           do nr=1,nar
              ja=ja+lar
              jb=ja+larh
              zx = x(ja)-x(jb)
              zy = y(ja)-y(jb)
              x(ja) = x(ja)+x(jb)
              y(ja) = y(ja)+y(jb)
              x(jb) = zx*excos-zy*exsin
              y(jb) = zx*exsin+zy*excos
           enddo
        enddo
     enddo
!
!     normalization
!
     if (iabs(is) == 3) gn=dsqrt(gn)
     if (iabs(is) == 2 .or. iabs(is) == 3) then
        y = y*gn
        x = x*gn
     endif
!
!     re-ordering by "bitreversed" indizes
!
     do j=1,ng
        js=j-1
        k=1
        nny=n+1
        do ny=1,n
           nny=nny-1
           if (js.lt.zh(nny)) cycle
           js=js-zh(nny)
           k=k+zh(ny)
        enddo
        if (j-k < 0) then
           zx = x(j)
           zy = y(j)
           x(j) = x(k)
           y(j) = y(k)
           x(k) = zx
           y(k) = zy
        else
           cycle
        endif
     enddo
   end subroutine doubleFastFourierTransform
!
 end module fourierTransform
