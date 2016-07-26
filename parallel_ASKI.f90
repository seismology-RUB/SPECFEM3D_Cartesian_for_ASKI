!
!----
!

  subroutine send_c_t(sendbuf, sendcount, dest)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  integer :: dest,sendcount,ier
  integer :: tag = 100
  complex, dimension(sendcount):: sendbuf

  call MPI_SEND(sendbuf,sendcount,MPI_COMPLEX,dest,tag, &
       MPI_COMM_WORLD,ier)

  end subroutine send_c_t

!
!----
!

  subroutine recv_c_t(recvbuf, recvcount, source)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  integer :: source,recvcount,ier
  integer :: tag = 100
  complex, dimension(recvcount):: recvbuf
  integer msg_status(MPI_STATUS_SIZE)

  call MPI_RECV(recvbuf,recvcount,MPI_COMPLEX,source,tag, &
       MPI_COMM_WORLD,msg_status,ier)

  end subroutine recv_c_t

