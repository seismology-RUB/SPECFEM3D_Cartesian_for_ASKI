!
!-------------------------------------------------------------------------------------------------
!

  subroutine send_c_t(sendbuf, sendcount, dest)

  use my_mpi

  implicit none

  integer :: dest,sendcount,ier
  integer :: tag = 100
  complex, dimension(sendcount):: sendbuf

  call MPI_SEND(sendbuf,sendcount,MPI_COMPLEX,dest,tag, &
       my_local_mpi_comm_world,ier)

  end subroutine send_c_t

!
!-------------------------------------------------------------------------------------------------
!

  subroutine recv_c_t(recvbuf, recvcount, source)

  use my_mpi

  implicit none

  integer :: source,recvcount,ier
  integer :: tag = 100
  complex, dimension(recvcount):: recvbuf

  call MPI_RECV(recvbuf,recvcount,MPI_COMPLEX,source,tag, &
       my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)

  end subroutine recv_c_t

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_l(buffer, countval)

  use my_mpi

  implicit none

  integer countval
  logical, dimension(countval) :: buffer

  integer ier

  call MPI_BCAST(buffer,countval,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_l
