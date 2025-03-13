
      subroutine stopcm1()
#ifdef MPI
      use mpi
#endif
      implicit none

#ifdef MPI
      integer :: errcode,ierr

      call mpi_abort( mpi_comm_world, errcode , ierr )
#endif

      stop

      end subroutine stopcm1

