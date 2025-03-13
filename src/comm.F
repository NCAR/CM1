  MODULE comm_module

  implicit none

  CONTAINS

!-----------------------------------------------------------------------
!  message passing routines
!-----------------------------------------------------------------------


      integer function nabor(i,j,nx,ny)
      implicit none
      integer i,j,nx,ny
      integer newi,newj

      newi=i
      newj=j

      if ( newi .lt.  1 ) newi = nx
      if ( newi .gt.  nx) newi = 1

      if ( newj .lt.  1 ) newj = ny
      if ( newj .gt.  ny) newj = 1

      nabor = (newi-1) + (newj-1)*nx

      end function nabor

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

#ifdef MPI
      subroutine getcorner(s,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      use input
      use mpi
      implicit none
 
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: s
      real, intent(inout), dimension(nk) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
 
      integer k,nn,nr,nrb,index
      integer :: index_nw,index_sw,index_ne,index_se
      integer reqs(8)
      integer tag1,tag2,tag3,tag4
      integer, dimension(mpi_status_size,8) :: status1

!-----

      nr = 0
      index_nw = -1
      index_sw = -1
      index_ne = -1
      index_se = -1

!------------------------------------------------------------------

      tag1=5001

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        nr = nr + 1
        call mpi_irecv(nw2,nk,MPI_REAL,mynw,tag1,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
        index_nw = nr
      endif

!-----

      tag2=5002

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        nr = nr + 1
        call mpi_irecv(sw2,nk,MPI_REAL,mysw,tag2,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
        index_sw = nr
      endif

!-----

      tag3=5003

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        nr = nr + 1
        call mpi_irecv(ne2,nk,MPI_REAL,myne,tag3,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
        index_ne = nr
      endif

!-----

      tag4=5004

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        nr = nr + 1
        call mpi_irecv(se2,nk,MPI_REAL,myse,tag4,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
        index_se = nr
      endif

!------------------------------------------------------------------

      nrb = 4

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          se1(k)=s(ni,1,k)
        enddo
        nrb = nrb + 1
        call mpi_isend(se1,nk,MPI_REAL,myse,tag1,MPI_COMM_WORLD,   &
                       reqs(nrb),ierr)
      endif
 
!-----

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          ne1(k)=s(ni,nj,k)
        enddo
        nrb = nrb + 1
        call mpi_isend(ne1,nk,MPI_REAL,myne,tag2,MPI_COMM_WORLD,   &
                       reqs(nrb),ierr)
      endif
 
!-----

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          sw1(k)=s(1,1,k)
        enddo
        nrb = nrb + 1
        call mpi_isend(sw1,nk,MPI_REAL,mysw,tag3,MPI_COMM_WORLD,   &
                       reqs(nrb),ierr)
      endif
 
!-----

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          nw1(k)=s(1,nj,k)
        enddo
        nrb = nrb + 1
        call mpi_isend(nw1,nk,MPI_REAL,mynw,tag4,MPI_COMM_WORLD,   &
                       reqs(nrb),ierr)
      endif
 
!-----

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(1:nr),index,MPI_STATUS_IGNORE,ierr)
        nn = nn + 1

      if(index.eq.index_nw)then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          s(0,nj+1,k)=nw2(k)
        enddo
      elseif(index.eq.index_sw)then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          s(0,0,k)=sw2(k)
        enddo
      elseif(index.eq.index_ne)then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          s(ni+1,nj+1,k)=ne2(k)
        enddo
      elseif(index.eq.index_se)then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          s(ni+1,0,k)=se2(k)
        enddo
      endif

      enddo

!-----

    nrb = nrb-4

    if( nrb.ge.1 )then
      call MPI_WAITALL(nrb,reqs(5:5+nrb-1),status1(1:mpi_status_size,5:5+nrb-1),ierr)
    endif

!-----

      if(timestats.ge.1) time_mptk1=time_mptk1+mytime()

      end subroutine getcorner


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getcornert(t,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      use input
      use mpi
      implicit none

      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: t
      real, intent(inout), dimension(nk+1) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2

      integer k,nn,nr,nrb,index
      integer :: index_nw,index_sw,index_ne,index_se
      integer reqs(8)
      integer tag1,tag2,tag3,tag4
      integer, dimension(mpi_status_size,8) :: status1

!-----

      nr = 0
      index_nw = -1
      index_sw = -1
      index_ne = -1
      index_se = -1

!-------------------------------------------------------------

      tag1=5001

      if(ibw.eq.0 .and. ibn.eq.0)then
        nr = nr + 1
        call mpi_irecv(nw2,nkp1,MPI_REAL,mynw,tag1,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
        index_nw = nr
      endif

!-----

      tag2=5002

      if(ibw.eq.0 .and. ibs.eq.0)then
        nr = nr + 1
        call mpi_irecv(sw2,nkp1,MPI_REAL,mysw,tag2,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
        index_sw = nr
      endif

!-----

      tag3=5003

      if(ibe.eq.0 .and. ibn.eq.0)then
        nr = nr + 1
        call mpi_irecv(ne2,nkp1,MPI_REAL,myne,tag3,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
        index_ne = nr
      endif

!-----

      tag4=5004

      if(ibe.eq.0 .and. ibs.eq.0)then
        nr = nr + 1
        call mpi_irecv(se2,nkp1,MPI_REAL,myse,tag4,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
        index_se = nr
      endif

!-------------------------------------------------------------

      nrb = 4

      if(ibe.eq.0 .and. ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          se1(k)=t(ni,1,k)
        enddo
        nrb = nrb + 1
        call mpi_isend(se1,nkp1,MPI_REAL,myse,tag1,MPI_COMM_WORLD,   &
                       reqs(nrb),ierr)
      endif

!-----

      if(ibe.eq.0 .and. ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          ne1(k)=t(ni,nj,k)
        enddo
        nrb = nrb + 1
        call mpi_isend(ne1,nkp1,MPI_REAL,myne,tag2,MPI_COMM_WORLD,   &
                       reqs(nrb),ierr)
      endif

!-----

      if(ibw.eq.0 .and. ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          sw1(k)=t(1,1,k)
        enddo
        nrb = nrb + 1
        call mpi_isend(sw1,nkp1,MPI_REAL,mysw,tag3,MPI_COMM_WORLD,   &
                       reqs(nrb),ierr)
      endif

!-----

      if(ibw.eq.0 .and. ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          nw1(k)=t(1,nj,k)
        enddo
        nrb = nrb + 1
        call mpi_isend(nw1,nkp1,MPI_REAL,mynw,tag4,MPI_COMM_WORLD,   &
                       reqs(nrb),ierr)
      endif

!-----

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(1:nr),index,MPI_STATUS_IGNORE,ierr)
        nn = nn + 1

      if(index.eq.index_nw)then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          t(0,nj+1,k)=nw2(k)
        enddo
      elseif(index.eq.index_sw)then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          t(0,0,k)=sw2(k)
        enddo
      elseif(index.eq.index_ne)then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          t(ni+1,nj+1,k)=ne2(k)
        enddo
      elseif(index.eq.index_se)then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          t(ni+1,0,k)=se2(k)
        enddo
      endif

      enddo

!-----

    nrb = nrb-4

    if( nrb.ge.1 )then
      call MPI_WAITALL(nrb,reqs(5:5+nrb-1),status1(1:mpi_status_size,5:5+nrb-1),ierr)
    endif

!-----

      if(timestats.ge.1) time_mptk1=time_mptk1+mytime()

      end subroutine getcornert


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getcorneru(u,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      use input
      use mpi
      implicit none

      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: u
      real, intent(inout), dimension(nk) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2

      integer k,nn,nr,nrb,index
      integer :: index_nw,index_sw,index_ne,index_se
      integer reqs(8)
      integer tag,count
      integer, dimension(mpi_status_size,8) :: status1

!-----

      nr = 0
      index_nw = -1
      index_sw = -1
      index_ne = -1
      index_se = -1

!-----------------------------------------------------------------------

      tag=5001
      count=nk

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        nr = nr + 1
        call mpi_irecv(nw2,count,MPI_REAL,mynw,tag,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
        index_nw = nr
      endif

!-----

      tag=5002
      count=nk

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        nr = nr + 1
        call mpi_irecv(sw2,count,MPI_REAL,mysw,tag,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
        index_sw = nr
      endif

!-----

      tag=5003
      count=nk

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        nr = nr + 1
        call mpi_irecv(ne2,count,MPI_REAL,myne,tag,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
        index_ne = nr
      endif

!-----

      tag=5004
      count=nk

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        nr = nr + 1
        call mpi_irecv(se2,count,MPI_REAL,myse,tag,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
        index_se = nr
      endif

!-----------------------------------------------------------------------

      nrb = 4

      tag=5001
      count=nk

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          se1(k)=u(ni,1,k)
        enddo
        nrb = nrb + 1
        call mpi_isend(se1,count,MPI_REAL,myse,tag,MPI_COMM_WORLD,   &
                       reqs(nrb),ierr)
      endif

!-----

      tag=5002
      count=nk

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          ne1(k)=u(ni,nj,k)
        enddo
        nrb = nrb + 1
        call mpi_isend(ne1,count,MPI_REAL,myne,tag,MPI_COMM_WORLD,   &
                       reqs(nrb),ierr)
      endif

!-----

      tag=5003
      count=nk

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          sw1(k)=u(2,1,k)
        enddo
        nrb = nrb + 1
        call mpi_isend(sw1,count,MPI_REAL,mysw,tag,MPI_COMM_WORLD,   &
                       reqs(nrb),ierr)
      endif

!-----

      tag=5004
      count=nk

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          nw1(k)=u(2,nj,k)
        enddo
        nrb = nrb + 1
        call mpi_isend(nw1,count,MPI_REAL,mynw,tag,MPI_COMM_WORLD,   &
                       reqs(nrb),ierr)
      endif

!-----------------------------------------------------------------------

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(1:nr),index,MPI_STATUS_IGNORE,ierr)
        nn = nn + 1

      if(index.eq.index_nw)then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          u(0,nj+1,k)=nw2(k)
        enddo
      elseif(index.eq.index_sw)then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          u(0,0,k)=sw2(k)
        enddo
      elseif(index.eq.index_ne)then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          u(ni+2,nj+1,k)=ne2(k)
        enddo
      elseif(index.eq.index_se)then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          u(ni+2,0,k)=se2(k)
        enddo
      endif

      enddo

!-----

    nrb = nrb-4

    if( nrb.ge.1 )then
      call MPI_WAITALL(nrb,reqs(5:5+nrb-1),status1(1:mpi_status_size,5:5+nrb-1),ierr)
    endif

!-----

      if(timestats.ge.1) time_mptk1=time_mptk1+mytime()

      end subroutine getcorneru


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getcornerv(v,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      use input
      use mpi
      implicit none

      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: v
      real, intent(inout), dimension(nk) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2

      integer k,nn,nr,nrb,index
      integer :: index_nw,index_sw,index_ne,index_se
      integer reqs(8)
      integer tag,count
      integer, dimension(mpi_status_size,8) :: status1

!-----

      nr = 0
      index_nw = -1
      index_sw = -1
      index_ne = -1
      index_se = -1

!-----------------------------------------------------------------------

      tag=5011
      count=nk

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        nr = nr + 1
        call mpi_irecv(nw2,count,MPI_REAL,mynw,tag,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
        index_nw = nr
      endif

!-----

      tag=5012
      count=nk

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        nr = nr + 1
        call mpi_irecv(sw2,count,MPI_REAL,mysw,tag,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
        index_sw = nr
      endif

!-----

      tag=5013
      count=nk

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        nr = nr + 1
        call mpi_irecv(ne2,count,MPI_REAL,myne,tag,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
        index_ne = nr
      endif

!-----

      tag=5014
      count=nk

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        nr = nr + 1
        call mpi_irecv(se2,count,MPI_REAL,myse,tag,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
        index_se = nr
      endif

!-----------------------------------------------------------------------

      nrb = 4

      tag=5011
      count=nk

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          se1(k)=v(ni,2,k)
        enddo
        nrb = nrb + 1
        call mpi_isend(se1,count,MPI_REAL,myse,tag,MPI_COMM_WORLD,   &
                       reqs(nrb),ierr)
      endif

!-----

      tag=5012
      count=nk

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          ne1(k)=v(ni,nj,k)
        enddo
        nrb = nrb + 1
        call mpi_isend(ne1,count,MPI_REAL,myne,tag,MPI_COMM_WORLD,   &
                       reqs(nrb),ierr)
      endif

!-----

      tag=5013
      count=nk

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          sw1(k)=v(1,2,k)
        enddo
        nrb = nrb + 1
        call mpi_isend(sw1,count,MPI_REAL,mysw,tag,MPI_COMM_WORLD,   &
                       reqs(nrb),ierr)
      endif

!-----

      tag=5014
      count=nk

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          nw1(k)=v(1,nj,k)
        enddo
        nrb = nrb + 1
        call mpi_isend(nw1,count,MPI_REAL,mynw,tag,MPI_COMM_WORLD,   &
                       reqs(nrb),ierr)
      endif

!-----------------------------------------------------------------------

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(1:nr),index,MPI_STATUS_IGNORE,ierr)
        nn = nn + 1

      if(index.eq.index_nw)then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          v(0,nj+2,k)=nw2(k)
        enddo
      elseif(index.eq.index_sw)then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          v(0,0,k)=sw2(k)
        enddo
      elseif(index.eq.index_ne)then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          v(ni+1,nj+2,k)=ne2(k)
        enddo
      elseif(index.eq.index_se)then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          v(ni+1,0,k)=se2(k)
        enddo
      endif

      enddo

!-----

    nrb = nrb-4

    if( nrb.ge.1 )then
      call MPI_WAITALL(nrb,reqs(5:5+nrb-1),status1(1:mpi_status_size,5:5+nrb-1),ierr)
    endif

!-----

      if(timestats.ge.1) time_mptk1=time_mptk1+mytime()

      end subroutine getcornerv


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getcornerw(w,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      use input
      use mpi
      implicit none

      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: w
      real, intent(inout), dimension(nk) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2

      integer k,nn,nr,nrb,index
      integer :: index_nw,index_sw,index_ne,index_se
      integer reqs(8)
      integer tag,count
      integer, dimension(mpi_status_size,8) :: status1

!-----

      nr = 0
      index_nw = -1
      index_sw = -1
      index_ne = -1
      index_se = -1

!-----------------------------------------------------------------------

      tag=5021
      count=nk

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        nr = nr + 1
        call mpi_irecv(nw2,count,MPI_REAL,mynw,tag,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
        index_nw = nr
      endif

!-----

      tag=5022
      count=nk

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        nr = nr + 1
        call mpi_irecv(sw2,count,MPI_REAL,mysw,tag,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
        index_sw = nr
      endif

!-----

      tag=5023
      count=nk

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        nr = nr + 1
        call mpi_irecv(ne2,count,MPI_REAL,myne,tag,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
        index_ne = nr
      endif

!-----

      tag=5024
      count=nk

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        nr = nr + 1
        call mpi_irecv(se2,count,MPI_REAL,myse,tag,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
        index_se = nr
      endif

!-----------------------------------------------------------------------

      nrb = 4

      tag=5021
      count=nk

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          se1(k)=w(ni,1,k)
        enddo
        nrb = nrb + 1
        call mpi_isend(se1,count,MPI_REAL,myse,tag,MPI_COMM_WORLD,   &
                       reqs(nrb),ierr)
      endif

!-----

      tag=5022
      count=nk

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          ne1(k)=w(ni,nj,k)
        enddo
        nrb = nrb + 1
        call mpi_isend(ne1,count,MPI_REAL,myne,tag,MPI_COMM_WORLD,   &
                       reqs(nrb),ierr)
      endif

!-----

      tag=5023
      count=nk

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          sw1(k)=w(1,1,k)
        enddo
        nrb = nrb + 1
        call mpi_isend(sw1,count,MPI_REAL,mysw,tag,MPI_COMM_WORLD,   &
                       reqs(nrb),ierr)
      endif

!-----

      tag=5024
      count=nk

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          nw1(k)=w(1,nj,k)
        enddo
        nrb = nrb + 1
        call mpi_isend(nw1,count,MPI_REAL,mynw,tag,MPI_COMM_WORLD,   &
                       reqs(nrb),ierr)
      endif

!-----------------------------------------------------------------------

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(1:nr),index,MPI_STATUS_IGNORE,ierr)
        nn = nn + 1

      if(index.eq.index_nw)then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          w(0,nj+1,k)=nw2(k)
        enddo
      elseif(index.eq.index_sw)then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          w(0,0,k)=sw2(k)
        enddo
      elseif(index.eq.index_ne)then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          w(ni+1,nj+1,k)=ne2(k)
        enddo
      elseif(index.eq.index_se)then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          w(ni+1,0,k)=se2(k)
        enddo
      endif

      enddo

!-----

    nrb = nrb-4

    if( nrb.ge.1 )then
      call MPI_WAITALL(nrb,reqs(5:5+nrb-1),status1(1:mpi_status_size,5:5+nrb-1),ierr)
    endif

!-----

      if(timestats.ge.1) time_mptk1=time_mptk1+mytime()

      end subroutine getcornerw


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine getcorner3(s,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      use input
      use mpi
      implicit none
 
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: s
      real, intent(inout), dimension(cmp,cmp,nk) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
 
      integer :: i,j,k,nn,nr,index
      integer :: index_nw,index_sw,index_ne,index_se
      integer reqs(8)
      integer tag1,tag2,tag3,tag4

!-----

      nr = 0
      index_nw = -1
      index_sw = -1
      index_ne = -1
      index_se = -1

!------------------------------------------------------------------

      tag1=5001

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        nr = nr + 1
        call mpi_irecv(nw2,cmp*cmp*nk,MPI_REAL,mynw,tag1,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
        index_nw = nr
      endif

!-----

      tag2=5002

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        nr = nr + 1
        call mpi_irecv(sw2,cmp*cmp*nk,MPI_REAL,mysw,tag2,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
        index_sw = nr
      endif

!-----

      tag3=5003

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        nr = nr + 1
        call mpi_irecv(ne2,cmp*cmp*nk,MPI_REAL,myne,tag3,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
        index_ne = nr
      endif

!-----

      tag4=5004

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        nr = nr + 1
        call mpi_irecv(se2,cmp*cmp*nk,MPI_REAL,myse,tag4,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
        index_se = nr
      endif

!------------------------------------------------------------------

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        do k=1,nk
        do j=1,cmp
        do i=1,cmp
!!!          se1(i,j,k)=s(ni,1,k)
          se1(i,j,k)=s(ni-cmp+i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(se1,cmp*cmp*nk,MPI_REAL,myse,tag1,MPI_COMM_WORLD,   &
                       reqs(5),ierr)
      endif
 
!-----

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        do k=1,nk
        do j=1,cmp
        do i=1,cmp
!!!          ne1(i,j,k)=s(ni,nj,k)
          ne1(i,j,k)=s(ni-cmp+i,nj-cmp+j,k)
        enddo
        enddo
        enddo
        call mpi_isend(ne1,cmp*cmp*nk,MPI_REAL,myne,tag2,MPI_COMM_WORLD,   &
                       reqs(6),ierr)
      endif
 
!-----

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        do k=1,nk
        do j=1,cmp
        do i=1,cmp
!!!          sw1(i,j,k)=s(1,1,k)
          sw1(i,j,k)=s(i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(sw1,cmp*cmp*nk,MPI_REAL,mysw,tag3,MPI_COMM_WORLD,   &
                       reqs(7),ierr)
      endif
 
!-----

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        do k=1,nk
        do j=1,cmp
        do i=1,cmp
!!!          nw1(i,j,k)=s(1,nj,k)
          nw1(i,j,k)=s(i,nj-cmp+j,k)
        enddo
        enddo
        enddo
        call mpi_isend(nw1,cmp*cmp*nk,MPI_REAL,mynw,tag4,MPI_COMM_WORLD,   &
                       reqs(8),ierr)
      endif
 
!-----

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(1:nr),index,MPI_STATUS_IGNORE,ierr)
        nn = nn + 1

      if(index.eq.index_nw)then
        do k=1,nk
        do j=1,cmp
        do i=1,cmp
!!!          s(0,nj+1,k)=nw2(i,j,k)
          s(-cmp+i,nj+j,k)=nw2(i,j,k)
        enddo
        enddo
        enddo
      elseif(index.eq.index_sw)then
        do k=1,nk
        do j=1,cmp
        do i=1,cmp
!!!          s(0,0,k)=sw2(i,j,k)
          s(-cmp+i,-cmp+j,k)=sw2(i,j,k)
        enddo
        enddo
        enddo
      elseif(index.eq.index_ne)then
        do k=1,nk
        do j=1,cmp
        do i=1,cmp
!!!          s(ni+1,nj+1,k)=ne2(i,j,k)
          s(ni+i,nj+j,k)=ne2(i,j,k)
        enddo
        enddo
        enddo
      elseif(index.eq.index_se)then
        do k=1,nk
        do j=1,cmp
        do i=1,cmp
!!!          s(ni+1,0,k)=se2(i,j,k)
          s(ni+i,-cmp+j,k)=se2(i,j,k)
        enddo
        enddo
        enddo
      endif

      enddo

!-----

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(5),MPI_STATUS_IGNORE,ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(6),MPI_STATUS_IGNORE,ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(7),MPI_STATUS_IGNORE,ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(8),MPI_STATUS_IGNORE,ierr)
      endif

!-----

      if(timestats.ge.1) time_mptk1=time_mptk1+mytime()

      end subroutine getcorner3


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_2we_start(s,west,newwest,east,neweast,reqs)
      use input
      use mpi
      implicit none

      real s(ib:ie,jb:je,kb:ke)
      real west(2,nj,nk),newwest(2,nj,nk)
      real east(2,nj,nk),neweast(2,nj,nk)
      integer reqs(4)

      integer i,j,k,nr
      integer tag1,tag2

!------------------------------------------------

      nr = 0

      nf=nf+1
      tag1=nf

      ! receive east
      if(ibe.eq.0)then
        nr = nr + 1
        call mpi_irecv(neweast,cs2we,MPI_REAL,myeast,tag1,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      nf=nf+1
      tag2=nf

      ! receive west
      if(ibw.eq.0)then
        nr = nr + 1
        call mpi_irecv(newwest,cs2we,MPI_REAL,mywest,tag2,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!------------------------------------------------

      nr = 2

      ! send west
      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,2
          west(i,j,k)=s(i,j,k)
        enddo
        enddo
        enddo
        nr = nr+1
        call mpi_isend(west,cs2we,MPI_REAL,mywest,tag1,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      ! send east
      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,2
          east(i,j,k)=s(ni-2+i,j,k)
        enddo
        enddo
        enddo
        nr = nr+1
        call mpi_isend(east,cs2we,MPI_REAL,myeast,tag2,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      if(timestats.ge.1) time_mps3=time_mps3+mytime()

      end subroutine comm_2we_start


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_2we_end(s,west,newwest,east,neweast,reqs)
      use input
      use mpi
      implicit none

      real s(ib:ie,jb:je,kb:ke)
      real west(2,nj,nk),newwest(2,nj,nk)
      real east(2,nj,nk),neweast(2,nj,nk)
      integer reqs(4)

      integer i,j,k,nn,nr,index
      integer :: index_east,index_west
      integer, dimension(mpi_status_size,8) :: status1

!-------------------------------------------------------------------

      index_east = -1
      index_west = -1

      nr = 0
      if(ibe.eq.0)then
        nr = nr + 1
        index_east = nr
      endif
      if(ibw.eq.0)then
        nr = nr + 1
        index_west = nr
      endif

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(1:nr),index,MPI_STATUS_IGNORE,ierr)
        nn = nn + 1

      if(index.eq.index_east)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,2
          s(ni+i,j,k)=neweast(i,j,k)
        enddo
        enddo
        enddo
      elseif(index.eq.index_west)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,2
          s(i-2,j,k)=newwest(i,j,k)
        enddo
        enddo
        enddo
      endif

      enddo

!----------

      nr = 0

      if(ibw.eq.0)then
        nr = nr+1
      endif
      if(ibe.eq.0)then
        nr = nr+1
      endif

    if( nr.ge.1 )then
      call MPI_WAITALL(nr,reqs(3:3+nr-1),status1(1:mpi_status_size,3:3+nr-1),ierr)
    endif

      if(timestats.ge.1) time_mps4=time_mps4+mytime()

      end subroutine comm_2we_end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_2sn_start(s,south,newsouth,north,newnorth,reqs)
      use input
      use mpi
      implicit none

      real s(ib:ie,jb:je,kb:ke)
      real south(ni,2,nk),newsouth(ni,2,nk)
      real north(ni,2,nk),newnorth(ni,2,nk)
      integer reqs(4)

      integer i,j,k,nr
      integer tag3,tag4

!----------

      nr = 0

      nf=nf+1
      tag3=nf

      ! receive south
      if(ibs.eq.0)then
        nr = nr + 1
        call mpi_irecv(newsouth,cs2sn,MPI_REAL,mysouth,tag3,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      nf=nf+1
      tag4=nf

      ! receive north
      if(ibn.eq.0)then
        nr = nr + 1
        call mpi_irecv(newnorth,cs2sn,MPI_REAL,mynorth,tag4,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      nr = 2

      ! send north
      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,2
        do i=1,ni
          north(i,j,k)=s(i,nj-2+j,k)
        enddo
        enddo
        enddo
        nr = nr+1
        call mpi_isend(north,cs2sn,MPI_REAL,mynorth,tag3,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      ! send south
      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,2
        do i=1,ni
          south(i,j,k)=s(i,j,k)
        enddo
        enddo
        enddo
        nr = nr+1
        call mpi_isend(south,cs2sn,MPI_REAL,mysouth,tag4,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      if(timestats.ge.1) time_mps3=time_mps3+mytime()

      end subroutine comm_2sn_start


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_2sn_end(s,south,newsouth,north,newnorth,reqs)
      use input
      use mpi
      implicit none

      real s(ib:ie,jb:je,kb:ke)
      real south(ni,2,nk),newsouth(ni,2,nk)
      real north(ni,2,nk),newnorth(ni,2,nk)
      integer reqs(4)

      integer i,j,k,nn,nr,index
      integer :: index_south,index_north
      integer, dimension(mpi_status_size,8) :: status1

!----------

      index_south = -1
      index_north = -1

      nr = 0
      if(ibs.eq.0)then
        nr = nr + 1
        index_south = nr
      endif
      if(ibn.eq.0)then
        nr = nr + 1
        index_north = nr
      endif

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(1:nr),index,MPI_STATUS_IGNORE,ierr)
        nn = nn + 1

      if(index.eq.index_south)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,2
        do i=1,ni
          s(i,j-2,k)=newsouth(i,j,k)
        enddo
        enddo
        enddo
      elseif(index.eq.index_north)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,2
        do i=1,ni
          s(i,nj+j,k)=newnorth(i,j,k)
        enddo
        enddo
        enddo
      endif

      enddo

!----------

      nr = 0

      if(ibn.eq.0)then
        nr = nr+1
      endif
      if(ibs.eq.0)then
        nr = nr+1
      endif

    if( nr.ge.1 )then
      call MPI_WAITALL(nr,reqs(3:3+nr-1),status1(1:mpi_status_size,3:3+nr-1),ierr)
    endif

      if(timestats.ge.1) time_mps4=time_mps4+mytime()

      end subroutine comm_2sn_end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_3s_start(s,west,newwest,east,neweast,   &
                                 south,newsouth,north,newnorth,reqs)
      use input
      use mpi
      implicit none
 
      real s(ib:ie,jb:je,kb:ke)
      real west(cmp,nj,nk),newwest(cmp,nj,nk)
      real east(cmp,nj,nk),neweast(cmp,nj,nk)
      real south(ni,cmp,nk),newsouth(ni,cmp,nk)
      real north(ni,cmp,nk),newnorth(ni,cmp,nk)
      integer reqs(8)
 
      integer i,j,k,nr
      integer tag1,tag2,tag3,tag4
 
!------------------------------------------------

      nr = 0

      nf=nf+1
      tag1=nf

      ! receive east
      if(ibe.eq.0)then
        nr = nr + 1
        call mpi_irecv(neweast,cs3we,MPI_REAL,myeast,tag1,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------
 
      nf=nf+1
      tag2=nf

      ! receive west
      if(ibw.eq.0)then
        nr = nr + 1
        call mpi_irecv(newwest,cs3we,MPI_REAL,mywest,tag2,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------
 
      nf=nf+1
      tag3=nf

      ! receive south
      if(ibs.eq.0)then
        nr = nr + 1
        call mpi_irecv(newsouth,cs3sn,MPI_REAL,mysouth,tag3,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------
 
      nf=nf+1
      tag4=nf

      ! receive north
      if(ibn.eq.0)then
        nr = nr + 1
        call mpi_irecv(newnorth,cs3sn,MPI_REAL,mynorth,tag4,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!------------------------------------------------
!------------------------------------------------
!------------------------------------------------

      nr = 4

      ! send west
      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,cmp
          west(i,j,k)=s(i,j,k)
        enddo
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(west,cs3we,MPI_REAL,mywest,tag1,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------

 
      ! send east
      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,cmp
          east(i,j,k)=s(ni-cmp+i,j,k)
        enddo
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(east,cs3we,MPI_REAL,myeast,tag2,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

 
      ! send north
      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,cmp
        do i=1,ni
          north(i,j,k)=s(i,nj-cmp+j,k)
        enddo
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(north,cs3sn,MPI_REAL,mynorth,tag3,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

 
      ! send south
      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,cmp
        do i=1,ni
          south(i,j,k)=s(i,j,k)
        enddo
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(south,cs3sn,MPI_REAL,mysouth,tag4,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------
 
      if(timestats.ge.1) time_mps1=time_mps1+mytime()

      end subroutine comm_3s_start


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_3s_end(s,west,newwest,east,neweast,   &
                               south,newsouth,north,newnorth,reqs)
      use input
      use mpi
      implicit none
 
      real s(ib:ie,jb:je,kb:ke)
      real west(cmp,nj,nk),newwest(cmp,nj,nk)
      real east(cmp,nj,nk),neweast(cmp,nj,nk)
      real south(ni,cmp,nk),newsouth(ni,cmp,nk)
      real north(ni,cmp,nk),newnorth(ni,cmp,nk)
      integer reqs(8)
 
      integer i,j,k,nn,nr,index
      integer :: index_east,index_west,index_south,index_north
      integer, dimension(mpi_status_size,8) :: status1

!-------------------------------------------------------------------

      index_east = -1
      index_west = -1
      index_south = -1
      index_north = -1

      nr = 0
      if(ibe.eq.0)then
        nr = nr + 1
        index_east = nr
      endif
      if(ibw.eq.0)then
        nr = nr + 1
        index_west = nr
      endif
      if(ibs.eq.0)then
        nr = nr + 1
        index_south = nr
      endif
      if(ibn.eq.0)then
        nr = nr + 1
        index_north = nr
      endif

    nn = 1
    do while( nn .le. nr )
      call MPI_WAITANY(nr,reqs(1:nr),index,MPI_STATUS_IGNORE,ierr)
      nn = nn + 1
 
      if(index.eq.index_east)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,cmp
          s(ni+i,j,k)=neweast(i,j,k)
        enddo
        enddo
        enddo
      elseif(index.eq.index_west)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,cmp
          s(i-cmp,j,k)=newwest(i,j,k)
        enddo
        enddo
        enddo
      elseif(index.eq.index_south)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,cmp
        do i=1,ni
          s(i,j-cmp,k)=newsouth(i,j,k)
        enddo
        enddo
        enddo
      elseif(index.eq.index_north)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,cmp
        do i=1,ni
          s(i,nj+j,k)=newnorth(i,j,k)
        enddo
        enddo
        enddo
      endif

    enddo

      if(timestats.ge.1) time_mps2=time_mps2+mytime()
 
!----------
!  patch for corner
 
      if( (ebc.eq.2.or.wbc.eq.2).and.(sbc.eq.1.or.nbc.eq.1) )then
 
        if(ibw.eq.1)then

          if(p2tchsww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(0,0,k)=s(1,0,k)
            enddo
          endif

          if(p2tchnww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(0,nj+1,k)=s(1,nj+1,k)
            enddo
          endif

        endif

        if(ibe.eq.1)then
 
          if(p2tchsee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(ni+1,0,k)=s(ni,0,k)
            enddo
          endif
 
          if(p2tchnee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(ni+1,nj+1,k)=s(ni,nj+1,k)
            enddo
          endif
 
        endif
 
      endif

      if( (ebc.eq.1.or.wbc.eq.1).and.(sbc.eq.2.or.nbc.eq.2) )then
 
        if(ibs.eq.1)then
 
          if(p2tchsws)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(0,0,k)=s(0,1,k)
            enddo
          endif
 
          if(p2tchses)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(ni+1,0,k)=s(ni+1,1,k)
            enddo
          endif
 
        endif
 
        if(ibn.eq.1)then
 
          if(p2tchnwn)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(0,nj+1,k)=s(0,nj,k)
            enddo
          endif
 
          if(p2tchnen)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(ni+1,nj+1,k)=s(ni+1,nj,k)
            enddo
          endif
 
        endif
 
      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

!----------

      nr = 0

      if(ibw.eq.0)then
        nr = nr+1
      endif
      if(ibe.eq.0)then
        nr = nr+1
      endif
      if(ibn.eq.0)then
        nr = nr+1
      endif
      if(ibs.eq.0)then
        nr = nr+1
      endif

    if( nr.ge.1 )then
      call MPI_WAITALL(nr,reqs(5:5+nr-1),status1(1:mpi_status_size,5:5+nr-1),ierr)
    endif

      if(timestats.ge.1) time_mps2=time_mps2+mytime()

!-----------------------------------------------------------
 
      end subroutine comm_3s_end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      ! tk1
      subroutine comm_3t_start(t,west,newwest,east,neweast,   &
                                 south,newsouth,north,newnorth,reqs)
      use input
      use mpi
      implicit none

      real t(ib:ie,jb:je,kb:ke+1)
      real west(cmp,nj,nk+1),newwest(cmp,nj,nk+1)
      real east(cmp,nj,nk+1),neweast(cmp,nj,nk+1)
      real south(ni,cmp,nk+1),newsouth(ni,cmp,nk+1)
      real north(ni,cmp,nk+1),newnorth(ni,cmp,nk+1)
      integer reqs(8)

      integer i,j,k,nr
      integer tag1,tag2,tag3,tag4

!------------------------------------------------

      nr = 0

      nf=nf+1
      tag1=nf

      ! receive east
      if(ibe.eq.0)then
        nr = nr + 1
        call mpi_irecv(neweast,ct3we,MPI_REAL,myeast,tag1,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      nf=nf+1
      tag2=nf

      ! receive west
      if(ibw.eq.0)then
        nr = nr + 1
        call mpi_irecv(newwest,ct3we,MPI_REAL,mywest,tag2,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      nf=nf+1
      tag3=nf

      ! receive south
      if(ibs.eq.0)then
        nr = nr + 1
        call mpi_irecv(newsouth,ct3sn,MPI_REAL,mysouth,tag3,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      nf=nf+1
      tag4=nf

      ! receive north
      if(ibn.eq.0)then
        nr = nr + 1
        call mpi_irecv(newnorth,ct3sn,MPI_REAL,mynorth,tag4,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!------------------------------------------------
!------------------------------------------------
!------------------------------------------------

      nr = 4

      ! send west
      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nkp1
        do j=1,nj
        do i=1,cmp
          west(i,j,k)=t(i,j,k)
        enddo
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(west,ct3we,MPI_REAL,mywest,tag1,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      ! send east
      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nkp1
        do j=1,nj
        do i=1,cmp
          east(i,j,k)=t(ni-cmp+i,j,k)
        enddo
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(east,ct3we,MPI_REAL,myeast,tag2,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      ! send north
      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nkp1
        do j=1,cmp
        do i=1,ni
          north(i,j,k)=t(i,nj-cmp+j,k)
        enddo
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(north,ct3sn,MPI_REAL,mynorth,tag3,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      ! send south
      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nkp1
        do j=1,cmp
        do i=1,ni
          south(i,j,k)=t(i,j,k)
        enddo
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(south,ct3sn,MPI_REAL,mysouth,tag4,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      if(timestats.ge.1) time_mps1=time_mps1+mytime()

      end subroutine comm_3t_start


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_3t_end(t,west,newwest,east,neweast,   &
                               south,newsouth,north,newnorth,reqs)
      use input
      use mpi
      implicit none

      real t(ib:ie,jb:je,kb:ke+1)
      real west(cmp,nj,nk+1),newwest(cmp,nj,nk+1)
      real east(cmp,nj,nk+1),neweast(cmp,nj,nk+1)
      real south(ni,cmp,nk+1),newsouth(ni,cmp,nk+1)
      real north(ni,cmp,nk+1),newnorth(ni,cmp,nk+1)
      integer reqs(8)

      integer i,j,k,nn,nr,index
      integer :: index_east,index_west,index_south,index_north
      integer, dimension(mpi_status_size,8) :: status1

!-------------------------------------------------------------------

      index_east = -1
      index_west = -1
      index_south = -1
      index_north = -1

      nr = 0
      if(ibe.eq.0)then
        nr = nr + 1
        index_east = nr
      endif
      if(ibw.eq.0)then
        nr = nr + 1
        index_west = nr
      endif
      if(ibs.eq.0)then
        nr = nr + 1
        index_south = nr
      endif
      if(ibn.eq.0)then
        nr = nr + 1
        index_north = nr
      endif

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(1:nr),index,MPI_STATUS_IGNORE,ierr)
        nn = nn + 1

      if(index.eq.index_east)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nkp1
        do j=1,nj
        do i=1,cmp
          t(ni+i,j,k)=neweast(i,j,k)
        enddo
        enddo
        enddo
      elseif(index.eq.index_west)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nkp1
        do j=1,nj
        do i=1,cmp
          t(i-cmp,j,k)=newwest(i,j,k)
        enddo
        enddo
        enddo
      elseif(index.eq.index_south)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nkp1
        do j=1,cmp
        do i=1,ni
          t(i,j-cmp,k)=newsouth(i,j,k)
        enddo
        enddo
        enddo
      elseif(index.eq.index_north)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nkp1
        do j=1,cmp
        do i=1,ni
          t(i,nj+j,k)=newnorth(i,j,k)
        enddo
        enddo
        enddo
      endif

      enddo

      if(timestats.ge.1) time_mps2=time_mps2+mytime()

!----------
!  patch for corner

     if( (ebc.eq.2.or.wbc.eq.2).and.(sbc.eq.1.or.nbc.eq.1) )then

        if(ibw.eq.1)then

          if(p2tchsww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(0,0,k)=t(1,0,k)
            enddo
          endif

          if(p2tchnww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(0,nj+1,k)=t(1,nj+1,k)
            enddo
          endif

        endif

        if(ibe.eq.1)then

          if(p2tchsee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(ni+1,0,k)=t(ni,0,k)
            enddo
          endif

          if(p2tchnee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(ni+1,nj+1,k)=t(ni,nj+1,k)
            enddo
          endif

        endif

      endif

      if( (ebc.eq.1.or.wbc.eq.1).and.(sbc.eq.2.or.nbc.eq.2) )then

        if(ibs.eq.1)then

          if(p2tchsws)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(0,0,k)=t(0,1,k)
            enddo
          endif

          if(p2tchses)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(ni+1,0,k)=t(ni+1,1,k)
            enddo
          endif

        endif

        if(ibn.eq.1)then

          if(p2tchnwn)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(0,nj+1,k)=t(0,nj,k)
            enddo
          endif

          if(p2tchnen)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(ni+1,nj+1,k)=t(ni+1,nj,k)
            enddo
          endif

        endif

      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

!----------

      nr = 0

      if(ibw.eq.0)then
        nr = nr+1
      endif
      if(ibe.eq.0)then
        nr = nr+1
      endif
      if(ibn.eq.0)then
        nr = nr+1
      endif
      if(ibs.eq.0)then
        nr = nr+1
      endif

    if( nr.ge.1 )then
      call MPI_WAITALL(nr,reqs(5:5+nr-1),status1(1:mpi_status_size,5:5+nr-1),ierr)
    endif

      if(timestats.ge.1) time_mps2=time_mps2+mytime()

      end subroutine comm_3t_end
      ! tk2


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_3u_start(u,west,newwest,east,neweast,   &
                                 south,newsouth,north,newnorth,reqs)
      use input
      use mpi
      implicit none
 
      real u(ib:ie+1,jb:je,kb:ke)
      real west(cmp,nj,nk),newwest(cmp,nj,nk)
      real east(cmp,nj,nk),neweast(cmp,nj,nk)
      real south(ni+1,cmp,nk),newsouth(ni+1,cmp,nk)
      real north(ni+1,cmp,nk),newnorth(ni+1,cmp,nk)
      integer reqs(8)
 
      integer i,j,k,nr
      integer tag1,tag2,tag3,tag4
 
!------------------------------------------------

      nr = 0

      nu=nu+1
      tag1=1000+nu

      ! receive east
      if(ibe.eq.0)then
        nr = nr + 1
        call mpi_irecv(neweast,cs3we,MPI_REAL,myeast,tag1,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------
 
      nu=nu+1
      tag2=1000+nu
 
      ! receive west
      if(ibw.eq.0)then
        nr = nr + 1
        call mpi_irecv(newwest,cs3we,MPI_REAL,mywest,tag2,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------
 
      nu=nu+1
      tag3=1000+nu

      ! receive north
      if(ibn.eq.0)then
        nr = nr + 1
        call mpi_irecv(newnorth,cu3sn,MPI_REAL,mynorth,tag3,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------
 
      nu=nu+1
      tag4=1000+nu

      ! receive south
      if(ibs.eq.0)then
        nr = nr + 1
        call mpi_irecv(newsouth,cu3sn,MPI_REAL,mysouth,tag4,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!------------------------------------------------
!------------------------------------------------
!------------------------------------------------

      nr = 4
 
      ! send west
      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,cmp
          west(i,j,k)=u(i+1,j,k)
        enddo
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(west,cs3we,MPI_REAL,mywest,tag1,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------

      ! send east
      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,cmp
          east(i,j,k)=u(ni-cmp+i,j,k)
        enddo
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(east,cs3we,MPI_REAL,myeast,tag2,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------
 
      ! send south
      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,cmp
        do i=1,ni+1
          south(i,j,k)=u(i,j,k)
        enddo
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(south,cu3sn,MPI_REAL,mysouth,tag3,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------
 
      ! send north
      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,cmp
        do i=1,ni+1
          north(i,j,k)=u(i,nj-cmp+j,k)
        enddo
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(north,cu3sn,MPI_REAL,mynorth,tag4,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------

      if(timestats.ge.1) time_mpu1=time_mpu1+mytime()
 
      end subroutine comm_3u_start


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_3u_end(u,west,newwest,east,neweast,   &
                               south,newsouth,north,newnorth,reqs)
      use input
      use mpi
      implicit none
 
      real u(ib:ie+1,jb:je,kb:ke)
      real west(cmp,nj,nk),newwest(cmp,nj,nk)
      real east(cmp,nj,nk),neweast(cmp,nj,nk)
      real south(ni+1,cmp,nk),newsouth(ni+1,cmp,nk)
      real north(ni+1,cmp,nk),newnorth(ni+1,cmp,nk)
      integer reqs(8)
 
      integer i,j,k,nn,nr,index
      integer :: index_east,index_west,index_south,index_north
      integer, dimension(mpi_status_size,8) :: status1

!----------

      index_east = -1
      index_west = -1
      index_south = -1
      index_north = -1

      nr = 0
      if(ibe.eq.0)then
        nr = nr + 1
        index_east = nr
      endif
      if(ibw.eq.0)then
        nr = nr + 1
        index_west = nr
      endif
      if(ibn.eq.0)then
        nr = nr + 1
        index_north = nr
      endif
      if(ibs.eq.0)then
        nr = nr + 1
        index_south = nr
      endif

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(1:nr),index,MPI_STATUS_IGNORE,ierr)
        nn = nn + 1
 
      if(index.eq.index_east)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,cmp
          u(ni+1+i,j,k)=neweast(i,j,k)
        enddo
        enddo
        enddo
      elseif(index.eq.index_west)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,cmp
          u(i-cmp,j,k)=newwest(i,j,k)
        enddo
        enddo
        enddo
      elseif(index.eq.index_north)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,cmp
        do i=1,ni+1
          u(i,nj+j,k)=newnorth(i,j,k)
        enddo
        enddo
        enddo
      elseif(index.eq.index_south)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,cmp
        do i=1,ni+1
          u(i,j-cmp,k)=newsouth(i,j,k)
        enddo
        enddo
        enddo
      endif

      enddo

      if(timestats.ge.1) time_mpu2=time_mpu2+mytime()

!----------
!  patch for corner
 
      if( (ebc.eq.2.or.wbc.eq.2).and.(sbc.eq.1.or.nbc.eq.1) )then
 
        if(ibw.eq.1)then
 
          if(p2tchsww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(0,0,k)=u(1,0,k)
            enddo
          endif
 
          if(p2tchnww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(0,nj+1,k)=u(1,nj+1,k)
            enddo
          endif
 
        endif
 
        if(ibe.eq.1)then
 
          if(p2tchsee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(ni+2,0,k)=u(ni+1,0,k)
            enddo
          endif
 
          if(p2tchnee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(ni+2,nj+1,k)=u(ni+1,nj+1,k)
            enddo
          endif
 
        endif
 
      endif

      if( (ebc.eq.1.or.wbc.eq.1).and.(sbc.eq.2.or.nbc.eq.2) )then
 
        if(ibs.eq.1)then
 
          if(p2tchsws)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(0,0,k)=u(0,1,k)
            enddo
          endif
 
          if(p2tchses)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(ni+2,0,k)=u(ni+2,1,k)
            enddo
          endif
 
        endif
 
        if(ibn.eq.1)then
 
          if(p2tchnwn)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(0,nj+1,k)=u(0,nj,k)
            enddo
          endif
 
          if(p2tchnen)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(ni+2,nj+1,k)=u(ni+2,nj,k)
            enddo
          endif
 
        endif
 
      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

!----------

      nr = 0

      if(ibw.eq.0)then
        nr = nr+1
      endif
      if(ibe.eq.0)then
        nr = nr+1
      endif
      if(ibn.eq.0)then
        nr = nr+1
      endif
      if(ibs.eq.0)then
        nr = nr+1
      endif

    if( nr.ge.1 )then
      call MPI_WAITALL(nr,reqs(5:5+nr-1),status1(1:mpi_status_size,5:5+nr-1),ierr)
    endif

      if(timestats.ge.1) time_mpu2=time_mpu2+mytime()
 
!----------
 
      end subroutine comm_3u_end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_3v_start(v,west,newwest,east,neweast,   &
                                 south,newsouth,north,newnorth,reqs)
      use input
      use mpi
      implicit none
 
      real v(ib:ie,jb:je+1,kb:ke)
      real west(cmp,nj+1,nk),newwest(cmp,nj+1,nk)
      real east(cmp,nj+1,nk),neweast(cmp,nj+1,nk)
      real south(ni,cmp,nk),newsouth(ni,cmp,nk)
      real north(ni,cmp,nk),newnorth(ni,cmp,nk)
      integer reqs(8)
 
      integer i,j,k,nr
      integer tag1,tag2,tag3,tag4
 
!------------------------------------------------

      nr = 0

      nv=nv+1
      tag1=2000+nv

      ! receive east
      if(ibe.eq.0)then
        nr = nr + 1
        call mpi_irecv(neweast,cv3we,MPI_REAL,myeast,tag1,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------
 
      nv=nv+1
      tag2=2000+nv
 
      ! receive west
      if(ibw.eq.0)then
        nr = nr + 1
        call mpi_irecv(newwest,cv3we,MPI_REAL,mywest,tag2,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------
 
      nv=nv+1
      tag3=2000+nv
 
      ! receive north
      if(ibn.eq.0)then
        nr = nr + 1
        call mpi_irecv(newnorth,cs3sn,MPI_REAL,mynorth,tag3,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------
 
      nv=nv+1
      tag4=2000+nv
 
      ! receive south
      if(ibs.eq.0)then
        nr = nr + 1
        call mpi_irecv(newsouth,cs3sn,MPI_REAL,mysouth,tag4,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!------------------------------------------------
!------------------------------------------------
!------------------------------------------------

      nr = 4
 
      ! send west
      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,cmp
          west(i,j,k)=v(i,j,k)
        enddo
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(west,cv3we,MPI_REAL,mywest,tag1,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------

      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,cmp
          east(i,j,k)=v(ni-cmp+i,j,k)
        enddo
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(east,cv3we,MPI_REAL,myeast,tag2,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      ! send south
      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,cmp
        do i=1,ni
          south(i,j,k)=v(i,j+1,k)
        enddo
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(south,cs3sn,MPI_REAL,mysouth,tag3,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------

      ! send north
      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,cmp
        do i=1,ni
          north(i,j,k)=v(i,nj-cmp+j,k)
        enddo
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(north,cs3sn,MPI_REAL,mynorth,tag4,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------

      if(timestats.ge.1) time_mpv1=time_mpv1+mytime()
 
      end subroutine comm_3v_start

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_3v_end(v,west,newwest,east,neweast,   &
                               south,newsouth,north,newnorth,reqs)
      use input
      use mpi
      implicit none
 
      real v(ib:ie,jb:je+1,kb:ke)
      real west(cmp,nj+1,nk),newwest(cmp,nj+1,nk)
      real east(cmp,nj+1,nk),neweast(cmp,nj+1,nk)
      real south(ni,cmp,nk),newsouth(ni,cmp,nk)
      real north(ni,cmp,nk),newnorth(ni,cmp,nk)
      integer reqs(8)
 
      integer i,j,k,nn,nr,index
      integer :: index_east,index_west,index_south,index_north
      integer, dimension(mpi_status_size,8) :: status1
 
!--------

      index_east = -1
      index_west = -1
      index_south = -1
      index_north = -1

      nr = 0
      if(ibe.eq.0)then
        nr = nr + 1
        index_east = nr
      endif
      if(ibw.eq.0)then
        nr = nr + 1
        index_west = nr
      endif
      if(ibn.eq.0)then
        nr = nr + 1
        index_north = nr
      endif
      if(ibs.eq.0)then
        nr = nr + 1
        index_south = nr
      endif

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(1:nr),index,MPI_STATUS_IGNORE,ierr)
        nn = nn + 1
 
      if(index.eq.index_east)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,cmp
          v(ni+i,j,k)=neweast(i,j,k)
        enddo
        enddo
        enddo
      elseif(index.eq.index_west)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,cmp
          v(i-cmp,j,k)=newwest(i,j,k)
        enddo
        enddo
        enddo
      elseif(index.eq.index_north)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,cmp
        do i=1,ni
          v(i,nj+1+j,k)=newnorth(i,j,k)
        enddo
        enddo
        enddo
      elseif(index.eq.index_south)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,cmp
        do i=1,ni
          v(i,j-cmp,k)=newsouth(i,j,k)
        enddo
        enddo
        enddo
      endif

      enddo

      if(timestats.ge.1) time_mpv2=time_mpv2+mytime()

!----------
!  patch for corner
 
      if( (ebc.eq.2.or.wbc.eq.2).and.(sbc.eq.1.or.nbc.eq.1) )then
 
        if(ibw.eq.1)then
 
          if(p2tchsww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(0,0,k)=v(1,0,k)
            enddo
          endif
 
          if(p2tchnww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(0,nj+2,k)=v(1,nj+2,k)
            enddo
          endif
 
        endif
 
        if(ibe.eq.1)then
 
          if(p2tchsee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(ni+1,0,k)=v(ni,0,k)
            enddo
          endif
 
          if(p2tchnee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(ni+1,nj+2,k)=v(ni,nj+2,k)
            enddo
          endif
 
        endif
 
      endif

      if( (ebc.eq.1.or.wbc.eq.1).and.(sbc.eq.2.or.nbc.eq.2) )then
 
        if(ibs.eq.1)then
 
          if(p2tchsws)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(0,0,k)=v(0,1,k)
            enddo
          endif
 
          if(p2tchses)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(ni+1,0,k)=v(ni+1,1,k)
            enddo
          endif
 
        endif
 
        if(ibn.eq.1)then
 
          if(p2tchnwn)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(0,nj+2,k)=v(0,nj+1,k)
            enddo
          endif
 
          if(p2tchnen)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(ni+1,nj+2,k)=v(ni+1,nj+1,k)
            enddo
          endif
 
        endif
 
      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

!--------

      nr = 0

      if(ibw.eq.0)then
        nr = nr+1
      endif
      if(ibe.eq.0)then
        nr = nr+1
      endif
      if(ibn.eq.0)then
        nr = nr+1
      endif
      if(ibs.eq.0)then
        nr = nr+1
      endif

    if( nr.ge.1 )then
      call MPI_WAITALL(nr,reqs(5:5+nr-1),status1(1:mpi_status_size,5:5+nr-1),ierr)
    endif

      if(timestats.ge.1) time_mpv2=time_mpv2+mytime()
 
!----------
 
      end subroutine comm_3v_end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_3w_start(w,west,newwest,east,neweast,   &
                                 south,newsouth,north,newnorth,reqs)
      use input
      use mpi
      implicit none
 
      real w(ib:ie,jb:je,kb:ke+1)
      real west(cmp,nj,nk-1),newwest(cmp,nj,nk-1)
      real east(cmp,nj,nk-1),neweast(cmp,nj,nk-1)
      real south(ni,cmp,nk-1),newsouth(ni,cmp,nk-1)
      real north(ni,cmp,nk-1),newnorth(ni,cmp,nk-1)
      integer reqs(8)
 
      integer i,j,k,nr
      integer tag1,tag2,tag3,tag4
 
!------------------------------------------------

      nr = 0

      nw=nw+1
      tag1=3000+nw

      ! receive east
      if(ibe.eq.0)then
        nr = nr + 1
        call mpi_irecv(neweast,cw3we,MPI_REAL,myeast,tag1,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------
 
      nw=nw+1
      tag2=3000+nw
 
      ! receive west
      if(ibw.eq.0)then
        nr = nr + 1
        call mpi_irecv(newwest,cw3we,MPI_REAL,mywest,tag2,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------

      nw=nw+1
      tag3=3000+nw
 
      ! receive north
      if(ibn.eq.0)then
        nr = nr + 1
        call mpi_irecv(newnorth,cw3sn,MPI_REAL,mynorth,tag3,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------
 
      nw=nw+1
      tag4=3000+nw
 
      ! receive south
      if(ibs.eq.0)then
        nr = nr + 1
        call mpi_irecv(newsouth,cw3sn,MPI_REAL,mysouth,tag4,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!------------------------------------------------
!------------------------------------------------
!------------------------------------------------

      nr = 4
 
      ! send west
      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,cmp
          west(i,j,k-1)=w(i,j,k)
        enddo
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(west,cw3we,MPI_REAL,mywest,tag1,    &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------

      ! send east
      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,cmp
          east(i,j,k-1)=w(ni-cmp+i,j,k)
        enddo
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(east,cw3we,MPI_REAL,myeast,tag2,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      ! send south
      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,cmp
        do i=1,ni
          south(i,j,k-1)=w(i,j,k)
        enddo
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(south,cw3sn,MPI_REAL,mysouth,tag3,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------

      ! send north
      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,cmp
        do i=1,ni
          north(i,j,k-1)=w(i,nj-cmp+j,k)
        enddo
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(north,cw3sn,MPI_REAL,mynorth,tag4,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------

      if(timestats.ge.1) time_mpw1=time_mpw1+mytime()
 
      end subroutine comm_3w_start

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_3w_end(w,west,newwest,east,neweast,   &
                               south,newsouth,north,newnorth,reqs)
      use input
      use mpi
      implicit none
 
      real w(ib:ie,jb:je,kb:ke+1)
      real west(cmp,nj,nk-1),newwest(cmp,nj,nk-1)
      real east(cmp,nj,nk-1),neweast(cmp,nj,nk-1)
      real south(ni,cmp,nk-1),newsouth(ni,cmp,nk-1)
      real north(ni,cmp,nk-1),newnorth(ni,cmp,nk-1)
      integer reqs(8)
 
      integer i,j,k,nn,nr,index
      integer :: index_east,index_west,index_south,index_north
      integer, dimension(mpi_status_size,8) :: status1

!--------

      index_east = -1
      index_west = -1
      index_south = -1
      index_north = -1

      nr = 0
      if(ibe.eq.0)then
        nr = nr + 1
        index_east = nr
      endif
      if(ibw.eq.0)then
        nr = nr + 1
        index_west = nr
      endif
      if(ibn.eq.0)then
        nr = nr + 1
        index_north = nr
      endif
      if(ibs.eq.0)then
        nr = nr + 1
        index_south = nr
      endif

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(1:nr),index,MPI_STATUS_IGNORE,ierr)
        nn = nn + 1
 
      if(index.eq.index_east)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,cmp
          w(ni+i,j,k)=neweast(i,j,k-1)
        enddo
        enddo
        enddo
      elseif(index.eq.index_west)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,cmp
          w(i-cmp,j,k)=newwest(i,j,k-1)
        enddo
        enddo
        enddo
      elseif(index.eq.index_north)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,cmp
        do i=1,ni
          w(i,nj+j,k)=newnorth(i,j,k-1)
        enddo
        enddo
        enddo
      elseif(index.eq.index_south)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,cmp
        do i=1,ni
          w(i,j-cmp,k)=newsouth(i,j,k-1)
        enddo
        enddo
        enddo
      endif

      enddo
 
      if(timestats.ge.1) time_mpw2=time_mpw2+mytime()

!----------
!  patch for corner
 
      if( (ebc.eq.2.or.wbc.eq.2).and.(sbc.eq.1.or.nbc.eq.1) )then
 
        if(ibw.eq.1)then
 
          if(p2tchsww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(0,0,k)=w(1,0,k)
            enddo
          endif
 
          if(p2tchnww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(0,nj+1,k)=w(1,nj+1,k)
            enddo
          endif
 
        endif
 
        if(ibe.eq.1)then
 
          if(p2tchsee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(ni+1,0,k)=w(ni,0,k)
            enddo
          endif
 
          if(p2tchnee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(ni+1,nj+1,k)=w(ni,nj+1,k)
            enddo
          endif
 
        endif
 
      endif

      if( (ebc.eq.1.or.wbc.eq.1).and.(sbc.eq.2.or.nbc.eq.2) )then
 
        if(ibs.eq.1)then
 
          if(p2tchsws)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(0,0,k)=w(0,1,k)
            enddo
          endif
 
          if(p2tchses)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(ni+1,0,k)=w(ni+1,1,k)
            enddo
          endif
 
        endif
 
        if(ibn.eq.1)then
 
          if(p2tchnwn)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(0,nj+1,k)=w(0,nj,k)
            enddo
          endif
 
          if(p2tchnen)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(ni+1,nj+1,k)=w(ni+1,nj,k)
            enddo
          endif
 
        endif
 
      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

!--------

      nr = 0

      if(ibw.eq.0)then
        nr = nr+1
      endif
      if(ibe.eq.0)then
        nr = nr+1
      endif
      if(ibn.eq.0)then
        nr = nr+1
      endif
      if(ibs.eq.0)then
        nr = nr+1
      endif

    if( nr.ge.1 )then
      call MPI_WAITALL(nr,reqs(5:5+nr-1),status1(1:mpi_status_size,5:5+nr-1),ierr)
    endif

      if(timestats.ge.1) time_mpw2=time_mpw2+mytime()
 
!----------
 
      end subroutine comm_3w_end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_1s2d_start(s,west,newwest,east,neweast,   &
                                 south,newsouth,north,newnorth,reqs)
      use input
      use mpi
      implicit none
 
      real s(ib:ie,jb:je)
      real west(nj),newwest(nj)
      real east(nj),neweast(nj)
      real south(ni),newsouth(ni)
      real north(ni),newnorth(ni)
      integer reqs(8)
 
      integer i,j,nr
      integer tag1,tag2,tag3,tag4
 
!------------------------------------------------

      nr = 0

      nf=nf+1
      tag1=nf

      ! receive east
      if(ibe.eq.0)then
        nr = nr + 1
        call mpi_irecv(neweast,nj,MPI_REAL,myeast,tag1,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      nf=nf+1
      tag2=nf

      ! receive west
      if(ibw.eq.0)then
        nr = nr + 1
        call mpi_irecv(newwest,nj,MPI_REAL,mywest,tag2,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      nf=nf+1
      tag3=nf

      ! receive south
      if(ibs.eq.0)then
        nr = nr + 1
        call mpi_irecv(newsouth,ni,MPI_REAL,mysouth,tag3,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      nf=nf+1
      tag4=nf

      ! receive north
      if(ibn.eq.0)then
        nr = nr + 1
        call mpi_irecv(newnorth,ni,MPI_REAL,mynorth,tag4,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!------------------------------------------------
!------------------------------------------------
!------------------------------------------------

      nr = 4
 
      ! send west
      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(j)
        do j=1,nj
          west(j)=s(1,j)
        enddo
        nr = nr + 1
        call mpi_isend(west,nj,MPI_REAL,mywest,tag1,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------
 
      ! send east
      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(j)
        do j=1,nj
          east(j)=s(ni,j)
        enddo
        nr = nr + 1
        call mpi_isend(east,nj,MPI_REAL,myeast,tag2,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------
 
      ! send north
      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i)
        do i=1,ni
          north(i)=s(i,nj)
        enddo
        nr = nr + 1
        call mpi_isend(north,ni,MPI_REAL,mynorth,tag3,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------
 
      ! send south
      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i)
        do i=1,ni
          south(i)=s(i,1)
        enddo
        nr = nr + 1
        call mpi_isend(south,ni,MPI_REAL,mysouth,tag4,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------
 
      if(timestats.ge.1) time_mps1=time_mps1+mytime()
 
      end subroutine comm_1s2d_start

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_1s2d_end(s,west,newwest,east,neweast,   &
                               south,newsouth,north,newnorth,reqs)
      use input
      use mpi
      implicit none
 
      real s(ib:ie,jb:je)
      real west(nj),newwest(nj)
      real east(nj),neweast(nj)
      real south(ni),newsouth(ni)
      real north(ni),newnorth(ni)
      integer reqs(8)
 
      integer i,j,nn,nr,index
      integer :: index_east,index_west,index_south,index_north
      integer, dimension(mpi_status_size,8) :: status1

!---------------------------------------------------------------------

      index_east = -1
      index_west = -1
      index_south = -1
      index_north = -1

      nr = 0
      if(ibe.eq.0)then
        nr = nr + 1
        index_east = nr
      endif
      if(ibw.eq.0)then
        nr = nr + 1
        index_west = nr
      endif
      if(ibs.eq.0)then
        nr = nr + 1
        index_south = nr
      endif
      if(ibn.eq.0)then
        nr = nr + 1
        index_north = nr
      endif

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(1:nr),index,MPI_STATUS_IGNORE,ierr)
        nn = nn + 1
 
      if(index.eq.index_east)then
!$omp parallel do default(shared)   &
!$omp private(j)
        do j=1,nj
          s(ni+1,j)=neweast(j)
        enddo
      elseif(index.eq.index_west)then
!$omp parallel do default(shared)   &
!$omp private(j)
        do j=1,nj
          s(0,j)=newwest(j)
        enddo
      elseif(index.eq.index_south)then
!$omp parallel do default(shared)   &
!$omp private(i)
        do i=1,ni
          s(i,0)=newsouth(i)
        enddo
      elseif(index.eq.index_north)then
!$omp parallel do default(shared)   &
!$omp private(i)
        do i=1,ni
          s(i,nj+1)=newnorth(i)
        enddo
      endif

      enddo

      if(timestats.ge.1) time_mps2=time_mps2+mytime()

!----------
!  patch for corner
 
      if( (ebc.eq.2.or.wbc.eq.2).and.(sbc.eq.1.or.nbc.eq.1) )then
 
        if(ibw.eq.1)then
 
          if(p2tchsww)then
              s(0,0)=s(1,0)
          endif
 
          if(p2tchnww)then
              s(0,nj+1)=s(1,nj+1)
          endif
 
        endif
 
        if(ibe.eq.1)then
 
          if(p2tchsee)then
              s(ni+1,0)=s(ni,0)
          endif
 
          if(p2tchnee)then
              s(ni+1,nj+1)=s(ni,nj+1)
          endif
 
        endif
 
      endif

      if( (ebc.eq.1.or.wbc.eq.1).and.(sbc.eq.2.or.nbc.eq.2) )then
 
        if(ibs.eq.1)then
 
          if(p2tchsws)then
              s(0,0)=s(0,1)
          endif
 
          if(p2tchses)then
              s(ni+1,0)=s(ni+1,1)
          endif
 
        endif
 
        if(ibn.eq.1)then
 
          if(p2tchnwn)then
              s(0,nj+1)=s(0,nj)
          endif
 
          if(p2tchnen)then
              s(ni+1,nj+1)=s(ni+1,nj)
          endif
 
        endif
 
      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

!----------

      nr = 0

      if(ibw.eq.0)then
        nr = nr+1
      endif
      if(ibe.eq.0)then
        nr = nr+1
      endif
      if(ibn.eq.0)then
        nr = nr+1
      endif
      if(ibs.eq.0)then
        nr = nr+1
      endif

    if( nr.ge.1 )then
      call MPI_WAITALL(nr,reqs(5:5+nr-1),status1(1:mpi_status_size,5:5+nr-1),ierr)
    endif

      if(timestats.ge.1) time_mps2=time_mps2+mytime()
 
!----------
 
      end subroutine comm_1s2d_end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_1s_start(s,west,newwest,east,neweast,   &
                                 south,newsouth,north,newnorth,reqs)
      use input
      use mpi
      implicit none
 
      real s(ib:ie,jb:je,kb:ke)
      real west(nj,nk),newwest(nj,nk)
      real east(nj,nk),neweast(nj,nk)
      real south(ni,nk),newsouth(ni,nk)
      real north(ni,nk),newnorth(ni,nk)
      integer reqs(8)
 
      integer i,j,k,nr
      integer tag1,tag2,tag3,tag4
 
!------------------------------------------------

      nr = 0

      nf=nf+1
      tag1=nf

      ! receive east
      if(ibe.eq.0)then
        nr = nr + 1
        call mpi_irecv(neweast,cs1we,MPI_REAL,myeast,tag1,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      nf=nf+1
      tag2=nf

      ! receive west
      if(ibw.eq.0)then
        nr = nr + 1
        call mpi_irecv(newwest,cs1we,MPI_REAL,mywest,tag2,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      nf=nf+1
      tag3=nf

      ! receive south
      if(ibs.eq.0)then
        nr = nr + 1
        call mpi_irecv(newsouth,cs1sn,MPI_REAL,mysouth,tag3,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      nf=nf+1
      tag4=nf

      ! receive north
      if(ibn.eq.0)then
        nr = nr + 1
        call mpi_irecv(newnorth,cs1sn,MPI_REAL,mynorth,tag4,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!------------------------------------------------
!------------------------------------------------
!------------------------------------------------

      nr = 4
 
      ! send west
      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          west(j,k)=s(1,j,k)
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(west,cs1we,MPI_REAL,mywest,tag1,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------
 
      ! send east
      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          east(j,k)=s(ni,j,k)
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(east,cs1we,MPI_REAL,myeast,tag2,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------
 
      ! send north
      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          north(i,k)=s(i,nj,k)
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(north,cs1sn,MPI_REAL,mynorth,tag3,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------
 
      ! send south
      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          south(i,k)=s(i,1,k)
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(south,cs1sn,MPI_REAL,mysouth,tag4,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------
 
      if(timestats.ge.1) time_mps1=time_mps1+mytime()
 
      end subroutine comm_1s_start


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_1s_end(s,west,newwest,east,neweast,   &
                               south,newsouth,north,newnorth,reqs)
      use input
      use mpi
      implicit none
 
      real s(ib:ie,jb:je,kb:ke)
      real west(nj,nk),newwest(nj,nk)
      real east(nj,nk),neweast(nj,nk)
      real south(ni,nk),newsouth(ni,nk)
      real north(ni,nk),newnorth(ni,nk)
      integer reqs(8)
 
      integer i,j,k,nn,nr,index
      integer :: index_east,index_west,index_south,index_north
      integer, dimension(mpi_status_size,8) :: status1

!---------------------------------------------------------------------

      index_east = -1
      index_west = -1
      index_south = -1
      index_north = -1

      nr = 0
      if(ibe.eq.0)then
        nr = nr + 1
        index_east = nr
      endif
      if(ibw.eq.0)then
        nr = nr + 1
        index_west = nr
      endif
      if(ibs.eq.0)then
        nr = nr + 1
        index_south = nr
      endif
      if(ibn.eq.0)then
        nr = nr + 1
        index_north = nr
      endif

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(1:nr),index,MPI_STATUS_IGNORE,ierr)
        nn = nn + 1
 
      if(index.eq.index_east)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          s(ni+1,j,k)=neweast(j,k)
        enddo
        enddo
      elseif(index.eq.index_west)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          s(0,j,k)=newwest(j,k)
        enddo
        enddo
      elseif(index.eq.index_south)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          s(i,0,k)=newsouth(i,k)
        enddo
        enddo
      elseif(index.eq.index_north)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          s(i,nj+1,k)=newnorth(i,k)
        enddo
        enddo
      endif

      enddo

      if(timestats.ge.1) time_mps2=time_mps2+mytime()

!----------
!  patch for corner
 
      if( (ebc.eq.2.or.wbc.eq.2).and.(sbc.eq.1.or.nbc.eq.1) )then
 
        if(ibw.eq.1)then
 
          if(p2tchsww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(0,0,k)=s(1,0,k)
            enddo
          endif
 
          if(p2tchnww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(0,nj+1,k)=s(1,nj+1,k)
            enddo
          endif
 
        endif
 
        if(ibe.eq.1)then
 
          if(p2tchsee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(ni+1,0,k)=s(ni,0,k)
            enddo
          endif
 
          if(p2tchnee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(ni+1,nj+1,k)=s(ni,nj+1,k)
            enddo
          endif
 
        endif
 
      endif

      if( (ebc.eq.1.or.wbc.eq.1).and.(sbc.eq.2.or.nbc.eq.2) )then
 
        if(ibs.eq.1)then
 
          if(p2tchsws)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(0,0,k)=s(0,1,k)
            enddo
          endif
 
          if(p2tchses)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(ni+1,0,k)=s(ni+1,1,k)
            enddo
          endif
 
        endif
 
        if(ibn.eq.1)then
 
          if(p2tchnwn)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(0,nj+1,k)=s(0,nj,k)
            enddo
          endif
 
          if(p2tchnen)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(ni+1,nj+1,k)=s(ni+1,nj,k)
            enddo
          endif
 
        endif
 
      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

!----------

      nr = 0

      if(ibw.eq.0)then
        nr = nr+1
      endif
      if(ibe.eq.0)then
        nr = nr+1
      endif
      if(ibn.eq.0)then
        nr = nr+1
      endif
      if(ibs.eq.0)then
        nr = nr+1
      endif

    if( nr.ge.1 )then
      call MPI_WAITALL(nr,reqs(5:5+nr-1),status1(1:mpi_status_size,5:5+nr-1),ierr)
    endif

      if(timestats.ge.1) time_mps2=time_mps2+mytime()
 
!----------
 
      end subroutine comm_1s_end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_1p_start(s,west,newwest,east,neweast,   &
                                 south,newsouth,north,newnorth,reqs)
      use input
      use mpi
      implicit none
 
      real s(ib:ie,jb:je,kb:ke)
      real west(nj,nk),newwest(nj,nk)
      real east(nj,nk),neweast(nj,nk)
      real south(ni,nk),newsouth(ni,nk)
      real north(ni,nk),newnorth(ni,nk)
      integer reqs(8)
 
      integer i,j,k,nr
      integer tag1,tag2,tag3,tag4
 
!------------------------------------------------

      nr = 0

      nf=nf+1
      tag1=nf

      ! receive east
      if(ibe.eq.0)then
        nr = nr + 1
        call mpi_irecv(neweast,cs1we,MPI_REAL,myeast,tag1,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      nf=nf+1
      tag2=nf

      ! receive west
      if(ibw.eq.0)then
        nr = nr + 1
        call mpi_irecv(newwest,cs1we,MPI_REAL,mywest,tag2,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      nf=nf+1
      tag3=nf

      ! receive south
      if(ibs.eq.0)then
        nr = nr + 1
        call mpi_irecv(newsouth,cs1sn,MPI_REAL,mysouth,tag3,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      nf=nf+1
      tag4=nf

      ! receive north
      if(ibn.eq.0)then
        nr = nr + 1
        call mpi_irecv(newnorth,cs1sn,MPI_REAL,mynorth,tag4,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!------------------------------------------------
!------------------------------------------------
!------------------------------------------------

      nr = 4
 
      ! send west
      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          west(j,k)=s(1,j,k)
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(west,cs1we,MPI_REAL,mywest,tag1,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------
 
      ! send east
      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          east(j,k)=s(ni,j,k)
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(east,cs1we,MPI_REAL,myeast,tag2,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------
 
      ! send north
      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          north(i,k)=s(i,nj,k)
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(north,cs1sn,MPI_REAL,mynorth,tag3,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------
 
      ! send south
      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          south(i,k)=s(i,1,k)
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(south,cs1sn,MPI_REAL,mysouth,tag4,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------
 
      if(timestats.ge.1) time_mpp1=time_mpp1+mytime()
 
      end subroutine comm_1p_start


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_1p_end1(s,west,newwest,east,neweast,   &
                               south,newsouth,north,newnorth,reqs)
      use input
      use mpi
      implicit none
 
      real s(ib:ie,jb:je,kb:ke)
      real west(nj,nk),newwest(nj,nk)
      real east(nj,nk),neweast(nj,nk)
      real south(ni,nk),newsouth(ni,nk)
      real north(ni,nk),newnorth(ni,nk)
      integer reqs(8)
 
      integer i,j,k,nn,nr,index
      integer :: index_east,index_west,index_south,index_north

!---------------------------------------------------------------------

      index_east = -1
      index_west = -1
      index_south = -1
      index_north = -1

      nr = 0
      if(ibe.eq.0)then
        nr = nr + 1
        index_east = nr
      endif
      if(ibw.eq.0)then
        nr = nr + 1
        index_west = nr
      endif
      if(ibs.eq.0)then
        nr = nr + 1
        index_south = nr
      endif
      if(ibn.eq.0)then
        nr = nr + 1
        index_north = nr
      endif

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(1:nr),index,MPI_STATUS_IGNORE,ierr)
        nn = nn + 1
 
      if(index.eq.index_east)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          s(ni+1,j,k)=neweast(j,k)
        enddo
        enddo
      elseif(index.eq.index_west)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          s(0,j,k)=newwest(j,k)
        enddo
        enddo
      elseif(index.eq.index_south)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          s(i,0,k)=newsouth(i,k)
        enddo
        enddo
      elseif(index.eq.index_north)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          s(i,nj+1,k)=newnorth(i,k)
        enddo
        enddo
      endif

      enddo

      if(timestats.ge.1) time_mpp2=time_mpp2+mytime()

!----------
 
      end subroutine comm_1p_end1


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_1p_end2(reqs)
      use input, only: ibw,ibe,ibs,ibn,timestats,time_mpp2,mytime
      use mpi
      implicit none
 
      integer reqs(8)
 
      integer :: nr,ierr
      integer, dimension(mpi_status_size,8) :: status1

!---------------------------------------------------------------------

      nr = 0

      if(ibw.eq.0)then
        nr = nr+1
      endif
      if(ibe.eq.0)then
        nr = nr+1
      endif
      if(ibn.eq.0)then
        nr = nr+1
      endif
      if(ibs.eq.0)then
        nr = nr+1
      endif

    if( nr.ge.1 )then
      call MPI_WAITALL(nr,reqs(5:5+nr-1),status1(1:mpi_status_size,5:5+nr-1),ierr)
    endif

      if(timestats.ge.1) time_mpp2=time_mpp2+mytime()
 
!----------
 
      end subroutine comm_1p_end2


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_1t_start(t,west,newwest,east,neweast,   &
                                 south,newsouth,north,newnorth,reqs)
      use input
      use mpi
      implicit none

      real t(ib:ie,jb:je,kb:ke+1)
      real west(nj,nk+1),newwest(nj,nk+1)
      real east(nj,nk+1),neweast(nj,nk+1)
      real south(ni,nk+1),newsouth(ni,nk+1)
      real north(ni,nk+1),newnorth(ni,nk+1)
      integer reqs(8)

      integer i,j,k,nr
      integer tag1,tag2,tag3,tag4

!------------------------------------------------

      nr = 0

      nf=nf+1
      tag1=nf

      ! receive east
      if(ibe.eq.0)then
        nr = nr + 1
        call mpi_irecv(neweast,ct1we,MPI_REAL,myeast,tag1,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      nf=nf+1
      tag2=nf

      ! receive west
      if(ibw.eq.0)then
        nr = nr + 1
        call mpi_irecv(newwest,ct1we,MPI_REAL,mywest,tag2,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      nf=nf+1
      tag3=nf

      ! receive south
      if(ibs.eq.0)then
        nr = nr + 1
        call mpi_irecv(newsouth,ct1sn,MPI_REAL,mysouth,tag3,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      nf=nf+1
      tag4=nf

      ! receive north
      if(ibn.eq.0)then
        nr = nr + 1
        call mpi_irecv(newnorth,ct1sn,MPI_REAL,mynorth,tag4,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!------------------------------------------------
!------------------------------------------------
!------------------------------------------------

      nr = 4

      ! send west
      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nkp1
        do j=1,nj
          west(j,k)=t(1,j,k)
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(west,ct1we,MPI_REAL,mywest,tag1,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      ! send east
      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nkp1
        do j=1,nj
          east(j,k)=t(ni,j,k)
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(east,ct1we,MPI_REAL,myeast,tag2,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      ! send north
      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nkp1
        do i=1,ni
          north(i,k)=t(i,nj,k)
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(north,ct1sn,MPI_REAL,mynorth,tag3,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      ! send south
      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nkp1
        do i=1,ni
          south(i,k)=t(i,1,k)
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(south,ct1sn,MPI_REAL,mysouth,tag4,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      if(timestats.ge.1) time_mps1=time_mps1+mytime()

      end subroutine comm_1t_start

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_1t_end(t,west,newwest,east,neweast,   &
                               south,newsouth,north,newnorth,reqs)
      use input
      use mpi
      implicit none

      real t(ib:ie,jb:je,kb:ke+1)
      real west(nj,nk+1),newwest(nj,nk+1)
      real east(nj,nk+1),neweast(nj,nk+1)
      real south(ni,nk+1),newsouth(ni,nk+1)
      real north(ni,nk+1),newnorth(ni,nk+1)
      integer reqs(8)

      integer i,j,k,nn,nr,index
      integer :: index_east,index_west,index_south,index_north
      integer, dimension(mpi_status_size,8) :: status1

!---------------------------------------------------------------------

      index_east = -1
      index_west = -1
      index_south = -1
      index_north = -1

      nr = 0
      if(ibe.eq.0)then
        nr = nr + 1
        index_east = nr
      endif
      if(ibw.eq.0)then
        nr = nr + 1
        index_west = nr
      endif
      if(ibs.eq.0)then
        nr = nr + 1
        index_south = nr
      endif
      if(ibn.eq.0)then
        nr = nr + 1
        index_north = nr
      endif

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(1:nr),index,MPI_STATUS_IGNORE,ierr)
        nn = nn + 1

      if(index.eq.index_east)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nkp1
        do j=1,nj
          t(ni+1,j,k)=neweast(j,k)
        enddo
        enddo
      elseif(index.eq.index_west)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nkp1
        do j=1,nj
          t(0,j,k)=newwest(j,k)
        enddo
        enddo
      elseif(index.eq.index_south)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nkp1
        do i=1,ni
          t(i,0,k)=newsouth(i,k)
        enddo
        enddo
      elseif(index.eq.index_north)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nkp1
        do i=1,ni
          t(i,nj+1,k)=newnorth(i,k)
        enddo
        enddo
      endif

      enddo

      if(timestats.ge.1) time_mps2=time_mps2+mytime()

!----------
!  patch for corner

      if( (ebc.eq.2.or.wbc.eq.2).and.(sbc.eq.1.or.nbc.eq.1) )then

        if(ibw.eq.1)then

          if(p2tchsww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(0,0,k)=t(1,0,k)
            enddo
          endif

          if(p2tchnww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(0,nj+1,k)=t(1,nj+1,k)
            enddo
          endif

        endif

        if(ibe.eq.1)then

          if(p2tchsee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(ni+1,0,k)=t(ni,0,k)
            enddo
          endif

          if(p2tchnee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(ni+1,nj+1,k)=t(ni,nj+1,k)
            enddo
          endif

        endif

      endif

      if( (ebc.eq.1.or.wbc.eq.1).and.(sbc.eq.2.or.nbc.eq.2) )then

        if(ibs.eq.1)then

          if(p2tchsws)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(0,0,k)=t(0,1,k)
            enddo
          endif

          if(p2tchses)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(ni+1,0,k)=t(ni+1,1,k)
            enddo
          endif

        endif

        if(ibn.eq.1)then

          if(p2tchnwn)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(0,nj+1,k)=t(0,nj,k)
            enddo
          endif

          if(p2tchnen)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(ni+1,nj+1,k)=t(ni+1,nj,k)
            enddo
          endif

        endif

      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

!----------

      nr = 0

      if(ibw.eq.0)then
        nr = nr+1
      endif
      if(ibe.eq.0)then
        nr = nr+1
      endif
      if(ibn.eq.0)then
        nr = nr+1
      endif
      if(ibs.eq.0)then
        nr = nr+1
      endif

    if( nr.ge.1 )then
      call MPI_WAITALL(nr,reqs(5:5+nr-1),status1(1:mpi_status_size,5:5+nr-1),ierr)
    endif

      if(timestats.ge.1) time_mps2=time_mps2+mytime()

!----------

      end subroutine comm_1t_end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_1u_start(u,west,newwest,east,neweast,   &
                                 south,newsouth,north,newnorth,reqs)
      use input
      use mpi
      implicit none
 
      real, intent(in   ), dimension(ib:ie+1,jb:je,kb:ke) :: u
      real west(nj,nk),newwest(nj,nk)
      real east(nj,nk),neweast(nj,nk)
      real south(ni+1,nk),newsouth(ni+1,nk)
      real north(ni+1,nk),newnorth(ni+1,nk)
      integer reqs(8)
 
      integer i,j,k,nr
      integer tag1,tag2,tag3,tag4
 
!------------------------------------------------

      nr = 0

      nu=nu+1
      tag1=1000+nu

      ! receive east
      if(ibe.eq.0)then
        nr = nr + 1
        call mpi_irecv(neweast,cs1we,MPI_REAL,myeast,tag1,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      nu=nu+1
      tag2=1000+nu

      ! receive west
      if(ibw.eq.0)then
        nr = nr + 1
        call mpi_irecv(newwest,cs1we,MPI_REAL,mywest,tag2,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      nu=nu+1
      tag3=1000+nu

      ! receive south
      if(ibs.eq.0)then
        nr = nr + 1
        call mpi_irecv(newsouth,cu1sn,MPI_REAL,mysouth,tag3,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      nu=nu+1
      tag4=1000+nu

      ! receive north
      if(ibn.eq.0)then
        nr = nr + 1
        call mpi_irecv(newnorth,cu1sn,MPI_REAL,mynorth,tag4,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!------------------------------------------------
!------------------------------------------------
!------------------------------------------------

      nr = 4
 
      ! send west
      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          west(j,k)=u(2,j,k)
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(west,cs1we,MPI_REAL,mywest,tag1,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------
 
      ! send east
      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          east(j,k)=u(ni,j,k)
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(east,cs1we,MPI_REAL,myeast,tag2,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------
 
      ! send north
      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni+1
          north(i,k)=u(i,nj,k)
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(north,cu1sn,MPI_REAL,mynorth,tag3,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------
 
      ! send south
      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni+1
          south(i,k)=u(i,1,k)
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(south,cu1sn,MPI_REAL,mysouth,tag4,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------
 
      if(timestats.ge.1) time_mpu1=time_mpu1+mytime()
 
      end subroutine comm_1u_start

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_1u_end(u,west,newwest,east,neweast,   &
                               south,newsouth,north,newnorth,reqs)
      use input
      use mpi
      implicit none
 
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: u
      real west(nj,nk),newwest(nj,nk)
      real east(nj,nk),neweast(nj,nk)
      real south(ni+1,nk),newsouth(ni+1,nk)
      real north(ni+1,nk),newnorth(ni+1,nk)
      integer reqs(8)
 
      integer i,j,k,nn,nr,index
      integer :: index_east,index_west,index_south,index_north
      integer, dimension(mpi_status_size,8) :: status1

!---------------------------------------------------------------------

      index_east = -1
      index_west = -1
      index_south = -1
      index_north = -1

      nr = 0
      if(ibe.eq.0)then
        nr = nr + 1
        index_east = nr
      endif
      if(ibw.eq.0)then
        nr = nr + 1
        index_west = nr
      endif
      if(ibs.eq.0)then
        nr = nr + 1
        index_south = nr
      endif
      if(ibn.eq.0)then
        nr = nr + 1
        index_north = nr
      endif

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(1:nr),index,MPI_STATUS_IGNORE,ierr)
        nn = nn + 1
 
      if(index.eq.index_east)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          u(ni+2,j,k)=neweast(j,k)
        enddo
        enddo
      elseif(index.eq.index_west)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          u(0,j,k)=newwest(j,k)
        enddo
        enddo
      elseif(index.eq.index_south)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni+1
          u(i,0,k)=newsouth(i,k)
        enddo
        enddo
      elseif(index.eq.index_north)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni+1
          u(i,nj+1,k)=newnorth(i,k)
        enddo
        enddo
      endif

      enddo

      if(timestats.ge.1) time_mpu2=time_mpu2+mytime()

!----------
!  patch for corner
 
      if( (ebc.eq.2.or.wbc.eq.2).and.(sbc.eq.1.or.nbc.eq.1) )then
 
        if(ibw.eq.1)then
 
          if(p2tchsww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(0,0,k)=u(1,0,k)
            enddo
          endif
 
          if(p2tchnww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(0,nj+1,k)=u(1,nj+1,k)
            enddo
          endif
 
        endif
 
        if(ibe.eq.1)then
 
          if(p2tchsee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(ni+2,0,k)=u(ni+1,0,k)
            enddo
          endif
 
          if(p2tchnee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(ni+2,nj+1,k)=u(ni+1,nj+1,k)
            enddo
          endif
 
        endif
 
      endif

 
      if( (ebc.eq.1.or.wbc.eq.1).and.(sbc.eq.2.or.nbc.eq.2) )then
 
        if(ibs.eq.1)then
 
          if(p2tchsws)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(0,0,k)=u(0,1,k)
            enddo
          endif
 
          if(p2tchses)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(ni+2,0,k)=u(ni+2,1,k)
            enddo
          endif
 
        endif
 
        if(ibn.eq.1)then
 
          if(p2tchnwn)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(0,nj+1,k)=u(0,nj,k)
            enddo
          endif
 
          if(p2tchnen)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(ni+2,nj+1,k)=u(ni+2,nj,k)
            enddo
          endif
 
        endif
 
      endif
 
      if(timestats.ge.1) time_bc=time_bc+mytime()

!----------

      nr = 0

      if(ibw.eq.0)then
        nr = nr+1
      endif
      if(ibe.eq.0)then
        nr = nr+1
      endif
      if(ibn.eq.0)then
        nr = nr+1
      endif
      if(ibs.eq.0)then
        nr = nr+1
      endif

    if( nr.ge.1 )then
      call MPI_WAITALL(nr,reqs(5:5+nr-1),status1(1:mpi_status_size,5:5+nr-1),ierr)
    endif

      if(timestats.ge.1) time_mpu2=time_mpu2+mytime()
 
!----------
 
      end subroutine comm_1u_end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_1v_start(v,west,newwest,east,neweast,   &
                                 south,newsouth,north,newnorth,reqs)
      use input
      use mpi
      implicit none
 
      real, intent(in   ), dimension(ib:ie,jb:je+1,kb:ke) :: v
      real west(nj+1,nk),newwest(nj+1,nk)
      real east(nj+1,nk),neweast(nj+1,nk)
      real south(ni,nk),newsouth(ni,nk)
      real north(ni,nk),newnorth(ni,nk)
      integer reqs(8)
 
      integer i,j,k,nr
      integer tag1,tag2,tag3,tag4
 
!------------------------------------------------

      nr = 0

      nv=nv+1
      tag1=2000+nv

      ! receive east
      if(ibe.eq.0)then
        nr = nr + 1
        call mpi_irecv(neweast,cv1we,MPI_REAL,myeast,tag1,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      nv=nv+1
      tag2=2000+nv

      ! receive west
      if(ibw.eq.0)then
        nr = nr + 1
        call mpi_irecv(newwest,cv1we,MPI_REAL,mywest,tag2,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      nv=nv+1
      tag3=2000+nv

      ! receive south
      if(ibs.eq.0)then
        nr = nr + 1
        call mpi_irecv(newsouth,cs1sn,MPI_REAL,mysouth,tag3,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      nv=nv+1
      tag4=2000+nv

      ! receive north
      if(ibn.eq.0)then
        nr = nr + 1
        call mpi_irecv(newnorth,cs1sn,MPI_REAL,mynorth,tag4,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!------------------------------------------------
!------------------------------------------------
!------------------------------------------------

      nr = 4
 
      ! send west
      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj+1
          west(j,k)=v(1,j,k)
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(west,cv1we,MPI_REAL,mywest,tag1,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------
 
      ! send east
      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj+1
          east(j,k)=v(ni,j,k)
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(east,cv1we,MPI_REAL,myeast,tag2,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------
 
      ! send north
      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          north(i,k)=v(i,nj,k)
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(north,cs1sn,MPI_REAL,mynorth,tag3,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
!----------
 
      ! send south
      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          south(i,k)=v(i,2,k)
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(south,cs1sn,MPI_REAL,mysouth,tag4,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------
 
      if(timestats.ge.1) time_mps1=time_mps1+mytime()
 
      end subroutine comm_1v_start

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_1v_end(v,west,newwest,east,neweast,   &
                               south,newsouth,north,newnorth,reqs)
      use input
      use mpi
      implicit none
 
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: v
      real west(nj+1,nk),newwest(nj+1,nk)
      real east(nj+1,nk),neweast(nj+1,nk)
      real south(ni,nk),newsouth(ni,nk)
      real north(ni,nk),newnorth(ni,nk)
      integer reqs(8)
 
      integer i,j,k,nn,nr,index
      integer :: index_east,index_west,index_south,index_north
      integer, dimension(mpi_status_size,8) :: status1

!---------------------------------------------------------------------

      index_east = -1
      index_west = -1
      index_south = -1
      index_north = -1

      nr = 0
      if(ibe.eq.0)then
        nr = nr + 1
        index_east = nr
      endif
      if(ibw.eq.0)then
        nr = nr + 1
        index_west = nr
      endif
      if(ibs.eq.0)then
        nr = nr + 1
        index_south = nr
      endif
      if(ibn.eq.0)then
        nr = nr + 1
        index_north = nr
      endif

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(1:nr),index,MPI_STATUS_IGNORE,ierr)
        nn = nn + 1
 
      if(index.eq.index_east)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj+1
          v(ni+1,j,k)=neweast(j,k)
        enddo
        enddo
      elseif(index.eq.index_west)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj+1
          v(0,j,k)=newwest(j,k)
        enddo
        enddo
      elseif(index.eq.index_south)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          v(i,0,k)=newsouth(i,k)
        enddo
        enddo
      elseif(index.eq.index_north)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          v(i,nj+2,k)=newnorth(i,k)
        enddo
        enddo
      endif

      enddo

      if(timestats.ge.1) time_mps2=time_mps2+mytime()

!----------
!  patch for corner
 
      if( (ebc.eq.2.or.wbc.eq.2).and.(sbc.eq.1.or.nbc.eq.1) )then
 
        if(ibw.eq.1)then
 
          if(p2tchsww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(0,0,k)=v(1,0,k)
            enddo
          endif
 
          if(p2tchnww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(0,nj+2,k)=v(1,nj+2,k)
            enddo
          endif
 
        endif
 
        if(ibe.eq.1)then
 
          if(p2tchsee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(ni+1,0,k)=v(ni,0,k)
            enddo
          endif
 
          if(p2tchnee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(ni+1,nj+2,k)=v(ni,nj+2,k)
            enddo
          endif
 
        endif
 
      endif

      if( (ebc.eq.1.or.wbc.eq.1).and.(sbc.eq.2.or.nbc.eq.2) )then
 
        if(ibs.eq.1)then
 
          if(p2tchsws)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(0,0,k)=v(0,1,k)
            enddo
          endif
 
          if(p2tchses)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(ni+1,0,k)=v(ni+1,1,k)
            enddo
          endif
 
        endif
 
        if(ibn.eq.1)then
 
          if(p2tchnwn)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(0,nj+2,k)=v(0,nj+1,k)
            enddo
          endif
 
          if(p2tchnen)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(ni+1,nj+2,k)=v(ni+1,nj+1,k)
            enddo
          endif
 
        endif
 
      endif
 
      if(timestats.ge.1) time_bc=time_bc+mytime()

!----------

      nr = 0

      if(ibw.eq.0)then
        nr = nr+1
      endif
      if(ibe.eq.0)then
        nr = nr+1
      endif
      if(ibn.eq.0)then
        nr = nr+1
      endif
      if(ibs.eq.0)then
        nr = nr+1
      endif

    if( nr.ge.1 )then
      call MPI_WAITALL(nr,reqs(5:5+nr-1),status1(1:mpi_status_size,5:5+nr-1),ierr)
    endif

      if(timestats.ge.1) time_mps2=time_mps2+mytime()
 
!----------
 
      end subroutine comm_1v_end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_1w_start(w,ww1,ww2,we1,we2,   &
                                 ws1,ws2,wn1,wn2,reqs)
      use input
      use mpi
      implicit none
 
      real w(ib:ie,jb:je,kb:ke+1)
      real ww1(nj,nk-1),ww2(nj,nk-1)
      real we1(nj,nk-1),we2(nj,nk-1)
      real ws1(ni,nk-1),ws2(ni,nk-1)
      real wn1(ni,nk-1),wn2(ni,nk-1)
      integer reqs(8)
 
      integer i,j,k,nr
      integer tag1,tag2,tag3,tag4

      nr = 0

!-----

      nw=nw+1
      tag1=3000+nw
 
        ! receive east
        if(ibe.eq.0)then
          nr = nr+1
          call mpi_irecv(we2,cw1we,MPI_REAL,myeast,tag1,   &
                        MPI_COMM_WORLD,reqs(nr),ierr)
        endif

!-----

      nw=nw+1
      tag2=3000+nw
 
        ! receive west
        if(ibw.eq.0)then
          nr = nr+1
          call mpi_irecv(ww2,cw1we,MPI_REAL,mywest,tag2,   &
                        MPI_COMM_WORLD,reqs(nr),ierr)
        endif

!-----

      nw=nw+1
      tag3=3000+nw
 
        ! receive north
        if(ibn.eq.0)then
          nr = nr+1
          call mpi_irecv(wn2,cw1sn,MPI_REAL,mynorth,tag3,   &
                        MPI_COMM_WORLD,reqs(nr),ierr)
        endif

!-----

      nw=nw+1
      tag4=3000+nw
 
        ! receive south
        if(ibs.eq.0)then
          nr = nr+1
          call mpi_irecv(ws2,cw1sn,MPI_REAL,mysouth,tag4,   &
                        MPI_COMM_WORLD,reqs(nr),ierr)
        endif

!------------------------------------------------
!------------------------------------------------
!------------------------------------------------

      nr = 4

        ! send west
        if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=2,nk
          do j=1,nj
            ww1(j,k-1)=w(1,j,k)
          enddo
          enddo
          nr = nr+1
          call mpi_isend(ww1,cw1we,MPI_REAL,mywest,tag1,   &
                         MPI_COMM_WORLD,reqs(nr),ierr)
        endif

        ! send east
        if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=2,nk
          do j=1,nj
            we1(j,k-1)=w(ni,j,k)
          enddo
          enddo
          nr = nr+1
          call mpi_isend(we1,cw1we,MPI_REAL,myeast,tag2,   &
                         MPI_COMM_WORLD,reqs(nr),ierr)
        endif

        ! send south
        if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=2,nk
          do i=1,ni
            ws1(i,k-1)=w(i,1,k)
          enddo
          enddo
          nr = nr+1
          call mpi_isend(ws1,cw1sn,MPI_REAL,mysouth,tag3,   &
                         MPI_COMM_WORLD,reqs(nr),ierr)
        endif

        ! send north
        if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=2,nk
          do i=1,ni
            wn1(i,k-1)=w(i,nj,k)
          enddo
          enddo
          nr = nr+1
          call mpi_isend(wn1,cw1sn,MPI_REAL,mynorth,tag4,   &
                         MPI_COMM_WORLD,reqs(nr),ierr)
        endif

!-----

      if(timestats.ge.1) time_mpw1=time_mpw1+mytime()
 
      end subroutine comm_1w_start

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_1w_end(w,ww1,ww2,we1,we2,   &
                               ws1,ws2,wn1,wn2,reqs)
      use input
      use mpi
      implicit none
 
      real w(ib:ie,jb:je,kb:ke+1)
      real ww1(nj,nk-1),ww2(nj,nk-1)
      real we1(nj,nk-1),we2(nj,nk-1)
      real ws1(ni,nk-1),ws2(ni,nk-1)
      real wn1(ni,nk-1),wn2(ni,nk-1)
      integer reqs(8)
 
      integer i,j,k,nn,nr,index
      integer :: index_east,index_west,index_south,index_north
      integer, dimension(mpi_status_size,8) :: status1
 
!-----

      index_east = -1
      index_west = -1
      index_south = -1
      index_north = -1

      nr = 0
      if(ibe.eq.0)then
        nr = nr + 1
        index_east = nr
      endif
      if(ibw.eq.0)then
        nr = nr + 1
        index_west = nr
      endif
      if(ibn.eq.0)then
        nr = nr + 1
        index_north = nr
      endif
      if(ibs.eq.0)then
        nr = nr + 1
        index_south = nr
      endif

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(1:nr),index,MPI_STATUS_IGNORE,ierr)
        nn = nn + 1

      if(index.eq.index_east)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=2,nk
        do j=1,nj
          w(ni+1,j,k)=we2(j,k-1)
        enddo
        enddo
      elseif(index.eq.index_west)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=2,nk
        do j=1,nj
          w(0,j,k)=ww2(j,k-1)
        enddo
        enddo
      elseif(index.eq.index_north)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=2,nk
        do i=1,ni
          w(i,nj+1,k)=wn2(i,k-1)
        enddo
        enddo
      elseif(index.eq.index_south)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=2,nk
        do i=1,ni
          w(i,0,k)=ws2(i,k-1)
        enddo
        enddo
      endif

      enddo

!-----

      nr = 0

      if(ibw.eq.0)then
        nr = nr+1
      endif
      if(ibe.eq.0)then
        nr = nr+1
      endif
      if(ibs.eq.0)then
        nr = nr+1
      endif
      if(ibn.eq.0)then
        nr = nr+1
      endif

    if( nr.ge.1 )then
      call MPI_WAITALL(nr,reqs(5:5+nr-1),status1(1:mpi_status_size,5:5+nr-1),ierr)
    endif

      if(timestats.ge.1) time_mpw2=time_mpw2+mytime()

!----------
!  patch for corner
 
      if( (ebc.eq.2.or.wbc.eq.2).and.(sbc.eq.1.or.nbc.eq.1) )then
 
        if(ibw.eq.1)then
 
          if(p2tchsww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(0,0,k)=w(1,0,k)
            enddo
          endif
 
          if(p2tchnww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(0,nj+1,k)=w(1,nj+1,k)
            enddo
          endif
 
        endif
 
        if(ibe.eq.1)then
 
          if(p2tchsee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(ni+1,0,k)=w(ni,0,k)
            enddo
          endif
 
          if(p2tchnee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(ni+1,nj+1,k)=w(ni,nj+1,k)
            enddo
          endif
 
        endif
 
      endif

      if( (ebc.eq.1.or.wbc.eq.1).and.(sbc.eq.2.or.nbc.eq.2) )then
 
        if(ibs.eq.1)then
 
          if(p2tchsws)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(0,0,k)=w(0,1,k)
            enddo
          endif
 
          if(p2tchses)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(ni+1,0,k)=w(ni+1,1,k)
            enddo
          endif
 
        endif
 
        if(ibn.eq.1)then
 
          if(p2tchnwn)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(0,nj+1,k)=w(0,nj,k)
            enddo
          endif
 
          if(p2tchnen)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(ni+1,nj+1,k)=w(ni+1,nj,k)
            enddo
          endif
 
        endif
 
      endif

!-----------------------------------------------------------
!  Mirror b.c. patch
!
!      if(ibn.eq.1)then
!!$omp parallel do default(shared)   &
!!$omp private(i,k)
!        do k=2,nk
!        do i=0,ni+1
!          w(i,nj+1,k)=w(i,nj  ,k)
!          w(i,nj+2,k)=w(i,nj-1,k)
!          w(i,nj+3,k)=w(i,nj-2,k)
!        enddo
!        enddo
!      endif
!
!      if(ibw.eq.1)then
!!$omp parallel do default(shared)   &
!!$omp private(i,k)
!        do k=2,nk
!        do j=0,nj+1
!          w(-2,j,k)=w(3,j,k)
!          w(-1,j,k)=w(2,j,k)
!          w( 0,j,k)=w(1,j,k)
!        enddo
!        enddo
!      endif
!
!-----------------------------------------------------------
 
      if(timestats.ge.1) time_bc=time_bc+mytime()
 
!----------

      end subroutine comm_1w_end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_2d_start(s,west,newwest,east,neweast,   &
                                 south,newsouth,north,newnorth,reqs)
      use input
      use mpi
      implicit none

      real, dimension(ib:ie,jb:je) :: s
      real, dimension(cmp,nj) :: west,newwest,east,neweast
      real, dimension(ni,cmp) :: south,newsouth,north,newnorth
      integer reqs(8)

      integer i,j,nr
      integer tag1,tag2,tag3,tag4

!------------------------------------------------

      nr = 0

      nf=nf+1
      tag1=nf

      ! receive east
      if(ibe.eq.0)then
        nr = nr + 1
        call mpi_irecv(neweast,cmp*nj,MPI_REAL,myeast,tag1,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      nf=nf+1
      tag2=nf

      ! receive west
      if(ibw.eq.0)then
        nr = nr + 1
        call mpi_irecv(newwest,cmp*nj,MPI_REAL,mywest,tag2,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      nf=nf+1
      tag3=nf

      ! receive south
      if(ibs.eq.0)then
        nr = nr + 1
        call mpi_irecv(newsouth,ni*cmp,MPI_REAL,mysouth,tag3,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      nf=nf+1
      tag4=nf

      ! receive north
      if(ibn.eq.0)then
        nr = nr + 1
        call mpi_irecv(newnorth,ni*cmp,MPI_REAL,mynorth,tag4,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!------------------------------------------------
!------------------------------------------------
!------------------------------------------------

      nr = 4

      ! send west
      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,cmp
          west(i,j)=s(i,j)
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(west,cmp*nj,MPI_REAL,mywest,tag1,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      ! send east
      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,cmp
          east(i,j)=s(ni-cmp+i,j)
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(east,cmp*nj,MPI_REAL,myeast,tag2,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      ! send north
      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,cmp
        do i=1,ni
          north(i,j)=s(i,nj-cmp+j)
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(north,ni*cmp,MPI_REAL,mynorth,tag3,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      ! send south
      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,cmp
        do i=1,ni
          south(i,j)=s(i,j)
        enddo
        enddo
        nr = nr + 1
        call mpi_isend(south,ni*cmp,MPI_REAL,mysouth,tag4,   &
                       MPI_COMM_WORLD,reqs(nr),ierr)
      endif

!----------

      if(timestats.ge.1) time_mpq1=time_mpq1+mytime()

      end subroutine comm_2d_start


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_2dew_end(s,west,newwest,east,neweast,reqs)
      use input
      use mpi
      implicit none

      real, dimension(ib:ie,jb:je) :: s
      real, dimension(cmp,nj) :: west,newwest,east,neweast
      integer reqs(8)

      integer i,j,nn,nr,index

!-------------------------------------------------------------------

      nr = 0
      if(ibe.eq.0)then
        nr = nr + 1
      endif
      if(ibw.eq.0)then
        nr = nr + 1
      endif

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(1:nr),index,MPI_STATUS_IGNORE,ierr)
        nn = nn + 1
      enddo

      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,cmp
          s(ni+i,j)=neweast(i,j)
        enddo
        enddo
      endif

      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,cmp
          s(i-cmp,j)=newwest(i,j)
        enddo
        enddo
      endif

!----------

      if(timestats.ge.1) time_mpq2=time_mpq2+mytime()

      end subroutine comm_2dew_end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_2dns_end(s,south,newsouth,north,newnorth,reqs)
      use input
      use mpi
      implicit none

      real, dimension(ib:ie,jb:je) :: s
      real, dimension(ni,cmp) :: south,newsouth,north,newnorth
      integer reqs(8)

      integer i,j,nn,nr1,nr,index
      integer, dimension(mpi_status_size,8) :: status1

!-------------------------------------------------------------------

      nr1 = 0
      if(ibe.eq.0)then
        nr1 = nr1 + 1
      endif
      if(ibw.eq.0)then
        nr1 = nr1 + 1
      endif

      nr = 0
      if(ibs.eq.0)then
        nr = nr + 1
      endif
      if(ibn.eq.0)then
        nr = nr + 1
      endif

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(nr1+1:nr1+nr),index,MPI_STATUS_IGNORE,ierr)
        nn = nn + 1
      enddo

      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,cmp
        do i=1,ni
          s(i,j-cmp)=newsouth(i,j)
        enddo
        enddo
      endif

      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,cmp
        do i=1,ni
          s(i,nj+j)=newnorth(i,j)
        enddo
        enddo
      endif

!----------

      nr = 0

      if(ibw.eq.0)then
        nr = nr+1
      endif

      if(ibe.eq.0)then
        nr = nr+1
      endif

      if(ibn.eq.0)then
        nr = nr+1
      endif

      if(ibs.eq.0)then
        nr = nr+1
      endif

    if( nr.ge.1 )then
      call MPI_WAITALL(nr,reqs(5:5+nr-1),status1(1:mpi_status_size,5:5+nr-1),ierr)
    endif

!----------

      if(timestats.ge.1) time_mpq2=time_mpq2+mytime()

      end subroutine comm_2dns_end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_2d_corner(s)
      use input
      use mpi
      implicit none

      real, dimension(ib:ie,jb:je) :: s

      integer reqs(8)
      integer :: tag1,tag2,tag3,tag4,nr,nrb
      integer, dimension(mpi_status_size,8) :: status1

!------------------------------------------------

      nr = 0

!-----

      tag1=5061

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        nr = nr + 1
        call mpi_irecv(s(0,nj+1),1,MPI_REAL,mynw,tag1,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
      endif

!-----

      tag2=5062

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        nr = nr + 1
        call mpi_irecv(s(0,0),1,MPI_REAL,mysw,tag2,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
      endif

!-----

      tag3=5063

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        nr = nr + 1
        call mpi_irecv(s(ni+1,nj+1),1,MPI_REAL,myne,tag3,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
      endif

!-----

      tag4=5064

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        nr = nr + 1
        call mpi_irecv(s(ni+1,0),1,MPI_REAL,myse,tag4,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
      endif

!------------------------------------------------
!------------------------------------------------
!------------------------------------------------

      nrb = 4

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        nrb = nrb + 1
        call mpi_isend(s(ni,1),1,MPI_REAL,myse,tag1,MPI_COMM_WORLD,   &
                       reqs(nrb),ierr)
      endif

!-----

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        nrb = nrb + 1
        call mpi_isend(s(ni,nj),1,MPI_REAL,myne,tag2,MPI_COMM_WORLD,   &
                       reqs(nrb),ierr)
      endif

!-----

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        nrb = nrb + 1
        call mpi_isend(s(1,1),1,MPI_REAL,mysw,tag3,MPI_COMM_WORLD,   &
                       reqs(nrb),ierr)
      endif

!-----

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        nrb = nrb + 1
        call mpi_isend(s(1,nj),1,MPI_REAL,mynw,tag4,MPI_COMM_WORLD,   &
                       reqs(nrb),ierr)
      endif

!------------------------------------------------
!------------------------------------------------
!------------------------------------------------

      if( nr.ge.1 )then
        call MPI_WAITALL(nr,reqs(1:nr),status1(1:mpi_status_size,1:nr),ierr)
      endif

!-----

    nrb = nrb-4

    if( nrb.ge.1 )then
      call MPI_WAITALL(nrb,reqs(5:5+nrb-1),status1(1:mpi_status_size,5:5+nrb-1),ierr)
    endif

!-----

      if(timestats.ge.1) time_mptk1=time_mptk1+mytime()

      end subroutine comm_2d_corner

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine getcorner3_2d(s)
      use input
      use mpi
      implicit none
 
      real, intent(inout), dimension(ib:ie,jb:je) :: s
 
      real, dimension(cmp,cmp) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      integer :: i,j,nn,nr,index
      integer :: index_nw,index_sw,index_ne,index_se
      integer reqs(8)
      integer tag1,tag2,tag3,tag4

!-----

      nr = 0
      index_nw = -1
      index_sw = -1
      index_ne = -1
      index_se = -1

!------------------------------------------------------------------

      tag1=5001

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        nr = nr + 1
        call mpi_irecv(nw2,cmp*cmp,MPI_REAL,mynw,tag1,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
        index_nw = nr
      endif

!-----

      tag2=5002

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        nr = nr + 1
        call mpi_irecv(sw2,cmp*cmp,MPI_REAL,mysw,tag2,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
        index_sw = nr
      endif

!-----

      tag3=5003

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        nr = nr + 1
        call mpi_irecv(ne2,cmp*cmp,MPI_REAL,myne,tag3,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
        index_ne = nr
      endif

!-----

      tag4=5004

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        nr = nr + 1
        call mpi_irecv(se2,cmp*cmp,MPI_REAL,myse,tag4,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
        index_se = nr
      endif

!------------------------------------------------------------------

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        do j=1,cmp
        do i=1,cmp
          se1(i,j)=s(ni-cmp+i,j)
        enddo
        enddo
        call mpi_isend(se1,cmp*cmp,MPI_REAL,myse,tag1,MPI_COMM_WORLD,   &
                       reqs(5),ierr)
      endif
 
!-----

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        do j=1,cmp
        do i=1,cmp
          ne1(i,j)=s(ni-cmp+i,nj-cmp+j)
        enddo
        enddo
        call mpi_isend(ne1,cmp*cmp,MPI_REAL,myne,tag2,MPI_COMM_WORLD,   &
                       reqs(6),ierr)
      endif
 
!-----

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        do j=1,cmp
        do i=1,cmp
          sw1(i,j)=s(i,j)
        enddo
        enddo
        call mpi_isend(sw1,cmp*cmp,MPI_REAL,mysw,tag3,MPI_COMM_WORLD,   &
                       reqs(7),ierr)
      endif
 
!-----

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        do j=1,cmp
        do i=1,cmp
          nw1(i,j)=s(i,nj-cmp+j)
        enddo
        enddo
        call mpi_isend(nw1,cmp*cmp,MPI_REAL,mynw,tag4,MPI_COMM_WORLD,   &
                       reqs(8),ierr)
      endif
 
!-----

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(1:nr),index,MPI_STATUS_IGNORE,ierr)
        nn = nn + 1

      if(index.eq.index_nw)then
        do j=1,cmp
        do i=1,cmp
          s(-cmp+i,nj+j)=nw2(i,j)
        enddo
        enddo
      elseif(index.eq.index_sw)then
        do j=1,cmp
        do i=1,cmp
          s(-cmp+i,-cmp+j)=sw2(i,j)
        enddo
        enddo
      elseif(index.eq.index_ne)then
        do j=1,cmp
        do i=1,cmp
          s(ni+i,nj+j)=ne2(i,j)
        enddo
        enddo
      elseif(index.eq.index_se)then
        do j=1,cmp
        do i=1,cmp
          s(ni+i,-cmp+j)=se2(i,j)
        enddo
        enddo
      endif

      enddo

!-----

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(5),MPI_STATUS_IGNORE,ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(6),MPI_STATUS_IGNORE,ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(7),MPI_STATUS_IGNORE,ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(8),MPI_STATUS_IGNORE,ierr)
      endif

!-----

      if(timestats.ge.1) time_mptk1=time_mptk1+mytime()

      end subroutine getcorner3_2d


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_all_s(s,sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,  &
                              n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2,reqs_s)

      use input
      use bc_module
      implicit none

      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: s
      real, intent(inout), dimension(cmp,jmp,kmp)   :: sw31,sw32,se31,se32
      real, intent(inout), dimension(imp,cmp,kmp)   :: ss31,ss32,sn31,sn32
      real, intent(inout), dimension(cmp,cmp,kmt+1) :: n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2
      integer, intent(inout), dimension(rmp) :: reqs_s

      call comm_3s_start(s,sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,reqs_s)
      call comm_3s_end(  s,sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,reqs_s)
      call getcorner3(s,n3w1(1,1,1),n3w2(1,1,1),n3e1(1,1,1),n3e2(1,1,1),   &
                        s3w1(1,1,1),s3w2(1,1,1),s3e1(1,1,1),s3e2(1,1,1))
      call bcs2(s)

      end subroutine comm_all_s


#endif
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine prepcorners(s,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,  &
                               pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,comm)
      use input
      use bc_module
      implicit none

      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: s
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      real, intent(inout), dimension(jmp,kmp) :: pw1,pw2,pe1,pe2
      real, intent(inout), dimension(imp,kmp) :: ps1,ps2,pn1,pn2
      integer, intent(inout), dimension(rmp) :: reqs_p
      integer, intent(in) :: comm

      integer :: i,j

!--------------------------------------------!
!  This subroutine is ONLY for parcel_interp !
!--------------------------------------------!

      IF( comm.eq.1 )THEN
        call bcs(s)
      ENDIF
#ifdef MPI
      IF( comm.eq.1 )THEN
        call comm_1s_start(s,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p)
        call comm_1s_end(  s,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p)
      ENDIF
      call getcorner(s,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
      call bcs2(s)
#endif

      IF( bbc.eq.1 .or. bbc.eq.2 .or. bbc.eq.3 )THEN
        ! extrapolate:
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=0,nj+1
        do i=0,ni+1
          s(i,j,0) = cgs1*s(i,j,1)+cgs2*s(i,j,2)+cgs3*s(i,j,3)
        enddo
        enddo
      ENDIF

      IF( tbc.eq.1 .or. tbc.eq.2 )THEN
        ! extrapolate:
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=0,nj+1
        do i=0,ni+1
          s(i,j,nk+1) = cgt1*s(i,j,nk)+cgt2*s(i,j,nk-1)+cgt3*s(i,j,nk-2)
        enddo
        enddo
      ENDIF

      end subroutine prepcorners


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine prepcornert(t,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,  &
                               tkw1,tkw2,tke1,tke2,tks1,tks2,tkn1,tkn2,reqs_p,comm)
      use input
      use bc_module
      implicit none

      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: t
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      real, intent(inout), dimension(cmp,jmp,kmt) :: tkw1,tkw2,tke1,tke2
      real, intent(inout), dimension(imp,cmp,kmt) :: tks1,tks2,tkn1,tkn2
      integer, intent(inout), dimension(rmp) :: reqs_p
      integer, intent(in) :: comm

      integer :: i,j

!--------------------------------------------!
!  This subroutine is ONLY for parcel_interp !
!--------------------------------------------!

      IF( comm.eq.1 )THEN
        call bcw(t,0)
      ENDIF
#ifdef MPI
      IF( comm.eq.1 )THEN
        call comm_1t_start(t,tkw1(1,1,1),tkw2(1,1,1),tke1(1,1,1),tke2(1,1,1),tks1(1,1,1),tks2(1,1,1),tkn1(1,1,1),tkn2(1,1,1),reqs_p)
        call comm_1t_end(  t,tkw1(1,1,1),tkw2(1,1,1),tke1(1,1,1),tke2(1,1,1),tks1(1,1,1),tks2(1,1,1),tkn1(1,1,1),tkn2(1,1,1),reqs_p)
      ENDIF
      call getcornert(t,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      call bct2(t)
#endif

      end subroutine prepcornert


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine prepcorner2d(s,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p)
      use input
      use bc_module
      implicit none

      real, intent(inout), dimension(ib:ie,jb:je) :: s
      real, intent(inout), dimension(jmp,kmp) :: pw1,pw2,pe1,pe2
      real, intent(inout), dimension(imp,kmp) :: ps1,ps2,pn1,pn2
      integer, intent(inout), dimension(rmp) :: reqs_p

      integer :: i,j

      call bc2d(s)
#ifdef MPI
      call comm_1s2d_start(s,pw1(1,1),pw2(1,1),pe1(1,1),pe2(1,1),  &
                             ps1(1,1),ps2(1,1),pn1(1,1),pn2(1,1),reqs_p)
      call comm_1s2d_end(  s,pw1(1,1),pw2(1,1),pe1(1,1),pe2(1,1),  &
                             ps1(1,1),ps2(1,1),pn1(1,1),pn2(1,1),reqs_p)
      call bcs2_2d(s)
      call comm_2d_corner(s)
#endif

      end subroutine prepcorner2d


#ifdef MPI
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getcorneru3(u,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      use input
      use mpi
      implicit none

      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: u
      real, intent(inout), dimension(cmp,cmp,nk) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2

      integer :: i,j,k,nn,index
      integer :: index_nw,index_sw,index_ne,index_se
      integer reqs(8)
      integer tag1,tag2,tag3,tag4,count,nr1,nr2
      integer, dimension(mpi_status_size,8) :: status1

      count=cmp*cmp*nk
      nr1 = 0
      index_nw = -1
      index_sw = -1
      index_ne = -1
      index_se = -1

!-----

      tag1=5031

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        nr1 = nr1+1
        call mpi_irecv(nw2,count,MPI_REAL,mynw,tag1,MPI_COMM_WORLD,   &
                       reqs(nr1),ierr)
        index_nw = nr1
      endif

!-----

      tag2=5032

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        nr1 = nr1+1
        call mpi_irecv(sw2,count,MPI_REAL,mysw,tag2,MPI_COMM_WORLD,   &
                       reqs(nr1),ierr)
        index_sw = nr1
      endif

!-----

      tag3=5033

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        nr1 = nr1+1
        call mpi_irecv(ne2,count,MPI_REAL,myne,tag3,MPI_COMM_WORLD,   &
                       reqs(nr1),ierr)
        index_ne = nr1
      endif

!-----

      tag4=5034

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        nr1 = nr1+1
        call mpi_irecv(se2,count,MPI_REAL,myse,tag4,MPI_COMM_WORLD,   &
                       reqs(nr1),ierr)
        index_se = nr1
      endif

!------------------------------------------------
!------------------------------------------------
!------------------------------------------------

      nr2 = 0

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        do k=1,nk
        do j=1,cmp
        do i=1,cmp
          se1(i,j,k)=u(ni-cmp+i,j,k)
        enddo
        enddo
        enddo
        nr2 = nr2+1
        call mpi_isend(se1,count,MPI_REAL,myse,tag1,MPI_COMM_WORLD,   &
                       reqs(4+nr2),ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        do k=1,nk
        do j=1,cmp
        do i=1,cmp
          ne1(i,j,k)=u(ni-cmp+i,nj-cmp+j,k)
        enddo
        enddo
        enddo
        nr2 = nr2+1
        call mpi_isend(ne1,count,MPI_REAL,myne,tag2,MPI_COMM_WORLD,   &
                       reqs(4+nr2),ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        do k=1,nk
        do j=1,cmp
        do i=1,cmp
          sw1(i,j,k)=u(1+i,j,k)
        enddo
        enddo
        enddo
        nr2 = nr2+1
        call mpi_isend(sw1,count,MPI_REAL,mysw,tag3,MPI_COMM_WORLD,   &
                       reqs(4+nr2),ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        do k=1,nk
        do j=1,cmp
        do i=1,cmp
          nw1(i,j,k)=u(1+i,nj-cmp+j,k)
        enddo
        enddo
        enddo
        nr2 = nr2+1
        call mpi_isend(nw1,count,MPI_REAL,mynw,tag4,MPI_COMM_WORLD,   &
                       reqs(4+nr2),ierr)
      endif

!------------------------------------------------
!------------------------------------------------
!------------------------------------------------

      nn = 1
      do while( nn .le. nr1 )

        call MPI_WAITANY(nr1,reqs(1:nr1),index,MPI_STATUS_IGNORE,ierr)
        nn = nn + 1

      if( index.eq.index_nw )then
        do k=1,nk
        do j=1,cmp
        do i=1,cmp
          u(-cmp+i,nj+j,k)=nw2(i,j,k)
        enddo
        enddo
        enddo
      elseif( index.eq.index_sw )then
        do k=1,nk
        do j=1,cmp
        do i=1,cmp
          u(-cmp+i,-cmp+j,k)=sw2(i,j,k)
        enddo
        enddo
        enddo
      elseif( index.eq.index_ne )then
        do k=1,nk
        do j=1,cmp
        do i=1,cmp
          u(ni+1+i,nj+j,k)=ne2(i,j,k)
        enddo
        enddo
        enddo
      elseif( index.eq.index_se )then
        do k=1,nk
        do j=1,cmp
        do i=1,cmp
          u(ni+1+i,-cmp+j,k)=se2(i,j,k)
        enddo
        enddo
        enddo
      endif

      enddo

!------------------------------------------------
!------------------------------------------------
!------------------------------------------------

    if( nr2.ge.1 )then
      call MPI_WAITALL(nr2,reqs(5:5+nr2-1),status1(1:mpi_status_size,5:5+nr2-1),ierr)
    endif

!-----

      if(timestats.ge.1) time_mptk1=time_mptk1+mytime()

      end subroutine getcorneru3


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getcornerv3(v,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      use input
      use mpi
      implicit none

      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: v
      real, intent(inout), dimension(cmp,cmp,nk) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2

      integer :: i,j,k,nn,index
      integer :: index_nw,index_sw,index_ne,index_se
      integer reqs(8)
      integer tag1,tag2,tag3,tag4,count,nr1,nr2
      integer, dimension(mpi_status_size,8) :: status1

      count=cmp*cmp*nk
      nr1 = 0
      index_nw = -1
      index_sw = -1
      index_ne = -1
      index_se = -1

!-----

      tag1=5041

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        nr1 = nr1+1
        call mpi_irecv(nw2,count,MPI_REAL,mynw,tag1,MPI_COMM_WORLD,   &
                       reqs(nr1),ierr)
        index_nw = nr1
      endif

!-----

      tag2=5042

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        nr1 = nr1+1
        call mpi_irecv(sw2,count,MPI_REAL,mysw,tag2,MPI_COMM_WORLD,   &
                       reqs(nr1),ierr)
        index_sw = nr1
      endif

!-----

      tag3=5043

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        nr1 = nr1+1
        call mpi_irecv(ne2,count,MPI_REAL,myne,tag3,MPI_COMM_WORLD,   &
                       reqs(nr1),ierr)
        index_ne = nr1
      endif

!-----

      tag4=5044

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        nr1 = nr1+1
        call mpi_irecv(se2,count,MPI_REAL,myse,tag4,MPI_COMM_WORLD,   &
                       reqs(nr1),ierr)
        index_se = nr1
      endif

!------------------------------------------------
!------------------------------------------------
!------------------------------------------------

      nr2 = 0

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        do k=1,nk
        do j=1,cmp
        do i=1,cmp
          se1(i,j,k)=v(ni-cmp+i,1+j,k)
        enddo
        enddo
        enddo
        nr2 = nr2+1
        call mpi_isend(se1,count,MPI_REAL,myse,tag1,MPI_COMM_WORLD,   &
                       reqs(4+nr2),ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        do k=1,nk
        do j=1,cmp
        do i=1,cmp
          ne1(i,j,k)=v(ni-cmp+i,nj-cmp+j,k)
        enddo
        enddo
        enddo
        nr2 = nr2+1
        call mpi_isend(ne1,count,MPI_REAL,myne,tag2,MPI_COMM_WORLD,   &
                       reqs(4+nr2),ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        do k=1,nk
        do j=1,cmp
        do i=1,cmp
          sw1(i,j,k)=v(i,1+j,k)
        enddo
        enddo
        enddo
        nr2 = nr2+1
        call mpi_isend(sw1,count,MPI_REAL,mysw,tag3,MPI_COMM_WORLD,   &
                       reqs(4+nr2),ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        do k=1,nk
        do j=1,cmp
        do i=1,cmp
          nw1(i,j,k)=v(i,nj-cmp+j,k)
        enddo
        enddo
        enddo
        nr2 = nr2+1
        call mpi_isend(nw1,count,MPI_REAL,mynw,tag4,MPI_COMM_WORLD,   &
                       reqs(4+nr2),ierr)
      endif

!------------------------------------------------
!------------------------------------------------
!------------------------------------------------

      nn = 1
      do while( nn .le. nr1 )

        call MPI_WAITANY(nr1,reqs(1:nr1),index,MPI_STATUS_IGNORE,ierr)
        nn = nn + 1

      if( index.eq.index_nw )then
        do k=1,nk
        do j=1,cmp
        do i=1,cmp
          v(-cmp+i,nj+1+j,k)=nw2(i,j,k)
        enddo
        enddo
        enddo
      elseif( index.eq.index_sw )then
        do k=1,nk
        do j=1,cmp
        do i=1,cmp
          v(-cmp+i,-cmp+j,k)=sw2(i,j,k)
        enddo
        enddo
        enddo
      elseif( index.eq.index_ne )then
        do k=1,nk
        do j=1,cmp
        do i=1,cmp
          v(ni+i,nj+1+j,k)=ne2(i,j,k)
        enddo
        enddo
        enddo
      elseif( index.eq.index_se )then
        do k=1,nk
        do j=1,cmp
        do i=1,cmp
          v(ni+i,-cmp+j,k)=se2(i,j,k)
        enddo
        enddo
        enddo
      endif

      enddo

!------------------------------------------------
!------------------------------------------------
!------------------------------------------------

    if( nr2.ge.1 )then
      call MPI_WAITALL(nr2,reqs(5:5+nr2-1),status1(1:mpi_status_size,5:5+nr2-1),ierr)
    endif

!-----

      if(timestats.ge.1) time_mptk1=time_mptk1+mytime()

      end subroutine getcornerv3


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getcornerw3(w,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      use input
      use mpi
      implicit none

      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: w
      real, intent(inout), dimension(cmp,cmp,nk+1) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2

      integer :: i,j,k,nn,index
      integer :: index_nw,index_sw,index_ne,index_se
      integer reqs(8)
      integer tag1,tag2,tag3,tag4,count,nr1,nr2
      integer, dimension(mpi_status_size,8) :: status1

      count=cmp*cmp*nkp1
      nr1 = 0
      index_nw = -1
      index_sw = -1
      index_ne = -1
      index_se = -1

!-----

      tag1=5051

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        nr1 = nr1+1
        call mpi_irecv(nw2,count,MPI_REAL,mynw,tag1,MPI_COMM_WORLD,   &
                       reqs(nr1),ierr)
        index_nw = nr1
      endif

!-----

      tag2=5052

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        nr1 = nr1+1
        call mpi_irecv(sw2,count,MPI_REAL,mysw,tag2,MPI_COMM_WORLD,   &
                       reqs(nr1),ierr)
        index_sw = nr1
      endif

!-----

      tag3=5053

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        nr1 = nr1+1
        call mpi_irecv(ne2,count,MPI_REAL,myne,tag3,MPI_COMM_WORLD,   &
                       reqs(nr1),ierr)
        index_ne = nr1
      endif

!-----

      tag4=5054

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        nr1 = nr1+1
        call mpi_irecv(se2,count,MPI_REAL,myse,tag4,MPI_COMM_WORLD,   &
                       reqs(nr1),ierr)
        index_se = nr1
      endif

!------------------------------------------------
!------------------------------------------------
!------------------------------------------------

      nr2 = 0

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        do k=1,nkp1
        do j=1,cmp
        do i=1,cmp
          se1(i,j,k)=w(ni-cmp+i,j,k)
        enddo
        enddo
        enddo
        nr2 = nr2+1
        call mpi_isend(se1,count,MPI_REAL,myse,tag1,MPI_COMM_WORLD,   &
                       reqs(4+nr2),ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        do k=1,nkp1
        do j=1,cmp
        do i=1,cmp
          ne1(i,j,k)=w(ni-cmp+i,nj-cmp+j,k)
        enddo
        enddo
        enddo
        nr2 = nr2+1
        call mpi_isend(ne1,count,MPI_REAL,myne,tag2,MPI_COMM_WORLD,   &
                       reqs(4+nr2),ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        do k=1,nkp1
        do j=1,cmp
        do i=1,cmp
          sw1(i,j,k)=w(i,j,k)
        enddo
        enddo
        enddo
        nr2 = nr2+1
        call mpi_isend(sw1,count,MPI_REAL,mysw,tag3,MPI_COMM_WORLD,   &
                       reqs(4+nr2),ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        do k=1,nkp1
        do j=1,cmp
        do i=1,cmp
          nw1(i,j,k)=w(i,nj-cmp+j,k)
        enddo
        enddo
        enddo
        nr2 = nr2+1
        call mpi_isend(nw1,count,MPI_REAL,mynw,tag4,MPI_COMM_WORLD,   &
                       reqs(4+nr2),ierr)
      endif

!------------------------------------------------
!------------------------------------------------
!------------------------------------------------

      nn = 1
      do while( nn .le. nr1 )

        call MPI_WAITANY(nr1,reqs(1:nr1),index,MPI_STATUS_IGNORE,ierr)
        nn = nn + 1

      if( index.eq.index_nw )then
        do k=1,nkp1
        do j=1,cmp
        do i=1,cmp
          w(-cmp+i,nj+j,k)=nw2(i,j,k)
        enddo
        enddo
        enddo
      elseif( index.eq.index_sw )then
        do k=1,nkp1
        do j=1,cmp
        do i=1,cmp
          w(-cmp+i,-cmp+j,k)=sw2(i,j,k)
        enddo
        enddo
        enddo
      elseif( index.eq.index_ne )then
        do k=1,nkp1
        do j=1,cmp
        do i=1,cmp
          w(ni+i,nj+j,k)=ne2(i,j,k)
        enddo
        enddo
        enddo
      elseif( index.eq.index_se )then
        do k=1,nkp1
        do j=1,cmp
        do i=1,cmp
          w(ni+i,-cmp+j,k)=se2(i,j,k)
        enddo
        enddo
        enddo
      endif

      enddo

!------------------------------------------------
!------------------------------------------------
!------------------------------------------------

    if( nr2.ge.1 )then
      call MPI_WAITALL(nr2,reqs(5:5+nr2-1),status1(1:mpi_status_size,5:5+nr2-1),ierr)
    endif

!-----

      if(timestats.ge.1) time_mptk1=time_mptk1+mytime()

      end subroutine getcornerw3


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#endif

  END MODULE comm_module
