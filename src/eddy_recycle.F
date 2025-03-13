  module eddy_recycle
  implicit none

  public

    logical :: do_recycle
    real :: recy_width,recy_depth,recy_cap_w,recy_cap_s,recy_cap_e,recy_cap_n
    real :: recy_inj_w,recy_inj_s,recy_inj_e,recy_inj_n
    integer :: irecywe,jrecywe,irecysn,jrecysn,krecy
    integer :: nrecy,urecy,vrecy,wrecy,trecy,qrecy,erecy,xrecy
    real :: tscale0

    integer :: wirc1,wirc2,wirc
    integer :: wircs1,wircs2,wircs,wirb
    integer :: wuircs1,wuircs2,wuircs,wuirb
    integer :: wiri1,wiri2,wiri
    integer :: wiris1,wiris2,wiris,wirisb
    integer :: wuiris1,wuiris2,wuiris,wuirisb

    integer :: eirc1,eirc2,eirc
    integer :: eircs1,eircs2,eircs,eirb
    integer :: euircs1,euircs2,euircs,euirb
    integer :: eiri1,eiri2,eiri
    integer :: eiris1,eiris2,eiris,eirisb
    integer :: euiris1,euiris2,euiris,euirisb

    integer :: sjrc1,sjrc2,sjrc
    integer :: sjrcs1,sjrcs2,sjrcs,sjrb
    integer :: svjrcs1,svjrcs2,svjrcs,svjrb
    integer :: sjri1,sjri2,sjri
    integer :: sjris1,sjris2,sjris,sjrisb
    integer :: svjris1,svjris2,svjris,svjrisb

    integer :: njrc1,njrc2,njrc
    integer :: njrcs1,njrcs2,njrcs,njrb
    integer :: nvjrcs1,nvjrcs2,nvjrcs,nvjrb
    integer :: njri1,njri2,njri
    integer :: njris1,njris2,njris,njrisb
    integer :: nvjris1,nvjris2,nvjris,nvjrisb

    ! 200719: shift data a tad, repeat excessive recycling
    ! displace data?
    integer, parameter :: idisp = 0
    integer, parameter :: jdisp = 0

    integer, parameter :: tscalefac = 10.0

  CONTAINS

!-----------------------------------------------------------------------

    subroutine eddy_recycle_setup


    ! cm1r21.0:  Eddy recycling parameters are now specified in 
    !            namelist.input (see param18 section)


    end subroutine eddy_recycle_setup

!-----------------------------------------------------------------------
!   Eddy recycler on west side of domain:

    subroutine do_eddy_recyw(dt,xh,xf,yh,yf,zh,zf,u0,v0,dum1,dum2,cm0,  &
                             u3d,v3d,w3d,th3d,q3d,tke3d,                &
                             uten1,vten1,wten1,thten1,qten,tketen,      &
                             recywe,out3d,recy_cap,recy_inj,timavg,adtlast)

    use input
    use constants
#ifdef MPI
    use mpi
#endif

      implicit none

      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: xh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: yh
      real, intent(in), dimension(jb:je+1) :: yf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2
      real, intent(in), dimension(ib:ie,jb:je) :: cm0
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0,u3d
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0,v3d
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: w3d
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: th3d
      real, intent(in), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: q3d
      real, intent(in), dimension(ibt:iet,jbt:jet,kbt:ket) :: tke3d
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: uten1
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: vten1
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: wten1
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: thten1
      real, intent(inout), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qten
      real, intent(inout), dimension(ibt:iet,jbt:jet,kbt:ket) :: tketen
      real, intent(inout), dimension(irecywe,jrecywe,krecy,nrecy) :: recywe
      real, intent(inout) , dimension(ib3d:ie3d,jb3d:je3d,kb3d:ke3d,nout3d) :: out3d
      integer, intent(inout), dimension(ib:ie,jb:je) :: recy_cap,recy_inj
      real, intent(in), dimension(ibta:ieta,jbta:jeta,kbta:keta,ntavr) :: timavg
      double precision, intent(in) :: adtlast

      integer :: i,j,k,n,nn,ii,proc,ircst,irbt,irist,irisbt
      real :: tscale
#ifdef MPI
      integer, dimension(mpi_status_size) :: status
      integer :: ndat,irbp
#endif

      if( adapt_dt.eq.1 )then
!!!        tscale = max( tscale0 , 4.0*adtlast )
        tscale = tscalefac*adtlast
      else
        tscale = tscalefac*dt
      endif
      if( myid.eq.0 ) print *,'  adtlast,tscale = ',adtlast,tscale

    !---------------------------------------------------------------------------
    !  step 1: populate recywe arrays (data to be recycled)

#ifndef MPI

          ! single-proc, shared memory:

          do k=1,krecy
          do j=1,nj
          do i=wuircs1,wuircs2
            ii = i-wuircs1+1
            recywe(ii,j,k,urecy) = u3d(i,j+jdisp,k)-timavg(i,j+jdisp,k,utav)
          enddo
          enddo
          enddo

          do k=1,krecy
          do j=1,nj+1
          do i=wirc1,wirc2
            ii = i-wirc1+1
            recywe(ii,j,k,vrecy) = v3d(i,j+jdisp,k)-timavg(i,j+jdisp,k,vtav)
          enddo
          enddo
          enddo

          do k=1,krecy
          do j=1,nj
          do i=wirc1,wirc2
            ii = i-wirc1+1
            recywe(ii,j,k,wrecy) = w3d(i,j+jdisp,k)-timavg(i,j+jdisp,k,wtav)
            recywe(ii,j,k,trecy) = th3d(i,j+jdisp,k)-timavg(i,j+jdisp,k,ttav)
          enddo
          enddo
          enddo

          if( imoist.eq.1 )then
            do k=1,krecy
            do j=1,nj
            do i=wirc1,wirc2
              ii = i-wirc1+1
              recywe(ii,j,k,qrecy) = q3d(i,j+jdisp,k,nqv)-timavg(i,j+jdisp,k,qtav)
            enddo
            enddo
            enddo
          endif

          if( etav.ge.1 )then
            do k=2,krecy
            do j=1,nj
            do i=wirc1,wirc2
              ii = i-wirc1+1
              recywe(ii,j,k,erecy) = tke3d(i,j+jdisp,k)-timavg(i,j+jdisp,k,etav)
            enddo
            enddo
            enddo
          endif

#else
        ! MPI version:

        IF( myi.eq.1 )THEN

          if( wuirb.ge.1 )then
            do k=1,krecy
            do j=1,nj
            do i=wuircs1,wuircs2
              ii = wuirb + (i-wuircs1)
              recywe(ii,j,k,urecy) = u3d(i,j+jdisp,k)-timavg(i,j+jdisp,k,utav)
            enddo
            enddo
            enddo
          endif

          if( wirb.ge.1 )then
            do k=1,krecy
            do j=1,nj+1
            do i=wircs1,wircs2
              ii = wirb + (i-wircs1)
              recywe(ii,j,k,vrecy) = v3d(i,j+jdisp,k)-timavg(i,j+jdisp,k,vtav)
            enddo
            enddo
            enddo

            do k=1,krecy
            do j=1,nj
            do i=wircs1,wircs2
              ii = wirb + (i-wircs1)
              recywe(ii,j,k,wrecy) = w3d(i,j+jdisp,k)-timavg(i,j+jdisp,k,wtav)
              recywe(ii,j,k,trecy) = th3d(i,j+jdisp,k)-timavg(i,j+jdisp,k,ttav)
            enddo
            enddo
            enddo

          if( cm1setup.eq.4 )then
            do j=1,nj
            do i=wircs1,wircs2
              if( cm0(i,j).gt.cmemin ) recy_cap(i,j) = 1.0
            enddo
            enddo
          else
            do j=1,nj
            do i=wircs1,wircs2
              recy_cap(i,j) = 1.0
            enddo
            enddo
          endif

            if( imoist.eq.1 )then
              do k=1,krecy
              do j=1,nj
              do i=wircs1,wircs2
                ii = wirb + (i-wircs1)
                recywe(ii,j,k,qrecy) = q3d(i,j+jdisp,k,nqv)-timavg(i,j+jdisp,k,qtav)
              enddo
              enddo
              enddo
            endif

            if( etav.ge.1 )then
              do k=1,krecy
              do j=1,nj
              do i=wircs1,wircs2
                ii = wirb + (i-wircs1)
                recywe(ii,j,k,erecy) = tke3d(i,j+jdisp,k)-timavg(i,j+jdisp,k,etav)
              enddo
              enddo
              enddo
            endif

          endif

        ENDIF

    !---------------------------------------------------------------------------
    !  step 2: communicate data (MPI only)

        IF( myi.eq.1 )THEN

          ! left-most processor receives data, collects, then sends
          ! I am left-most proc in this row:
          ! loop through all procs in this row:
          ! receive data, add to recywe arrays:
          do n=2,nodex
            proc = (myj-1)*nodex + n - 1
            ! ircst = width of data to receive
            call mpi_recv(ircst,1,MPI_INTEGER,proc,8001,MPI_COMM_WORLD,status,ierr)
            if( ircst.gt.0 )then
              ! irbt = starting point of data in the recycle array
              call mpi_recv(irbt,1,MPI_INTEGER,proc,8002,MPI_COMM_WORLD,status,ierr)
              call recv_recywe(ircst,nj  ,krecy  ,recywe(1,1,1,urecy),irbt,dum1(ib,jb,kb),8011,proc)
            endif
            call mpi_recv(ircst,1,MPI_INTEGER,proc,8003,MPI_COMM_WORLD,status,ierr)
            if( ircst.gt.0 )then
              call mpi_recv(irbt,1,MPI_INTEGER,proc,8004,MPI_COMM_WORLD,status,ierr)
              call recv_recywe(ircst,nj+1,krecy  ,recywe(1,1,1,vrecy),irbt,dum1(ib,jb,kb),8012,proc)
              call recv_recywe(ircst,nj  ,krecy  ,recywe(1,1,1,wrecy),irbt,dum1(ib,jb,kb),8013,proc)
              call recv_recywe(ircst,nj  ,krecy  ,recywe(1,1,1,trecy),irbt,dum1(ib,jb,kb),8014,proc)
                if( etav.ge.1 )  &
              call recv_recywe(ircst,nj  ,krecy  ,recywe(1,1,1,erecy),irbt,dum1(ib,jb,kb),8015,proc)
                if( imoist.eq.1 )  &
              call recv_recywe(ircst,nj  ,krecy  ,recywe(1,1,1,qrecy),irbt,dum1(ib,jb,kb),8016,proc)
            endif
          enddo

        ELSE

          ! I am NOT left-most proc in this row:
          ! send my data to left-most proc (if applicable):
          proc = (myj-1)*nodex
          ! ircs = width of data to send
          call mpi_send(wuircs,1,MPI_INTEGER,proc,8001,MPI_COMM_WORLD,ierr)
          if( wuircs1.gt.0 .and. wuircs2.gt.0 )then
            ! irb = starting point of data in the recycle array
            call mpi_send(wuirb,1,MPI_INTEGER,proc,8002,MPI_COMM_WORLD,ierr)
            call send_recywe(wuircs,nj  ,krecy  ,recywe(1,1,1,urecy),1 ,0 ,0 ,u3d,  timavg(ibta,jbta,kbta,utav),dum1(ib,jb,kb),proc,8011,recy_cap,cm0,wuircs1)
          endif
          call mpi_send(wircs,1,MPI_INTEGER,proc,8003,MPI_COMM_WORLD,ierr)
          if( wircs1.gt.0 .and. wircs2.gt.0 )then
            call mpi_send(wirb,1,MPI_INTEGER,proc,8004,MPI_COMM_WORLD,ierr)
            call send_recywe(wircs,nj+1,krecy  ,recywe(1,1,1,vrecy),0 ,1 ,0 ,v3d,  timavg(ibta,jbta,kbta,vtav),dum1(ib,jb,kb),proc,8012,recy_cap,cm0,wircs1)
            call send_recywe(wircs,nj  ,krecy  ,recywe(1,1,1,wrecy),0 ,0 ,1 ,w3d,  timavg(ibta,jbta,kbta,wtav),dum1(ib,jb,kb),proc,8013,recy_cap,cm0,wircs1)
            call send_recywe(wircs,nj  ,krecy  ,recywe(1,1,1,trecy),0 ,0 ,0 ,th3d, timavg(ibta,jbta,kbta,ttav),dum1(ib,jb,kb),proc,8014,recy_cap,cm0,wircs1)
              if( etav.ge.1 )  &
            call send_recywe(wircs,nj  ,krecy  ,recywe(1,1,1,erecy),0 ,0 ,1 ,tke3d,timavg(ibta,jbta,kbta,etav),dum1(ib,jb,kb),proc,8015,recy_cap,cm0,wircs1)
              if( imoist.eq.1 )  &
            call send_recywe(wircs,nj  ,krecy  ,recywe(1,1,1,qrecy),0 ,0 ,0 ,q3d(ib,jb,kb,nqv),timavg(ibta,jbta,kbta,qtav),dum1(ib,jb,kb),proc,8016,recy_cap,cm0,wircs1)
          endif

        ENDIF

      !-------------------------------------------------------------------------

        IF( myi.eq.1 )THEN

          ! loop through all procs in this row:
          ! send data to procs in injection region:
          do n=2,nodex
            proc = (myj-1)*nodex + n - 1
            ! irist = width of data to receive
            call mpi_recv(irist,1,MPI_INTEGER,proc,8021,MPI_COMM_WORLD,status,ierr)
            if( irist.gt.0 )then
              ! irisbt = starting point of data in the recycle array
              call mpi_recv(irisbt,1,MPI_INTEGER,proc,8022,MPI_COMM_WORLD,status,ierr)
              call send_injwe(irist,nj  ,krecy  ,recywe(1,1,1,urecy),irisbt,dum1(ib,jb,kb),8031,proc)
            endif 
            call mpi_recv(irist,1,MPI_INTEGER,proc,8023,MPI_COMM_WORLD,status,ierr)
            if( irist.gt.0 )then
              call mpi_recv(irisbt,1,MPI_INTEGER,proc,8024,MPI_COMM_WORLD,status,ierr)
              call send_injwe(irist,nj+1,krecy  ,recywe(1,1,1,vrecy),irisbt,dum1(ib,jb,kb),8032,proc)
              call send_injwe(irist,nj  ,krecy  ,recywe(1,1,1,wrecy),irisbt,dum1(ib,jb,kb),8033,proc)
              call send_injwe(irist,nj  ,krecy  ,recywe(1,1,1,trecy),irisbt,dum1(ib,jb,kb),8034,proc)
                if( etav.ge.1 )  &
              call send_injwe(irist,nj  ,krecy  ,recywe(1,1,1,erecy),irisbt,dum1(ib,jb,kb),8035,proc)
                if( imoist.eq.1 )  &
              call send_injwe(irist,nj  ,krecy  ,recywe(1,1,1,qrecy),irisbt,dum1(ib,jb,kb),8036,proc)
            endif
          enddo


        ELSE

          ! receive from left-most proc (if applicable):
          proc = (myj-1)*nodex
          ! iris = width of data to send
          call mpi_send(wuiris,1,MPI_INTEGER,proc,8021,MPI_COMM_WORLD,ierr)
          if( wuiris1.gt.0 .and. wuiris2.gt.0 )then
            ! irisb = starting point of data in the recycle array
            call mpi_send(wuirisb,1,MPI_INTEGER,proc,8022,MPI_COMM_WORLD,ierr)
            call recv_injwe(wuiris,nj  ,krecy  ,recywe(1,1,1,urecy),proc,8031,dum1(ib,jb,kb),wuiris1,wuirisb)
          endif
          call mpi_send(wiris,1,MPI_INTEGER,proc,8023,MPI_COMM_WORLD,ierr)
          if( wiris1.gt.0 .and. wiris2.gt.0 )then
            call mpi_send(wirisb,1,MPI_INTEGER,proc,8024,MPI_COMM_WORLD,ierr)
            call recv_injwe(wiris,nj+1,krecy  ,recywe(1,1,1,vrecy),proc,8032,dum1(ib,jb,kb),wiris1,wirisb)
            call recv_injwe(wiris,nj  ,krecy  ,recywe(1,1,1,wrecy),proc,8033,dum1(ib,jb,kb),wiris1,wirisb)
            call recv_injwe(wiris,nj  ,krecy  ,recywe(1,1,1,trecy),proc,8034,dum1(ib,jb,kb),wiris1,wirisb)
              if( etav.ge.1 )  &
            call recv_injwe(wiris,nj  ,krecy  ,recywe(1,1,1,erecy),proc,8035,dum1(ib,jb,kb),wiris1,wirisb)
              if( imoist.eq.1 )  &
            call recv_injwe(wiris,nj  ,krecy  ,recywe(1,1,1,qrecy),proc,8036,dum1(ib,jb,kb),wiris1,wirisb)
          endif

        ENDIF

      !-------------------------------------------------------------------------

#endif

      ! get recycle tendencies:

      call     recy_tendencywe(recywe(1,1,1,urecy),1 ,0 ,0 ,u3d  ,timavg(ibta,jbta,kbta,utav),uten1 ,recy_inj,cm0,out3d, 1, 7,13,dt,adtlast,zh,wuiris1,wuiris2,wuirisb)
      call     recy_tendencywe(recywe(1,1,1,vrecy),0 ,1 ,0 ,v3d  ,timavg(ibta,jbta,kbta,vtav),vten1 ,recy_inj,cm0,out3d, 2, 8,14,dt,adtlast,zh,wiris1,wiris2,wirisb)
      call     recy_tendencywe(recywe(1,1,1,wrecy),0 ,0 ,1 ,w3d  ,timavg(ibta,jbta,kbta,wtav),wten1 ,recy_inj,cm0,out3d, 3, 9,15,dt,adtlast,zf,wiris1,wiris2,wirisb)
      call     recy_tendencywe(recywe(1,1,1,trecy),0 ,0 ,0 ,th3d ,timavg(ibta,jbta,kbta,ttav),thten1,recy_inj,cm0,out3d, 4,10,16,dt,adtlast,zh,wiris1,wiris2,wirisb)
      if( etav.ge.1 )  &
      call     recy_tendencywe(recywe(1,1,1,erecy),0 ,0 ,1 ,tke3d,timavg(ibta,jbta,kbta,etav),tketen,recy_inj,cm0,out3d, 5,11,17,dt,adtlast,zf,wiris1,wiris2,wirisb)
      if( imoist.eq.1 )  &
      call     recy_tendencywe(recywe(1,1,1,qrecy),0 ,0 ,0 ,q3d(ib,jb,kb,nqv),timavg(ibta,jbta,kbta,qtav),qten(ib,jb,kb,nqv),recy_inj,cm0,out3d, 6,12,18,dt,adtlast,zh,wiris1,wiris2,wirisb)

    end subroutine do_eddy_recyw

!-----------------------------------------------------------------------
!   Eddy recycler on east side of domain:

    subroutine do_eddy_recye(dt,xh,xf,yh,yf,zh,zf,u0,v0,dum1,dum2,cm0,  &
                             u3d,v3d,w3d,th3d,q3d,tke3d,                &
                             uten1,vten1,wten1,thten1,qten,tketen,      &
                             recywe,out3d,recy_cap,recy_inj,timavg,adtlast)

    use input
    use constants
#ifdef MPI
    use mpi
#endif

      implicit none

      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: xh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: yh
      real, intent(in), dimension(jb:je+1) :: yf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2
      real, intent(in), dimension(ib:ie,jb:je) :: cm0
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0,u3d
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0,v3d
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: w3d
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: th3d
      real, intent(in), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: q3d
      real, intent(in), dimension(ibt:iet,jbt:jet,kbt:ket) :: tke3d
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: uten1
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: vten1
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: wten1
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: thten1
      real, intent(inout), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qten
      real, intent(inout), dimension(ibt:iet,jbt:jet,kbt:ket) :: tketen
      real, intent(inout), dimension(irecywe,jrecywe,krecy,nrecy) :: recywe
      real, intent(inout) , dimension(ib3d:ie3d,jb3d:je3d,kb3d:ke3d,nout3d) :: out3d
      integer, intent(inout), dimension(ib:ie,jb:je) :: recy_cap,recy_inj
      real, intent(in), dimension(ibta:ieta,jbta:jeta,kbta:keta,ntavr) :: timavg
      double precision, intent(in) :: adtlast

      integer :: i,j,k,n,nn,ii,proc,ircst,irbt,irist,irisbt
#ifdef MPI
      integer, dimension(mpi_status_size) :: status
      integer :: ndat,irbp
#endif


    !---------------------------------------------------------------------------
    !  step 1: populate recywe arrays (data to be recycled)

#ifndef MPI

          ! single-proc, shared memory:

          do k=1,krecy
          do j=1,nj
          do i=euircs1,euircs2
            ii = i-euircs1+1
            recywe(ii,j,k,urecy) = u3d(i,j+jdisp,k)-timavg(i,j+jdisp,k,utav)
          enddo
          enddo
          enddo

          do k=1,krecy
          do j=1,nj+1
          do i=eirc1,eirc2
            ii = i-eirc1+1
            recywe(ii,j,k,vrecy) = v3d(i,j+jdisp,k)-timavg(i,j+jdisp,k,vtav)
          enddo
          enddo
          enddo

          do k=1,krecy
          do j=1,nj
          do i=eirc1,eirc2
            ii = i-eirc1+1
            recywe(ii,j,k,wrecy) = w3d(i,j+jdisp,k)-timavg(i,j+jdisp,k,wtav)
            recywe(ii,j,k,trecy) = th3d(i,j+jdisp,k)-timavg(i,j+jdisp,k,ttav)
          enddo
          enddo
          enddo

          if( imoist.eq.1 )then
            do k=1,krecy
            do j=1,nj
            do i=eirc1,eirc2
              ii = i-eirc1+1
              recywe(ii,j,k,qrecy) = q3d(i,j+jdisp,k,nqv)-timavg(i,j+jdisp,k,qtav)
            enddo
            enddo
            enddo
          endif

          if( etav.ge.1 )then
            do k=2,krecy
            do j=1,nj
            do i=eirc1,eirc2
              ii = i-eirc1+1
              recywe(ii,j,k,erecy) = tke3d(i,j+jdisp,k)-timavg(i,j+jdisp,k,etav)
            enddo
            enddo
            enddo
          endif

#else
        ! MPI version:

        IF( myi.eq.1 )THEN

          if( euirb.ge.1 )then
            do k=1,krecy
            do j=1,nj
            do i=euircs1,euircs2
              ii = euirb + (i-euircs1)
              recywe(ii,j,k,urecy) = u3d(i,j+jdisp,k)-timavg(i,j+jdisp,k,utav)
            enddo
            enddo
            enddo
          endif

          if( eirb.ge.1 )then
            do k=1,krecy
            do j=1,nj+1
            do i=eircs1,eircs2
              ii = eirb + (i-eircs1)
              recywe(ii,j,k,vrecy) = v3d(i,j+jdisp,k)-timavg(i,j+jdisp,k,vtav)
            enddo
            enddo
            enddo

            do k=1,krecy
            do j=1,nj
            do i=eircs1,eircs2
              ii = eirb + (i-eircs1)
              recywe(ii,j,k,wrecy) = w3d(i,j+jdisp,k)-timavg(i,j+jdisp,k,wtav)
              recywe(ii,j,k,trecy) = th3d(i,j+jdisp,k)-timavg(i,j+jdisp,k,ttav)
            enddo
            enddo
            enddo

          if( cm1setup.eq.4 )then
            do j=1,nj
            do i=eircs1,eircs2
              if( cm0(i,j).gt.cmemin ) recy_cap(i,j) = 1.0
            enddo
            enddo
          else
            do j=1,nj
            do i=eircs1,eircs2
              recy_cap(i,j) = 1.0
            enddo
            enddo
          endif

            if( imoist.eq.1 )then
              do k=1,krecy
              do j=1,nj
              do i=eircs1,eircs2
                ii = eirb + (i-eircs1)
                recywe(ii,j,k,qrecy) = q3d(i,j+jdisp,k,nqv)-timavg(i,j+jdisp,k,qtav)
              enddo
              enddo
              enddo
            endif

            if( etav.ge.1 )then
              do k=1,krecy
              do j=1,nj
              do i=eircs1,eircs2
                ii = eirb + (i-eircs1)
                recywe(ii,j,k,erecy) = tke3d(i,j+jdisp,k)-timavg(i,j+jdisp,k,etav)
              enddo
              enddo
              enddo
            endif

          endif

        ENDIF

    !---------------------------------------------------------------------------
    !  step 2: communicate data (MPI only)

        IF( myi.eq.1 )THEN

          ! left-most processor receives data, collects, then sends
          ! I am left-most proc in this row:
          ! loop through all procs in this row:
          ! receive data, add to recywe arrays:
          do n=2,nodex
            proc = (myj-1)*nodex + n - 1

            call mpi_recv(ircst,1,MPI_INTEGER,proc,8101,MPI_COMM_WORLD,status,ierr)
            if( ircst.gt.0 )then
              ! irbt = starting point of data in the recycle array
              call mpi_recv(irbt,1,MPI_INTEGER,proc,8102,MPI_COMM_WORLD,status,ierr)
              call recv_recywe(ircst,nj  ,krecy  ,recywe(1,1,1,urecy),irbt,dum1(ib,jb,kb),8111,proc)
            endif
            call mpi_recv(ircst,1,MPI_INTEGER,proc,8103,MPI_COMM_WORLD,status,ierr)
            if( ircst.gt.0 )then
              call mpi_recv(irbt,1,MPI_INTEGER,proc,8104,MPI_COMM_WORLD,status,ierr)
              call recv_recywe(ircst,nj+1,krecy  ,recywe(1,1,1,vrecy),irbt,dum1(ib,jb,kb),8112,proc)
              call recv_recywe(ircst,nj  ,krecy  ,recywe(1,1,1,wrecy),irbt,dum1(ib,jb,kb),8113,proc)
              call recv_recywe(ircst,nj  ,krecy  ,recywe(1,1,1,trecy),irbt,dum1(ib,jb,kb),8114,proc)
                if( etav.ge.1 )  &
              call recv_recywe(ircst,nj  ,krecy  ,recywe(1,1,1,erecy),irbt,dum1(ib,jb,kb),8115,proc)
                if( imoist.eq.1 )  &
              call recv_recywe(ircst,nj  ,krecy  ,recywe(1,1,1,qrecy),irbt,dum1(ib,jb,kb),8116,proc)
            endif
          enddo

        ELSE

          ! I am NOT left-most proc in this row:
          ! send my data to left-most proc (if applicable):
          proc = (myj-1)*nodex

          call mpi_send(euircs,1,MPI_INTEGER,proc,8101,MPI_COMM_WORLD,ierr)
          if( euircs1.gt.0 .and. euircs2.gt.0 )then
            ! irb = starting point of data in the recycle array
            call mpi_send(euirb,1,MPI_INTEGER,proc,8102,MPI_COMM_WORLD,ierr)
            call send_recywe(euircs,nj  ,krecy  ,recywe(1,1,1,urecy),1 ,0 ,0 ,u3d,  timavg(ibta,jbta,kbta,utav),dum1(ib,jb,kb),proc,8111,recy_cap,cm0,euircs1)
          endif
          call mpi_send(eircs,1,MPI_INTEGER,proc,8103,MPI_COMM_WORLD,ierr)
          if( eircs1.gt.0 .and. eircs2.gt.0 )then
            call mpi_send(eirb,1,MPI_INTEGER,proc,8104,MPI_COMM_WORLD,ierr)
            call send_recywe(eircs,nj+1,krecy  ,recywe(1,1,1,vrecy),0 ,1 ,0 ,v3d,  timavg(ibta,jbta,kbta,vtav),dum1(ib,jb,kb),proc,8112,recy_cap,cm0,eircs1)
            call send_recywe(eircs,nj  ,krecy  ,recywe(1,1,1,wrecy),0 ,0 ,1 ,w3d,  timavg(ibta,jbta,kbta,wtav),dum1(ib,jb,kb),proc,8113,recy_cap,cm0,eircs1)
            call send_recywe(eircs,nj  ,krecy  ,recywe(1,1,1,trecy),0 ,0 ,0 ,th3d, timavg(ibta,jbta,kbta,ttav),dum1(ib,jb,kb),proc,8114,recy_cap,cm0,eircs1)
              if( etav.ge.1 )  &
            call send_recywe(eircs,nj  ,krecy  ,recywe(1,1,1,erecy),0 ,0 ,1 ,tke3d,timavg(ibta,jbta,kbta,etav),dum1(ib,jb,kb),proc,8115,recy_cap,cm0,eircs1)
              if( imoist.eq.1 )  &
            call send_recywe(eircs,nj  ,krecy  ,recywe(1,1,1,qrecy),0 ,0 ,0 ,q3d(ib,jb,kb,nqv),timavg(ibta,jbta,kbta,qtav),dum1(ib,jb,kb),proc,8116,recy_cap,cm0,eircs1)
          endif

        ENDIF

      !-------------------------------------------------------------------------

        IF( myi.eq.1 )THEN

          ! loop through all procs in this row:
          ! send data to procs in injection region:
          do n=2,nodex
            proc = (myj-1)*nodex + n - 1
            ! irist = width of data to receive
            call mpi_recv(irist,1,MPI_INTEGER,proc,8121,MPI_COMM_WORLD,status,ierr)
            if( irist.gt.0 )then
              ! irisbt = starting point of data in the recycle array
              call mpi_recv(irisbt,1,MPI_INTEGER,proc,8122,MPI_COMM_WORLD,status,ierr)
              call send_injwe(irist,nj  ,krecy  ,recywe(1,1,1,urecy),irisbt,dum1(ib,jb,kb),8131,proc)
            endif 
            call mpi_recv(irist,1,MPI_INTEGER,proc,8123,MPI_COMM_WORLD,status,ierr)
            if( irist.gt.0 )then
              call mpi_recv(irisbt,1,MPI_INTEGER,proc,8124,MPI_COMM_WORLD,status,ierr)
              call send_injwe(irist,nj+1,krecy  ,recywe(1,1,1,vrecy),irisbt,dum1(ib,jb,kb),8132,proc)
              call send_injwe(irist,nj  ,krecy  ,recywe(1,1,1,wrecy),irisbt,dum1(ib,jb,kb),8133,proc)
              call send_injwe(irist,nj  ,krecy  ,recywe(1,1,1,trecy),irisbt,dum1(ib,jb,kb),8134,proc)
                if( etav.ge.1 )  &
              call send_injwe(irist,nj  ,krecy  ,recywe(1,1,1,erecy),irisbt,dum1(ib,jb,kb),8135,proc)
                if( imoist.eq.1 )  &
              call send_injwe(irist,nj  ,krecy  ,recywe(1,1,1,qrecy),irisbt,dum1(ib,jb,kb),8136,proc)
            endif
          enddo


        ELSE

          ! receive from left-most proc (if applicable):
          proc = (myj-1)*nodex
          ! iris = width of data to send
          call mpi_send(euiris,1,MPI_INTEGER,proc,8121,MPI_COMM_WORLD,ierr)
          if( euiris1.gt.0 .and. euiris2.gt.0 )then
            ! irisb = starting point of data in the recycle array
            call mpi_send(euirisb,1,MPI_INTEGER,proc,8122,MPI_COMM_WORLD,ierr)
            call recv_injwe(euiris,nj  ,krecy  ,recywe(1,1,1,urecy),proc,8131,dum1(ib,jb,kb),euiris1,euirisb)
          endif
          call mpi_send(eiris,1,MPI_INTEGER,proc,8123,MPI_COMM_WORLD,ierr)
          if( eiris1.gt.0 .and. eiris2.gt.0 )then
            call mpi_send(eirisb,1,MPI_INTEGER,proc,8124,MPI_COMM_WORLD,ierr)
            call recv_injwe(eiris,nj+1,krecy  ,recywe(1,1,1,vrecy),proc,8132,dum1(ib,jb,kb),eiris1,eirisb)
            call recv_injwe(eiris,nj  ,krecy  ,recywe(1,1,1,wrecy),proc,8133,dum1(ib,jb,kb),eiris1,eirisb)
            call recv_injwe(eiris,nj  ,krecy  ,recywe(1,1,1,trecy),proc,8134,dum1(ib,jb,kb),eiris1,eirisb)
              if( etav.ge.1 )  &
            call recv_injwe(eiris,nj  ,krecy  ,recywe(1,1,1,erecy),proc,8135,dum1(ib,jb,kb),eiris1,eirisb)
              if( imoist.eq.1 )  &
            call recv_injwe(eiris,nj  ,krecy  ,recywe(1,1,1,qrecy),proc,8136,dum1(ib,jb,kb),eiris1,eirisb)
          endif

        ENDIF

      !-------------------------------------------------------------------------

#endif

      ! get recycle tendencies:

      call     recy_tendencywe(recywe(1,1,1,urecy),1 ,0 ,0 ,u3d  ,timavg(ibta,jbta,kbta,utav),uten1 ,recy_inj,cm0,out3d, 1, 7,13,dt,adtlast,zh,euiris1,euiris2,euirisb)
      call     recy_tendencywe(recywe(1,1,1,vrecy),0 ,1 ,0 ,v3d  ,timavg(ibta,jbta,kbta,vtav),vten1 ,recy_inj,cm0,out3d, 2, 8,14,dt,adtlast,zh,eiris1,eiris2,eirisb)
      call     recy_tendencywe(recywe(1,1,1,wrecy),0 ,0 ,1 ,w3d  ,timavg(ibta,jbta,kbta,wtav),wten1 ,recy_inj,cm0,out3d, 3, 9,15,dt,adtlast,zf,eiris1,eiris2,eirisb)
      call     recy_tendencywe(recywe(1,1,1,trecy),0 ,0 ,0 ,th3d ,timavg(ibta,jbta,kbta,ttav),thten1,recy_inj,cm0,out3d, 4,10,16,dt,adtlast,zh,eiris1,eiris2,eirisb)
      if( etav.ge.1 )  &
      call     recy_tendencywe(recywe(1,1,1,erecy),0 ,0 ,1 ,tke3d,timavg(ibta,jbta,kbta,etav),tketen,recy_inj,cm0,out3d, 5,11,17,dt,adtlast,zf,eiris1,eiris2,eirisb)
      if( imoist.eq.1 )  &
      call     recy_tendencywe(recywe(1,1,1,qrecy),0 ,0 ,0 ,q3d(ib,jb,kb,nqv),timavg(ibta,jbta,kbta,qtav),qten(ib,jb,kb,nqv),recy_inj,cm0,out3d, 6,12,18,dt,adtlast,zh,eiris1,eiris2,eirisb)

    end subroutine do_eddy_recye

!-----------------------------------------------------------------------

    subroutine recy_tendencywe(recywe ,is,js,ks,var  ,timavg,varten,recy_inj,cm0,out3d,o1,o2,o3,dt,adtlast,zz,tiris1,tiris2,tirisb)
    use input
    use constants, only : cmemin
    implicit none

    real, intent(in), dimension(irecywe,jrecywe,krecy) :: recywe
    integer, intent(in) :: is,js,ks,o1,o2,o3,tiris1,tiris2,tirisb
    real, intent(in   ), dimension(ib:ie+is,jb:je+js,kb:ke+ks) :: var
    real, intent(inout), dimension(ib:ie+is,jb:je+js,kb:ke+ks) :: varten
    real, intent(in),    dimension(ibta:ieta,jbta:jeta,kbta:keta) :: timavg
    integer, intent(inout), dimension(ib:ie,jb:je) :: recy_inj
    real, intent(in),    dimension(ib:ie,jb:je) :: cm0
    real, intent(inout) , dimension(ib3d:ie3d,jb3d:je3d,kb3d:ke3d,nout3d) :: out3d
    real, intent(in) :: dt
    double precision, intent(in) :: adtlast
    real, intent(in), dimension(ib:ie,jb:je,kb:ke+ks) :: zz

    integer :: i,j,k,ii
    real :: tem,tscale,rts

    IF( tiris1.ge.1 )THEN

      if( adapt_dt.eq.1 )then
!!!        tscale = max( tscale0 , 4.0*adtlast )
        tscale = tscalefac*adtlast
      else
        tscale = tscalefac*dt
      endif
      rts = 1.0/tscale

    IF( cm1setup.eq.4 )THEN
      do k=1+ks,krecy
        tem = max( 0.0 , min( 1.0 , 1.0-(zz(1,1,k)-0.9*recy_depth)/(0.1*recy_depth)  )  )
        do j=1,nj+js
        do i=tiris1,tiris2
        if( cm0(i,j).gt.cmemin )then
          ii = tirisb-1 + i - tiris1+1
!!!          varten(i,j,k) = varten(i,j,k)-tem*rts*( var(i,j,k)-recywe(ii,j,k) )
          varten(i,j,k) = varten(i,j,k)-tem*rts*( (var(i,j,k)-timavg(i,j,k))-recywe(ii,j,k) )
        endif
        enddo
        enddo
      enddo
    ELSE
      do k=1+ks,krecy
        tem = max( 0.0 , min( 1.0 , 1.0-(zz(1,1,k)-0.9*recy_depth)/(0.1*recy_depth)  )  )
        do j=1,nj+js
        do i=tiris1,tiris2
          ii = tirisb-1 + i - tiris1+1
!!!          varten(i,j,k) = varten(i,j,k)-tem*rts*( var(i,j,k)-recywe(ii,j,k) )
          varten(i,j,k) = varten(i,j,k)-tem*rts*( (var(i,j,k)-timavg(i,j,k))-recywe(ii,j,k) )
        enddo
        enddo
      enddo
    ENDIF

    if( is.eq.0 .and. js.eq.0 .and. ks.eq.0 )then
      if( cm1setup.eq.4 )then
        do j=1,nj
        do i=tiris1,tiris2
          if( cm0(i,j).gt.cmemin ) recy_inj(i,j) = 1.0
        enddo
        enddo
      else
        do j=1,nj
        do i=tiris1,tiris2
          recy_inj(i,j) = 1.0
        enddo
        enddo
      endif
    endif

    ENDIF

    end subroutine recy_tendencywe

!-----------------------------------------------------------------------
#ifdef MPI

    ! send data from capture zone:
    subroutine   send_recywe(numi  ,numj,numk   , recywe,is,js,ks,var,timavg,dum,proc,tag,recy_cap,cm0,tircs1)
    use input
    use constants, only : cmemin
    use mpi
    implicit none

    integer, intent(in) :: numi,numj,numk,is,js,ks,proc,tag,tircs1
    real, intent(inout), dimension(irecywe,jrecywe,krecy) :: recywe
    real, intent(in), dimension(ib:ie+is,jb:je+js,kb:ke+ks) :: var
    real, intent(in), dimension(ibta:ieta,jbta:jeta,kbta:keta) :: timavg
    real, intent(inout), dimension(numi,numj,numk) :: dum
    integer, intent(inout), dimension(ib:ie,jb:je) :: recy_cap
    real, intent(in),    dimension(ib:ie,jb:je) :: cm0

    integer :: i,j,k

    do k=1,numk
    do j=1,numj
    do i=1,numi
      dum(i,j,k) = var(tircs1-1+i,j+jdisp,k)-timavg(tircs1-1+i,j+jdisp,k)
    enddo
    enddo
    enddo

  if( is.eq.0 .and. js.eq.0 .and. ks.eq.0 )then
    if( cm1setup.eq.4 )then
      do j=1,numj
      do i=1,numi
        if( cm0(tircs1-1+i,j).gt.cmemin ) recy_cap(tircs1-1+i,j) = 1.0
      enddo
      enddo
    else
      do j=1,numj
      do i=1,numi
        recy_cap(tircs1-1+i,j) = 1.0
      enddo
      enddo
    endif
  endif

    call mpi_send(dum,numi*numj*numk,MPI_REAL,proc,tag,MPI_COMM_WORLD,ierr)

    end subroutine send_recywe

!-----------------------------------------------------------------------

    ! p0, recv data from capture zone:
    subroutine     recv_recywe(numi,numj,numk   , recywe,irbt,dum           ,tag,proc)
    use input
    use mpi
    implicit none

    integer, intent(in) :: numi,numj,numk,irbt,tag,proc
    real, intent(inout), dimension(irecywe,jrecywe,krecy) :: recywe
    real, intent(inout), dimension(numi,numj,numk) :: dum

    integer :: i,j,k
    integer, dimension(mpi_status_size) :: status

    call mpi_recv(dum,numi*numj*numk,MPI_REAL,proc,tag,MPI_COMM_WORLD,status,ierr)

    do k=1,numk
    do j=1,numj
    do i=1,numi
      recywe(irbt-1+i,j,k) = dum(i,j,k)
    enddo
    enddo
    enddo

    end subroutine recv_recywe

!-----------------------------------------------------------------------

    ! recv data for injection zone:
    subroutine   recv_injwe(numi  ,numj,numk   , recywe,proc,tag,dum,tiris1,tirisb)
    use input
    use mpi
    implicit none

    integer, intent(in) :: numi,numj,numk,proc,tag,tiris1,tirisb
    real, intent(inout), dimension(irecywe,jrecywe,krecy) :: recywe
    real, intent(inout), dimension(numi,numj,numk) :: dum

    integer :: i,j,k
    integer, dimension(mpi_status_size) :: status

    call mpi_recv(dum,numi*numj*numk,MPI_REAL,proc,tag,MPI_COMM_WORLD,status,ierr)

    do k=1,numk
    do j=1,numj
    do i=1,numi
      recywe(tirisb-1+i,j,k) = dum(i,j,k)
    enddo
    enddo
    enddo

    end subroutine recv_injwe

!-----------------------------------------------------------------------

    ! p0, send data to injection zone:
    subroutine     send_injwe(numi   ,numj,numk   , recywe,irisbt,dum,tag,proc)
    use input
    use mpi
    implicit none

    integer, intent(in) :: numi,numj,numk,irisbt,tag,proc
    real, intent(inout), dimension(irecywe,jrecywe,krecy) :: recywe
    real, intent(inout), dimension(numi,numj,numk) :: dum

    integer :: i,j,k


    do k=1,numk
    do j=1,numj
    do i=1,numi
      dum(i,j,k) = recywe(irisbt-1+i,j,k)
    enddo
    enddo
    enddo

    call mpi_send(dum,numi*numj*numk,MPI_REAL,proc,tag,MPI_COMM_WORLD,ierr)

    end subroutine send_injwe

#endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   Eddy recycler on south side of domain:

    subroutine do_eddy_recys(dt,xh,xf,yh,yf,zh,zf,u0,v0,dum1,dum2,cm0,  &
                             u3d,v3d,w3d,th3d,q3d,tke3d,                &
                             uten1,vten1,wten1,thten1,qten,tketen,      &
                             recysn,out3d,recy_cap,recy_inj,timavg,adtlast)

    use input
    use constants
#ifdef MPI
    use mpi
#endif

      implicit none

      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: xh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: yh
      real, intent(in), dimension(jb:je+1) :: yf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2
      real, intent(in), dimension(ib:ie,jb:je) :: cm0
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0,u3d
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0,v3d
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: w3d
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: th3d
      real, intent(in), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: q3d
      real, intent(in), dimension(ibt:iet,jbt:jet,kbt:ket) :: tke3d
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: uten1
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: vten1
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: wten1
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: thten1
      real, intent(inout), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qten
      real, intent(inout), dimension(ibt:iet,jbt:jet,kbt:ket) :: tketen
      real, intent(inout), dimension(irecysn,jrecysn,krecy,nrecy) :: recysn
      real, intent(inout) , dimension(ib3d:ie3d,jb3d:je3d,kb3d:ke3d,nout3d) :: out3d
      integer, intent(inout), dimension(ib:ie,jb:je) :: recy_cap,recy_inj
      real, intent(in), dimension(ibta:ieta,jbta:jeta,kbta:keta,ntavr) :: timavg
      double precision, intent(in) :: adtlast

      integer :: i,j,k,n,nn,jj,proc,jrcst,jrbt,jrist,jrisbt
#ifdef MPI
      integer, dimension(mpi_status_size) :: status
      integer :: ndat,irbp
#endif


    !---------------------------------------------------------------------------
    !  step 1: populate recysn arrays (data to be recycled)

#ifndef MPI

          ! single-proc, shared memory:

          do k=1,krecy
          do j=sjrc1,sjrc2
            jj = j-sjrc1+1
          do i=1,ni+1
            recysn(i,jj,k,urecy) = u3d(i+idisp,j,k)-timavg(i+idisp,j,k,utav)
          enddo
          enddo
          enddo

          do k=1,krecy
          do j=svjrcs1,svjrcs2
            jj = j-svjrcs1+1
          do i=1,ni
            recysn(i,jj,k,vrecy) = v3d(i+idisp,j,k)-timavg(i+idisp,j,k,vtav)
          enddo
          enddo
          enddo

          do k=1,krecy
          do j=sjrc1,sjrc2
            jj = j-sjrc1+1
          do i=1,ni
            recysn(i,jj,k,wrecy) = w3d(i+idisp,j,k)-timavg(i+idisp,j,k,wtav)
            recysn(i,jj,k,trecy) = th3d(i+idisp,j,k)-timavg(i+idisp,j,k,ttav)
          enddo
          enddo
          enddo

          if( imoist.eq.1 )then
            do k=1,krecy
            do j=sjrc1,sjrc2
              jj = j-sjrc1+1
            do i=1,ni
              recysn(i,jj,k,qrecy) = q3d(i+idisp,j,k,nqv)-timavg(i+idisp,j,k,qtav)
            enddo
            enddo
            enddo
          endif

          if( etav.ge.1 )then
            do k=2,krecy
            do j=sjrc1,sjrc2
              jj = j-sjrc1+1
            do i=1,ni
              recysn(i,jj,k,erecy) = tke3d(i+idisp,j,k)-timavg(i+idisp,j,k,etav)
            enddo
            enddo
            enddo
          endif

#else
        ! MPI version:

        IF( myj.eq.1 )THEN

          if( svjrb.ge.1 )then
            do k=1,krecy
            do j=sjrcs1,sjrcs2
              jj = svjrb + (j-sjrcs1)
            do i=1,ni+1
              recysn(i,jj,k,urecy) = u3d(i+idisp,j,k)-timavg(i+idisp,j,k,utav)
            enddo
            enddo
            enddo
          endif

          if( sjrb.ge.1 )then
            do k=1,krecy
            do j=svjrcs1,svjrcs2
              jj = sjrb + (j-svjrcs1)
            do i=1,ni
              recysn(i,jj,k,vrecy) = v3d(i+idisp,j,k)-timavg(i+idisp,j,k,vtav)
            enddo
            enddo
            enddo

            do k=1,krecy
            do j=sjrcs1,sjrcs2
              jj = sjrb + (j-sjrcs1)
            do i=1,ni
              recysn(i,jj,k,wrecy) = w3d(i+idisp,j,k)-timavg(i+idisp,j,k,wtav)
              recysn(i,jj,k,trecy) = th3d(i+idisp,j,k)-timavg(i+idisp,j,k,ttav)
            enddo
            enddo
            enddo

          if( cm1setup.eq.4 )then
            do j=sjrcs1,sjrcs2
            do i=1,ni
              if( cm0(i,j).gt.cmemin ) recy_cap(i,j) = 1.0
            enddo
            enddo
          else
            do j=sjrcs1,sjrcs2
            do i=1,ni
              recy_cap(i,j) = 1.0
            enddo
            enddo
          endif

            if( imoist.eq.1 )then
              do k=1,krecy
              do j=sjrcs1,sjrcs2
                jj = sjrb + (j-sjrcs1)
              do i=1,ni
                recysn(i,jj,k,qrecy) = q3d(i+idisp,j,k,nqv)-timavg(i+idisp,j,k,qtav)
              enddo
              enddo
              enddo
            endif

            if( etav.ge.1 )then
              do k=1,krecy
              do j=sjrcs1,sjrcs2
                jj = sjrb + (j-sjrcs1)
              do i=1,ni
                recysn(i,jj,k,erecy) = tke3d(i+idisp,j,k)-timavg(i+idisp,j,k,etav)
              enddo
              enddo
              enddo
            endif

          endif

        ENDIF

    !---------------------------------------------------------------------------
    !  step 2: communicate data (MPI only)

        IF( myj.eq.1 )THEN

          ! south-most processor receives data, collects, then sends
          ! I am south-most proc in this column:
          ! loop through all procs in this column:
          ! receive data, add to recysn arrays:
          do n=2,nodey
            proc = myid + (n-1)*nodex
            ! jrcst = width of data to receive
            call mpi_recv(jrcst,1,MPI_INTEGER,proc,8201,MPI_COMM_WORLD,status,ierr)
            if( jrcst.gt.0 )then
              ! jrbt = starting point of data in the recycle array
              call mpi_recv(jrbt,1,MPI_INTEGER,proc,8202,MPI_COMM_WORLD,status,ierr)
              call recv_recysn(ni  ,jrcst,krecy  ,recysn(1,1,1,vrecy),jrbt,dum1(ib,jb,kb),8211,proc)
            endif
            call mpi_recv(jrcst,1,MPI_INTEGER,proc,8203,MPI_COMM_WORLD,status,ierr)
            if( jrcst.gt.0 )then
              call mpi_recv(jrbt,1,MPI_INTEGER,proc,8204,MPI_COMM_WORLD,status,ierr)
              call recv_recysn(ni+1,jrcst,krecy  ,recysn(1,1,1,urecy),jrbt,dum1(ib,jb,kb),8212,proc)
              call recv_recysn(ni  ,jrcst,krecy  ,recysn(1,1,1,wrecy),jrbt,dum1(ib,jb,kb),8213,proc)
              call recv_recysn(ni  ,jrcst,krecy  ,recysn(1,1,1,trecy),jrbt,dum1(ib,jb,kb),8214,proc)
                if( etav.ge.1 )  &
              call recv_recysn(ni  ,jrcst,krecy  ,recysn(1,1,1,erecy),jrbt,dum1(ib,jb,kb),8215,proc)
                if( imoist.eq.1 )  &
              call recv_recysn(ni  ,jrcst,krecy  ,recysn(1,1,1,qrecy),jrbt,dum1(ib,jb,kb),8216,proc)
            endif
          enddo

        ELSE

          ! I am NOT south-most proc in this column: 
          ! send my data to south-most proc (if applicable):
          proc = (myi-1)
          ! jrcs = width of data to send
          call mpi_send(svjrcs,1,MPI_INTEGER,proc,8201,MPI_COMM_WORLD,ierr)
          if( svjrcs1.gt.0 .and. svjrcs2.gt.0 )then
            ! jrb = starting point of data in the recycle array
            call mpi_send(svjrb,1,MPI_INTEGER,proc,8202,MPI_COMM_WORLD,ierr)
            call send_recysn(ni  ,svjrcs,krecy  ,recysn(1,1,1,vrecy),0 ,1 ,0 ,v3d,  timavg(ibta,jbta,kbta,vtav),dum1(ib,jb,kb),proc,8211,recy_cap,cm0,svjrcs1)
          endif
          call mpi_send(sjrcs,1,MPI_INTEGER,proc,8203,MPI_COMM_WORLD,ierr)
          if( sjrcs1.gt.0 .and. sjrcs2.gt.0 )then
            call mpi_send(sjrb,1,MPI_INTEGER,proc,8204,MPI_COMM_WORLD,ierr)
            call send_recysn(ni+1,sjrcs,krecy  ,recysn(1,1,1,urecy),1 ,0 ,0 ,u3d,  timavg(ibta,jbta,kbta,utav),dum1(ib,jb,kb),proc,8212,recy_cap,cm0,sjrcs1)
            call send_recysn(ni  ,sjrcs,krecy  ,recysn(1,1,1,wrecy),0 ,0 ,1 ,w3d,  timavg(ibta,jbta,kbta,wtav),dum1(ib,jb,kb),proc,8213,recy_cap,cm0,sjrcs1)
            call send_recysn(ni  ,sjrcs,krecy  ,recysn(1,1,1,trecy),0 ,0 ,0 ,th3d, timavg(ibta,jbta,kbta,ttav),dum1(ib,jb,kb),proc,8214,recy_cap,cm0,sjrcs1)
              if( etav.ge.1 )  &
            call send_recysn(ni  ,sjrcs,krecy  ,recysn(1,1,1,erecy),0 ,0 ,1 ,tke3d,timavg(ibta,jbta,kbta,etav),dum1(ib,jb,kb),proc,8215,recy_cap,cm0,sjrcs1)
              if( imoist.eq.1 )  &
            call send_recysn(ni  ,sjrcs,krecy  ,recysn(1,1,1,qrecy),0 ,0 ,0 ,q3d(ib,jb,kb,nqv),timavg(ibta,jbta,kbta,qtav),dum1(ib,jb,kb),proc,8216,recy_cap,cm0,sjrcs1)
          endif

        ENDIF

      !-------------------------------------------------------------------------

        IF( myj.eq.1 )THEN

          ! loop through all procs in this column:
          ! send data to procs in injection region:
          do n=2,nodey
            proc = myid + (n-1)*nodex
            ! jrist = width of data to receive
            call mpi_recv(jrist,1,MPI_INTEGER,proc,8221,MPI_COMM_WORLD,status,ierr)
            if( jrist.gt.0 )then
              ! jrisbt = starting point of data in the recycle array
              call mpi_recv(jrisbt,1,MPI_INTEGER,proc,8222,MPI_COMM_WORLD,status,ierr)
              call send_injsn(ni  ,jrist,krecy  ,recysn(1,1,1,vrecy),jrisbt,dum1(ib,jb,kb),8231,proc)
            endif 
            call mpi_recv(jrist,1,MPI_INTEGER,proc,8223,MPI_COMM_WORLD,status,ierr)
            if( jrist.gt.0 )then
              call mpi_recv(jrisbt,1,MPI_INTEGER,proc,8224,MPI_COMM_WORLD,status,ierr)
              call send_injsn(ni+1,jrist,krecy  ,recysn(1,1,1,urecy),jrisbt,dum1(ib,jb,kb),8232,proc)
              call send_injsn(ni  ,jrist,krecy  ,recysn(1,1,1,wrecy),jrisbt,dum1(ib,jb,kb),8233,proc)
              call send_injsn(ni  ,jrist,krecy  ,recysn(1,1,1,trecy),jrisbt,dum1(ib,jb,kb),8234,proc)
                if( etav.ge.1 )  &
              call send_injsn(ni  ,jrist,krecy  ,recysn(1,1,1,erecy),jrisbt,dum1(ib,jb,kb),8235,proc)
                if( imoist.eq.1 )  &
              call send_injsn(ni  ,jrist,krecy  ,recysn(1,1,1,qrecy),jrisbt,dum1(ib,jb,kb),8236,proc)
            endif
          enddo


        ELSE

          ! receive from south-most proc (if applicable):
          proc = (myi-1)
          ! iris = width of data to send
          call mpi_send(svjris,1,MPI_INTEGER,proc,8221,MPI_COMM_WORLD,ierr)
          if( svjris1.gt.0 .and. svjris2.gt.0 )then
            ! irisb = starting point of data in the recycle array
            call mpi_send(svjrisb,1,MPI_INTEGER,proc,8222,MPI_COMM_WORLD,ierr)
            call recv_injsn(ni  ,svjris,krecy  ,recysn(1,1,1,vrecy),proc,8231,dum1(ib,jb,kb),svjris1,svjrisb)
          endif
          call mpi_send(sjris,1,MPI_INTEGER,proc,8223,MPI_COMM_WORLD,ierr)
          if( sjris1.gt.0 .and. sjris2.gt.0 )then
            call mpi_send(sjrisb,1,MPI_INTEGER,proc,8224,MPI_COMM_WORLD,ierr)
            call recv_injsn(ni+1,sjris,krecy  ,recysn(1,1,1,urecy),proc,8232,dum1(ib,jb,kb),sjris1,sjrisb)
            call recv_injsn(ni  ,sjris,krecy  ,recysn(1,1,1,wrecy),proc,8233,dum1(ib,jb,kb),sjris1,sjrisb)
            call recv_injsn(ni  ,sjris,krecy  ,recysn(1,1,1,trecy),proc,8234,dum1(ib,jb,kb),sjris1,sjrisb)
              if( etav.ge.1 )  &
            call recv_injsn(ni  ,sjris,krecy  ,recysn(1,1,1,erecy),proc,8235,dum1(ib,jb,kb),sjris1,sjrisb)
              if( imoist.eq.1 )  &
            call recv_injsn(ni  ,sjris,krecy  ,recysn(1,1,1,qrecy),proc,8236,dum1(ib,jb,kb),sjris1,sjrisb)
          endif

        ENDIF

      !-------------------------------------------------------------------------

#endif

      ! get recycle tendencies:

      call     recy_tendencysn(recysn(1,1,1,urecy),1 ,0 ,0 ,u3d  ,timavg(ibta,jbta,kbta,utav),uten1 ,recy_inj,cm0,out3d, 1, 7,13,dt,adtlast,zh,sjris1,sjris2,sjrisb)
      call     recy_tendencysn(recysn(1,1,1,vrecy),0 ,1 ,0 ,v3d  ,timavg(ibta,jbta,kbta,vtav),vten1 ,recy_inj,cm0,out3d, 2, 8,14,dt,adtlast,zh,svjris1,svjris2,svjrisb)
      call     recy_tendencysn(recysn(1,1,1,wrecy),0 ,0 ,1 ,w3d  ,timavg(ibta,jbta,kbta,wtav),wten1 ,recy_inj,cm0,out3d, 3, 9,15,dt,adtlast,zf,sjris1,sjris2,sjrisb)
      call     recy_tendencysn(recysn(1,1,1,trecy),0 ,0 ,0 ,th3d ,timavg(ibta,jbta,kbta,ttav),thten1,recy_inj,cm0,out3d, 4,10,16,dt,adtlast,zh,sjris1,sjris2,sjrisb)
      if( etav.ge.1 )  &
      call     recy_tendencysn(recysn(1,1,1,erecy),0 ,0 ,1 ,tke3d,timavg(ibta,jbta,kbta,etav),tketen,recy_inj,cm0,out3d, 5,11,17,dt,adtlast,zf,sjris1,sjris2,sjrisb)
      if( imoist.eq.1 )  &
      call     recy_tendencysn(recysn(1,1,1,qrecy),0 ,0 ,0 ,q3d(ib,jb,kb,nqv),timavg(ibta,jbta,kbta,qtav),qten(ib,jb,kb,nqv),recy_inj,cm0,out3d, 6,12,18,dt,adtlast,zh,sjris1,sjris2,sjrisb)

    end subroutine do_eddy_recys

!-----------------------------------------------------------------------
!   Eddy recycler on north side of domain:

    subroutine do_eddy_recyn(dt,xh,xf,yh,yf,zh,zf,u0,v0,dum1,dum2,cm0,  &
                             u3d,v3d,w3d,th3d,q3d,tke3d,                &
                             uten1,vten1,wten1,thten1,qten,tketen,      &
                             recysn,out3d,recy_cap,recy_inj,timavg,adtlast)

    use input
    use constants
#ifdef MPI
    use mpi
#endif

      implicit none

      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: xh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: yh
      real, intent(in), dimension(jb:je+1) :: yf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2
      real, intent(in), dimension(ib:ie,jb:je) :: cm0
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0,u3d
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0,v3d
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: w3d
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: th3d
      real, intent(in), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: q3d
      real, intent(in), dimension(ibt:iet,jbt:jet,kbt:ket) :: tke3d
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: uten1
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: vten1
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: wten1
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: thten1
      real, intent(inout), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qten
      real, intent(inout), dimension(ibt:iet,jbt:jet,kbt:ket) :: tketen
      real, intent(inout), dimension(irecysn,jrecysn,krecy,nrecy) :: recysn
      real, intent(inout) , dimension(ib3d:ie3d,jb3d:je3d,kb3d:ke3d,nout3d) :: out3d
      integer, intent(inout), dimension(ib:ie,jb:je) :: recy_cap,recy_inj
      real, intent(in), dimension(ibta:ieta,jbta:jeta,kbta:keta,ntavr) :: timavg
      double precision, intent(in) :: adtlast

      integer :: i,j,k,n,nn,jj,proc,jrcst,jrbt,jrist,jrisbt
#ifdef MPI
      integer, dimension(mpi_status_size) :: status
      integer :: ndat,irbp
#endif


    !---------------------------------------------------------------------------
    !  step 1: populate recysn arrays (data to be recycled)

#ifndef MPI

          ! single-proc, shared memory:

          do k=1,krecy
          do j=njrc1,njrc2
            jj = j-njrc1+1
          do i=1,ni+1
            recysn(i,jj,k,urecy) = u3d(i+idisp,j,k)-timavg(i+idisp,j,k,utav)
          enddo
          enddo
          enddo

          do k=1,krecy
          do j=nvjrcs1,nvjrcs2
            jj = j-nvjrcs1+1
          do i=1,ni
            recysn(i,jj,k,vrecy) = v3d(i+idisp,j,k)-timavg(i+idisp,j,k,vtav)
          enddo
          enddo
          enddo

          do k=1,krecy
          do j=njrc1,njrc2
            jj = j-njrc1+1
          do i=1,ni
            recysn(i,jj,k,wrecy) = w3d(i+idisp,j,k)-timavg(i+idisp,j,k,wtav)
            recysn(i,jj,k,trecy) = th3d(i+idisp,j,k)-timavg(i+idisp,j,k,ttav)
          enddo
          enddo
          enddo

          if( imoist.eq.1 )then
            do k=1,krecy
            do j=njrc1,njrc2
              jj = j-njrc1+1
            do i=1,ni
              recysn(i,jj,k,qrecy) = q3d(i+idisp,j,k,nqv)-timavg(i+idisp,j,k,qtav)
            enddo
            enddo
            enddo
          endif

          if( etav.ge.1 )then
            do k=2,krecy
            do j=njrc1,njrc2
              jj = j-njrc1+1
            do i=1,ni
              recysn(i,jj,k,erecy) = tke3d(i+idisp,j,k)-timavg(i+idisp,j,k,etav)
            enddo
            enddo
            enddo
          endif

#else
        ! MPI version:

        IF( myj.eq.1 )THEN

          if( nvjrb.ge.1 )then
            do k=1,krecy
            do j=njrcs1,njrcs2
              jj = nvjrb + (j-njrcs1)
            do i=1,ni+1
              recysn(i,jj,k,urecy) = u3d(i+idisp,j,k)-timavg(i+idisp,j,k,utav)
            enddo
            enddo
            enddo
          endif

          if( njrb.ge.1 )then
            do k=1,krecy
            do j=nvjrcs1,nvjrcs2
              jj = njrb + (j-nvjrcs1)
            do i=1,ni
              recysn(i,jj,k,vrecy) = v3d(i+idisp,j,k)-timavg(i+idisp,j,k,vtav)
            enddo
            enddo
            enddo

            do k=1,krecy
            do j=njrcs1,njrcs2
              jj = njrb + (j-njrcs1)
            do i=1,ni
              recysn(i,jj,k,wrecy) = w3d(i+idisp,j,k)-timavg(i+idisp,j,k,wtav)
              recysn(i,jj,k,trecy) = th3d(i+idisp,j,k)-timavg(i+idisp,j,k,ttav)
            enddo
            enddo
            enddo

          if( cm1setup.eq.4 )then
            do j=njrcs1,njrcs2
            do i=1,ni
              if( cm0(i,j).gt.cmemin ) recy_cap(i,j) = 1.0
            enddo
            enddo
          else
            do j=njrcs1,njrcs2
            do i=1,ni
              recy_cap(i,j) = 1.0
            enddo
            enddo
          endif

            if( imoist.eq.1 )then
              do k=1,krecy
              do j=njrcs1,njrcs2
                jj = njrb + (j-njrcs1)
              do i=1,ni
                recysn(i,jj,k,qrecy) = q3d(i+idisp,j,k,nqv)-timavg(i+idisp,j,k,qtav)
              enddo
              enddo
              enddo
            endif

            if( etav.ge.1 )then
              do k=1,krecy
              do j=njrcs1,njrcs2
                jj = njrb + (j-njrcs1)
              do i=1,ni
                recysn(i,jj,k,erecy) = tke3d(i+idisp,j,k)-timavg(i+idisp,j,k,etav)
              enddo
              enddo
              enddo
            endif

          endif

        ENDIF

    !---------------------------------------------------------------------------
    !  step 2: communicate data (MPI only)

        IF( myj.eq.1 )THEN

          ! south-most processor receives data, collects, then sends
          ! I am south-most proc in this column:
          ! loop through all procs in this column:
          ! receive data, add to recysn arrays:
          do n=2,nodey
            proc = myid + (n-1)*nodex
            ! jrcst = width of data to receive
            call mpi_recv(jrcst,1,MPI_INTEGER,proc,8301,MPI_COMM_WORLD,status,ierr)
            if( jrcst.gt.0 )then
              ! jrbt = starting point of data in the recycle array
              call mpi_recv(jrbt,1,MPI_INTEGER,proc,8302,MPI_COMM_WORLD,status,ierr)
              call recv_recysn(ni  ,jrcst,krecy  ,recysn(1,1,1,vrecy),jrbt,dum1(ib,jb,kb),8311,proc)
            endif
            call mpi_recv(jrcst,1,MPI_INTEGER,proc,8303,MPI_COMM_WORLD,status,ierr)
            if( jrcst.gt.0 )then
              call mpi_recv(jrbt,1,MPI_INTEGER,proc,8304,MPI_COMM_WORLD,status,ierr)
              call recv_recysn(ni+1,jrcst,krecy  ,recysn(1,1,1,urecy),jrbt,dum1(ib,jb,kb),8312,proc)
              call recv_recysn(ni  ,jrcst,krecy  ,recysn(1,1,1,wrecy),jrbt,dum1(ib,jb,kb),8313,proc)
              call recv_recysn(ni  ,jrcst,krecy  ,recysn(1,1,1,trecy),jrbt,dum1(ib,jb,kb),8314,proc)
                if( etav.ge.1 )  &
              call recv_recysn(ni  ,jrcst,krecy  ,recysn(1,1,1,erecy),jrbt,dum1(ib,jb,kb),8315,proc)
                if( imoist.eq.1 )  &
              call recv_recysn(ni  ,jrcst,krecy  ,recysn(1,1,1,qrecy),jrbt,dum1(ib,jb,kb),8316,proc)
            endif
          enddo

        ELSE

          ! I am NOT south-most proc in this column: 
          ! send my data to south-most proc (if applicable):
          proc = (myi-1)
          ! jrcs = width of data to send
          call mpi_send(nvjrcs,1,MPI_INTEGER,proc,8301,MPI_COMM_WORLD,ierr)
          if( nvjrcs1.gt.0 .and. nvjrcs2.gt.0 )then
            ! jrb = starting point of data in the recycle array
            call mpi_send(nvjrb,1,MPI_INTEGER,proc,8302,MPI_COMM_WORLD,ierr)
            call send_recysn(ni  ,nvjrcs,krecy  ,recysn(1,1,1,vrecy),0 ,1 ,0 ,v3d,  timavg(ibta,jbta,kbta,vtav),dum1(ib,jb,kb),proc,8311,recy_cap,cm0,nvjrcs1)
          endif
          call mpi_send(njrcs,1,MPI_INTEGER,proc,8303,MPI_COMM_WORLD,ierr)
          if( njrcs1.gt.0 .and. njrcs2.gt.0 )then
            call mpi_send(njrb,1,MPI_INTEGER,proc,8304,MPI_COMM_WORLD,ierr)
            call send_recysn(ni+1,njrcs,krecy  ,recysn(1,1,1,urecy),1 ,0 ,0 ,u3d,  timavg(ibta,jbta,kbta,utav),dum1(ib,jb,kb),proc,8312,recy_cap,cm0,njrcs1)
            call send_recysn(ni  ,njrcs,krecy  ,recysn(1,1,1,wrecy),0 ,0 ,1 ,w3d,  timavg(ibta,jbta,kbta,wtav),dum1(ib,jb,kb),proc,8313,recy_cap,cm0,njrcs1)
            call send_recysn(ni  ,njrcs,krecy  ,recysn(1,1,1,trecy),0 ,0 ,0 ,th3d, timavg(ibta,jbta,kbta,ttav),dum1(ib,jb,kb),proc,8314,recy_cap,cm0,njrcs1)
              if( etav.ge.1 )  &
            call send_recysn(ni  ,njrcs,krecy  ,recysn(1,1,1,erecy),0 ,0 ,1 ,tke3d,timavg(ibta,jbta,kbta,etav),dum1(ib,jb,kb),proc,8315,recy_cap,cm0,njrcs1)
              if( imoist.eq.1 )  &
            call send_recysn(ni  ,njrcs,krecy  ,recysn(1,1,1,qrecy),0 ,0 ,0 ,q3d(ib,jb,kb,nqv),timavg(ibta,jbta,kbta,qtav),dum1(ib,jb,kb),proc,8316,recy_cap,cm0,njrcs1)
          endif

        ENDIF

      !-------------------------------------------------------------------------

        IF( myj.eq.1 )THEN

          ! loop through all procs in this column:
          ! send data to procs in injection region:
          do n=2,nodey
            proc = myid + (n-1)*nodex
            ! jrist = width of data to receive
            call mpi_recv(jrist,1,MPI_INTEGER,proc,8321,MPI_COMM_WORLD,status,ierr)
            if( jrist.gt.0 )then
              ! jrisbt = starting point of data in the recycle array
              call mpi_recv(jrisbt,1,MPI_INTEGER,proc,8322,MPI_COMM_WORLD,status,ierr)
              call send_injsn(ni  ,jrist,krecy  ,recysn(1,1,1,vrecy),jrisbt,dum1(ib,jb,kb),8331,proc)
            endif 
            call mpi_recv(jrist,1,MPI_INTEGER,proc,8323,MPI_COMM_WORLD,status,ierr)
            if( jrist.gt.0 )then
              call mpi_recv(jrisbt,1,MPI_INTEGER,proc,8324,MPI_COMM_WORLD,status,ierr)
              call send_injsn(ni+1,jrist,krecy  ,recysn(1,1,1,urecy),jrisbt,dum1(ib,jb,kb),8332,proc)
              call send_injsn(ni  ,jrist,krecy  ,recysn(1,1,1,wrecy),jrisbt,dum1(ib,jb,kb),8333,proc)
              call send_injsn(ni  ,jrist,krecy  ,recysn(1,1,1,trecy),jrisbt,dum1(ib,jb,kb),8334,proc)
                if( etav.ge.1 )  &
              call send_injsn(ni  ,jrist,krecy  ,recysn(1,1,1,erecy),jrisbt,dum1(ib,jb,kb),8335,proc)
                if( imoist.eq.1 )  &
              call send_injsn(ni  ,jrist,krecy  ,recysn(1,1,1,qrecy),jrisbt,dum1(ib,jb,kb),8336,proc)
            endif
          enddo


        ELSE

          ! receive from south-most proc (if applicable):
          proc = (myi-1)
          ! iris = width of data to send
          call mpi_send(nvjris,1,MPI_INTEGER,proc,8321,MPI_COMM_WORLD,ierr)
          if( nvjris1.gt.0 .and. nvjris2.gt.0 )then
            ! irisb = starting point of data in the recycle array
            call mpi_send(nvjrisb,1,MPI_INTEGER,proc,8322,MPI_COMM_WORLD,ierr)
            call recv_injsn(ni  ,nvjris,krecy  ,recysn(1,1,1,vrecy),proc,8331,dum1(ib,jb,kb),nvjris1,nvjrisb)
          endif
          call mpi_send(njris,1,MPI_INTEGER,proc,8323,MPI_COMM_WORLD,ierr)
          if( njris1.gt.0 .and. njris2.gt.0 )then
            call mpi_send(njrisb,1,MPI_INTEGER,proc,8324,MPI_COMM_WORLD,ierr)
            call recv_injsn(ni+1,njris,krecy  ,recysn(1,1,1,urecy),proc,8332,dum1(ib,jb,kb),njris1,njrisb)
            call recv_injsn(ni  ,njris,krecy  ,recysn(1,1,1,wrecy),proc,8333,dum1(ib,jb,kb),njris1,njrisb)
            call recv_injsn(ni  ,njris,krecy  ,recysn(1,1,1,trecy),proc,8334,dum1(ib,jb,kb),njris1,njrisb)
              if( etav.ge.1 )  &
            call recv_injsn(ni  ,njris,krecy  ,recysn(1,1,1,erecy),proc,8335,dum1(ib,jb,kb),njris1,njrisb)
              if( imoist.eq.1 )  &
            call recv_injsn(ni  ,njris,krecy  ,recysn(1,1,1,qrecy),proc,8336,dum1(ib,jb,kb),njris1,njrisb)
          endif

        ENDIF

      !-------------------------------------------------------------------------

#endif

      ! get recycle tendencies:

      call     recy_tendencysn(recysn(1,1,1,urecy),1 ,0 ,0 ,u3d  ,timavg(ibta,jbta,kbta,utav),uten1 ,recy_inj,cm0,out3d, 1, 7,13,dt,adtlast,zh,njris1,njris2,njrisb)
      call     recy_tendencysn(recysn(1,1,1,vrecy),0 ,1 ,0 ,v3d  ,timavg(ibta,jbta,kbta,vtav),vten1 ,recy_inj,cm0,out3d, 2, 8,14,dt,adtlast,zh,nvjris1,nvjris2,nvjrisb)
      call     recy_tendencysn(recysn(1,1,1,wrecy),0 ,0 ,1 ,w3d  ,timavg(ibta,jbta,kbta,wtav),wten1 ,recy_inj,cm0,out3d, 3, 9,15,dt,adtlast,zf,njris1,njris2,njrisb)
      call     recy_tendencysn(recysn(1,1,1,trecy),0 ,0 ,0 ,th3d ,timavg(ibta,jbta,kbta,ttav),thten1,recy_inj,cm0,out3d, 4,10,16,dt,adtlast,zh,njris1,njris2,njrisb)
      if( etav.ge.1 )  &
      call     recy_tendencysn(recysn(1,1,1,erecy),0 ,0 ,1 ,tke3d,timavg(ibta,jbta,kbta,etav),tketen,recy_inj,cm0,out3d, 5,11,17,dt,adtlast,zf,njris1,njris2,njrisb)
      if( imoist.eq.1 )  &
      call     recy_tendencysn(recysn(1,1,1,qrecy),0 ,0 ,0 ,q3d(ib,jb,kb,nqv),timavg(ibta,jbta,kbta,qtav),qten(ib,jb,kb,nqv),recy_inj,cm0,out3d, 6,12,18,dt,adtlast,zh,njris1,njris2,njrisb)

    end subroutine do_eddy_recyn

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine recy_tendencysn(recysn,is,js,ks,var  ,timavg,varten,recy_inj,cm0,out3d,o1,o2,o3,dt,adtlast,zz,tjris1,tjris2,tjrisb)
    use input
    use constants, only : cmemin
    implicit none

    real, intent(in), dimension(irecysn,jrecysn,krecy) :: recysn
    integer, intent(in) :: is,js,ks,o1,o2,o3,tjris1,tjris2,tjrisb
    real, intent(in   ), dimension(ib:ie+is,jb:je+js,kb:ke+ks) :: var
    real, intent(inout), dimension(ib:ie+is,jb:je+js,kb:ke+ks) :: varten
    real, intent(in),    dimension(ibta:ieta,jbta:jeta,kbta:keta) :: timavg
    integer, intent(inout), dimension(ib:ie,jb:je) :: recy_inj
    real, intent(in),    dimension(ib:ie,jb:je) :: cm0
    real, intent(inout) , dimension(ib3d:ie3d,jb3d:je3d,kb3d:ke3d,nout3d) :: out3d
    real, intent(in) :: dt
    double precision, intent(in) :: adtlast
    real, intent(in), dimension(ib:ie,jb:je,kb:ke+ks) :: zz

    integer :: i,j,k,jj
    real :: tem,tscale,rts

    IF( tjris1.ge.1 )THEN

      if( adapt_dt.eq.1 )then
!!!        tscale = max( tscale0 , 4.0*adtlast )
        tscale = tscalefac*adtlast
      else
        tscale = tscalefac*dt
      endif
      rts = 1.0/tscale

    IF( cm1setup.eq.4 )THEN
      do k=1+ks,krecy
        tem = max( 0.0 , min( 1.0 , 1.0-(zz(1,1,k)-0.9*recy_depth)/(0.1*recy_depth)  )  )
        do j=tjris1,tjris2
          jj = tjrisb-1 + j - tjris1+1
        do i=1,ni+is
        if( cm0(i,j).gt.cmemin )then
!!!          varten(i,j,k) = varten(i,j,k)-tem*rts*( var(i,j,k)-recysn(i,jj,k) )
          varten(i,j,k) = varten(i,j,k)-tem*rts*( (var(i,j,k)-timavg(i,j,k))-recysn(i,jj,k) )
        endif
        enddo
        enddo
      enddo
    ELSE
      do k=1+ks,krecy
        tem = max( 0.0 , min( 1.0 , 1.0-(zz(1,1,k)-0.9*recy_depth)/(0.1*recy_depth)  )  )
        do j=tjris1,tjris2
          jj = tjrisb-1 + j - tjris1+1
        do i=1,ni+is
!!!          varten(i,j,k) = varten(i,j,k)-tem*rts*( var(i,j,k)-recysn(i,jj,k) )
          varten(i,j,k) = varten(i,j,k)-tem*rts*( (var(i,j,k)-timavg(i,j,k))-recysn(i,jj,k) )
        enddo
        enddo
      enddo
    ENDIF

    if( is.eq.0 .and. js.eq.0 .and. ks.eq.0 )then
      if( cm1setup.eq.4 )then
        do j=tjris1,tjris2
        do i=1,ni
          if( cm0(i,j).gt.cmemin ) recy_inj(i,j) = 1.0
        enddo
        enddo
      else
        do j=tjris1,tjris2
        do i=1,ni
          recy_inj(i,j) = 1.0
        enddo
        enddo
      endif
    endif

    ENDIF

    end subroutine recy_tendencysn

!-----------------------------------------------------------------------
#ifdef MPI

    ! send data from capture zone:
    subroutine   send_recysn(numi  ,numj,numk   , recysn,is,js,ks,var,timavg,dum,proc,tag,recy_cap,cm0,tjrcs1)
    use input
    use constants, only : cmemin
    use mpi
    implicit none

    integer, intent(in) :: numi,numj,numk,is,js,ks,proc,tag,tjrcs1
    real, intent(inout), dimension(irecysn,jrecysn,krecy) :: recysn
    real, intent(in), dimension(ib:ie+is,jb:je+js,kb:ke+ks) :: var
    real, intent(in), dimension(ibta:ieta,jbta:jeta,kbta:keta) :: timavg
    real, intent(inout), dimension(numi,numj,numk) :: dum
    integer, intent(inout), dimension(ib:ie,jb:je) :: recy_cap
    real, intent(in),    dimension(ib:ie,jb:je) :: cm0

    integer :: i,j,k

    do k=1,numk
    do j=1,numj
    do i=1,numi
      dum(i,j,k) = var(i+idisp,tjrcs1-1+j,k)-timavg(i+idisp,tjrcs1-1+j,k)
    enddo
    enddo
    enddo

  if( is.eq.0 .and. js.eq.0 .and. ks.eq.0 )then
    if( cm1setup.eq.4 )then
      do j=1,numj
      do i=1,numi
        if( cm0(i,tjrcs1-1+j).gt.cmemin ) recy_cap(i,tjrcs1-1+j) = 1.0
      enddo
      enddo
    else
      do j=1,numj
      do i=1,numi
        recy_cap(i,tjrcs1-1+j) = 1.0
      enddo
      enddo
    endif
  endif

    call mpi_send(dum,numi*numj*numk,MPI_REAL,proc,tag,MPI_COMM_WORLD,ierr)

    end subroutine send_recysn

!-----------------------------------------------------------------------

    ! p0, recv data from capture zone:
    subroutine     recv_recysn(numi,numj,numk   , recysn,jrbt,dum           ,tag,proc)
    use input
    use mpi
    implicit none

    integer, intent(in) :: numi,numj,numk,jrbt,tag,proc
    real, intent(inout), dimension(irecysn,jrecysn,krecy) :: recysn
    real, intent(inout), dimension(numi,numj,numk) :: dum

    integer :: i,j,k
    integer, dimension(mpi_status_size) :: status

    call mpi_recv(dum,numi*numj*numk,MPI_REAL,proc,tag,MPI_COMM_WORLD,status,ierr)

    do k=1,numk
    do j=1,numj
    do i=1,numi
      recysn(i,jrbt-1+j,k) = dum(i,j,k)
    enddo
    enddo
    enddo

    end subroutine recv_recysn

!-----------------------------------------------------------------------

    ! recv data for injection zone:
    subroutine   recv_injsn(numi  ,numj,numk   , recysn,proc,tag,dum,tjris1,tjrisb)
    use input
    use mpi
    implicit none

    integer, intent(in) :: numi,numj,numk,proc,tag,tjris1,tjrisb
    real, intent(inout), dimension(irecysn,jrecysn,krecy) :: recysn
    real, intent(inout), dimension(numi,numj,numk) :: dum

    integer :: i,j,k
    integer, dimension(mpi_status_size) :: status

    call mpi_recv(dum,numi*numj*numk,MPI_REAL,proc,tag,MPI_COMM_WORLD,status,ierr)

    do k=1,numk
    do j=1,numj
    do i=1,numi
      recysn(i,tjrisb-1+j,k) = dum(i,j,k)
    enddo
    enddo
    enddo

    end subroutine recv_injsn

!-----------------------------------------------------------------------

    ! p0, send data to injection zone:
    subroutine     send_injsn(numi   ,numj,numk   , recysn,jrisbt,dum,tag,proc)
    use input
    use mpi
    implicit none

    integer, intent(in) :: numi,numj,numk,jrisbt,tag,proc
    real, intent(inout), dimension(irecysn,jrecysn,krecy) :: recysn
    real, intent(inout), dimension(numi,numj,numk) :: dum

    integer :: i,j,k


    do k=1,numk
    do j=1,numj
    do i=1,numi
      dum(i,j,k) = recysn(i,jrisbt-1+j,k)
    enddo
    enddo
    enddo

    call mpi_send(dum,numi*numj*numk,MPI_REAL,proc,tag,MPI_COMM_WORLD,ierr)

    end subroutine send_injsn

!-----------------------------------------------------------------------
#endif

    subroutine eddy_recycling_indices(xh,xf,yh,yf,xhref,xfref,yhref,yfref)
    use input
#ifdef MPI
    use mpi
#endif
    implicit none

    real, intent(in), dimension(ib:ie) :: xh
    real, intent(in), dimension(ib:ie+1) :: xf
    real, intent(in), dimension(jb:je) :: yh
    real, intent(in), dimension(jb:je+1) :: yf
    real, intent(in), dimension(1-ngxy:nx+ngxy+1) :: xfref,xhref
    real, intent(in), dimension(1-ngxy:ny+ngxy+1) :: yfref,yhref

    integer :: i,j,n
    real :: ds
    real, dimension(:), allocatable :: foor
    integer, dimension(:), allocatable :: fooi
#ifdef MPI
    integer, dimension(MPI_STATUS_SIZE) :: status
#endif

      !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
      !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
      ! west&east sides of domain:


      dorecyew:  &
      IF( do_recycle_w .or. do_recycle_e )THEN

      ! ........  capture zone info  ........ !

        ! begin/end points of "capture zone" on total domain:
        wirc1 = 0
        wirc2 = 0
        eirc1 = 0
        eirc2 = 0
       ! loop through entire domain:
        do i=1,nx
          ds = 0.1*(xfref(i+1)-xfref(i))
          if( xhref(i).ge.(recy_cap_w-ds) .and. xhref(i).le.(recy_cap_w+recy_width+ds) )then
            if( wirc1.eq.0 ) wirc1 = i
            wirc2 = i
          endif
          if( xhref(i).le.(recy_cap_e+ds) .and. xhref(i).ge.(recy_cap_e-recy_width-ds) )then
            if( eirc1.eq.0 ) eirc1 = i
            eirc2 = i
          endif
        enddo

        wirc  = 0
        ! total width of capture zone
        if( wirc2.gt.0 .and. wirc1.gt.0 )then
          wirc  = wirc2-wirc1+1
        endif
        eirc  = 0
        ! total width of capture zone
        if( eirc2.gt.0 .and. eirc1.gt.0 )then
          eirc  = eirc2-eirc1+1
        endif

        ! begin/end points of "capture zone" on my subdomain (if any):
        ! irb denotes beginning i value
        ! note: my subdomain has no "capture" points if irb=0
        wircs1 = 0
        wircs2 = 0
        wirb = 0
        eircs1 = 0
        eircs2 = 0
        eirb = 0
        ! loop through my subdomain:
        do i=1,ni
          ds = 0.1*(xf(i+1)-xf(i))
          if( xh(i).ge.(recy_cap_w-ds) .and. xh(i).le.(recy_cap_w+recy_width+ds) )then
            if( wircs1.eq.0 )then
              wircs1 = i
 !                    i for total grid    minus beginning pt    + 1
              wirb =      (myi1-1+i)       - wirc1                + 1
            endif
            wircs2 = i
          endif
          if( xh(i).le.(recy_cap_e+ds) .and. xh(i).ge.(recy_cap_e-recy_width-ds) )then
            if( eircs1.eq.0 )then
              eircs1 = i
 !                    i for total grid    minus beginning pt    + 1
              eirb =      (myi1-1+i)       - eirc1               + 1
            endif
            eircs2 = i
          endif
        enddo

        wircs  = 0
        ! total width of capture zone on my subdomain
        if( wircs2.gt.0 .and. wircs1.gt.0 )then
          wircs = wircs2-wircs1+1
        endif
        eircs  = 0
        ! total width of capture zone on my subdomain
        if( eircs2.gt.0 .and. eircs1.gt.0 )then
          eircs = eircs2-eircs1+1
        endif

        wuircs1 = 0
        wuircs2 = 0
        wuirb = 0
        euircs1 = 0
        euircs2 = 0
        euirb = 0
        ! loop through my subdomain:
        do i=1,ni+1
          ds = 0.1*(xh(i)-xh(i-1))
          if( xf(i).ge.(recy_cap_w-ds) .and. xf(i).le.(recy_cap_w+recy_width+ds) )then
            if( wuircs1.eq.0 )then
              wuircs1 = i
 !                    i for total grid    minus beginning pt    + 1
              wuirb =      (myi1-1+i)       - wirc1                + 1
            endif
            wuircs2 = i
          endif
          if( xf(i).le.(recy_cap_e+ds) .and. xf(i).ge.(recy_cap_e-recy_width-ds) )then
            if( euircs1.eq.0 )then
              euircs1 = i
 !                    i for total grid    minus beginning pt    + 1
              euirb =      (myi1-1+i)       - eirc1                + 1
            endif
            euircs2 = i
          endif
        enddo

        wuircs = 0
        if( wuircs2.gt.0 .and. wuircs1.gt.0 )then
          wuircs = wuircs2-wuircs1+1
        endif
        euircs = 0
        if( euircs2.gt.0 .and. euircs1.gt.0 )then
          euircs = euircs2-euircs1+1
        endif

      ! ........  injection zone info  ........ !

        ! begin/end points of "injection zone" on total domain:
        wiri1 = 0
        wiri2 = 0
        eiri1 = 0
        eiri2 = 0
        ! loop through entire domain:
        do i=1,nx
          ds = 0.1*(xfref(i+1)-xfref(i))
          if( xhref(i).ge.(recy_inj_w-ds) .and. xhref(i).le.(recy_inj_w+recy_width+ds) )then
            if( wiri1.eq.0 ) wiri1 = i
            wiri2 = i
          endif
          if( xhref(i).le.(recy_inj_e+ds) .and. xhref(i).ge.(recy_inj_e-recy_width-ds) )then
            if( eiri1.eq.0 ) eiri1 = i
            eiri2 = i
          endif
        enddo

        wiri  = 0
        ! total width of injection zone
        if( wiri2.gt.0 .and. wiri1.gt.0 )then
          wiri  = wiri2-wiri1+1
        endif
        eiri  = 0
        ! total width of injection zone
        if( eiri2.gt.0 .and. eiri1.gt.0 )then
          eiri  = eiri2-eiri1+1
        endif

        ! begin/end points of "injection zone" on my subdomain (if any):

        wiris1 = 0
        wiris2 = 0
        wirisb = 0   ! (my beginning value of recycle grid)
        eiris1 = 0
        eiris2 = 0
        eirisb = 0   ! (my beginning value of recycle grid)

       ! loop through my subdomain:
        do i=1,ni
          ds = 0.1*(xf(i+1)-xf(i))
          if( xh(i).ge.(recy_inj_w-ds) .and. xh(i).le.(recy_inj_w+recy_width+ds) )then
            if( wiris1.eq.0 )then
              wiris1 = i
 !                    i for total grid    minus beginning pt    + 1
              wirisb =     (myi1-1+i)       - wiri1                + 1
            endif
            wiris2 = i
          endif
          if( xh(i).le.(recy_inj_e+ds) .and. xh(i).ge.(recy_inj_e-recy_width-ds) )then
            if( eiris1.eq.0 )then
              eiris1 = i
 !                    i for total grid    minus beginning pt    + 1
              eirisb =     (myi1-1+i)       - eiri1                + 1
            endif
            eiris2 = i
          endif
        enddo

        wiris  = 0
        ! total width of injection zone on my subdomain
        if( wiris2.gt.0 .and. wiris1.gt.0 )then
          wiris = wiris2-wiris1+1
        endif
        eiris  = 0
        ! total width of injection zone on my subdomain
        if( eiris2.gt.0 .and. eiris1.gt.0 )then
          eiris = eiris2-eiris1+1
        endif

        wuiris1 = 0
        wuiris2 = 0
        wuirisb = 0   ! (my beginning value of recycle grid)
        euiris1 = 0
        euiris2 = 0
        euirisb = 0   ! (my beginning value of recycle grid)

       ! loop through my subdomain:
        do i=1,ni+1
          ds = 0.1*(xh(i)-xh(i-1))
          if( xf(i).ge.(recy_inj_w-ds) .and. xf(i).le.(recy_inj_w+recy_width+ds) )then
            if( wuiris1.eq.0 )then
              wuiris1 = i
 !                    i for total grid    minus beginning pt    + 1
              wuirisb =     (myi1-1+i)       - wiri1                + 1
            endif
            wuiris2 = i
          endif
          if( xf(i).le.(recy_inj_e+ds) .and. xf(i).ge.(recy_inj_e-recy_width-ds) )then
            if( euiris1.eq.0 )then
              euiris1 = i
 !                    i for total grid    minus beginning pt    + 1
              euirisb =     (myi1-1+i)       - eiri1                + 1
            endif
            euiris2 = i
          endif
        enddo

        wuiris = 0
        if( wuiris2.gt.0 .and. wuiris1.gt.0 )then
          wuiris = wuiris2-wuiris1+1
        endif
        euiris = 0
        if( euiris2.gt.0 .and. euiris1.gt.0 )then
          euiris = euiris2-euiris1+1
        endif

      !.............

        if( wiri .ne. wirc )then
          print *,'  wiri,wirc = ',wiri,wirc
          stop 13987
        endif
        if( eiri .ne. eirc )then
          print *,'  eiri,eirc = ',eiri,eirc
          stop 13988
        endif
        if( wiri .ne. eiri )then
          print *,'  wiri,eiri = ',wiri,eiri
          stop 13989
        endif

        ! this code assumes a S-N slice over whole domain:
        ! add 1 to account for staggered vars
        irecywe = wiri  + 1
        jrecywe =  nj   + 1

      ENDIF  dorecyew


      !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
      !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
      ! south&north sides of domain:


      dorecysn:  &
      IF( do_recycle_s .or. do_recycle_n )THEN

      ! ........  capture zone info  ........ !

        ! begin/end points of "capture zone" on total domain:
        sjrc1 = 0
        sjrc2 = 0
        njrc1 = 0
        njrc2 = 0
       ! loop through entire domain:
        do j=1,ny
          ds = 0.1*(yfref(j+1)-yfref(j))
          if( yhref(j).ge.(recy_cap_s-ds) .and. yhref(j).le.(recy_cap_s+recy_width+ds) )then
            if( sjrc1.eq.0 ) sjrc1 = j
            sjrc2 = j
          endif
          if( yhref(j).le.(recy_cap_n+ds) .and. yhref(j).ge.(recy_cap_n-recy_width-ds) )then
            if( njrc1.eq.0 ) njrc1 = j
            njrc2 = j
          endif
        enddo

        sjrc  = 0
        ! total sjdth of capture zone
        if( sjrc2.gt.0 .and. sjrc1.gt.0 )then
          sjrc  = sjrc2-sjrc1+1
        endif
        njrc  = 0
        ! total width of capture zone
        if( njrc2.gt.0 .and. njrc1.gt.0 )then
          njrc  = njrc2-njrc1+1
        endif

        ! begin/end points of "capture zone" on my subdomain (if any):
        ! irb denotes beginning i value
        ! note: my subdomain has no "capture" points if irb=0
        sjrcs1 = 0
        sjrcs2 = 0
        sjrb = 0
        njrcs1 = 0
        njrcs2 = 0
        njrb = 0
        ! loop through my subdomain:
        do j=1,nj
          ds = 0.1*(yf(j+1)-yf(j))
          if( yh(j).ge.(recy_cap_s-ds) .and. yh(j).le.(recy_cap_s+recy_width+ds) )then
            if( sjrcs1.eq.0 )then
              sjrcs1 = j
 !                    j for total grid    minus beginning pt    + 1
              sjrb =      (myj1-1+j)       - sjrc1                + 1
            endif
            sjrcs2 = j
          endif
          if( yh(j).le.(recy_cap_n+ds) .and. yh(j).ge.(recy_cap_n-recy_width-ds) )then
            if( njrcs1.eq.0 )then
              njrcs1 = j
 !                    j for total grid    minus beginning pt    + 1
              njrb =      (myj1-1+j)       - njrc1               + 1
            endif
            njrcs2 = j
          endif
        enddo

        sjrcs  = 0
        ! total width of capture zone on my subdomain
        if( sjrcs2.gt.0 .and. sjrcs1.gt.0 )then
          sjrcs = sjrcs2-sjrcs1+1
        endif
        njrcs  = 0
        ! total width of capture zone on my subdomain
        if( njrcs2.gt.0 .and. njrcs1.gt.0 )then
          njrcs = njrcs2-njrcs1+1
        endif

        svjrcs1 = 0
        svjrcs2 = 0
        svjrb = 0
        nvjrcs1 = 0
        nvjrcs2 = 0
        nvjrb = 0
        ! loop through my subdomain:
        do j=1,nj+1
          ds = 0.1*(yh(j)-yh(j-1))
          if( yf(j).ge.(recy_cap_s-ds) .and. yf(j).le.(recy_cap_s+recy_width+ds) )then
            if( svjrcs1.eq.0 )then
              svjrcs1 = j
 !                    j for total grid    minus beginning pt    + 1
              svjrb =      (myj1-1+j)       - sjrc1                + 1
            endif
            svjrcs2 = j
          endif
          if( yf(j).le.(recy_cap_n+ds) .and. yf(j).ge.(recy_cap_n-recy_width-ds) )then
            if( nvjrcs1.eq.0 )then
              nvjrcs1 = j
 !                    j for total grid    minus beginning pt    + 1
              nvjrb =      (myj1-1+j)       - njrc1                + 1
            endif
            nvjrcs2 = j
          endif
        enddo

        svjrcs = 0
        if( svjrcs2.gt.0 .and. svjrcs1.gt.0 )then
          svjrcs = svjrcs2-svjrcs1+1
        endif
        nvjrcs = 0
        if( nvjrcs2.gt.0 .and. nvjrcs1.gt.0 )then
          nvjrcs = nvjrcs2-nvjrcs1+1
        endif

      ! ........  injection zone info  ........ !

        ! begin/end points of "injection zone" on total domain:
        sjri1 = 0
        sjri2 = 0
        njri1 = 0
        njri2 = 0
        ! loop through entire domain:
        do j=1,ny
          ds = 0.1*(yfref(j+1)-yfref(j))
          if( yhref(j).ge.(recy_inj_s-ds) .and. yhref(j).le.(recy_inj_s+recy_width+ds) )then
            if( sjri1.eq.0 ) sjri1 = j
            sjri2 = j
          endif
          if( yhref(j).le.(recy_inj_n+ds) .and. yhref(j).ge.(recy_inj_n-recy_width-ds) )then
            if( njri1.eq.0 ) njri1 = j
            njri2 = j
          endif
        enddo

        sjri  = 0
        ! total width of injection zone
        if( sjri2.gt.0 .and. sjri1.gt.0 )then
          sjri  = sjri2-sjri1+1
        endif
        njri  = 0
        ! total width of injection zone
        if( njri2.gt.0 .and. njri1.gt.0 )then
          njri  = njri2-njri1+1
        endif

        ! begin/end points of "injection zone" on my subdomain (if any):

        sjris1 = 0
        sjris2 = 0
        sjrisb = 0   ! (my beginning value of recycle grid)
        njris1 = 0
        njris2 = 0
        njrisb = 0   ! (my beginning value of recycle grid)

       ! loop through my subdomain:
        do j=1,nj
          ds = 0.1*(yf(j+1)-yf(j))
          if( yh(j).ge.(recy_inj_s-ds) .and. yh(j).le.(recy_inj_s+recy_width+ds) )then
            if( sjris1.eq.0 )then
              sjris1 = j
 !                    j for total grid    minus beginning pt    + 1
              sjrisb =     (myj1-1+j)       - sjri1                + 1
            endif
            sjris2 = j
          endif
          if( yh(j).le.(recy_inj_n+ds) .and. yh(j).ge.(recy_inj_n-recy_width-ds) )then
            if( njris1.eq.0 )then
              njris1 = j
 !                    j for total grid    minus beginning pt    + 1
              njrisb =     (myj1-1+j)       - njri1                + 1
            endif
            njris2 = j
          endif
        enddo

        sjris  = 0
        ! total width of injection zone on my subdomain
        if( sjris2.gt.0 .and. sjris1.gt.0 )then
          sjris = sjris2-sjris1+1
        endif
        njris  = 0
        ! total width of injection zone on my subdomain
        if( njris2.gt.0 .and. njris1.gt.0 )then
          njris = njris2-njris1+1
        endif

        svjris1 = 0
        svjris2 = 0
        svjrisb = 0   ! (my beginning value of recycle grid)
        nvjris1 = 0
        nvjris2 = 0
        nvjrisb = 0   ! (my beginning value of recycle grid)

       ! loop through my subdomain:
        do j=1,nj+1
          ds = 0.1*(yh(j)-yh(j-1))
          if( yf(j).ge.(recy_inj_s-ds) .and. yf(j).le.(recy_inj_s+recy_width+ds) )then
            if( svjris1.eq.0 )then
              svjris1 = j
 !                    j for total grid    minus beginning pt    + 1
              svjrisb =     (myj1-1+j)       - sjri1                + 1
            endif
            svjris2 = j
          endif
          if( yf(j).le.(recy_inj_n+ds) .and. yf(j).ge.(recy_inj_n-recy_width-ds) )then
            if( nvjris1.eq.0 )then
              nvjris1 = j
 !                    j for total grid    minus beginning pt    + 1
              nvjrisb =     (myj1-1+j)       - njri1                + 1
            endif
            nvjris2 = j
          endif
        enddo

        svjris = 0
        if( svjris2.gt.0 .and. svjris1.gt.0 )then
          svjris = svjris2-svjris1+1
        endif
        nvjris = 0
        if( nvjris2.gt.0 .and. nvjris1.gt.0 )then
          nvjris = nvjris2-nvjris1+1
        endif

      !.............

        if( sjri .ne. sjrc )then
          print *,'  sjri,sjrc = ',sjri,sjrc
          stop 13987
        endif
        if( njri .ne. njrc )then
          print *,'  njri,njrc = ',njri,njrc
          stop 13988
        endif
        if( sjri .ne. njri )then
          print *,'  sjri,njri = ',sjri,njri
          stop 13989
        endif

        ! this code assumes a S-N slice over whole domain:
        ! add 1 to account for staggered vars
        irecysn =  ni   + 1
        jrecysn = sjri  + 1

      ENDIF  dorecysn


      !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
      !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c


          ! print info to file:

            allocate( foor(10) )
            foor = 0.0
            allocate( fooi(125) )
            fooi = 0

            foor(1)  = recy_width
            foor(2)  = recy_depth
            foor(3)  = recy_cap_w
            foor(4)  = recy_cap_e
            foor(5)  = recy_cap_s
            foor(6)  = recy_cap_n

            fooi(  1) = wirc1
            fooi(  2) = wirc2
            fooi(  3) = wirc
            fooi(  4) = wircs1
            fooi(  5) = wircs2
            fooi(  6) = wircs
            fooi(  7) = wirb
            fooi(  8) = wuircs1
            fooi(  9) = wuircs2
            fooi( 10) = wuircs
            fooi( 11) = wuirb
            fooi( 12) = wiri1
            fooi( 13) = wiri2
            fooi( 14) = wiri
            fooi( 15) = wiris1
            fooi( 16) = wiris2
            fooi( 17) = wiris
            fooi( 18) = wirisb
            fooi( 19) = wuiris1
            fooi( 20) = wuiris2
            fooi( 21) = wuiris
            fooi( 22) = wuirisb

            fooi( 31) = eirc1
            fooi( 32) = eirc2
            fooi( 33) = eirc
            fooi( 34) = eircs1
            fooi( 35) = eircs2
            fooi( 36) = eircs
            fooi( 37) = eirb
            fooi( 38) = euircs1
            fooi( 39) = euircs2
            fooi( 40) = euircs
            fooi( 41) = euirb
            fooi( 42) = eiri1
            fooi( 43) = eiri2
            fooi( 44) = eiri
            fooi( 45) = eiris1
            fooi( 46) = eiris2
            fooi( 47) = eiris
            fooi( 48) = eirisb
            fooi( 49) = euiris1
            fooi( 50) = euiris2
            fooi( 51) = euiris
            fooi( 52) = euirisb

            fooi( 61) = sjrc1
            fooi( 62) = sjrc2
            fooi( 63) = sjrc
            fooi( 64) = sjrcs1
            fooi( 65) = sjrcs2
            fooi( 66) = sjrcs
            fooi( 67) = sjrb
            fooi( 68) = svjrcs1
            fooi( 69) = svjrcs2
            fooi( 70) = svjrcs
            fooi( 71) = svjrb
            fooi( 72) = sjri1
            fooi( 73) = sjri2
            fooi( 74) = sjri
            fooi( 75) = sjris1
            fooi( 76) = sjris2
            fooi( 77) = sjris
            fooi( 78) = sjrisb
            fooi( 79) = svjris1
            fooi( 80) = svjris2
            fooi( 81) = svjris
            fooi( 82) = svjrisb

            fooi( 91) = sjrc1
            fooi( 92) = sjrc2
            fooi( 93) = sjrc
            fooi( 94) = sjrcs1
            fooi( 95) = sjrcs2
            fooi( 96) = sjrcs
            fooi( 97) = sjrb
            fooi( 98) = svjrcs1
            fooi( 99) = svjrcs2
            fooi(100) = svjrcs
            fooi(101) = svjrb
            fooi(102) = sjri1
            fooi(103) = sjri2
            fooi(104) = sjri
            fooi(105) = sjris1
            fooi(106) = sjris2
            fooi(107) = sjris
            fooi(108) = sjrisb
            fooi(109) = svjris1
            fooi(110) = svjris2
            fooi(111) = svjris
            fooi(112) = svjrisb

            fooi(121) = ni
            fooi(122) = nj
            fooi(123) = myi
            fooi(124) = myj

        IF( myid.eq.0 )THEN

          open(unit=10,file='recycle.info',status='unknown')

        do n = 0 , (numprocs-1)

#ifdef MPI
          if( n.ne.0 )then
            call mpi_recv(foor  ,10 ,MPI_REAL   ,n,8001,MPI_COMM_WORLD,status,ierr)
            call mpi_recv(fooi  ,125,MPI_INTEGER,n,8003,MPI_COMM_WORLD,status,ierr)
          endif
#endif

          write(10,*) ' --------------------------- '
          write(10,*) '  proc,ni,nj = ',n,fooi(121),fooi(122)
          write(10,*) '  myi,myj    = ',fooi(123),fooi(124)
          write(10,*) 
          write(10,*) '  recy_width,recy_depth = ',foor(1),foor(2)
          write(10,*) 
          write(10,*) '  recy_cap_w         = ',foor(3)
          write(10,*) '  wirc1,wirc2,wirc       = ',fooi( 1),fooi( 2),fooi( 3)
          write(10,*) '  wircs1,wircs2,wircs    = ',fooi( 4),fooi( 5),fooi( 6)
          write(10,*) '  wuircs1,wuircs2,wuircs = ',fooi( 8),fooi( 9),fooi(10)
          write(10,*) '  wirb,wuirb             = ',fooi( 7),fooi(11)
          write(10,*) '  wiri1,wiri2,wiri       = ',fooi(12),fooi(13),fooi(14)
          write(10,*) '  wiris1,wiris2,wiris    = ',fooi(15),fooi(16),fooi(17)
          write(10,*) '  wuiris1,wuiris2,wuiris = ',fooi(19),fooi(20),fooi(21)
          write(10,*) '  wirisb,wuirisb         = ',fooi(18),fooi(22)
          write(10,*) 
          write(10,*) '  recy_cap_e          = ',foor(4)
          write(10,*) '  eirc1,eirc2,eirc       = ',fooi(31),fooi(32),fooi(33)
          write(10,*) '  eircs1,eircs2,eircs    = ',fooi(34),fooi(35),fooi(36)
          write(10,*) '  euircs1,euircs2,euircs = ',fooi(38),fooi(39),fooi(40)
          write(10,*) '  eirb,euirb             = ',fooi(37),fooi(41)
          write(10,*) '  eiri1,eiri2,eiri       = ',fooi(42),fooi(43),fooi(44)
          write(10,*) '  eiris1,eiris2,eiris    = ',fooi(45),fooi(46),fooi(47)
          write(10,*) '  euiris1,euiris2,euiris = ',fooi(49),fooi(50),fooi(51)
          write(10,*) '  eirisb,euirisb         = ',fooi(48),fooi(52)
          write(10,*) 
          write(10,*) '  recy_cap_s          = ',foor(5)
          write(10,*) '  sjrc1,sjrc2,sjrc       = ',fooi(61),fooi(62),fooi(63)
          write(10,*) '  sjrcs1,sjrcs2,sjrcs    = ',fooi(64),fooi(65),fooi(66)
          write(10,*) '  svjrcs1,svjrcs2,svjrcs = ',fooi(68),fooi(69),fooi(70)
          write(10,*) '  sjrb,svjrb             = ',fooi(67),fooi(71)
          write(10,*) '  sjri1,sjri2,sjri       = ',fooi(72),fooi(73),fooi(74)
          write(10,*) '  sjris1,sjris2,sjris    = ',fooi(75),fooi(76),fooi(77)
          write(10,*) '  svjris1,svjris2,svjris = ',fooi(79),fooi(80),fooi(81)
          write(10,*) '  sjrisb,svjrisb         = ',fooi(78),fooi(82)
          write(10,*) 
          write(10,*) '  recy_cap_n          = ',foor(6)
          write(10,*) '  sjrc1,sjrc2,sjrc       = ',fooi(91),fooi(92),fooi(93)
          write(10,*) '  sjrcs1,sjrcs2,sjrcs    = ',fooi(94),fooi(95),fooi(96)
          write(10,*) '  svjrcs1,svjrcs2,svjrcs = ',fooi(98),fooi(99),fooi(100)
          write(10,*) '  sjrb,svjrb             = ',fooi(97),fooi(101)
          write(10,*) '  sjri1,sjri2,sjri       = ',fooi(102),fooi(103),fooi(104)
          write(10,*) '  sjris1,sjris2,sjris    = ',fooi(105),fooi(106),fooi(107)
          write(10,*) '  svjris1,svjris2,svjris = ',fooi(109),fooi(110),fooi(111)
          write(10,*) '  sjrisb,svjrisb         = ',fooi(108),fooi(112)
          write(10,*) 
        enddo

          write(10,*) ' --------------------------- '

          close(unit=10)

#ifdef MPI
        ELSE
            call mpi_send(foor  ,10 ,MPI_REAL   ,0,8001,MPI_COMM_WORLD,ierr)
            call mpi_send(fooi  ,125,MPI_INTEGER,0,8003,MPI_COMM_WORLD,ierr)
#endif
        ENDIF

        deallocate( foor )
        deallocate( fooi )


      !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
      !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c

    end subroutine eddy_recycling_indices

!-----------------------------------------------------------------------

  end module eddy_recycle
