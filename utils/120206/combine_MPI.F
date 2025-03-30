    program newcombine
    implicit none

    character*60 fname
    character*6 type
    character*8 var
    character*30 description
    integer :: n,nx,ny,nz,nt,nv,nlev,num2d,num3d,nodex,nodey
    integer :: nstart,nend,ninc,mix,mjx,mkx,nfiles,irec
    integer, dimension(:), allocatable :: nrec
    real :: lon,lon1,lon2,lat,lat1,lat2,lev,lev1,lev2

!------

    write(6,*)
    write(6,*) '  Enter nodex .... '
    write(6,*)
    read(5,*) nodex

    write(6,*)
    write(6,*) '  Enter nodey .... '
    write(6,*)
    read(5,*) nodey

    write(6,*)
    write(6,*) '  Enter name of descriptor file (e.g., cm1out_s.ctl) .... '
    write(6,*)
    read(5,*) fname

!------

    write(6,*)
    write(6,*) '  Enter first time to process ...'
    write(6,*)
    read(5,*) nstart

    write(6,*)
    write(6,*) '  Enter last time to process ...'
    write(6,*)
    read(5,*) nend

    write(6,*)
    write(6,*) '  Enter time increment ...'
    write(6,*)
    read(5,*) ninc

!------

    write(6,*)
    write(6,*) ' nstart,nend,ninc=',nstart,nend,ninc
    write(6,*) ' fname = ',fname
    write(6,*)
    write(6,*) '  Processing ...'
    write(6,*)

    open(unit=20,file=fname,status='old')
    open(unit=30,file='cm1out_MPI.ctl')

    read(20,*)
    write(30,201)
201 format('dset ^cm1out_%y4.dat')

!!!    read(20,*)
    write(30,101)
101 format('options template')

    read(20,*)
    write(30,202)
202 format('title CM1 output')

    read(20,*)
    write(30,203)
203 format('undef -99999999.')

    read(20,104) nx,type
104 format(5x,i6,1x,a6)
    if(type.eq.'linear')then
      backspace(20)
      read(20,204) lon1,lon2
204   format(19x,f13.6,1x,f13.6)
      write(30,304) nx,lon1,lon2
304   format('xdef ',i6,' linear ',f13.6,1x,f13.6)
    elseif(type.eq.'levels')then
      write(30,314) nx
314   format('xdef ',i6,' levels ')
      do n=1,nx
        read(20,217) lon
        write(30,217) lon
217     format(2x,f13.6)
      enddo
    else
      print *,'  unknown type'
      stop 11111
    endif

    read(20,105) ny,type
105 format(5x,i6,1x,a6)
    if(type.eq.'linear')then
      backspace(20)
      read(20,205) lat1,lat2
205   format(19x,f13.6,1x,f13.6)
      write(30,305) ny,lat1,lat2
305   format('ydef ',i6,' linear ',f13.6,1x,f13.6)
    elseif(type.eq.'levels')then
      write(30,315) ny
315   format('ydef ',i6,' levels ')
      do n=1,ny
        read(20,217) lat
        write(30,217) lat
      enddo
    else
      print *,'  unknown type'
      stop 11111
    endif

    read(20,106) nz,type
106 format(5x,i6,1x,a6)
    if(type.eq.'linear')then
      backspace(20)
      read(20,206) lev1,lev2
206   format(19x,f13.6,1x,f13.6)
      write(30,306) nz,lev1,lev2
306   format('zdef ',i6,' linear ',f13.6,1x,f13.6)
    elseif(type.eq.'levels')then
      write(30,316) nz
316   format('zdef ',i6,' levels ')
      do n=1,nz
        read(20,217) lev
        write(30,217) lev
      enddo
    else
      print *,'  unknown type'
      stop 11111
    endif

    read(20,107) nt
107 format(5x,i10)
    write(30,207) nt
207 format('tdef ',i10,' linear 00Z01JAN0001 1YR')

    read(20,108) nv
108 format(5x,i4)
    write(30,208) nv
208 format('vars ',i4)

    num2d = 0
    num3d = 0

    DO n=1,nv

      read(20,209) var,nlev,description
      write(30,209) var,nlev,description
209   format(a8,2x,i6,'  99  ',a30)

      if(nlev.eq.0)then
        num2d = num2d + 1
      else
        num3d = num3d + 1
      endif

    ENDDO

    write(30,210)
210 format('endvars')

    close(unit=20)
    close(unit=30)

!-----------------------------------------------------------------------

    print *
    print *,'  nodex,nodey = ',nodex,nodey
    print *,'  nx,ny,nz    = ',nx,ny,nz
    print *,'  num2d,num3d = ',num2d,num3d
    print *

    mix=nx/nodex
    mjx=ny/nodey
    mkx=nz

    nfiles=nodex*nodey

    allocate( nrec(nfiles) )

!-----------------------------------------------------------------------

    DO nt=nstart,nend,ninc

      print *
      print *
      print *,'  nt=',nt

      irec=0

      do n=1,nfiles
        nrec(n)=0
      enddo

      call readwrite(mix,mjx,mkx,nx,ny,nz,nodex,nodey,nfiles,irec,nrec,   &
                     1,mkx,nt,num2d,num3d)

    ENDDO

!-----------------------------------------------------------------------

    stop 99999
    end program newcombine


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine readwrite(                                           &
                 mix,mjx,mkx,it,jt,kt,nodex,nodey,nfiles,irec,nrec,   &
                 k1,k2,nt,num2d,num3d)
      implicit none

      integer mix,mjx,mkx,it,jt,kt,nodex,nodey,nfiles,irec,nrec(nfiles)
      integer orec
      integer k1,k2,nt,n,myi,myj
      integer num2d,num3d

      integer i,j,k,nx,ny,np
      integer ifoo,i1,i2
      character*50 newname
      character*50 fname
      real var

      newname='cm1out_XXXX.dat'
      write(newname(8:11),101) nt
101   format(i4.4)
102   format(i3.3)

      open(unit=57,file=newname,                 &
           form='unformatted',access='direct',   &
           recl=4,status='unknown')

      do ny=1,nodey
      do nx=1,nodex
        np=nx+(ny-1)*nodex-1
        myj = np / nodex + 1
        myi = np - (myj-1)*nodex  + 1
        ifoo=(myi-1)*mix
          i1=1
          i2=mix
          print *,'    ifoo,i1,i2=',ifoo,i1,i2
          fname="cm1out_XXXX_YYYY_s.dat"
          write(fname(8:11),100) np
100       format(i4.4)
          write(fname(13:16),100) nt
          print *,'  np=',np
          open(unit=60,file=fname,                   &
               form='unformatted',access='direct',   &
               recl=4,status='old')
          irec=1
          do n=1,num2d
          do j=1,mjx
            irec=1+(mix)*(mjx)*(n-1)   &
                  +(j-1)*mix           &
                  +(i1-1)
            orec=1+it*(mjx*nodey)*(n-1)   &
                  +(myj-1)*it*mjx         &
                  +(j-1)*it               &
                  +max( (ifoo),0 )
            do i=i1,i2
              read(60,rec=irec) var
              irec=irec+1
              write(57,rec=orec) var
              orec=orec+1
            enddo
          enddo
          enddo
          do n=1,num3d
          do k=1,mkx
          do j=1,mjx
            irec=1+(mix)*(mjx)*(num2d+(n-1)*mkx+(k-1))   &
                  +(j-1)*mix                             &
                  +(i1-1)
            orec=1+it*(mjx*nodey)*(num2d+(n-1)*mkx+(k-1))   &
                  +(myj-1)*it*mjx                           &
                  +(j-1)*it                                 &
                  +max( (ifoo),0 )
            do i=i1,i2
              read(60,rec=irec) var
              irec=irec+1
              write(57,rec=orec) var
              orec=orec+1
            enddo
          enddo
          enddo
          enddo
      enddo
      enddo

      close(unit=57)

      return
      end subroutine readwrite
