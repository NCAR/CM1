    program newcombine
    implicit none

    character*60 fname
    character*6 type
    character*8 var
    character*30 description
    integer :: n,ni,nj,nx,ny,nz,nt,nv,nlev,num2d,num3d,nodex,nodey
    integer :: nstart,nend,ninc,nfiles
    real :: lon,lon1,lon2,lat,lat1,lat2,lev,lev1,lev2
    integer, dimension(:), allocatable :: nrec
    real, dimension(:,:), allocatable :: dat1,dat2

    integer i,j,k,np,fnum
    character*60 newname
    integer :: myi,myj,orec,fooi,fooj

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
    write(6,*) '  Enter name of descriptor file (e.g., cm1out_w.ctl) .... '
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
    open(unit=30,file='cm1out_MPI_w.ctl')

    read(20,*)
    write(30,201)
201 format('dset ^cm1out_00%y4_w.dat')

    read(20,*)
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

    ni=nx/nodex
    nj=ny/nodey

    print *
    print *,'  nodex,nodey = ',nodex,nodey
    print *,'  nx,ny,nz    = ',nx,ny,nz
    print *,'  ni,nj       = ',ni,nj
    print *,'  num2d,num3d = ',num2d,num3d
    print *

    nfiles=nodex*nodey

    allocate( nrec(nfiles) )
    allocate( dat1(ni,nj) )
    allocate( dat2(nx,ny) )

    dat1 = 0.0
    dat2 = 0.0

!-----------------------------------------------------------------------

    DO nt=nstart,nend,ninc

      print *
      print *
      print *,'  nt=',nt

      newname='cm1out_XXXXXX_w.dat'
      write(newname(8:13),401) nt
401   format(i6.6)

      open(unit=57,file=newname,                 &
           form='unformatted',access='direct',   &
           recl=4*nx*ny,status='unknown')
      orec = 1

      DO np=0,nodex*nodey-1
        fnum = 100+np
        fname="cm1out_XXXXXX_YYYYYY_w.dat"
        write(fname( 8:13),401) np
        write(fname(15:20),401) nt
!!!        print *,'  np,fnum,fname=',np,fnum,' ',fname
        open(unit=fnum,file=fname,                 &
             form='unformatted',access='direct',   &
             recl=4*ni*nj,status='old')
      ENDDO

      do n=1,num2d
      print *,'  2d = ',n
        do np=0,nodex*nodey-1
          fnum = 100+np
          read(fnum,rec=orec) ((dat1(i,j),i=1,ni),j=1,nj)
          fooj = np / nodex + 1
          fooi = np - (fooj-1)*nodex  + 1
          do j=1,nj
          do i=1,ni
            dat2((fooi-1)*ni+i,(fooj-1)*nj+j)=dat1(i,j)
          enddo
          enddo
        enddo
        write(57,rec=orec) ((dat2(i,j),i=1,nx),j=1,ny)
        orec = orec + 1
      enddo

      do n=1,num3d
      print *,'  3d = ',n
      do k=1,nz
        do np=0,nodex*nodey-1
          fnum = 100+np
          read(fnum,rec=orec) ((dat1(i,j),i=1,ni),j=1,nj)
          fooj = np / nodex + 1
          fooi = np - (fooj-1)*nodex  + 1
          do j=1,nj
          do i=1,ni
            dat2((fooi-1)*ni+i,(fooj-1)*nj+j)=dat1(i,j)
          enddo
          enddo
        enddo
        write(57,rec=orec) ((dat2(i,j),i=1,nx),j=1,ny)
        orec = orec + 1
      enddo
      enddo

      close(unit=57)

    ENDDO

!-----------------------------------------------------------------------

    stop 99999
    end program newcombine
