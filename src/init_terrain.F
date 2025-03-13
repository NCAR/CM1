  MODULE init_terrain_module

  implicit none

  CONTAINS

      subroutine init_terrain(xh,uh,xf,uf,yh,vh,yf,vf,rds,sigma,rdsf,sigmaf,  &
                              zh,zf,zs,gz,rgz,gzu,rgzu,gzv,rgzv,         &
                              dzdx,dzdy,gx,gxu,gy,gyv,                   &
                              reqs_u,reqs_v,reqs_s,reqs_p,               &
                              nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,           &
                              sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,   &
                              uw31,uw32,ue31,ue32,us31,us32,un31,un32,   &
                              vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,   &
                              west,newwest,east,neweast,                 &
                              south,newsouth,north,newnorth)

      use input
      use constants
      use bc_module
      use adv_routines , only : zsgrad
      use comm_module
#ifdef MPI
      use mpi
#endif
      implicit none

      real, intent(in), dimension(ib:ie) :: xh,uh
      real, intent(in), dimension(ib:ie+1) :: xf,uf
      real, intent(in), dimension(jb:je) :: yh,vh
      real, intent(in), dimension(jb:je+1) :: yf,vf
      real, intent(inout), dimension(kb:ke) :: rds,sigma
      real, intent(inout), dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: zf
      real, intent(inout), dimension(ib:ie,jb:je) :: zs
      real, intent(inout), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy
      real, intent(inout), dimension(itb:ite,jtb:jte,ktb:kte) :: gx,gxu,gy,gyv
      integer, intent(inout), dimension(rmp) :: reqs_u,reqs_v,reqs_s,reqs_p
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      real, intent(inout), dimension(cmp,jmp,kmp)   :: sw31,sw32,se31,se32
      real, intent(inout), dimension(imp,cmp,kmp)   :: ss31,ss32,sn31,sn32
      real, intent(inout), dimension(cmp,jmp,kmp)   :: uw31,uw32,ue31,ue32
      real, intent(inout), dimension(imp+1,cmp,kmp) :: us31,us32,un31,un32
      real, intent(inout), dimension(cmp,jmp+1,kmp) :: vw31,vw32,ve31,ve32
      real, intent(inout), dimension(imp,cmp,kmp)   :: vs31,vs32,vn31,vn32
      real, intent(inout), dimension(cmp,jmp) :: west,newwest,east,neweast
      real, intent(inout), dimension(imp,cmp) :: south,newsouth,north,newnorth

      integer :: i,j,k,irec,ii,jj
      real :: hh,aa,xval,xc,yc
      real :: tem1,tem2,rr,angle
      real, dimension(:,:), allocatable :: zs_all

#ifdef MPI
      character(len=50) :: fname
#endif

!-----------------------------------------------------------------------
!     SPECIFY TERRAIN HERE
!-----------------------------------------------------------------------

!----------------------------------------------------------
!  itern = 1
!  bell-shaped  (eg, Durran and Klemp, 1983, MWR, pg 2341)

        IF(itern.eq.1)THEN

          hh =      400.0              ! max. height (m)
          aa =     1000.0              ! half width (m)
          xc =        0.0 + 0.5*dx     ! x-location (m)

          do j=jb,je
          do i=ib,ie
            zs(i,j)=hh/( 1.0+( (xh(i)-xc)/aa )**2 )
          enddo
          enddo

!---------------
!  itern = 2
!  Schaer case  (eg, Klemp et al, 2003, MWR, pg 1229)

        ELSEIF(itern.eq.2)THEN

          do j=jb,je
          do i=ib,ie
            xval=dx*(i-ni/2)
            zs(i,j)=250.0*exp(-(xval/5000.0)**2)*(cos(pi*xval/4000.0)**2)
          enddo
          enddo

!---------------
!  itern = 3
!  Doernbrack et al, 2005, ASL, pg 118

        ELSEIF(itern.eq.3)THEN

          hh =      500.0     ! max. height (m)
          aa =    20000.0     ! half width (m)

          do j=jb,je
          do i=ib,ie
            xval = sqrt( (xh(i)-0.25*maxx)**2   &
                        +(yh(j)-0.25*maxy)**2   &
                                             )
            zs(i,j)=hh*( (1.0+(xval/aa)**2 )**(-1.5) )
          enddo
          enddo

!----------------------------------------------------------
!  itern = 4
!  read from GrADS file "perts.dat"

        ELSEIF(itern.eq.4)THEN

          ! 170313:  zs_all array stores terrain for the entire domain
          allocate( zs_all(nx,ny) )

          ! only processor 0 reads file:

          IF( myid.eq.0 )THEN

            open(unit=73,file='perts.dat',status='old',   &
                 form='unformatted',access='direct',recl=4*nx*ny)

            read(73,rec=1) ((zs_all(i,j),i=1,nx),j=1,ny)

            close(unit=73)

          ENDIF

#ifdef MPI
          ! communicate data to all processors:
          call MPI_BCAST(zs_all,nx*ny,MPI_REAL,0,MPI_COMM_WORLD,ierr)
#endif

          ! 170313:
          do jj=1,ny
          do ii=1,nx
            i = ii - myi1 + 1
            j = jj - myj1 + 1
            ! check to see if this processor has this gridpoint:
            IF( i.ge.ib .and. i.le.ie .and. j.ge.jb .and. j.le.je )THEN
              zs(i,j) = zs_all(ii,jj)
            ENDIF
          enddo
          enddo

          ! zs_all no longer needed ... deallocate:
          deallocate( zs_all )

!----------------------------------------------------------

        ENDIF

!--------------------------------------------------------------
!  DO NOT CHANGE ANYTHING BELOW HERE !
!--------------------------------------------------------------

        call bc2d(zs)
#ifdef MPI
        nf=0
        nu=0
        nv=0
        nw=0
        call comm_2d_start(zs,west,newwest,east,neweast,   &
                              south,newsouth,north,newnorth,reqs_p)
        call comm_2dew_end(zs,west,newwest,east,neweast,reqs_p)
        call comm_2dns_end(zs,south,newsouth,north,newnorth,reqs_p)
        call bcs2_2d(zs)
        call bc2d(zs)
        call getcorner3_2d(zs)
#endif

        zt = maxz
        rzt = 1.0/maxz

        if(dowr) write(outfile,*)
        do k=1,nk+1
          if(dowr) write(outfile,*) '  sigmaf:',k,sigmaf(k)
        enddo
        if(dowr) write(outfile,*)

        do k=1,nk
        do j=jb,je
        do i=ib,ie
          zh(i,j,k)=zs(i,j)+sigma(k)*(zt-zs(i,j))/zt
        enddo
        enddo
        enddo

        do k=kb,ke+1
        do j=jb,je
        do i=ib,ie
          zf(i,j,k)=zs(i,j)+sigmaf(k)*(zt-zs(i,j))/zt
        enddo
        enddo
        enddo

        call zsgrad(zs,gz,dzdx,dzdy,uh,vh)

!--------------------------------
!  set boundary points

        call bc2d(dzdx)
        call bc2d(dzdy)

#ifdef MPI
        call comm_2d_start(dzdx,west,newwest,east,neweast,   &
                                south,newsouth,north,newnorth,reqs_p)
        call comm_2dew_end(dzdx,west,newwest,east,neweast,reqs_p)
        call comm_2dns_end(dzdx,south,newsouth,north,newnorth,reqs_p)
        call bcs2_2d(dzdx)
        call bc2d(dzdx)
        call getcorner3_2d(dzdx)

        call comm_2d_start(dzdy,west,newwest,east,neweast,   &
                                south,newsouth,north,newnorth,reqs_p)
        call comm_2dew_end(dzdy,west,newwest,east,neweast,reqs_p)
        call comm_2dns_end(dzdy,south,newsouth,north,newnorth,reqs_p)
        call bcs2_2d(dzdy)
        call bc2d(dzdy)
        call getcorner3_2d(dzdy)
#endif

!--------------------------------
!  the metric terms:

        do j=jb,je
        do i=ib,ie
           gz(i,j)=zt/(zt-zs(i,j))
          rgz(i,j)=(zt-zs(i,j))/zt
        enddo
        enddo

        call bc2d(rgz)
#ifdef MPI
        call comm_2d_start(rgz,west,newwest,east,neweast,   &
                               south,newsouth,north,newnorth,reqs_p)
        call comm_2dew_end(rgz,west,newwest,east,neweast,reqs_p)
        call comm_2dns_end(rgz,south,newsouth,north,newnorth,reqs_p)
        call bcs2_2d(rgz)
        call bc2d(rgz)
        call getcorner3_2d(rgz)
#endif

        do j=jb+1,je
        do i=ib+1,ie
           gzu(i,j)=zt/(zt-0.5*(zs(i-1,j)+zs(i,j)))
          rgzu(i,j)=(zt-0.5*(zs(i-1,j)+zs(i,j)))/zt
           gzv(i,j)=zt/(zt-0.5*(zs(i,j-1)+zs(i,j)))
          rgzv(i,j)=(zt-0.5*(zs(i,j-1)+zs(i,j)))/zt
        enddo
        enddo

        do k=1,nk+1
        do j=jb+1,je-1
        do i=ib+1,ie-1
          gx(i,j,k)=(zt-sigmaf(k))*gz(i,j)*(rgzu(i+1,j)-rgzu(i,j))*rdx*uh(i)
          gxu(i,j,k)=(zt-sigmaf(k))*gzu(i,j)*(rgz(i,j)-rgz(i-1,j))*rdx*uf(i)
          gy(i,j,k)=(zt-sigmaf(k))*gz(i,j)*(rgzv(i,j+1)-rgzv(i,j))*rdy*vh(j)
          gyv(i,j,k)=(zt-sigmaf(k))*gzv(i,j)*(rgz(i,j)-rgz(i,j-1))*rdy*vf(j)
        enddo
        enddo
        enddo

!--------------------------------

        do j=jb,je
        do i=ib,ie
          zf(i,j,0)=zf(i,j,1)-(zf(i,j,2)-zf(i,j,1))
          zf(i,j,nk+2)=zf(i,j,nk+1)+(zf(i,j,nk+1)-zf(i,j,nk))
          zh(i,j,0)=0.5*(zf(i,j,0)+zf(i,j,1))
          zh(i,j,nk+1)=0.5*(zf(i,j,nk+1)+zf(i,j,nk+2))
        enddo
        enddo

        if(dowr) write(outfile,*)
        do i=ib,ie
          if(dowr) write(outfile,*) '  zs at nj/2:',i,zs(i,nj/2)
        enddo
        if(dowr) write(outfile,*)

        if(dowr) write(outfile,*)
        do j=jb,je
          if(dowr) write(outfile,*) '  zs at ni/2:',j,zs(ni/2,j)
        enddo
        if(dowr) write(outfile,*)

!---------------------------------------

      do k=2,nk
        rds(k) = 1.0/(sigma(k)-sigma(k-1))
      enddo

      do k=1,nk
        rdsf(k) = 1.0/(sigmaf(k+1)-sigmaf(k))
      enddo

!-----------------------------------------------------------------------

      end subroutine init_terrain


  END MODULE init_terrain_module
