  module init_physics_module

  implicit none

  contains

      subroutine init_physics(prs0,rf0,dum1,dum2,dum3,u0,ua,v0,va,o30,   &
                             lu_index,xland,emiss,thc,albd,znt,mavail,tsk,u1,v1,s1, &
                             zh,u10,v10,wspd,lowlyr,ivgtyp,ust,sice)
      use input
      use constants
      use module_sf_slab
      use module_sf_sfclay
      use radtrns3d_module
      use module_ra_rrtmg_lw , only : rrtmg_lwinit
      use module_ra_rrtmg_sw , only : rrtmg_swinit
      use module_sf_mynn , only : mynn_sf_init_driver
      use module_sf_myjsfc , only : myjsfcinit
      use sf_sfclayrev , only : sf_sfclayrev_init
      use cu_ntiedtke, only : cu_ntiedtke_init
#ifdef MPI
      use mpi
#endif
      implicit none

      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: prs0,rf0
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0,ua
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0,va
      real, intent(inout), dimension(ibr:ier,jbr:jer,kbr:ker) :: o30
      integer, intent(in), dimension(ibl:iel,jbl:jel) :: lu_index
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: emiss,thc,albd,mavail
      real, intent(inout), dimension(ib:ie,jb:je) :: tsk,znt,u1,v1,s1,xland
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: u10,v10,wspd
      integer, intent(inout), dimension(ibmyj:iemyj,jbmyj:jemyj) :: lowlyr,ivgtyp
      real, intent(inout), dimension(ib:ie,jb:je) :: ust
      real, intent(inout), dimension(ibmyj:iemyj,jbmyj:jemyj) :: sice

      real :: foo1,foo2,foo3,foo4,foo5,foo6
      real :: p_top

      character(len=80) :: errmsg
      integer :: errflg

!-----------------------------------------------------------------------
!-----  USERS SHOULD NOT NEED TO CHANGE ANYTHING IN THIS FILE ----------
!-----  (unless you really, really know what you are doing -------------
!-----------------------------------------------------------------------

      if(radopt.ge.1)then
        ! initialize radiation code:
        p_top = cgt1*prs0(1,1,nk)+cgt2*prs0(1,1,nk-1)+cgt3*prs0(1,1,nk-2)
        if(myid.eq.0) print *,'  p_top = ',p_top
        call setradwrk(nir,njr,nkr)
        call julday( year, month, day, jday )
        if(dowr) write(outfile,*) '  jday = ',jday
        if( radopt.eq.1 )  &
        call initrad(myid,year,month,day,hour,minute,second,jday,nir,njr,nkr)
        o30 = 1.0e-12
        if( radopt.eq.1 )  &
        call fito3(nir,njr,1,1,nkr,dum1(1,1,1),dum2(1,1,1),prs0,o30,ib,ie,jb,je,kb,ke,nk)
        ! Settings from Goddard scheme:
        call getgoddardvars(foo1,foo2,foo3,foo4,foo5,foo6,ptype)
        roqr = foo1
        tnw  = foo2
        roqs = foo3
        tns  = foo4
        roqg = foo5
        tng  = foo6
      endif

      if( radopt.eq.2 )then

!!!             aclwalloc = .false.
!!!             acswalloc = .false.

#ifdef MPI
             call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
             if(myid.eq.0) print *,'  rrtmg_lwinit '
#ifdef MPI
             call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
             call rrtmg_lwinit(                                             &
                  p_top=p_top,                                              &
                  allowed_to_read = .true.         ,                        &
                   myid = myid ,                                            &
                  ids=1  ,ide=ni+1 , jds= 1 ,jde=nj+1 , kds=1  ,kde=nk+1 ,  &
                  ims=ib ,ime=ie   , jms=jb ,jme=je   , kms=kb ,kme=ke ,    &
                  its=1  ,ite=ni   , jts=1  ,jte=nj   , kts=1  ,kte=nk )

!!!             aclwalloc = .true.

#ifdef MPI
             call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
             if(myid.eq.0) print *,'  rrtmg_swinit '
#ifdef MPI
             call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
             call rrtmg_swinit(                                              &
                   allowed_to_read = .true.         ,                        &
                   myid = myid ,                                             &
                  ids=1  ,ide=ni+1 , jds= 1 ,jde=nj+1 , kds=1  ,kde=nk+1 ,   &
                  ims=ib ,ime=ie   , jms=jb ,jme=je   , kms=kb ,kme=ke ,     &
                  its=1  ,ite=ni   , jts=1  ,jte=nj   , kts=1  ,kte=nk )
#ifdef MPI
             call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
             if(myid.eq.0) print *,'  done rrtmg init '
#ifdef MPI
             call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif

!!!             acswalloc = .true.

      endif

      if( isfcflx.eq.1 )then
        if( sfcmodel.eq.2 )then
          call sfclayinit
        endif
        if( sfcmodel.eq.3 )then
          if(myid.eq.0) print *,'  sf_sfclayrev_init '
          call sf_sfclayrev_init(errmsg=errmsg,errflg=errflg)
        endif
        if( sfcmodel.eq.6 )then
          call mynn_sf_init_driver( .true. )
        endif
        if( sfcmodel.eq.7 )then
          call   MYJSFCINIT(LOWLYR=lowlyr,USTAR=ust,Z0=znt                            &
                           ,SEAMASK=xland,XICE=sice,IVGTYP=ivgtyp,RESTART=.false.                &
                           ,ALLOWED_TO_READ=.false.                    &
                 ,ids=1  ,ide=ni+1 , jds= 1 ,jde=nj+1 , kds=1  ,kde=nk+1 ,                     &
                  ims=ib ,ime=ie   , jms=jb ,jme=je   , kms=kb ,kme=ke ,                       &
                  its=1  ,ite=ni   , jts=1  ,jte=nj   , kts=1  ,kte=nk  )
        endif
      endif

      if( cuparam.eq.1 )then
      if(myid.eq.0) print *
      if(myid.eq.0) print *,'  cu_ntiedtke_init '
      if(myid.eq.0) print *
      call  cu_ntiedtke_init(con_cp=cp,      &
                             con_rd=rd,      &
                             con_rv=rv,      &
                             con_xlv=xlv,    &
                             con_xls=xls,    &
                             con_xlf=xlf,    &
                             con_grav=g,     &
                             errmsg=errmsg,  &
                             errflg=errflg)
      endif

!-----------------------------------------------------------------------

      end subroutine init_physics


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getgoddardvars(foo1,foo2,foo3,foo4,foo5,foo6,ptype)
      use goddard_module, only : ROQR,ROQG,ROQS,TNW,TNG,TNSS
      implicit none

      real, intent(inout) :: foo1,foo2,foo3,foo4,foo5,foo6
      integer, intent(in) :: ptype

    if(ptype.eq.2)then
      foo1 = roqr
      foo2 = tnw
      foo3 = roqs
      foo4 = tnss
      foo5 = roqg
      foo6 = tng
    else
      ! 130903: set to reasonable values to prevent divide-by-zeros
      ! if goodard microphysics scheme is not being used
      foo1 = 1.
      foo2 = .08
      foo3 = .1
      foo4 = 1.
      foo5 = .4
      foo6 = .04
    endif

      end subroutine getgoddardvars


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine initrad(myid,year,month,day,hour,minute,second,jday,   &
                         nir,njr,nkr)
      use radtrns3d_module
      use irrad3d_module
      implicit none

      integer, intent(in) :: myid,year,month,day,hour,minute,second,jday,   &
                             nir,njr,nkr

      integer :: ip,iw,it
      logical :: high

      real, parameter :: pi  = 3.14159265358979323

!----------------------------------------------------------------------

  if ( rlwopt == 0 ) then
    high = .false.
  else
    high = .true.
  end if

  if(myid.eq.0) print *,'  high = ',high

!----------------------------------------------------------------------
!  from zenangl:

    pi2 = 2.0 * pi
    deg2rad = pi/180.0
    rad2deg = 1./deg2rad

    hour0 = FLOAT(hour)                                                 &
          + FLOAT(minute)/60.0                                          &
          + FLOAT(second)/3600.0

    if ( MOD(year, 4) == 0 ) then
      yrday = 366.
    else
      yrday = 365.
    end if

!!! not using arps MPI code:  GHB, 100720
!!! hard-wire these in, just in case:
    nxmid = 1
    nymid = 1
    source = 0

!----------------------------------------------------------------------
!  from irrad:

!-----tables co2 and h2o are only used with 'high' option

    if (high) then

      do iw=1,nh
        do ip=1,nx
          h11(ip,iw,1)=1.0-h11(ip,iw,1)
          h21(ip,iw,1)=1.0-h21(ip,iw,1)
          h71(ip,iw,1)=1.0-h71(ip,iw,1)
        end do
      end do

      do iw=1,nc
        do ip=1,nx
          c1(ip,iw,1)=1.0-c1(ip,iw,1)
        end do
      end do

!-----tables are replicated to avoid memory bank conflicts

      do it=2,nt
        do iw=1,nc
          do ip=1,nx
            c1 (ip,iw,it)= c1(ip,iw,1)
            c2 (ip,iw,it)= c2(ip,iw,1)
            c3 (ip,iw,it)= c3(ip,iw,1)
          end do
        end do
        do iw=1,nh
          do ip=1,nx
            h11(ip,iw,it)=h11(ip,iw,1)
            h12(ip,iw,it)=h12(ip,iw,1)
            h13(ip,iw,it)=h13(ip,iw,1)
            h21(ip,iw,it)=h21(ip,iw,1)
            h22(ip,iw,it)=h22(ip,iw,1)
            h23(ip,iw,it)=h23(ip,iw,1)
            h71(ip,iw,it)=h71(ip,iw,1)
            h72(ip,iw,it)=h72(ip,iw,1)
            h73(ip,iw,it)=h73(ip,iw,1)
          end do
        end do
      end do

    end if

!-----always use table look-up for ozone transmittance

    do iw=1,no
      do ip=1,nx
        o1(ip,iw,1)=1.0-o1(ip,iw,1)
      end do
    end do

    do it=2,nt
      do iw=1,no
        do ip=1,nx
          o1 (ip,iw,it)= o1(ip,iw,1)
          o2 (ip,iw,it)= o2(ip,iw,1)
          o3 (ip,iw,it)= o3(ip,iw,1)
        end do
      end do
    end do

      return
      end subroutine initrad

  end module init_physics_module
