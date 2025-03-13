  MODULE turbnba_module

  implicit none

  private
  public :: turbnba,turbnba2

      real, parameter :: tke_min         =  1.0e-10
      real, parameter :: nm_min          =  1.0e-6
      real, parameter :: small_len_frac  =  0.001

      real, parameter :: odt = 1.0/3.0

      real, parameter :: nba_cb = 0.36
      real, parameter :: nba_sk = 0.5

  CONTAINS


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine turbnba(nstep,uh,ruh,uf,ruf,vh,rvh,vf,rvf,mh,rmh,mf,rmf,zf,c1,c2,rho,rf,zntmp,ust,cm0,  &
                         s11 ,s12 ,s13 ,s22 ,s23 ,s33 ,            &
                         r12 ,r13 ,r23  ,lenscl,grdscl,tk  ,cme,   &
                         m11,m12,m13,m22,m23,m33,                  &
                         ua ,va ,wa ,tkea ,nm,                     &
                         kw1,kw2,ke1,ke2,ks1,ks2,kn1,kn2,reqs_s,   &
                         nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      use input
      use constants
      use bc_module
      use comm_module
      implicit none

!-----------------------------------------------------------------------

      integer, intent(in) :: nstep
      real, intent(in), dimension(ib:ie) :: uh,ruh
      real, intent(in), dimension(ib:ie+1) :: uf,ruf
      real, intent(in), dimension(jb:je) :: vh,rvh
      real, intent(in), dimension(jb:je+1) :: vf,rvf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh,rmh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf,rmf,zf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2,rho,rf
      real, intent(in), dimension(ib:ie,jb:je) :: zntmp,ust,cm0
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: r12,r13,r23
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: lenscl,grdscl,tk
      real, intent(in),    dimension(ibc:iec,jbc:jec,kbc:kec) :: cme
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: s11,s12,s13,s22,s23,s33
      real, intent(inout), dimension(ibnba:ienba,jbnba:jenba,kbnba:kenba) :: m11,m12,m13,m22,m23,m33
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(inout), dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: nm
      real, intent(inout), dimension(jmp,kmt,4)     :: kw1,kw2,ke1,ke2
      real, intent(inout), dimension(imp,kmt,4)     :: ks1,ks2,kn1,kn2
      integer, intent(inout), dimension(rmp) :: reqs_s
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2

!-----------------------------------------------------------------------

      integer :: i,j,k
      real :: tem,tem1,tkk
      real :: nba_cs,nba_c1,nba_c2,nba_c3,nba_ce,nba_ck
      real :: nba_a,nba_b,delta,tmp11,tmp22,tmp33,tmp12,tmp13,tmp23,tmr12,tmr13,tmr23,smn

      real, dimension(cmp,jmp) :: west,newwest,east,neweast
      real, dimension(imp,cmp) :: south,newsouth,north,newnorth

!-----------------------------------------------------------------------

        nba_ck = (2.0/3.0)*( (pi*c_m)**(-2.0/3.0) )

        nba_cs = sqrt( ( 8.0*( 1.0+nba_cb ) )/( 27.0*pi*pi ) )  &
                 * cme(1,1,2) / c_m
      if( myid.eq.0 ) print *,'  tem = ',cme(1,1,2) / c_m

        nba_c1 = ( sqrt(960.0)*nba_cb )/( 7.0*( 1.0+nba_cb )*nba_sk )
        nba_c2 = nba_c1
        nba_ce = ( ( 8.0*pi/27.0 )**(1.0/3.0) )*( nba_cs**(4.0/3.0) )
        nba_c3 = ( ( 27.0/( 8.0*pi ) )**(1.0/3.0) )*( nba_cs**(2.0/3.0) )

        if( myid.eq.0 ) print *,'  cs,ce,ck = ',nba_cs,nba_ce,nba_ck
        if( myid.eq.0 ) print *,'  c1,c2,c3 = ',nba_c1,nba_c2,nba_c3

!------------------------------------------------------------------

    DO k=2,nk

    !  get grid scale
    IF(tconfig.eq.1)THEN
      ! single length scale:  appropriate if dx,dy are nearly the same as dz
      do j=1,nj
      do i=1,ni
        grdscl(i,j,k)=( ((dx*ruh(i))*(dy*rvh(j)))*(dz*rmf(i,j,k)) )**0.33333333
        ! cm1r17:  wall condition near surface
        ! cm1r20.1: revisit near-sfc grid scale at a later date
!!!        grdscl(i,j,k) = sqrt(1.0/( 1.0/(grdscl(i,j,k)**2)                                  &
!!!                                  +1.0/((karman*((zf(i,j,k)-zf(i,j,1))+zntmp(i,j))*rcs)**2)  &
!!!                               ) )
      enddo
      enddo
    ELSEIF(tconfig.eq.2)THEN
      ! two length scales:  one for horizontal, one for vertical
      stop 15151
    ENDIF

      ! Get turbulence length scale
      do j=1,nj
      do i=1,ni
        if( tkea(i,j,k).le.tke_min )then
          ! 170718:
          tk(i,j,k) = tke_min
          lenscl(i,j,k) = small_len_frac*grdscl(i,j,k)
        else
          tk(i,j,k)=tkea(i,j,k)
          lenscl(i,j,k)=grdscl(i,j,k)
          if(nm(i,j,k).gt.nm_min)then
            lenscl(i,j,k)=c_l*sqrt(tk(i,j,k)/nm(i,j,k))
            lenscl(i,j,k)=min(lenscl(i,j,k),grdscl(i,j,k))
            lenscl(i,j,k)=max(lenscl(i,j,k),small_len_frac*grdscl(i,j,k))
          endif 
        endif
      enddo
      enddo

    ENDDO

    if( nstep.eq.0 .and. myid.eq.0 )then
      print *
      print *,'  zf,grdscl:'
      i = 1
      j = 1
      do k=2,nk
        print *,k,(zf(i,j,k)-zf(i,j,1)),grdscl(i,j,k)
      enddo
    endif

#ifdef MPI
      call comm_1s_start(lenscl,kw1(1,1,1),kw2(1,1,1),ke1(1,1,1),ke2(1,1,1),  &
                                ks1(1,1,1),ks2(1,1,1),kn1(1,1,1),kn2(1,1,1),reqs_s)
      call comm_1s_end(lenscl,kw1(1,1,1),kw2(1,1,1),ke1(1,1,1),ke2(1,1,1),  &
                              ks1(1,1,1),ks2(1,1,1),kn1(1,1,1),kn2(1,1,1),reqs_s)
      call bcs2(lenscl)
      call getcorner(lenscl,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
#endif

      do j=0,nj+1
      do i=0,ni+1
        lenscl(i,j,1) = karman*zntmp(i,j)*rcs
        lenscl(i,j,nk+1) = lenscl(i,j,nk)
      enddo
      enddo

      if( cm1setup.eq.4 )then
        ! zero-out lenscl in mesoscale-model part of domain:

        DO k=1,nk
        do j=0,nj+1
        do i=0,ni+1
          if( cm0(i,j).le.cmemin ) lenscl(i,j,k) = 0.0
        enddo
        enddo
        ENDDO

      endif

!------------------------------------------------------------------

        DO k=1,nk

          do j=0,nj+2
          do i=0,ni+2
            s11(i,j,k)=(ua(i+1,j,k)-ua(i,j,k))*rdx*uh(i)
            s22(i,j,k)=(va(i,j+1,k)-va(i,j,k))*rdy*vh(j)
            s33(i,j,k)=(wa(i,j,k+1)-wa(i,j,k))*rdz*mh(i,j,k)
            s12(i,j,k)=0.5*( (ua(i,j,k)-ua(i,j-1,k))*rdy*vf(j)   &
                            +(va(i,j,k)-va(i-1,j,k))*rdx*uf(i) )
            s13(i,j,k)=0.5*( (wa(i,j,k)-wa(i-1,j,k))*rdx*uf(i)   &
                            +(ua(i,j,k)-ua(i,j,k-1))*rdz*mf(1,1,k) )
            s23(i,j,k)=0.5*( (wa(i,j,k)-wa(i,j-1,k))*rdy*vf(j)   &
                            +(va(i,j,k)-va(i,j,k-1))*rdz*mf(1,1,k) )
            r12(i,j,k) = 0.5*( (ua(i,j,k)-ua(i,j-1,k))*rdy*vf(j) &
                              -(va(i,j,k)-va(i-1,j,k))*rdx*uf(i) )
            r13(i,j,k) = 0.5*(-(wa(i,j,k)-wa(i-1,j,k))*rdx*uf(i) &
                              +(ua(i,j,k)-ua(i,j,k-1))*rdz*mf(1,1,k) )
            r23(i,j,k) = 0.5*(-(wa(i,j,k)-wa(i,j-1,k))*rdy*vf(j) &
                              +(va(i,j,k)-va(i,j,k-1))*rdz*mf(1,1,k) )
          enddo
          enddo

        ENDDO

      do j=-1,nj+2
      do i=-1,ni+2
        grdscl(i,j,1) = ust(i,j)/(karman*zntmp(i,j))
        grdscl(i,j,2) = sqrt( (0.5*(ua(i,j,1)+ua(i+1,j,1)))**2 &
                             +(0.5*(va(i,j,1)+va(i,j+1,1)))**2 )
      enddo
      enddo

      do k=1,2
        call bc2d(grdscl(ib,jb,k))
#ifdef MPI
        call comm_2d_start(grdscl(ib,jb,k),west,newwest,east,neweast,   &
                                           south,newsouth,north,newnorth,reqs_s)
        call comm_2dew_end(grdscl(ib,jb,k),west,newwest,east,neweast,reqs_s)
        call comm_2dns_end(grdscl(ib,jb,k),south,newsouth,north,newnorth,reqs_s)
        call bcs2_2d(grdscl(ib,jb,k))
        call comm_2d_corner(grdscl(ib,jb,k))
#endif
      enddo

      tem1 = 1.0/(zf(1,1,2)-zf(1,1,1))

      do j=0,nj+2
      do i=0,ni+2
!!!        s13(i,j,1) = 0.5*(dgs3*ua(i,j,3)+dgs2*ua(i,j,2)+dgs1*ua(i,j,1))*tem1
!!!        s23(i,j,1) = 0.5*(dgs3*va(i,j,3)+dgs2*va(i,j,2)+dgs1*va(i,j,1))*tem1
!!!        s13(i,j,1) = (grdscl(i-1,j,1)+grdscl(i,j,1))*ua(i,j,1)/max(0.01,grdscl(i-1,j,2)+grdscl(i,j,2))
!!!        s23(i,j,1) = (grdscl(i,j-1,1)+grdscl(i,j,1))*va(i,j,1)/max(0.01,grdscl(i,j-1,2)+grdscl(i,j,2))
        s13(i,j,1) = 0.0
        s23(i,j,1) = 0.0
        r13(i,j,1) = s13(i,j,1)
        r23(i,j,1) = s23(i,j,1)
        s13(i,j,nk+1) = 0.0
        s23(i,j,nk+1) = 0.0
        r13(i,j,nk+1) = s13(i,j,nk+1)
        r23(i,j,nk+1) = s23(i,j,nk+1)
      enddo
      enddo

#ifdef MPI
      call bct2(tkea)
      call getcornert(tkea,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
#endif

  kloop:  DO k=1,nk

    do j=0,nj+1
    do i=0,ni+1

      delta = 0.5*( lenscl(i,j,k  ) &
                   +lenscl(i,j,k+1) )

      nba_a = nba_ce
      nba_b = nba_ce*delta*nba_c3*delta

      tmp11 = s11(i,j,k)**2
      tmp22 = s22(i,j,k)**2
      tmp33 = s33(i,j,k)**2

      tmp12 = 0.25*( ( s12(i,j  ,k)**2 + s12(i+1,j+1,k)**2 ) &
                   + ( s12(i,j+1,k)**2 + s12(i+1,j  ,k)**2 ) )

      tmp13 = 0.25*( ( s13(i,j,k  )**2 + s13(i+1,j,k  )**2 ) &
                   + ( s13(i,j,k+1)**2 + s13(i+1,j,k+1)**2 ) )

      tmp23 = 0.25*( ( s23(i,j,k  )**2 + s23(i,j+1,k  )**2 ) &
                   + ( s23(i,j,k+1)**2 + s23(i,j+1,k+1)**2 ) )

      smn = ( tmp11 + tmp22 ) + tmp33 + 2.0*( tmp12 + ( tmp13 + tmp23 ) )

      tmr12 = 0.25*( ( s12(i,j  ,k)*r12(i,j  ,k) + s12(i+1,j+1,k)*r12(i+1,j+1,k) ) &
                   + ( s12(i,j+1,k)*r12(i,j+1,k) + s12(i+1,j  ,k)*r12(i+1,j  ,k) ) )

      tmr13 = 0.25*( ( s13(i,j,k  )*r13(i,j,k  ) + s13(i+1,j,k  )*r13(i+1,j,k  ) ) &
                   + ( s13(i,j,k+1)*r13(i,j,k+1) + s13(i+1,j,k+1)*r13(i+1,j,k+1) ) )

      tmr23 = 0.25*( ( s23(i,j,k  )*r23(i,j,k  ) + s23(i,j+1,k  )*r23(i,j+1,k  ) ) &
                   + ( s23(i,j,k+1)*r23(i,j,k+1) + s23(i,j+1,k+1)*r23(i,j+1,k+1) ) )

      tkk = 0.5*( lenscl(i,j,k)*sqrt(tkea(i,j,k)) + lenscl(i,j,k+1)*sqrt(tkea(i,j,k+1)) )

      m11(i,j,k) =  rho(i,j,k)*( nba_ce*2.0*tkk*s11(i,j,k)                       &
                   + nba_b*( nba_c1*( ( tmp11 + tmp12 ) + ( tmp13 - smn*odt ) )  &
                           + nba_c2*( -2.0*( tmr12 + tmr13 ) )                   &
                         ) )
      m22(i,j,k) =  rho(i,j,k)*( nba_ce*2.0*tkk*s22(i,j,k)                       &
                   + nba_b*( nba_c1*( ( tmp22 + tmp12 ) + ( tmp23 - smn*odt ) )  &
                           + nba_c2*(  2.0*( tmr12 - tmr23 ) )                   &
                         ) )

      m33(i,j,k) =  rho(i,j,k)*( nba_ce*2.0*tkk*s33(i,j,k)                       &
                   + nba_b*( nba_c1*( ( tmp13 + tmp23 ) + ( tmp33 - smn*odt ) )  &
                           + nba_c2*(  2.0*( tmr13 + tmr23 ) )                   &
                         ) )

!!!      t11(i,j,k) = rho(i,j,k)*nba_a*(   2.0*tkk*s11(i,j,k) )
!!!      t22(i,j,k) = rho(i,j,k)*nba_a*(   2.0*tkk*s22(i,j,k) )
!!!      t33(i,j,k) = rho(i,j,k)*nba_a*(   2.0*tkk*s33(i,j,k) )

    enddo
    enddo

    do j=1,nj+1
    do i=1,ni+1

      delta = 0.125*( ( ( lenscl(i  ,j  ,k  )     &
                         +lenscl(i-1,j-1,k  ) )   &
                       +( lenscl(i-1,j  ,k  )     &
                         +lenscl(i  ,j-1,k  ) ) ) &
                     +( ( lenscl(i  ,j  ,k+1)     &
                         +lenscl(i-1,j-1,k+1) )   &
                       +( lenscl(i-1,j  ,k+1)     &
                         +lenscl(i  ,j-1,k+1) ) ) )

      nba_a = nba_ce
      nba_b = nba_ce*delta*nba_c3*delta

      tmp11 = 0.25*(  (s11(i  ,j,k)+s11(i-1,j-1,k)) &
                     +(s11(i-1,j,k)+s11(i  ,j-1,k)) )
      tmp22 = 0.25*(  (s22(i  ,j,k)+s22(i-1,j-1,k)) &
                     +(s22(i-1,j,k)+s22(i  ,j-1,k)) )
      tmp12 = s12(i,j,k)
      tmp13 = 0.25*( (s13(i,j,k  )+s13(i,j-1,k  )) &
                    +(s13(i,j,k+1)+s13(i,j-1,k+1)) )
      tmp23 = 0.25*( (s23(i,j,k  )+s23(i-1,j,k  )) &
                    +(s23(i,j,k+1)+s23(i-1,j,k+1)) )

      tkk = 0.125*( ( (lenscl(i  ,j,k  )*sqrt(tkea(i  ,j,k  ))+lenscl(i-1,j-1,k  )*sqrt(tkea(i-1,j-1,k  )))   &
                     +(lenscl(i-1,j,k  )*sqrt(tkea(i-1,j,k  ))+lenscl(i  ,j-1,k  )*sqrt(tkea(i  ,j-1,k  ))) ) &
                   +( (lenscl(i  ,j,k+1)*sqrt(tkea(i  ,j,k+1))+lenscl(i-1,j-1,k+1)*sqrt(tkea(i-1,j-1,k+1)))   &
                     +(lenscl(i-1,j,k+1)*sqrt(tkea(i-1,j,k+1))+lenscl(i  ,j-1,k+1)*sqrt(tkea(i  ,j-1,k+1))) ) )

      tmr12 = r12(i,j,k)
      tmr13 = 0.25*( (r13(i,j,k  )+r13(i,j-1,k  )) &
                    +(r13(i,j,k+1)+r13(i,j-1,k+1)) )
      tmr23 = 0.25*( (r23(i,j,k  )+r23(i-1,j,k  )) &
                    +(r23(i,j,k+1)+r23(i-1,j,k+1)) )

      m12(i,j,k) = nba_ce*2.0*tkk*s12(i,j,k)  &
                   + nba_b*( nba_c1*( ( tmp11*tmp12 + tmp12*tmp22 ) + tmp13*tmp23 )  &
                           + nba_c2*( (   tmp11*tmr12 - tmp22*tmr12 )                &
                                     +( - tmp13*tmr23 - tmp23*tmr13 ) )              &
                           )
      m12(i,j,k) = m12(i,j,k)*0.25*( (rho(i-1,j-1,k)+rho(i,j,k))+(rho(i-1,j,k)+rho(i,j-1,k)) )

!!!      t12(i,j,k) = nba_a*(   2.0*tkk*s12(i,j,k) )
!!!      t12(i,j,k) = t12(i,j,k)*0.25*( (rho(i-1,j-1,k)+rho(i,j,k))+(rho(i-1,j,k)+rho(i,j-1,k)) )

    enddo
    enddo

  kcheck:  IF( k.ge.2 )THEN

    do j=1,nj
    do i=1,ni+1

      delta = 0.5*( lenscl(i-1,j,k) &
                   +lenscl(i  ,j,k) )
      nba_a = nba_ce
      nba_b = nba_ce*delta*nba_c3*delta

      tmp11 = 0.5*( c1(1,1,k)*(s11(i-1,j,k-1)+s11(i,j,k-1)) &
                   +c2(1,1,k)*(s11(i-1,j,k  )+s11(i,j,k  )) )
      tmp13 = s13(i,j,k)
      tmp12 = 0.5*( c1(1,1,k)*(s12(i,j+1,k-1)+s12(i,j,k-1)) &
                   +c2(1,1,k)*(s12(i,j+1,k  )+s12(i,j,k  )) )
      tmp23 = 0.25*( (s23(i-1,j+1,k)+s23(i,j  ,k)) &
                    +(s23(i-1,j  ,k)+s23(i,j+1,k)) )
      tmp33 = 0.5*( c1(1,1,k)*(s33(i-1,j,k-1)+s33(i,j,k-1)) &
                   +c2(1,1,k)*(s33(i-1,j,k  )+s33(i,j,k  )) )

      tkk = 0.5*( lenscl(i,j,k)*sqrt(tkea(i,j,k)) + lenscl(i-1,j,k)*sqrt(tkea(i-1,j,k)) )

      tmr13 = r13(i,j,k)
      tmr23 = 0.25*( (r23(i-1,j+1,k)+r23(i,j  ,k)) &
                    +(r23(i-1,j  ,k)+r23(i,j+1,k)) )
      tmr12 = 0.5*( c1(1,1,k)*(r12(i,j+1,k-1)+r12(i,j,k-1)) &
                   +c2(1,1,k)*(r12(i,j+1,k  )+r12(i,j,k  )) )

      m13(i,j,k) = nba_ce*2.0*tkk*s13(i,j,k)  &
                   + nba_b*( nba_c1*( ( tmp11*tmp13 + tmp13*tmp33 ) + tmp12*tmp23 )  &
                           + nba_c2*(   ( tmp11*tmr13 - tmp33*tmr13 )                &
                                      + ( tmp12*tmr23 - tmp23*tmr12 )  )             &
                           )
      m13(i,j,k) = m13(i,j,k)*0.5*( rf(i-1,j,k)+rf(i,j,k) )

!!!      t13(i,j,k) = nba_a*(   2.0*tkk*s13(i,j,k) )
!!!      t13(i,j,k) = t13(i,j,k)*0.5*( rf(i-1,j,k)+rf(i,j,k) )

    enddo
    enddo

    do j=1,nj+1
    do i=1,ni

      delta = 0.5*( lenscl(i,j-1,k) &
                   +lenscl(i,j  ,k) )
      nba_a = nba_ce
      nba_b = nba_ce*delta*nba_c3*delta

      tmp23 = s23(i,j,k)
      tmp12 = 0.5*( c1(1,1,k)*(s12(i,j,k-1)+s12(i+1,j,k-1)) &
                   +c2(1,1,k)*(s12(i,j,k  )+s12(i+1,j,k  )) )
      tmp13 = 0.25*( (s13(i+1,j-1,k)+s13(i,j  ,k)) &
                    +(s13(i+1,j  ,k)+s13(i,j-1,k)) )
      tmp33 = 0.5*( c1(1,1,k)*(s33(i,j,k-1)+s33(i,j-1,k-1)) &
                   +c2(1,1,k)*(s33(i,j,k  )+s33(i,j-1,k  )) )
      tmp22 = 0.5*( c1(1,1,k)*(s22(i,j,k-1)+s22(i,j-1,k-1)) &
                   +c2(1,1,k)*(s22(i,j,k  )+s22(i,j-1,k  )) )

      tkk = 0.5*( lenscl(i,j,k)*sqrt(tkea(i,j,k)) + lenscl(i,j-1,k)*sqrt(tkea(i,j-1,k)) )

      tmr23 = r23(i,j,k)
      tmr13 = 0.25*( (r13(i+1,j-1,k)+r13(i,j  ,k)) &
                    +(r13(i+1,j  ,k)+r13(i,j-1,k)) )
      tmr12 = 0.5*( c1(1,1,k)*(r12(i,j,k-1)+r12(i+1,j,k-1)) &
                   +c2(1,1,k)*(r12(i,j,k  )+r12(i+1,j,k  )) )

      m23(i,j,k) = nba_ce*2.0*tkk*s23(i,j,k)  &
                   + nba_b*( nba_c1*( ( tmp22*tmp23 + tmp23*tmp33 ) + tmp12*tmp13 )  &
                           + nba_c2*(   ( tmp22*tmr23 - tmp33*tmr23 )                &
                                      + ( tmp12*tmr13 + tmp13*tmr12 )  )             &
                           )
      m23(i,j,k) = m23(i,j,k)*0.5*( rf(i,j-1,k)+rf(i,j,k) )

!!!      t23(i,j,k) = nba_a*(   2.0*tkk*s23(i,j,k) )
!!!      t23(i,j,k) = t23(i,j,k)*0.5*( rf(i,j-1,k)+rf(i,j,k) )

    enddo
    enddo

  ENDIF  kcheck

  ENDDO  kloop

!------------------------------------------------------------------

      IF(bbc.eq.1)THEN
        ! free slip:
        do j=1,nj+1
        do i=1,ni+1
          m13(i,j,1) = m13(i,j,2)
          m23(i,j,1) = m23(i,j,2)
        enddo
        enddo
      ELSEIF(bbc.eq.2.or.bbc.eq.3)THEN
        IF(axisymm.eq.0)THEN
!!!          do j=1,nj+1
!!!          do i=1,ni+1
!!!            m13(i,j,1) = t13(i,j,1)
!!!            m23(i,j,1) = t23(i,j,1)
!!!          enddo
!!!          enddo
        ELSE
          stop 1212
        ENDIF
      ELSE
        stop 2222
      ENDIF

!------------------------------------------------------------------

      end subroutine turbnba


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine turbnba2(nstep,uh,ruh,uf,ruf,vh,rvh,vf,rvf,mh,rmh,mf,rmf,zf,c1,c2,rho,rf,zntmp,ust,cm0,  &
                         s11 ,s12 ,s13 ,s22 ,s23 ,s33 ,            &
                         r12 ,r13 ,r23  ,lenscl,grdscl,tk  ,cme,   &
                         m11,m12,m13,m22,m23,m33,                  &
                         ua ,va ,wa ,tkea ,nm,                     &
                         kw1,kw2,ke1,ke2,ks1,ks2,kn1,kn2,reqs_s,   &
                         nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      use input
      use constants
      use bc_module
      use comm_module
      implicit none

!-----------------------------------------------------------------------

      integer, intent(in) :: nstep
      real, intent(in), dimension(ib:ie) :: uh,ruh
      real, intent(in), dimension(ib:ie+1) :: uf,ruf
      real, intent(in), dimension(jb:je) :: vh,rvh
      real, intent(in), dimension(jb:je+1) :: vf,rvf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh,rmh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf,rmf,zf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2,rho,rf
      real, intent(in), dimension(ib:ie,jb:je) :: zntmp,ust,cm0
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: r12,r13,r23
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: lenscl,grdscl,tk
      real, intent(in),    dimension(ibc:iec,jbc:jec,kbc:kec) :: cme
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: s11,s12,s13,s22,s23,s33
      real, intent(inout), dimension(ibnba:ienba,jbnba:jenba,kbnba:kenba) :: m11,m12,m13,m22,m23,m33
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(inout), dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: nm
      real, intent(inout), dimension(jmp,kmt,4)     :: kw1,kw2,ke1,ke2
      real, intent(inout), dimension(imp,kmt,4)     :: ks1,ks2,kn1,kn2
      integer, intent(inout), dimension(rmp) :: reqs_s
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2

!-----------------------------------------------------------------------

      integer :: i,j,k
      real :: tem,tem1,tkk
      real :: nba_cs,nba_c1,nba_c2,nba_c3,nba_ce,nba_ck
      real :: nba_a,nba_b,delta,tmp11,tmp22,tmp33,tmp12,tmp13,tmp23,tmr12,tmr13,tmr23,smn

      real, dimension(cmp,jmp) :: west,newwest,east,neweast
      real, dimension(imp,cmp) :: south,newsouth,north,newnorth

!-----------------------------------------------------------------------

        nba_ck = (2.0/3.0)*( (pi*c_m)**(-2.0/3.0) )

        nba_cs = sqrt( ( 8.0*( 1.0+nba_cb ) )/( 27.0*pi*pi ) )  &
                 * cme(1,1,2) / c_m
      if( myid.eq.0 ) print *,'  tem = ',cme(1,1,2) / c_m

        nba_c1 = ( sqrt(960.0)*nba_cb )/( 7.0*( 1.0+nba_cb )*nba_sk )
        nba_c2 = nba_c1
        nba_ce = ( ( 8.0*pi/27.0 )**(1.0/3.0) )*( nba_cs**(4.0/3.0) )
        nba_c3 = ( ( 27.0/( 8.0*pi ) )**(1.0/3.0) )*( nba_cs**(2.0/3.0) )

        if( myid.eq.0 ) print *,'  cs,ce,ck = ',nba_cs,nba_ce,nba_ck
        if( myid.eq.0 ) print *,'  c1,c2,c3 = ',nba_c1,nba_c2,nba_c3

!------------------------------------------------------------------

    DO k=2,nk

    !  get grid scale
    IF(tconfig.eq.1)THEN
      ! single length scale:  appropriate if dx,dy are nearly the same as dz
      do j=1,nj
      do i=1,ni
        grdscl(i,j,k)=( ((dx*ruh(i))*(dy*rvh(j)))*(dz*rmf(i,j,k)) )**0.33333333
        ! cm1r17:  wall condition near surface
        ! cm1r20.1: revisit near-sfc grid scale at a later date
!!!        grdscl(i,j,k) = sqrt(1.0/( 1.0/(grdscl(i,j,k)**2)                                  &
!!!                                  +1.0/((karman*((zf(i,j,k)-zf(i,j,1))+zntmp(i,j))*rcs)**2)  &
!!!                               ) )
      enddo
      enddo
    ELSEIF(tconfig.eq.2)THEN
      ! two length scales:  one for horizontal, one for vertical
      stop 15151
    ENDIF

    ENDDO

    if( nstep.eq.0 .and. myid.eq.0 )then
      print *
      print *,'  zf,grdscl:'
      i = 1
      j = 1
      do k=2,nk
        print *,k,(zf(i,j,k)-zf(i,j,1)),grdscl(i,j,k)
      enddo
    endif

#ifdef MPI
      call comm_1s_start(grdscl,kw1(1,1,1),kw2(1,1,1),ke1(1,1,1),ke2(1,1,1),  &
                                ks1(1,1,1),ks2(1,1,1),kn1(1,1,1),kn2(1,1,1),reqs_s)
      call comm_1s_end(grdscl,kw1(1,1,1),kw2(1,1,1),ke1(1,1,1),ke2(1,1,1),  &
                              ks1(1,1,1),ks2(1,1,1),kn1(1,1,1),kn2(1,1,1),reqs_s)
      call bcs2(grdscl)
      call getcorner(grdscl,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
#endif

      do j=0,nj+1
      do i=0,ni+1
        grdscl(i,j,1) = grdscl(i,j,2)
        grdscl(i,j,nk+1) = grdscl(i,j,nk)
      enddo
      enddo

      if( cm1setup.eq.4 )then
        ! zero-out grdscl in mesoscale-model part of domain:

        DO k=1,nk
        do j=0,nj+1
        do i=0,ni+1
          if( cm0(i,j).le.cmemin ) grdscl(i,j,k) = 0.0
        enddo
        enddo
        ENDDO

      endif

!------------------------------------------------------------------

        DO k=1,nk

          do j=0,nj+2
          do i=0,ni+2
            s11(i,j,k)=(ua(i+1,j,k)-ua(i,j,k))*rdx*uh(i)
            s22(i,j,k)=(va(i,j+1,k)-va(i,j,k))*rdy*vh(j)
            s33(i,j,k)=(wa(i,j,k+1)-wa(i,j,k))*rdz*mh(i,j,k)
            s12(i,j,k)=0.5*( (ua(i,j,k)-ua(i,j-1,k))*rdy*vf(j)   &
                            +(va(i,j,k)-va(i-1,j,k))*rdx*uf(i) )
            s13(i,j,k)=0.5*( (wa(i,j,k)-wa(i-1,j,k))*rdx*uf(i)   &
                            +(ua(i,j,k)-ua(i,j,k-1))*rdz*mf(1,1,k) )
            s23(i,j,k)=0.5*( (wa(i,j,k)-wa(i,j-1,k))*rdy*vf(j)   &
                            +(va(i,j,k)-va(i,j,k-1))*rdz*mf(1,1,k) )
            r12(i,j,k) = 0.5*( (ua(i,j,k)-ua(i,j-1,k))*rdy*vf(j) &
                              -(va(i,j,k)-va(i-1,j,k))*rdx*uf(i) )
            r13(i,j,k) = 0.5*(-(wa(i,j,k)-wa(i-1,j,k))*rdx*uf(i) &
                              +(ua(i,j,k)-ua(i,j,k-1))*rdz*mf(1,1,k) )
            r23(i,j,k) = 0.5*(-(wa(i,j,k)-wa(i,j-1,k))*rdy*vf(j) &
                              +(va(i,j,k)-va(i,j,k-1))*rdz*mf(1,1,k) )
          enddo
          enddo

        ENDDO

      tem1 = 1.0/(zf(1,1,2)-zf(1,1,1))

      do j=0,nj+2
      do i=0,ni+2
!!!        s13(i,j,1) = 0.5*(dgs3*ua(i,j,3)+dgs2*ua(i,j,2)+dgs1*ua(i,j,1))*tem1
!!!        s23(i,j,1) = 0.5*(dgs3*va(i,j,3)+dgs2*va(i,j,2)+dgs1*va(i,j,1))*tem1
        s13(i,j,1) = 0.0
        s23(i,j,1) = 0.0
        r13(i,j,1) = s13(i,j,1)
        r23(i,j,1) = s23(i,j,1)
        s13(i,j,nk+1) = 0.0
        s23(i,j,nk+1) = 0.0
        r13(i,j,nk+1) = s13(i,j,nk+1)
        r23(i,j,nk+1) = s23(i,j,nk+1)
      enddo
      enddo

  kloop:  DO k=1,nk

    do j=0,nj+1
    do i=0,ni+1

      delta = 0.5*( grdscl(i,j,k  ) &
                   +grdscl(i,j,k+1) )

      nba_a = nba_cs*nba_cs*delta*delta

      tmp11 = s11(i,j,k)**2
      tmp22 = s22(i,j,k)**2
      tmp33 = s33(i,j,k)**2

      tmp12 = 0.25*( ( s12(i,j  ,k)**2 + s12(i+1,j+1,k)**2 ) &
                   + ( s12(i,j+1,k)**2 + s12(i+1,j  ,k)**2 ) )

      tmp13 = 0.25*( ( s13(i,j,k  )**2 + s13(i+1,j,k  )**2 ) &
                   + ( s13(i,j,k+1)**2 + s13(i+1,j,k+1)**2 ) )

      tmp23 = 0.25*( ( s23(i,j,k  )**2 + s23(i,j+1,k  )**2 ) &
                   + ( s23(i,j,k+1)**2 + s23(i,j+1,k+1)**2 ) )

      smn = ( tmp11 + tmp22 ) + tmp33 + 2.0*( tmp12 + ( tmp13 + tmp23 ) )

      tmr12 = 0.25*( ( s12(i,j  ,k)*r12(i,j  ,k) + s12(i+1,j+1,k)*r12(i+1,j+1,k) ) &
                   + ( s12(i,j+1,k)*r12(i,j+1,k) + s12(i+1,j  ,k)*r12(i+1,j  ,k) ) )

      tmr13 = 0.25*( ( s13(i,j,k  )*r13(i,j,k  ) + s13(i+1,j,k  )*r13(i+1,j,k  ) ) &
                   + ( s13(i,j,k+1)*r13(i,j,k+1) + s13(i+1,j,k+1)*r13(i+1,j,k+1) ) )

      tmr23 = 0.25*( ( s23(i,j,k  )*r23(i,j,k  ) + s23(i,j+1,k  )*r23(i,j+1,k  ) ) &
                   + ( s23(i,j,k+1)*r23(i,j,k+1) + s23(i,j+1,k+1)*r23(i,j+1,k+1) ) )

      m11(i,j,k) =  rho(i,j,k)*nba_a*( 2.0*sqrt( 2.0*smn )*s11(i,j,k)        &
                   +       ( nba_c1*( ( tmp11 + tmp12 ) + ( tmp13 - smn*odt ) )  &
                           + nba_c2*( -2.0*( tmr12 + tmr13 ) )                   &
                         ) )
      m22(i,j,k) =  rho(i,j,k)*nba_a*( 2.0*sqrt( 2.0*smn )*s22(i,j,k)        &
                   +       ( nba_c1*( ( tmp22 + tmp12 ) + ( tmp23 - smn*odt ) )  &
                           + nba_c2*(  2.0*( tmr12 - tmr23 ) )                   &
                         ) )

      m33(i,j,k) =  rho(i,j,k)*nba_a*( 2.0*sqrt( 2.0*smn )*s33(i,j,k)        &
                   +       ( nba_c1*( ( tmp13 + tmp23 ) + ( tmp33 - smn*odt ) )  &
                           + nba_c2*(  2.0*( tmr13 + tmr23 ) )                   &
                         ) )

!!!      t11(i,j,k) = rho(i,j,k)*nba_a*(   2.0*tkk*s11(i,j,k) )
!!!      t22(i,j,k) = rho(i,j,k)*nba_a*(   2.0*tkk*s22(i,j,k) )
!!!      t33(i,j,k) = rho(i,j,k)*nba_a*(   2.0*tkk*s33(i,j,k) )

    enddo
    enddo

    do j=1,nj+1
    do i=1,ni+1

      delta = 0.125*( ( ( grdscl(i  ,j  ,k  )     &
                         +grdscl(i-1,j-1,k  ) )   &
                       +( grdscl(i-1,j  ,k  )     &
                         +grdscl(i  ,j-1,k  ) ) ) &
                     +( ( grdscl(i  ,j  ,k+1)     &
                         +grdscl(i-1,j-1,k+1) )   &
                       +( grdscl(i-1,j  ,k+1)     &
                         +grdscl(i  ,j-1,k+1) ) ) )

      nba_a = nba_cs*nba_cs*delta*delta

      tmp11 = 0.25*(  (s11(i  ,j,k)+s11(i-1,j-1,k)) &
                     +(s11(i-1,j,k)+s11(i  ,j-1,k)) )
      tmp22 = 0.25*(  (s22(i  ,j,k)+s22(i-1,j-1,k)) &
                     +(s22(i-1,j,k)+s22(i  ,j-1,k)) )
      tmp33 = 0.25*(  (s33(i  ,j,k)+s33(i-1,j-1,k)) &
                     +(s33(i-1,j,k)+s33(i  ,j-1,k)) )
      tmp12 = s12(i,j,k)
      tmp13 = 0.25*( (s13(i,j,k  )+s13(i,j-1,k  )) &
                    +(s13(i,j,k+1)+s13(i,j-1,k+1)) )
      tmp23 = 0.25*( (s23(i,j,k  )+s23(i-1,j,k  )) &
                    +(s23(i,j,k+1)+s23(i-1,j,k+1)) )

      tmr12 = r12(i,j,k)
      tmr13 = 0.25*( (r13(i,j,k  )+r13(i,j-1,k  )) &
                    +(r13(i,j,k+1)+r13(i,j-1,k+1)) )
      tmr23 = 0.25*( (r23(i,j,k  )+r23(i-1,j,k  )) &
                    +(r23(i,j,k+1)+r23(i-1,j,k+1)) )

      smn = ( tmp11*tmp11 + tmp22*tmp22 ) + tmp33*tmp33 + 2.0*( tmp12*tmp12 + ( tmp13*tmp13 + tmp23*tmp23 ) )

      m12(i,j,k) = nba_a*( 2.0*sqrt( 2.0*smn )*s12(i,j,k)  &
                   +       ( nba_c1*( ( tmp11*tmp12 + tmp12*tmp22 ) + tmp13*tmp23 )  &
                           + nba_c2*( (   tmp11*tmr12 - tmp22*tmr12 )                &
                                     +( - tmp13*tmr23 - tmp23*tmr13 ) )              &
                         ) )
      m12(i,j,k) = m12(i,j,k)*0.25*( (rho(i-1,j-1,k)+rho(i,j,k))+(rho(i-1,j,k)+rho(i,j-1,k)) )

!!!      t12(i,j,k) = nba_a*(   2.0*tkk*s12(i,j,k) )
!!!      t12(i,j,k) = t12(i,j,k)*0.25*( (rho(i-1,j-1,k)+rho(i,j,k))+(rho(i-1,j,k)+rho(i,j-1,k)) )

    enddo
    enddo

  kcheck:  IF( k.ge.2 )THEN

    do j=1,nj
    do i=1,ni+1

      delta = 0.5*( grdscl(i-1,j,k) &
                   +grdscl(i  ,j,k) )

      nba_a = nba_cs*nba_cs*delta*delta

      tmp11 = 0.5*( c1(1,1,k)*(s11(i-1,j,k-1)+s11(i,j,k-1)) &
                   +c2(1,1,k)*(s11(i-1,j,k  )+s11(i,j,k  )) )
      tmp22 = 0.5*( c1(1,1,k)*(s22(i-1,j,k-1)+s22(i,j,k-1)) &
                   +c2(1,1,k)*(s22(i-1,j,k  )+s22(i,j,k  )) )
      tmp33 = 0.5*( c1(1,1,k)*(s33(i-1,j,k-1)+s33(i,j,k-1)) &
                   +c2(1,1,k)*(s33(i-1,j,k  )+s33(i,j,k  )) )
      tmp12 = 0.5*( c1(1,1,k)*(s12(i,j+1,k-1)+s12(i,j,k-1)) &
                   +c2(1,1,k)*(s12(i,j+1,k  )+s12(i,j,k  )) )
      tmp13 = s13(i,j,k)
      tmp23 = 0.25*( (s23(i-1,j+1,k)+s23(i,j  ,k)) &
                    +(s23(i-1,j  ,k)+s23(i,j+1,k)) )

      tmr13 = r13(i,j,k)
      tmr23 = 0.25*( (r23(i-1,j+1,k)+r23(i,j  ,k)) &
                    +(r23(i-1,j  ,k)+r23(i,j+1,k)) )
      tmr12 = 0.5*( c1(1,1,k)*(r12(i,j+1,k-1)+r12(i,j,k-1)) &
                   +c2(1,1,k)*(r12(i,j+1,k  )+r12(i,j,k  )) )

      smn = ( tmp11*tmp11 + tmp22*tmp22 ) + tmp33*tmp33 + 2.0*( tmp12*tmp12 + ( tmp13*tmp13 + tmp23*tmp23 ) )

      m13(i,j,k) = nba_a*( 2.0*sqrt( 2.0*smn )*s13(i,j,k)  &
                   +       ( nba_c1*( ( tmp11*tmp13 + tmp13*tmp33 ) + tmp12*tmp23 )  &
                           + nba_c2*(   ( tmp11*tmr13 - tmp33*tmr13 )                &
                                      + ( tmp12*tmr23 - tmp23*tmr12 )  )             &
                         ) )
      m13(i,j,k) = m13(i,j,k)*0.5*( rf(i-1,j,k)+rf(i,j,k) )

!!!      t13(i,j,k) = nba_a*(   2.0*tkk*s13(i,j,k) )
!!!      t13(i,j,k) = t13(i,j,k)*0.5*( rf(i-1,j,k)+rf(i,j,k) )

    enddo
    enddo

    do j=1,nj+1
    do i=1,ni

      delta = 0.5*( grdscl(i,j-1,k) &
                   +grdscl(i,j  ,k) )

      nba_a = nba_cs*nba_cs*delta*delta

      tmp11 = 0.5*( c1(1,1,k)*(s11(i,j,k-1)+s11(i,j-1,k-1)) &
                   +c2(1,1,k)*(s11(i,j,k  )+s11(i,j-1,k  )) )
      tmp22 = 0.5*( c1(1,1,k)*(s22(i,j,k-1)+s22(i,j-1,k-1)) &
                   +c2(1,1,k)*(s22(i,j,k  )+s22(i,j-1,k  )) )
      tmp33 = 0.5*( c1(1,1,k)*(s33(i,j,k-1)+s33(i,j-1,k-1)) &
                   +c2(1,1,k)*(s33(i,j,k  )+s33(i,j-1,k  )) )
      tmp12 = 0.5*( c1(1,1,k)*(s12(i,j,k-1)+s12(i+1,j,k-1)) &
                   +c2(1,1,k)*(s12(i,j,k  )+s12(i+1,j,k  )) )
      tmp13 = 0.25*( (s13(i+1,j-1,k)+s13(i,j  ,k)) &
                    +(s13(i+1,j  ,k)+s13(i,j-1,k)) )
      tmp23 = s23(i,j,k)

      tmr23 = r23(i,j,k)
      tmr13 = 0.25*( (r13(i+1,j-1,k)+r13(i,j  ,k)) &
                    +(r13(i+1,j  ,k)+r13(i,j-1,k)) )
      tmr12 = 0.5*( c1(1,1,k)*(r12(i,j,k-1)+r12(i+1,j,k-1)) &
                   +c2(1,1,k)*(r12(i,j,k  )+r12(i+1,j,k  )) )

      smn = ( tmp11*tmp11 + tmp22*tmp22 ) + tmp33*tmp33 + 2.0*( tmp12*tmp12 + ( tmp13*tmp13 + tmp23*tmp23 ) )

      m23(i,j,k) = nba_a*( 2.0*sqrt( 2.0*smn )*s23(i,j,k)  &
                   +       ( nba_c1*( ( tmp22*tmp23 + tmp23*tmp33 ) + tmp12*tmp13 )  &
                           + nba_c2*(   ( tmp22*tmr23 - tmp33*tmr23 )                &
                                      + ( tmp12*tmr13 + tmp13*tmr12 )  )             &
                         ) )
      m23(i,j,k) = m23(i,j,k)*0.5*( rf(i,j-1,k)+rf(i,j,k) )

!!!      t23(i,j,k) = nba_a*(   2.0*tkk*s23(i,j,k) )
!!!      t23(i,j,k) = t23(i,j,k)*0.5*( rf(i,j-1,k)+rf(i,j,k) )

    enddo
    enddo

  ENDIF  kcheck

  ENDDO  kloop

!------------------------------------------------------------------

      IF(bbc.eq.1)THEN
        ! free slip:
        do j=1,nj+1
        do i=1,ni+1
          m13(i,j,1) = m13(i,j,2)
          m23(i,j,1) = m23(i,j,2)
        enddo
        enddo
      ELSEIF(bbc.eq.2.or.bbc.eq.3)THEN
        IF(axisymm.eq.0)THEN
!!!          do j=1,nj+1
!!!          do i=1,ni+1
!!!            m13(i,j,1) = t13(i,j,1)
!!!            m23(i,j,1) = t23(i,j,1)
!!!          enddo
!!!          enddo
        ELSE
          stop 1212
        ENDIF
      ELSE
        stop 2222
      ENDIF

!------------------------------------------------------------------

      end subroutine turbnba2


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


  END MODULE turbnba_module
