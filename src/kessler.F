  MODULE kessler_module

  implicit none

  CONTAINS


      subroutine kessler(dt,tauto,taccr,tevar,ruh,rvh,rmh,pi0,th0,tmp,   &
                         rho,rr,pp3d,th3d,prs,                        &
                         qv3d,qc3d,qr3d)
      use input
      use constants
      implicit none

      real :: dt
      double precision :: tauto,taccr,tevar
      real, dimension(ib:ie) :: ruh
      real, dimension(jb:je) :: rvh
      real, dimension(ib:ie,jb:je,kb:ke) :: rmh,pi0,th0
      real, dimension(ib:ie,jb:je,kb:ke) :: tmp,rho,rr,pp3d,    &
                                            th3d,prs,qv3d,qc3d,qr3d

      integer :: i,j,k
      real :: qvnew,qcnew,qrnew
      real :: ar,cr,esl,qvs,er,term1,cpml,cvml,rm,tem
      double precision :: dum
      double precision, dimension(nk) :: bud1,bud2,bud3

!-------------------------------------------------------------------

!$omp parallel do default(shared)              &
!$omp private(i,j,k,qvnew,qcnew,qrnew,ar,cr,esl,qvs,er,term1,cpml,cvml,rm,tem)
      do k=1,nk
      bud1(k)=0.0d0
      bud2(k)=0.0d0
      bud3(k)=0.0d0
      do j=1,nj
      do i=1,ni

        tmp(i,j,k)=(th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k))

        qvnew=qv3d(i,j,k)
        qcnew=qc3d(i,j,k)
        qrnew=qr3d(i,j,k)

        ar=0.0
        cr=0.0
        er=0.0

        !! autoconversion of cloud to rain
        if(qcnew.gt.0.0)then
          ar=max(0.001*(qcnew-0.001),0.0)
          ar=ar*dt
        endif

        !! accretion of cloud by rain
        if(qcnew.gt.0.0 .and. qrnew.gt.0.0)then
          cr=2.2*qcnew*(qrnew**0.875)
          cr=cr*dt
        endif

        !! evap of rain to vapor
        if(qrnew.gt.0.0)then
          esl=611.2*exp( 17.67 * ( tmp(i,j,k) - 273.15 ) / ( tmp(i,j,k) - 29.65 ) )
          ! 171023 (fix for very cold temps):
          esl = min( esl , prs(i,j,k)*0.5 )
          qvs=eps*esl/(prs(i,j,k)-esl)
          if(qvnew.lt.qvs)then
            er=(1.6+30.3922*((rho(i,j,k)*qrnew)**0.2046))*          &
                (1.0-(qvnew/qvs))*                                  &
                ((rho(i,j,k)*qrnew)**0.525)/                        &
               ( (2.03e4+9.584e6/(qvs*prs(i,j,k))) * rho(i,j,k) )
            er=min(er*dt,qrnew)
            if( (qvnew+er).gt.qvs )then
              er=qvs-qvnew
            endif
          endif
        endif

        if((ar+cr).gt.qcnew)then
          term1=ar+cr
          ar=qcnew*ar/term1
          cr=qcnew*cr/term1
        endif

        qvnew=qvnew+er
        qcnew=qcnew-(ar+cr)
        qrnew=qrnew+(ar+cr-er)

      if(er.gt.1.0e-7)then
        if(eqtset.eq.2)then
          cpml=cp+cpv*qvnew+cpl*(qcnew+qrnew)
          cvml=cv+cvv*qvnew+cpl*(qcnew+qrnew)
          rm=rd+rv*qvnew
          th3d(i,j,k)=th3d(i,j,k)-er*(                                     &
                  (lv1-lv2*tmp(i,j,k))/(cpdcv*cvml*(pi0(i,j,k)+pp3d(i,j,k)))   &
                - (th0(i,j,k)+th3d(i,j,k))*(rv/cvml)*(1.0-rovcp*cpml/rm) )
          pp3d(i,j,k)=((rho(i,j,k)*rm*(th0(i,j,k)+th3d(i,j,k))*rp00)**rddcv)-pi0(i,j,k)
          prs(i,j,k)=p00*((pi0(i,j,k)+pp3d(i,j,k))**cpdrd)
        else
          th3d(i,j,k)=th3d(i,j,k)-er*( (lv1-lv2*tmp(i,j,k))         &
                                      /(cp*(pi0(i,j,k)+pp3d(i,j,k))) )
          rho(i,j,k)=prs(i,j,k)   &
             /(rd*(th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k))*(1.0+qvnew*reps))
        endif
      endif

        tem=ruh(i)*rvh(j)*rmh(i,j,k)

        bud1(k)=bud1(k)+rr(i,j,k)*ar*tem
        bud2(k)=bud2(k)+rr(i,j,k)*cr*tem
        bud3(k)=bud3(k)+rr(i,j,k)*er*tem

        qv3d(i,j,k)=qvnew
        qc3d(i,j,k)=qcnew
        qr3d(i,j,k)=qrnew

      enddo
      enddo
      enddo

      dum=dx*dy*dz

      do k=1,nk
        tauto=tauto+bud1(k)*dum
      enddo

      do k=1,nk
        taccr=taccr+bud2(k)*dum
      enddo

      do k=1,nk
        tevar=tevar+bud3(k)*dum
      enddo

      if(timestats.ge.1) time_microphy=time_microphy+mytime()
 
      end subroutine kessler


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine satadj(nrk,dt,tcond,tevac,ruh,rvh,rmh,pi0,th0,   &
                        rho,rr,pp3d,prs,th3d,q3d)
      use input
      use constants
#ifdef MPI
      use mpi
#endif
      implicit none

      integer nrk
      real, intent(in) :: dt
      double precision :: tcond,tevac
      real, dimension(ib:ie) :: ruh
      real, dimension(jb:je) :: rvh
      real, dimension(ib:ie,jb:je,kb:ke) :: rmh,pi0,th0
      real, dimension(ib:ie,jb:je,kb:ke) :: rho,rr,pp3d,prs,th3d
      real, dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: q3d

      integer :: i,j,k,n,nmax,omax,iflag
      real :: tnew,pnew
      real :: esl,qvs,qvnew,qcnew,cvml,rm,lhv,tlast,dqv
      real :: converge,t1,d1,tem,ql,qi,rdt
      double precision :: dum
      double precision, dimension(nk) :: bud1,bud2
      logical :: doit

!--------------------------------------------------------------------
!  iterative sat adj.

    nmax=0
    iflag=0

    IF(eqtset.eq.2)THEN

      if(nrk.eq.4)then
!!!        converge=0.0005
        converge=2.0*tsmall
      else
!!!        converge=0.01
        converge=20.0*tsmall
      endif

      rdt = 1.0/dt

!$omp parallel do default(shared)  &
!$omp private(i,j,k,n,tnew,pnew,esl,qvs,qvnew,qcnew,cvml,rm,lhv,   &
!$omp tlast,dqv,t1,d1,tem,ql,qi,doit)
      do k=1,nk
      bud1(k)=0.0d0
      bud2(k)=0.0d0
      do j=1,nj
      do i=1,ni

        tnew=(th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k))
        esl=611.2*exp( 17.67 * ( tnew - 273.15 ) / ( tnew - 29.65 ) )
        ! 171023 (fix for very cold temps):
        esl = min( esl , prs(i,j,k)*0.5 )
        qvs=eps*esl/(prs(i,j,k)-esl)

        IF(q3d(i,j,k,nqc).gt.1.0e-12 .or. q3d(i,j,k,nqv).gt.qvs)THEN

          qvnew=q3d(i,j,k,nqv)
          qcnew=q3d(i,j,k,nqc)
          ql=0.0
          qi=0.0
          do n=nql1,nql2
            ql=ql+q3d(i,j,k,n)
          enddo
          if(iice.eq.1)then
            do n=nqs1,nqs2
              qi=qi+q3d(i,j,k,n)
            enddo
          endif
          ql=max(0.0,ql)
          qi=max(0.0,qi)
          cvml=cv+cvv*qvnew+cpl*ql+cpi*qi
          lhv=lv1-lv2*tnew

          t1=(lhv-rv*tnew)/cvml
          d1=t1*17.67*243.5

          n=0
          tlast=tnew
          doit=.true.

          do while( doit )
            n=n+1
            dqv=(qvs-qvnew)/(1.0+d1*qvs/((tnew-29.65)**2) )
            dqv=min(dqv,qcnew)
            if(  (qvnew+dqv).lt.1.0e-20 ) dqv=1.0e-20-qvnew

            qvnew=qvnew+dqv
            qcnew=qcnew-dqv
            tnew=tnew-dqv*t1
            pnew=rho(i,j,k)*(rd+rv*qvnew)*tnew

            doit = .false.
            if( abs(tnew-tlast).gt.converge )then
              tlast=tnew
              esl=611.2*exp( 17.67 * ( tnew - 273.15 ) / ( tnew - 29.65 ) )
              ! 171023 (fix for very cold temps):
              esl = min( esl , prs(i,j,k)*0.5 )
              qvs=eps*esl/(pnew-esl)
              doit = .true.
            endif

            if(n.gt.50) print *,'  satadj:',myid,n,tnew,pnew
            if(n.eq.100)then
              print *,'  infinite loop!'
              print *,'  i,j,k=',i,j,k
              iflag=1
              doit=.false.
            endif

          enddo

          tem=ruh(i)*rvh(j)*rmh(i,j,k)

          bud1(k)=bud1(k)+rr(i,j,k)*max(qcnew-q3d(i,j,k,nqc),0.0)*tem
          bud2(k)=bud2(k)+rr(i,j,k)*max(qvnew-q3d(i,j,k,nqv),0.0)*tem

          prs(i,j,k) = pnew
          pp3d(i,j,k) = (pnew*rp00)**rovcp - pi0(i,j,k)
          th3d(i,j,k) = tnew/(pi0(i,j,k)+pp3d(i,j,k))-th0(i,j,k)
          q3d(i,j,k,nqc) = qcnew
          q3d(i,j,k,nqv) = qvnew

          nmax=max(n,nmax)

        ENDIF


      enddo
      enddo
      enddo

    ELSE

      nmax=1

!$omp parallel do default(shared)  &
!$omp private(i,j,k,qvnew,qcnew,tnew,esl,qvs,lhv,dqv,tem,rm)
      do k=1,nk
      bud1(k)=0.0d0
      bud2(k)=0.0d0
      do j=1,nj
      do i=1,ni

        qvnew=q3d(i,j,k,nqv)
        qcnew=q3d(i,j,k,nqc)
        tnew=(th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k))
        esl=611.2*exp( 17.67 * ( tnew - 273.15 ) / ( tnew - 29.65 ) )
        ! 171023 (fix for very cold temps):
        esl = min( esl , prs(i,j,k)*0.5 )
        qvs=eps*esl/(prs(i,j,k)-esl)
        lhv=lv1-lv2*tnew
!!!        dqv=(qvs-qvnew)/(1.0+lhv*qvs*4097.8531/(cp*((tnew-35.86)**2)))
        dqv=(qvs-qvnew)/(1.0+lhv*qvs*17.67*243.5/(cp*((tnew-29.65)**2)))
        if(  (qvnew+dqv).lt.1.0e-20 ) dqv=1.0e-20-qvnew
        dqv=min(dqv,max(0.0,qcnew))

        qvnew=qvnew+dqv
        qcnew=qcnew-dqv
        tem=ruh(i)*rvh(j)*rmh(i,j,k)
        bud1(k)=bud1(k)+rr(i,j,k)*max(qcnew-q3d(i,j,k,nqc),0.0)*tem
        bud2(k)=bud2(k)+rr(i,j,k)*max(qvnew-q3d(i,j,k,nqv),0.0)*tem

        th3d(i,j,k)=th3d(i,j,k)-dqv*( lhv/(cp*(pi0(i,j,k)+pp3d(i,j,k))) )
        q3d(i,j,k,nqc)=qcnew
        q3d(i,j,k,nqv)=qvnew
        rho(i,j,k)=prs(i,j,k)   &
             /(rd*(th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k))*(1.0+qvnew*reps))

      enddo
      enddo
      enddo

    ENDIF

    IF(nrk.ge.nrkmax)THEN
      dum=dx*dy*dz
      do k=1,nk
        tcond=tcond+bud1(k)*dum
      enddo

      do k=1,nk
        tevac=tevac+bud2(k)*dum
      enddo
    ENDIF

!!!#ifdef MPI
!!!      omax=0
!!!      call MPI_REDUCE(nmax,omax,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,ierr)
!!!      nmax=omax
!!!#endif
!!!      if(myid.eq.0) print *,'  nmax = ',nmax

      if(iflag.ne.0)then
        print *
        print *,' Convergence cannot be reached in satadj subroutine.'
        print *
        print *,' This may be a problem with the algorithm in satadj.'
        print *,' However, the model may have became unstable somewhere'
        print *,' else and the symptoms first appeared here.'
        print *
        print *,' Try decreasing the timestep.'
        print *
        print *,'  ... stopping cm1 ... '
        print *
        call stopcm1
      endif

      if(timestats.ge.1) time_satadj=time_satadj+mytime()

      end subroutine satadj


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine k_fallout(rho,qr3d,vr)
      use input
      use constants
      implicit none
 
      real, dimension(ib:ie,jb:je,kb:ke) :: rho,qr3d,vr
 
      integer i,j,k

!--------------------------------------------------------------------
!  Get fall velocities

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          vr(i,j,k)=0.0
          if(qr3d(i,j,k).gt.1.0e-12)then
            vr(i,j,k)=14.34*((rho(i,j,k)*qr3d(i,j,k))**0.1346)    &
                           *sqrt(1.15/rho(i,j,k))
          endif
        enddo
        enddo
        enddo

      if(timestats.ge.1) time_fall=time_fall+mytime()
 
      end subroutine k_fallout


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine fallout(dt,train,ruh,rvh,zh,mh,mf,rain,prate,rr,rho,   &
                         q3d,vq)
      use input
      use constants
      implicit none

      real :: dt
      double precision :: train
      real, dimension(ib:ie) :: ruh
      real, dimension(jb:je) :: rvh
      real, dimension(ib:ie,jb:je,kb:ke) :: zh,mh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, dimension(ib:ie,jb:je,nrain) :: rain
      real, intent(inout), dimension(ib:ie,jb:je) :: prate
      real, dimension(ib:ie,jb:je,kb:ke) :: rr,rho,q3d,vq

      integer :: i,j,k,n,nr,nrk
      real, dimension(-1:nk+3):: rq
      real, dimension(0:nk) :: ffk
      real, dimension(nk) :: qtmp

      integer :: nfall
      real :: crmax,dtfall

      double precision :: tem
      double precision, dimension(nj) :: bud

!--------------------------------------------------------------------

!$omp parallel do default(shared)            &
!$omp private(i,j,k,crmax,nfall,rq,n,dtfall,nr)
    do j=1,nj
    bud(j)=0.0d0
    do i=1,ni
      crmax = 0.0
      do k=1,nk
        crmax = max( crmax , vq(i,j,k)*dt*rdz*mh(i,j,k) )
      enddo
      nfall = max( 1 , int(crmax+1.0) )
      ! cm1r17:  following code is needed if nfall is large.
      !  - below edge of falling precip, set fallspeed to
      !    value at gridpoint above:
      do k=(nk-1),1,-1
        if( q3d(i,j,k).lt.1.0e-8 ) vq(i,j,k) = vq(i,j,k+1)
      enddo
      dtfall=dt/nfall
      do n=1,nfall
        do k=1,nk
          rq(k)=rho(i,j,k)*vq(i,j,k)*max(0.0,q3d(i,j,k))
        enddo
        rq(nk+1)=0.0
        do k=1,nk
          q3d(i,j,k)=q3d(i,j,k)+dtfall*(rq(k+1)-rq(k))   &
                                     *rdz*mh(i,j,k)/rho(i,j,k)
        enddo
        do nr=1,nrain
          rain(i,j,nr)=rain(i,j,nr)+0.1*dtfall*rq(1)
                                ! mult by 100 cm/m  (convert from m to cm)
                                ! div by 1000 kg/m^3  (density of water)
        enddo
        prate(i,j) = rq(1)
        bud(j)=bud(j)+dtfall*rr(i,j,1)*vq(i,j,1)*max(0.0,q3d(i,j,1))*ruh(i)*rvh(j)
      enddo
    enddo
    enddo

      tem=dx*dy
      do j=1,nj
        train=train+bud(j)*tem
      enddo

      if(timestats.ge.1) time_fall=time_fall+mytime()

      end subroutine fallout


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getefall(setzero,cpx,mf,t,cvm,tten,q3d,vr)
      use input
      use constants
      implicit none

      integer, intent(in) :: setzero
      real, intent(in) :: cpx
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: t,cvm
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: tten
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: q3d,vr

      integer :: i,j,k

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk-1
        if(setzero.eq.1)then
          do j=1,nj
          do i=1,ni
            tten(i,j,k)=0.0
          enddo
          enddo
        endif
        do j=1,nj
        do i=1,ni
          tten(i,j,k)=tten(i,j,k)+q3d(i,j,k)*vr(i,j,k)*(              &
                       cpx*(t(i,j,k+1)-t(i,j,k))*rdz*mf(i,j,k+1) + g  &
                                    )/cvm(i,j,k)
        enddo
        enddo
      enddo

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_fall=time_fall+mytime()

      end subroutine getefall


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine geterain(dt,cpx,lx1,erain,ruh,rvh,t,rr,q3d,vr)
      use input
      use constants
      implicit none

      real :: dt,cpx,lx1
      double precision :: erain
      real, dimension(ib:ie) :: ruh
      real, dimension(jb:je) :: rvh
      real, dimension(ib:ie,jb:je,kb:ke) :: t,rr
      real, dimension(ib:ie,jb:je,kb:ke) :: q3d,vr

      integer :: i,j
      double precision :: bud(nj)

!$omp parallel do default(shared)  &
!$omp private(i,j)
      do j=1,nj
      bud(j)=0.0d0
      do i=1,ni
        bud(j)=bud(j)+dt*vr(i,j,1)*rr(i,j,1)*q3d(i,j,1)*(cpx*t(i,j,1)-lx1)*ruh(i)*rvh(j)
      enddo
      enddo

      do j=1,nj
        erain=erain+bud(j)*dx*dy
      enddo

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_fall=time_fall+mytime()

      end subroutine geterain


  END MODULE kessler_module
