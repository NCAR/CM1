  MODULE goddard_module

  implicit none

  private
  public :: satadj_ice,goddard,consat,consat2,  &
            ROQR,ROQG,ROQS,TNW,TNG,TNSS,T0K,T00K,RT0

      real TNW,TNSS,TNG,ROQR,ROQS,ROQS_POTTER,ROQG

      real C76,C358,C172,C409,C218,C580,C380,C610,C149,C879,C141

      real ZRC,ZGC,ZSC,VRC,VGC,VSC

      real AGG,BGG,ASS,BSS,AWW,BWW,BGH,BGQ,BSH,BSQ,BWH,BWQ

      real ALV,ALF,ALS,T0K,T00K,AVC,AFC,ASC,RN1,BND1,RN2,BND2,             &
         RN3,RN4,RN5,RN6,RN7,RN8,RN9,RN10,RN101,RN10A,RN11,RN11A,          &
         RN12,RN14,RN15,RN15A,RN16,RN17,                                   &
         RN17A,RN17B,RN17C,RN18,RN18A,RN19,RN19A,RN19B,RN20,RN20A,RN20B,   &
         BND3,RN21,RN22,RN23,RN23A,RN23B,RN25,RN30A,RN30B,                 &
         RN30C,RN31,BETA,RN32,AMI40,AMI50,AMI100,RT0,BW3,BS3,BG3,BWH5,     &
         BSH5,BGH5,BW6,BS6,BETAH,R10T,R11AT,R19BT,R20T,R23T
      real, dimension(31) :: RN12A,RN12B,RN25A
      real, dimension(0:151) :: pwr81
      real, dimension(0:50001) :: pwr2

  CONTAINS

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine satadj_ice(nrk,dt,tcond,tevac,ruh,rvh,rmh,pi0,th0,   &
                            rho,rr,pp3d,prs,th3d,                  &
                            qv3d,qc3d,qr3d,                        &
                            qi3d,qs3d,qg3d)
      use input
      use constants
    use cm1libs , only : rslf,rsif
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
      real, dimension(ibm:iem,jbm:jem,kbm:kem) :: qv3d,qc3d,qr3d
      real, dimension(ibi:iei,jbi:jei,kbi:kei) :: qi3d,qs3d,qg3d

      integer :: i,j,k,n,nmax,omax,iflag
      real :: tnew,pnew,qvs,qvnew,qcnew,qinew,cvml,rm,lhv,lhs
      real :: fliq,fice,tem,tlast,dqv,qsw,qsi,cnd,dep,term1,term2
      real :: converge,t1,t2,dum,rdt
      double precision :: tem6
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
!$omp private(i,j,k,n,tnew,pnew,qvs,qvnew,qcnew,qinew,cvml,rm,   &
!$omp lhv,lhs,fliq,fice,tem,tlast,dqv,qsw,qsi,cnd,dep,term1,term2,      &
!$omp t1,t2,dum,doit)
      do k=1,nk
      bud1(k)=0.0d0
      bud2(k)=0.0d0
      do j=1,nj
      do i=1,ni

        tnew=(th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k))
        fliq=max(min((tnew-t00k)*rt0,1.0),0.0)
        fice=1.0-fliq
        qsw=0.0
        if(tnew.gt.t00k)then
          qsw=fliq*rslf(prs(i,j,k),tnew)
        endif
        qsi=0.0
        if(tnew.lt.t0k)then
          qsi=fice*rsif(prs(i,j,k),tnew)
        endif
        qvs=qsw+qsi

        IF(qc3d(i,j,k).gt.1.0e-12 .or. qi3d(i,j,k).gt.1.0e-12     &
           .or. qv3d(i,j,k).gt.qvs)THEN

          qvnew=qv3d(i,j,k)
          qcnew=qc3d(i,j,k)
          qinew=qi3d(i,j,k)

          tem=cpl*(qc3d(i,j,k)+qr3d(i,j,k))               &
             +cpi*(qi3d(i,j,k)+qs3d(i,j,k)+qg3d(i,j,k))
          cvml=cv+cvv*qv3d(i,j,k)+tem
          lhv=lv1-lv2*tnew
          lhs=ls1-ls2*tnew

          t1=1.0/cvml
          t2=rv*tnew/cvml

          n=0
          tlast=tnew
          doit=.true.

          do while( doit )
            n=n+1
            term1=0.0
            if(tnew.gt.t00k)then
              term1=qsw*(fliq*C409/((tnew-C358)**2))
            endif
            term2=0.0
            if(tnew.lt.t0k)then
              term2=qsi*(fice*C580/((tnew-C76)**2))
            endif
            dqv=(qvs-qvnew)/(1.0+((lhv*fliq+lhs*fice)*t1-t2)*(term1+term2) )
            dqv=min(dqv,qcnew+qinew)
            if( (dqv*qcnew.gt.0.0.and.qinew.eq.0.0) .or.   &
                (dqv*qinew.gt.0.0.and.qcnew.eq.0.0) )then
              fliq=qcnew/(qcnew+qinew)
              fice=1.0-fliq
            endif
            if(  (qvnew+dqv).lt.1.0e-20 ) dqv=1.0e-20-qvnew
            cnd=min(fliq*dqv,qcnew)
            dep=min(fice*dqv,qinew)
            dqv=cnd+dep

            qvnew=qvnew+dqv
            qcnew=max(qcnew-cnd,0.0)
            qinew=max(qinew-dep,0.0)

            tnew=tnew-( (lhv*cnd+lhs*dep)*t1 - dqv*t2 )
            pnew=rho(i,j,k)*(rd+rv*qvnew)*tnew

            doit = .false.
            if( abs(tnew-tlast).gt.converge )then
              tlast=tnew
              fliq=max(min((tnew-t00k)*rt0,1.0),0.0)
              fice=1.0-fliq
              qsw=0.0
              if(tnew.gt.t00k)then
                qsw=fliq*rslf(prs(i,j,k),tnew)
              endif
              qsi=0.0
              if(tnew.lt.t0k)then
                qsi=fice*rsif(prs(i,j,k),tnew)
              endif
              qvs=qsw+qsi
              doit = .true.
            endif

            if(n.gt.50) print *,'  satadj_ice:',myid,n,tnew,pnew
            if(n.eq.100)then
              print *,'  infinite loop!'
              print *,'  i,j,k=',i,j,k
              iflag=1
              doit=.false.
            endif

          enddo

          dum=ruh(i)*rvh(j)*rmh(i,j,k)

          bud1(k)=bud1(k)+rr(i,j,k)*max(qcnew-qc3d(i,j,k),0.0)*dum
          bud1(k)=bud1(k)+rr(i,j,k)*max(qinew-qi3d(i,j,k),0.0)*dum
          bud2(k)=bud2(k)-rr(i,j,k)*min(qcnew-qc3d(i,j,k),0.0)*dum
          bud2(k)=bud2(k)-rr(i,j,k)*min(qinew-qi3d(i,j,k),0.0)*dum
          
          prs(i,j,k) = pnew
          pp3d(i,j,k) = (pnew*rp00)**rovcp - pi0(i,j,k)
          th3d(i,j,k) = tnew/(pi0(i,j,k)+pp3d(i,j,k))-th0(i,j,k)
          qc3d(i,j,k)=qcnew
          qv3d(i,j,k)=qvnew
          qi3d(i,j,k)=qinew
          
          nmax=max(n,nmax)

        ENDIF

      enddo
      enddo
      enddo

    ELSE

      nmax=1

!$omp parallel do default(shared)  &
!$omp private(i,j,k,qvnew,qcnew,qinew,tnew,fliq,fice,qvs,lhv,lhs,dqv,   &
!$omp qsw,qsi,cnd,dep,term1,term2,dum,rm)
      do k=1,nk
      bud1(k)=0.0d0
      bud2(k)=0.0d0
      do j=1,nj
      do i=1,ni

        qvnew=qv3d(i,j,k)
        qcnew=qc3d(i,j,k)
        qinew=qi3d(i,j,k)
        tnew=(th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k))
        fliq=max(min((tnew-t00k)*rt0,1.0),0.0)
        fice=1.0-fliq
        qsw=0.0
        term1=0.0
        if(tnew.gt.t00k)then
          qsw=fliq*rslf(prs(i,j,k),tnew)
          term1=qsw*(fliq*C409/((tnew-C358)**2))
        endif
        qsi=0.0
        term2=0.0
        if(tnew.lt.t0k)then
          qsi=fice*rsif(prs(i,j,k),tnew)
          term2=qsi*(fice*C580/((tnew-C76)**2))
        endif
        qvs=qsw+qsi
        lhv=lv1-lv2*tnew
        lhs=ls1-ls2*tnew
        dqv=(qvs-qvnew)/(1.0+(lhv*fliq+lhs*fice)*(term1+term2)*rcp)
        dqv=min(dqv,qcnew+qinew)
        if( (dqv*qcnew.gt.0.0.and.qinew.eq.0.0) .or.   &
            (dqv*qinew.gt.0.0.and.qcnew.eq.0.0) )then
          fliq=qcnew/(qcnew+qinew)
          fice=qinew/(qcnew+qinew)
        endif
        if(  (qvnew+dqv).lt.1.0e-20 ) dqv=1.0e-20-qvnew
        cnd=min(fliq*dqv,qcnew)
        dep=min(fice*dqv,qinew)

        qvnew=qvnew+cnd+dep
        qcnew=qcnew-cnd
        qinew=qinew-dep

        dum=ruh(i)*rvh(j)*rmh(i,j,k)

        bud1(k)=bud1(k)+rr(i,j,k)*max(qcnew-qc3d(i,j,k),0.0)*dum
        bud1(k)=bud1(k)+rr(i,j,k)*max(qinew-qi3d(i,j,k),0.0)*dum
        bud2(k)=bud2(k)-rr(i,j,k)*min(qcnew-qc3d(i,j,k),0.0)*dum
        bud2(k)=bud2(k)-rr(i,j,k)*min(qinew-qi3d(i,j,k),0.0)*dum

        th3d(i,j,k)=th3d(i,j,k)-( (lhv*cnd+lhs*dep)   &
                                 /(cp*(pi0(i,j,k)+pp3d(i,j,k))) )
        qc3d(i,j,k)=qcnew
        qv3d(i,j,k)=qvnew
        qi3d(i,j,k)=qinew
        rho(i,j,k)=prs(i,j,k)   &
           /(rd*(th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k))*(1.0+qvnew*reps))

      enddo
      enddo
      enddo

    ENDIF

    IF(nrk.ge.nrkmax)THEN
      tem6=dx*dy*dz
      do k=1,nk
        tcond=tcond+bud1(k)*tem6
      enddo

      do k=1,nk
        tevac=tevac+bud2(k)*tem6
      enddo
    ENDIF

!!!#ifdef MPI
!!!      omax=0
!!!      call MPI_REDUCE(nmax,omax,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,ierr)
!!!      nmax=omax
!!!#endif
!!!

      if(iflag.ne.0)then
        print *
        print *,' Convergence cannot be reached in satadj_ice subroutine.'
        print *
        print *,' This may be a problem with the algorithm in satadj_ice.'
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

      END subroutine satadj_ice


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine goddard(dt,tauto,taccr,tevar,ruh,rvh,rmh,pi0,th0,   &
                         rhod,rr,prs,pp3d,th3d,                   &
                         qv3d,qc3d,qr3d,vtr,                      &
                         qi3d,vti,qs3d,vts,                       &
                         qg3d,vtg)
      use input
      use constants
      implicit none

      real :: dt
      double precision :: tauto,taccr,tevar
      real, dimension(ib:ie) :: ruh
      real, dimension(jb:je) :: rvh
      real, dimension(ib:ie,jb:je,kb:ke) :: rmh,pi0,th0
      real, dimension(ib:ie,jb:je,kb:ke) :: rhod,rr,prs,pp3d,th3d,qv3d,qc3d,   &
                                            qr3d,vtr,qi3d,vti,qs3d,vts,qg3d,vtg

!----------------------------------------------------------------------------
!                                                                       
!***********************************************************************
!     LIN ET AL (83) ICE PHASE MICROPHYSICAL PROCESSES                 *
!     MODIFIED AND CODED BY GODDARD CUMULUS ENSEMBLE MODELING GROUP    *
!     TAO, SIMPSON AND MCCUMBER`S SATURATION TECHNIQUE (MWR, 1989)     *
!     TAO AND SIMPSON (JAS, 1989; TAO, 1993),                          *
!     BRAUN AND TAO (MWR, 2000), BRAUN ET AL. (JAM, 2001, submitted)   *
!                                                                      *
!     D2T - LEAPFROG TIME STEP (S)                                     *
!     DPT - DEVIATION POTENTIAL TEMPERATURE FIELD (K)                  *
!     DQV - DEVIATION WATER VAPOR FIELD (G/G)                          *
!     QCL - CLOUD WATER FIELD (G/G)                                    *
!     QRN - RAIN FIELD (G/G)                                           *
!     QCI - CLOUD ICE FIELD (G/G)                                      *
!     QCS - SNOW FIELD (G/G)                                           *
!     QCG - HAIL FIELD (G/G)                                           *
!     RHO - AIR DENSITY (G/CM3)                                        *
!     TA1 - BASE AIR POTENTIAL TEMPERATURE AT THE LEVEL (K)            *
!     QA1 - BASE WATER VAPOR AT THE LEVEL (G/G)                        *
!     P0 - AIR PRESSURE (UB=1.E-3MB)                                   *
!     PI - EXNER FUNCTION (P/1000.E3)**R/CP                            *
!                                                                      *
!     MIX: DIMENSION IN X-DIRECTION                                    *
!     MKX: DIMENSION IN Z-DIRECTION                                    *
!     NOTE: PHYSICAL DOMAIN EXTENDS FROM K=1 TO K=MKX-1                *
!           AND I=2 TO I=MIX-1                                         *
!     C.G.S UNITS                                                      *
!                                                                      *
!     THIS ICE SCHEME HAS BEEN IMPLEMENTED INTO                        *
!     ARPS (CAPS - V. WONG, M. XUE, K. DROEGEMEIER) - 1993             *
!     GODDARD`S SCHLENSINGER MODEL (SIMPSON ET AL., 1991)              *
!     CHEN`S MODEL (C. CHEN) - 1992                                    *
!     GMASS (BIAK ET AL., 1992)                                        *
!     MM5V1 (LIU ET AL., 1994)                                         *
!                                                                      *
!     RECENT CHANGES (6.25.01) BY SCOTT BRAUN INCLUDE                  *
!       * OPTION TO CHOOSE HAIL OR GRAUPEL AS THIRD ICE CATEGORY       *
!	* OPTION OF SATURATION ADJUSTMENT SCHEMES                      *      
!	* NEW FORMULATION OF PSFI                                      *      
!	* REDUCED COLLECTION EFFICIENCIES (EGS) FOR COLLECTION OF SNOW *      
!	  BY GRAUPEL                                                   *      
!	* MODIFICATIONS TO PIDEP, PINT INCLUDING NEW COEFFICIENTS FOR  *      
!	  THE FLETCHER EQ., A LIMIT OF 1 CM**-3 FOR THE NUMBER         *      
!         CONCENTRATION OF ICE NUCLEI, INCLUSION OF A PISUB TERM       *
!	* ELIMINATED THE PGAUT TERM                                    *      
!***********************************************************************
!                                                                       
!     THE FOLLOWINGS ARE LOCALLY USED VARIABLES                         
!                                                                       
      real d2t,VFR,VFS,VFG
      real RHO,P0,PPI,rho0
      integer i,j,k
      real COL,DEP,RGMP,DD,DD1,QVS,DM,RSUB1,WGACR,CND,RQ,             &
       ERN,SCV,TCA,DWV,ZR,VR,ZS,VS,ZG,VG,EGS,ESI,                     &
       QSI,SSI,QSW,SSW,PIHOM,PIDW,PIMLT,PSAUT,PSACI,PSACW,QSACW,      &
       PRACI,PIACR,PRAUT,PRACW,PSFW,PSFI,DGACS,DGACW,DGACI,DGACR,     &
       PGACS,WGACS,QGACW,WGACI,QGACR,PGWET,PGAUT,PRACS,PSACR,QSACR,   &
       PGFR,PSMLT,PGMLT,PSDEP,PSSUB,PGSUB,PINT,PIDEP,PISUB,           &
       PT,QV,QC,QR,QI,QS,QG,TAIR,TAIRC,PR,PS,PG,PRN,                  &
       PSN,DLT1,DLT2,DLT3,RTAIR,DDA,DDB,Y1,Y2,Y3,Y4,Y5,FV

      integer index,IT
      real cmin,cmin1,vgcr,vgcf,r3456,        &
           a2,ee1,ee2,a1,r7r,r8r,del,r11rt,nci
      real tem1,tem2,tem3,tem4
      real rrho,rp0,rpi
      real ene1,ene2,lhv,lhs,lhf
      real tem5,rdt,temp,f1,f2,r1,dum1,cpml,cvml,rm,dum

      double precision :: tem6
      double precision, dimension(nk) :: bud1,bud2,bud3

      parameter(cmin=1.e-14)
      parameter(cmin1=1.e-12)
      parameter(r1=1.0/1.0e-9)

      tem5=.5/(SQRT(AMI100)-SQRT(AMI40))
      d2t=dt
      rdt=1.0/D2T

!     ******   TWO CLASS WATER AND THREE CLASS OF ICE-PHASE    *********


!$omp parallel do default(shared)  &
!$omp private(i,j,k,VFR,VFS,VFG,RHO,P0,PPI,rho0,COL,DEP,RGMP,DD,DD1,QVS,    &
!$omp DM,RSUB1,WGACR,CND,RQ,ERN,SCV,TCA,DWV,ZR,VR,ZS,VS,ZG,VG,EGS,ESI,QSI,  &
!$omp SSI,QSW,SSW,PIHOM,PIDW,PIMLT,PSAUT,PSACI,PSACW,QSACW,PRACI,PIACR,     &
!$omp PRAUT,PRACW,PSFW,PSFI,DGACS,DGACW,DGACI,DGACR,PGACS,WGACS,QGACW,      &
!$omp WGACI,QGACR,PGWET,PGAUT,PRACS,PSACR,QSACR,PGFR,PSMLT,PGMLT,PSDEP,     &
!$omp PSSUB,PGSUB,PINT,PIDEP,PISUB,PT,QV,QC,QR,QI,QS,QG,TAIR,TAIRC,PR,PS,   &
!$omp PG,PRN,PSN,DLT1,DLT2,DLT3,RTAIR,DDA,DDB,Y1,Y2,Y3,Y4,Y5,FV,index,IT,   &
!$omp vgcr,vgcf,r3456,a2,ee1,ee2,a1,r7r,r8r,del,r11rt,nci,tem1,tem2,tem3,   &
!$omp tem4,rrho,rp0,rpi,ene1,ene2,lhv,lhs,lhf,temp,f1,f2,dum1,cpml,cvml,    &
!$omp rm,dum)
      DO k=1,nk
      bud1(k)=0.0d0
      bud2(k)=0.0d0
      bud3(k)=0.0d0
      DO j=1,nj
      DO i=1,ni
          rho0=1.0e-3*rhod(i,j,1)

          vtr(i,j,k)=0.0
          vts(i,j,k)=0.0
          vtg(i,j,k)=0.0

        IT=1                                                           
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C    ******************************************************************
          fv=rho0                                               
          RHO=1.0e-3*rhod(i,j,k)
          PPI=pi0(i,j,k)+pp3d(i,j,k)
          P0=10.0*prs(i,j,k)
          PT=th3d(i,j,k)
          QV=QV3D(I,J,K)
          QC=QC3D(I,J,K)
          QR=QR3D(I,J,K)
          QI=QI3D(I,J,K)
          QS=QS3D(I,J,K)
          QG=QG3D(I,J,K)
          QV=max(QV,cmin1)
          IF(QC.LT.cmin)  QC=0.0                                  
          IF(QR.LT.cmin)  QR=0.0                                  
          IF(QI.LT.cmin)  QI=0.0                                  
          IF(QS.LT.cmin)  QS=0.0                                  
          IF(QG.LT.cmin)  QG=0.0                                  
          TAIR=(th0(i,j,k)+PT)*PPI                                     
          TAIRC=TAIR-T0K                                          
!----------------------------------------------------
!  variables added to improve performance
          tem1=1.0/(sqrt(sqrt(rho)))
          tem2=sqrt(FV/RHO)
          tem3=1.0/(sqrt(rho))
          tem4=sqrt(sqrt(FV/RHO))
          rrho=1.0/rho
          rp0=1.0/p0
          rpi=1.0/ppi
!----------------------------------------------------
!  variables for new equation set
          if(eqtset.eq.2)then
            dum1=cpl*(qc+qr)+cpi*(qi+qs+qg)
            cvml=cv+cvv*qv+dum1
            cpml=cp+cpv*qv+dum1
            rm=rd+rv*qv
            ene1=cv/(cp*cvml*ppi)
            ene2=(th0(i,j,k)+pt)*(rv/cvml)*(1.0-rovcp*cpml/rm)
          else
            ene1=1.0/(cp*ppi)
            ene2=0.0
          endif
          lhv=lv1-lv2*tair
          lhs=ls1-ls2*tair
          lhf=lhs-lhv
!----------------------------------------------------
!     ***   COMPUTE ZR,ZS,ZG,VR,VS,VG      *****************************
          ZR=1.E5*ZRC*tem1
          ZS=1.E5*ZSC*tem1
          ZG=1.E5*ZGC*tem1
          VR=0.0                                                     
          VS=0.0                                                     
          VG=0.0                                                     
          IF(QR.GT.cmin)THEN                                         
            DD=RHO*QR                                          
            Y1=sqrt(sqrt(DD))
            ZR=ZRC/Y1                                             
            temp=r1*min(50.0e-6,max(0.0,DD))
            index=int(temp)
            f1=pwr2(index)
            f2=pwr2(index+1)
            VR=VRC*tem2*(f1+(f2-f1)*(temp-index))
          ENDIF                                                         
          IF(QS.GT.cmin)THEN                                         
            DD=RHO*QS                                          
            Y1=sqrt(sqrt(DD))
            ZS=ZSC/Y1                                             
!!!            VS=VSC*tem2*sqrt(sqrt(y1))
            VS=VSC*tem2*DD**BSQ
          ENDIF                                                         
          if (qg .gt. cmin) then                                     
            dd=rho*qg                                         
            y1=sqrt(sqrt(dd))
            zg=zgc/y1                                             
            if (ihail .eq. 1) then                                      
              vgcr=vgc*tem3
              vg=vgcr*sqrt(y1)
            else                                                        
              vgcf=vgc*tem2
              vg=vgcf*exp(alog(dd)*bgq)
            endif                                                       
          endif                                                         
          VFR=VR                                                
          VFS=VS                                                
          VFG=VG                                                
!     ******************************************************************
!     ***   Y1 : DYNAMIC VISCOSITY OF AIR (U)                           
!     ***   DWV : DIFFUSIVITY OF WATER VAPOR IN AIR (PI)                
!     ***   TCA : THERMAL CONDUCTIVITY OF AIR (KA)                      
!     ***   Y2 : KINETIC VISCOSITY (V)                                  
          Y1=C149*TAIR*sqrt(TAIR)/(TAIR+120.)
          temp=min(150.0,max(0.0,TAIR-173.15))
          index=int(temp)
          f1=pwr81(index)
          f2=pwr81(index+1)
          DWV=C879*rp0*TAIR*(f1+(f2-f1)*(temp-index))
          TCA=C141*Y1                                             
          SCV=exp(alog(RHO/(Y1*DWV**2))*.1666667)
!*  1 * PSAUT : AUTOCONVERSION OF QI TO QS                        ***1**
!*  3 * PSACI : ACCRETION OF QI TO QS                             ***3**
!*  4 * PSACW : ACCRETION OF QC BY QS (RIMING) (QSACW FOR PSMLT)  ***4**
!*  5 * PRACI : ACCRETION OF QI BY QR                             ***5**
!*  6 * PIACR : ACCRETION OF QR OR QG BY QI                       ***6**
          PSAUT=0.0                                                  
          PSACI=0.0                                                  
          PRACI=0.0                                                  
          PIACR=0.0                                                  
          PSACW=0.0                                                  
          QSACW=0.0                                                  
          PRAUT=0.0                                                  
          PRACW=0.0                                                  
          r3456=tem2
!!!          DD=1.0/(sqrt(sqrt(ZS))*ZS**3)
          DD=1.0/ZS**BS3
          IF( (QI.GT.CMIN) .or. (QC.GT.CMIN) )THEN
            temp=exp(alog(ZR)*(-BW3))
          ENDIF
          IF(QI.GT.CMIN) THEN                                        
           IF(TAIR.LT.T0K)THEN                                       
             ESI=EXP(.025*TAIRC)                                  
             PSAUT=MAX(RN1*ESI*(QI-BND1),0.0)                
             PSACI=RN3*R3456*.1*QI*DD                          
             PRACI=RN5*R3456*QI*temp
             PIACR=RN6*R3456*QI*exp(alog(ZR)*(-BW6))
           ENDIF                                                        
          ENDIF                                                         
          IF(QC.GT.CMIN) THEN                                        
            IF(TAIR.LT.T0K)THEN                                      
              PSACW=RN4*R3456*QC*DD                            
            ELSE                                                        
              QSACW=RN4*R3456*QC*DD                            
            ENDIF                                                       
!* 21 * PRAUT : AUTOCONVERSION OF QC TO QR                        **21**
!* 22 * PRACW : ACCRETION OF QC BY QR                             **22**
           PRACW=RN22*tem2*QC*temp
           IF(iautoc.eq.1)THEN
             Y1=QC-BND3                                             
             IF(Y1.GT.cmin)THEN                                        
               PRAUT=RHO*Y1*Y1/(1.2E-4+RN21/Y1)            
             ENDIF                                                        
           ENDIF
          ENDIF                                                         
!* 12 * PSFW : BERGERON PROCESSES FOR QS (KOENING, 1971)          **12**
!* 13 * PSFI : BERGERON PROCESSES FOR QI                          **13**
!* 32 * PIDEP : DEPOSITION OF QI                                  **32**
          PSFW=0.0                                                   
          PSFI=0.0                                                   
          PIDEP=0.0                                                  
          PISUB=0.0                                                  
          IF(TAIR.LT.T0K.AND.QI.GT.cmin)THEN                      
            Y1=MAX(MIN(TAIRC,-1.),-31.)                       
            IT=INT(ABS(Y1))                                       
            Y1=RN12A(IT)                                          
            Y2=RN12B(IT)                                          
            PSFW=MAX(D2T*Y1*(RN12*RHO*QC)*QI,0.0)      
            RTAIR=1./(TAIR-C76)                                   
            Y2=EXP(C218-C580*RTAIR)                               
            QSI=C380*rp0*Y2                                   
            ESI=C610*Y2                                           
            SSI=QV/QSI-1.                                
            NCI=MIN(RN25*EXP(BETA*TAIRC),1.)                       
            DM=QV-QSI                                      
            RSUB1=C580*ASC*QSI*RTAIR*RTAIR                  
            Y3=1./TAIR                                            
            DD=Y3*(RN30A*Y3-RN30B)+RN30C*TAIR/ESI        
            Y1=206.18*SSI/DD                                   
            DEP=DM/(1.+RSUB1)*rdt
            PIDEP=Y1*SQRT(NCI*QI*rrho)                       
            IF(DM.GT.cmin) THEN                                      
               a2=1.                                                    
               if(pidep.gt.dep .and. pidep .gt. cmin) then     
                  a2=dep/pidep                                    
                  pidep=dep                                       
               endif                                                    
               PSFI=a2*QI*Y1*tem5
            ELSEIF(DM.LT.-cmin) THEN                                 
               PISUB=MAX(-QI*rdt,PIDEP)                      
               PISUB=MIN(-PISUB,-DM*rdt)                     
               PIDEP=0.                                              
               PSFI=0.                                               
            ELSE                                                        
               PISUB=0.                                              
               PIDEP=0.                                              
               PSFI=0.                                               
            ENDIF                                                       
          ENDIF                                                         
!TTT***** QG=QG+MIN(PGDRY,PGWET)                                        
!*  9 * PGACS : ACCRETION OF QS BY QG (DGACS,WGACS: DRY AND WET)  ***9**
!* 14 * DGACW : ACCRETION OF QC BY QG (QGACW FOR PGMLT)           **14**
!* 16 * DGACR : ACCRETION OF QR TO QG (QGACR FOR PGMLT)           **16**
          DGACS=0.0                                                  
          WGACS=0.0                                                  
          DGACW=0.0                                                  
          DGACR=0.0                                                  
          PGACS=0.0                                                  
          QGACW=0.0                                                  
          QGACR=0.0                                                  
          DGACI=0.0                                                  
          WGACI=0.0                                                  
          PGWET=0.0                                                  
        IF(QG.GT.CMIN)THEN                                           
                                                                        
          IF(QC+QR.LT.1.e-4) THEN
             EE1=.01
          ELSE
             EE1=1.
          ENDIF
          EE2=0.09                                                      
          EGS=EE1*EXP(EE2*TAIRC)                                  
          IF(TAIR.GE.T0K)EGS=1.0                                  
          Y1=ABS(VG-VS)                                        
          Y2=ZS*ZG                                             
          Y3=5./Y2                                                
          Y4=.08*Y3*Y3                                         
          Y5=.05*Y3*Y4                                         
          DD=Y1*(Y3/ZS**5+Y4/ZS**3+Y5/ZS)       
          PGACS=RN9*rrho*EGS*DD                              
          DGACS=PGACS                                             
          WGACS=RN9*rrho*DD                                     
          IF(IHAIL.EQ.1) THEN                                           
             Y1=1.0/(sqrt(ZG)*ZG**3)
             DGACW=MAX(RN14*QC*Y1*tem3,0.0)          
          ELSE                                                          
             Y1=exp(alog(ZG)*(-BG3))
             DGACW=MAX(RN14*QC*Y1*R3456,0.0)                 
          ENDIF                                                         
          QGACW=DGACW                                             
          Y1=ABS(VG-VR)                                        
          Y2=ZR*ZG                                             
          Y3=5./Y2                                                
          Y4=.08*Y3*Y3                                         
          Y5=.05*Y3*Y4                                         
          DD=Y1*(Y3/ZR**5+Y4/ZR**3+Y5/ZR)*RN16   &
                *rrho                                                  
          DGACR=MAX(DD,0.0)                                     
          QGACR=DGACR                                             
          IF(TAIR.GE.T0K)THEN                                        
            DGACS=0.0                                                
            WGACS=0.0                                                
            DGACW=0.0                                                
            DGACR=0.0                                                
          ELSE                                                          
            PGACS=0.0                                                
            QGACW=0.0                                                
            QGACR=0.0                                                
          ENDIF                                                         
!*******PGDRY : DGACW+DGACI+DGACR+DGACS                           ******
!* 15 * DGACI : ACCRETION OF QI BY QG (WGACI FOR WET GROWTH)      **15**
!* 17 * PGWET : WET GROWTH OF QG                                  **17**
          IF(TAIR.LT.T0K.AND.TAIR.GT.T0K-40.)THEN                 
            IF(IHAIL.EQ.1) THEN                                         
               Y1=QI/(sqrt(ZG)*ZG**3)
               DGACI=Y1*RN15*tem3
               WGACI=Y1*RN15A*tem3
            ELSE                                                        
               Y1=QI/ZG**BG3                                      
               DGACI=Y1*RN15*R3456                                
               WGACI=Y1*RN15A*R3456                               
            ENDIF                                                       
            Y1=1./(ALF+RN17C*TAIRC)                               
            IF(IHAIL.EQ.1) THEN                                         
               Y3=.78/ZG**2+SCV*RN17A*tem1*exp(alog(ZG)*(-BGH5))
            ELSE                                                        
               Y3=.78/ZG**2+SCV*RN17A*tem4*exp(alog(ZG)*(-BGH5))
            ENDIF                                                       
            Y4=RHO*ALV*DWV*(C380*rp0-QV)-TCA*   &
                  TAIRC                                              
            DD=Y1*(Y4*Y3*RN17*rrho+(WGACI+WGACS)*(   &
                  ALF+RN17B*TAIRC))                                  
            PGWET=MAX(DD,0.0)                                   
          ENDIF                                                         
        ENDIF                                                           
!********   HANDLING THE NEGATIVE CLOUD WATER (QC)    ******************
!********   HANDLING THE NEGATIVE CLOUD ICE (QI)      ******************
!********   DIFFERENT FROM LIN ET AL                  ******************
           y1=qc*rdt
          psacw=MIN(y1, psacw)                               
          praut=MIN(y1, praut)                               
          pracw=MIN(y1, pracw)                               
          psfw= MIN(y1, psfw)                                
          dgacw=MIN(y1, dgacw)                               
          qsacw=MIN(y1, qsacw)                               
          qgacw=MIN(y1, qgacw)                               
                                                                        
          Y1=(PSACW+PRAUT+PRACW+PSFW+DGACW+QSACW+   &
                QGACW)*D2T                                           
          QC=QC-Y1                                             
                                                                        
          IF(QC.LT.0.0) THEN                                         
             a1=1.                                                      
              if (y1 .ne. 0.0) A1=QC/Y1+1.                     
            PSACW=PSACW*A1                                        
            PRAUT=PRAUT*A1                                        
            PRACW=PRACW*A1                                        
            PSFW=PSFW*A1                                          
            DGACW=DGACW*A1                                        
            QSACW=QSACW*A1                                        
            QGACW=QGACW*A1                                        
            QC=0.0                                                   
          ENDIF                                                         

          dum=ruh(i)*rvh(j)*rmh(i,j,k)

          bud1(k)=bud1(k)+rr(i,j,k)*PRAUT*dum
          bud2(k)=bud2(k)+rr(i,j,k)*PRACW*dum

!                                                                       
!******** SHED PROCESS (WGACR=PGWET-DGACW-WGACI-WGACS)                  
!     CALCULATIONS OF THIS TERM HAS BEEN MOVED TO THIS LOCATION TO      
!     ACCOUNT FOR RESCALING OF THE DGACW TERM ABOVE. ALTHOUGH RESCALING 
!     OF THE PGWET, WGACI AND WGACS TERMS OCCURS BELOW, THESE CHANGES   
!     ARE EXPECTED TO BE LESS THAN THAT ASSOCIATED WITH DGACW. THIS CALC
!     IS NOT DONE AFTER RESCALING OF ALL TERMS SINCE WGACR IS NEEDED BEL
!     FOR RAIN AND WGACI AND DGACI ARE NEEDED FOR ICE.                  
          WGACR=PGWET-DGACW-WGACI-WGACS                  
          Y2=DGACW+DGACI+DGACR+DGACS                     
          IF(PGWET.GE.Y2)THEN                                     
            WGACR=0.0                                                
            WGACI=0.0                                                
            WGACS=0.0                                                
          ELSE                                                          
            DGACR=0.0                                                
            DGACI=0.0                                                
            DGACS=0.0                                                
          ENDIF                                                         
                                                                        
            y1=qi*rdt
           psaut=MIN(y1, psaut)                              
           psaci=MIN(y1, psaci)                              
           praci=MIN(y1, praci)                              
           psfi= MIN(y1, psfi)                               
           dgaci=MIN(y1, dgaci)                              
           wgaci=MIN(y1, wgaci)                              

          Y1=(PSAUT+PSACI+PRACI+PSFI+DGACI+WGACI   &
                +PISUB)*D2T                                          

           qi=qi-y1+PIDEP*D2T                               
           if(qi.lt.0.0) then                                        
               a2=1.                                                    
                if (y1 .ne. 0.0) a2=qi/y1+1.                   
            psaut=psaut*a2                                        
            psaci=psaci*a2                                        
            praci=praci*a2                                        
            psfi=psfi*a2                                          
            dgaci=dgaci*a2                                        
            wgaci=wgaci*a2                                        
            pisub=pisub*a2                                        
            qi=0.0                                                   
           endif                                                        
!                                                                       
          DLT3=0.0                                                   
          DLT2=0.0                                                   
          IF(TAIR.LT.T0K)THEN                                        
            IF(QR.LT.1.E-4)THEN                                      
              DLT3=1.0                                               
              DLT2=1.0                                               
            ENDIF                                                       
            IF(QS.GE.1.E-4)DLT2=0.0                               
          ENDIF                                                         
          PR=(QSACW+PRAUT+PRACW+QGACW)*D2T               
          PS=(PSAUT+PSACI+PSACW+PSFW+PSFI+DLT3*    &
                PRACI)*D2T                                           
          PG=((1.-DLT3)*PRACI+DGACI+WGACI+DGACW)*D2T  
!*  7 * PRACS : ACCRETION OF QS BY QR                             ***7**
!*  8 * PSACR : ACCRETION OF QR BY QS (QSACR FOR PSMLT)           ***8**
            PRACS=0.0                                                
            PSACR=0.0                                                
            QSACR=0.0                                                
            PGFR=0.0                                                 
            PGAUT=0.0                                                
        IF(QR.GT.CMIN)THEN                                           
          Y1=ABS(VR-VS)                                        
          Y2=ZR*ZS                                             
          Y3=5./Y2                                                
          Y4=.08*Y3*Y3                                         
          Y5=.05*Y3*Y4                                         
          R7R=RN7*rrho                                                
          PRACS=R7R*Y1*(Y3/ZS**5+Y4/ZS**3+Y5/ZS)
          R8R=RN8*rrho                                                
          PSACR=R8R*Y1*(Y3/ZR**5+Y4/ZR**3+Y5/ZR)
          QSACR=PSACR                                             
          IF(TAIR.GE.T0K)THEN                                        
            PRACS=0.0                                                
            PSACR=0.0                                                
          ELSE                                                          
            QSACR=0.0                                                
          ENDIF                                                         
!*  2 * PGAUT : AUTOCONVERSION OF QS TO QG                        ***2**
!* 18 * PGFR  : FREEZING OF QR TO QG                              **18**
          IF(TAIR.LT.T0K)THEN                                        
!           Y1=EXP(.09*TAIRC)                                     
!           if (ihail .eq. 1) PGAUT=MAX(RN2*Y1*(QS-BND2),0.0)
            Y2=EXP(RN18A*(T0K-TAIR))                              
            PGFR=MAX((Y2-1.)*RN18*rrho/ZR**7,0.0)       
          ENDIF                                                         
        ENDIF                                                           
!********   HANDLING THE NEGATIVE RAIN WATER (QR)    *******************
!********   HANDLING THE NEGATIVE SNOW (QS)          *******************
!********   DIFFERENT FROM LIN ET AL                 *******************
          y1=qr*rdt
          Y2=-QG*rdt
         piacr=MIN(y1, piacr)                                
         dgacr=MIN(y1, dgacr)                                
         wgacr=MIN(y1, wgacr)                                
         wgacr=MAX(y2, wgacr)                                
         psacr=MIN(y1, psacr)                                
         pgfr= MIN(y1, pgfr)                                 

         del=0.                                                         
         IF(WGACR.LT.0.) del=1.                                      
          Y1=(PIACR+DGACR+(1.-del)*WGACR+PSACR+          &
                  PGFR)*D2T                                          
          qr=qr+pr-y1-del*WGACR*D2T                      
                                                                        
          if(qr.lt.0.0) then                                         
            a1=1.0                                                      
             if(y1 .ne. 0.) a1=qr/y1+1.                        
            piacr=piacr*a1                                        
            dgacr=dgacr*a1                                        
            if(wgacr.gt.0) wgacr=wgacr*a1                      
            pgfr=pgfr*a1                                          
            psacr=psacr*a1                                        
            qr=0.0                                                   
          endif                                                         

          PRN=D2T*((1.-DLT3)*PIACR+DGACR+WGACR+(1.-      &
                 DLT2)*PSACR+PGFR)                             
          PS=PS+D2T*(DLT3*PIACR+DLT2*PSACR)           

           y1=qs*rdt
          pgacs=MIN(y1, pgacs)                               
          dgacs=MIN(y1, dgacs)                               
          wgacs=MIN(y1, wgacs)                               
          pgaut=MIN(y1, pgaut)                               
          pracs=MIN(y1, pracs)                               

          PRACS=(1.-DLT2)*PRACS                                
          PSN=D2T*(PGACS+DGACS+WGACS+PGAUT+PRACS)     
          QS=QS+PS-PSN                                      
          if(qs .lt. 0.0) then                                       
            a2=1.                                                       
              if(psn .ne. 0.) a2=qs/psn+1.                     
            pgacs=pgacs*a2                                        
            dgacs=dgacs*a2                                        
            wgacs=wgacs*a2                                        
            pgaut=pgaut*a2                                        
            pracs=pracs*a2                                        
            psn=psn*a2                                            
            qs=0.0                                                   
          endif                                                         
          Y2=D2T*(PSACW+PSFW+DGACW+PIACR+DGACR+       &
                WGACR+PSACR+PGFR)                              
          PT=PT+lhf*ene1*Y2                                   
          QG=QG+PG+PRN+PSN                               
!* 11 * PSMLT : MELTING OF QS                                     **11**
!* 19 * PGMLT : MELTING OF QG TO QR                               **19**
          PSMLT=0.0                                                  
          PGMLT=0.0                                                  
          TAIR=(th0(i,j,k)+PT)*PPI                                     
          r3456=tem2
          IF(TAIR.GE.T0K .AND. (QS+QG).GT.CMIN) THEN           
            TAIRC=TAIR-T0K                                        
            Y1=TCA*TAIRC-RHO*ALV*DWV    &
                       *(C380*rp0-QV)
            Y2=.78/ZS**2+RN101*tem4*SCV    &
                 *exp(alog(ZS)*(-BSH5))
            R11RT=RN11*D2T*rrho                                       
            DD=R11RT*Y1*Y2+R11AT*TAIRC*(QSACW+QSACR)  
            PSMLT=MAX(0.0,MIN(DD,QS))                      
            IF(IHAIL.EQ.1) THEN                                         
               Y3=.78/ZG**2+SCV*RN19A*tem1*exp(alog(ZG)*(-BGH5))
            ELSE                                                        
               Y3=.78/ZG**2+SCV*RN19A*tem4*exp(alog(ZG)*(-BGH5))
            ENDIF                                                       
            DD1=Y1*Y3*RN19*D2T*rrho+R19BT*TAIRC*(QGACW   &
                   +QGACR)                                            
            PGMLT=MAX(0.0,MIN(DD1,QG))                     
            PT=PT-lhf*ene1*(PSMLT+PGMLT)
            QR=QR+PSMLT+PGMLT                               
            QS=QS-PSMLT                                        
            QG=QG-PGMLT                                        
          ENDIF                                                         
!* 24 * PIHOM : HOMOGENEOUS FREEZING OF QC TO QI (T < T00K)    **24**   
!* 25 * PIDW : DEPOSITION GROWTH OF QC TO QI ( T0K < T <= T00K)**25**   
!* 26 * PIMLT : MELTING OF QI TO QC (T >= T0K)                 **26**   
          IF(QC.LE.cmin)QC=0.0                                    
          IF(QI.LE.cmin)QI=0.0                                    
          PIHOM=0.                                                   
          PIMLT=0.                                                   
          TAIR=(th0(i,j,k)+PT)*PPI                                     
          IF((TAIR-T00K).LE.0.) PIHOM=QC                       
          IF((TAIR-T0K) .GE.0.) PIMLT=QI                       
          PIDW=0.0                                                   
          IF(TAIR.LT.T0K.AND.TAIR.GT.T00K.AND.QC.GT.CMIN)THEN  
            TAIRC=TAIR-T0K                                        
            Y1=MAX(MIN(TAIRC,-1.),-31.)                       
            Y2=RN25A(INT(ABS(Y1)))                                
            Y3=MIN(RN25*EXP(BETA*TAIRC),1.)                     
            PIDW=MIN(Y2*Y3*D2T*rrho,QC)                 
          ENDIF                                                         
          Y1=PIHOM-PIMLT+PIDW                               
          PT=PT+lhf*ene1*Y1    &
                     +(lhs*ene1-ene2)*(PIDEP-PISUB)*D2T
          QV=QV-(PIDEP-PISUB)*D2T                           
          QC=QC-Y1                                             
          QI=QI+Y1                                             
!* 31 * PINT  : INITIATION OF QI                                  **31**
!****** DIFFERENT FROM LIN ET AL.                                 ******
          PINT=0.0                                                   
          TAIR=(th0(i,j,k)+PT)*PPI                                     
           if (qi .le. cmin) qi=0.0                               
          IF(TAIR.LT.T0K)THEN                                        
            TAIRC=TAIR-T0K                                        
            RTAIR=1./(TAIR-C76)                                   
            Y2=EXP(C218-C580*RTAIR)                               
            QSI=C380*rp0*Y2                                   
            ESI=C610*Y2                                           
            SSI=QV/QSI-1.                                
            NCI=MIN(RN25*EXP(BETA*TAIRC),1.)                       
!     Do not initiate the maximum number of nuclei. Some ice may already
!     present, so only initiate that needed to get up to the maximum num
!     Use ami50 g as mass of 50 micron size particle (assumed to be the 
!     size of ice particles), and use 1.e-9 g as initial mass of ice    
!     particles                                                         
            dd=MAX(1.e-9*nci*rrho-qi*1.e-9/ami50 , 0.)        
            DM=MAX((QV-QSI),0.)                          
            RSUB1=C580*ASC*QSI*RTAIR*RTAIR                  
            PINT=MIN(DD,DM)                                  
            DEP=DM/(1.+RSUB1)                                  
            PINT=MIN(PINT,DEP)                               
            if (pint .le. cmin) pint=0.0                          
            PT=PT+(lhs*ene1-ene2)*PINT
            QV=QV-PINT                                         
            QI=QI+PINT                                         
          ENDIF                                                         
!--------------------------------------
!  cut saturation .... do later
!--------------------------------------
!* 10 * PSDEP : DEPOSITION OR SUBLIMATION OF QS                   **10**
!* 20 * PGSUB : SUBLIMATION OF QG                                 **20**
          PSDEP=0.0                                                  
          PSSUB=0.0                                                  
          PGSUB=0.0                                                  
          r3456=tem2
          TAIR=(th0(i,j,k)+PT)*PPI                                     
          IF((QS+QG).GT.CMIN)THEN              
            IF(QS.LT.cmin) QS=0.0                                 
            IF(QG.LT.cmin) QG=0.0                                 
            RTAIR=1./(TAIR-C76)                                   
            QSI=C380*rp0*EXP(C218-C580*RTAIR)                 
            SSI=QV/QSI-1.                                
            Y1=RN10A*RHO/(TCA*TAIR**2)+1./(DWV*QSI)   
            Y2=.78/ZS**2+RN101*tem4*SCV    &
                 *exp(alog(ZS)*(-BSH5))
            PSDEP=R10T*SSI*Y2/Y1                            
            PSSUB=PSDEP                                           
            PSDEP=MAX(PSDEP,0.)                                 
            PSSUB=MAX(-QS,MIN(PSSUB,0.))                   
            IF(IHAIL.EQ.1) THEN                                         
               Y2=.78/ZG**2+SCV*RN20B*tem1*exp(alog(ZG)*(-BGH5))
            ELSE                                                        
               Y2=.78/ZG**2+SCV*RN20B*tem4*exp(alog(ZG)*(-BGH5))
            ENDIF                                                       
            PGSUB=R20T*SSI*Y2/Y1                            
            DM=QV-QSI                                      
            RSUB1=C580*ASC*QSI*RTAIR*RTAIR                  
!     ********   DEPOSITION OR SUBLIMATION OF QS  **********************
            Y1=DM/(1.+RSUB1)                                   
            PSDEP=MIN(PSDEP,MAX(Y1,0.))                    
            Y2=MIN(Y1,0.)                                       
            PSSUB=MAX(PSSUB,Y2)                              
!     ********   SUBLIMATION OF QG   ***********************************
            DD=MAX((-Y2-QS),0.)                              
            PGSUB=MIN(DD,QG,MAX(PGSUB,0.))              
!      DLT1=CVMGP(1.,0.,QC+QI-1.E-5)                           
            IF((QC+QI-1.E-5).GE.0.)THEN                           
              DLT1=1.                                                
            ELSE                                                        
              DLT1=0.                                                
            ENDIF                                                       
            PSDEP=DLT1*PSDEP                                   
            PSSUB=(1.-DLT1)*PSSUB                              
            PGSUB=(1.-DLT1)*PGSUB                              
            PT=PT+(lhs*ene1-ene2)*(PSDEP+PSSUB-PGSUB)
            QV=QV+PGSUB-PSSUB-PSDEP                      
            QS=QS+PSDEP+PSSUB                               
            QG=QG-PGSUB                                        
          ENDIF                                                         
!* 23 * ERN : EVAPORATION OF QR (SUBSATURATION)                   **23**
          ERN=0.0                                                    
          IF(QR.GT.CMIN)THEN                                         
            TAIR=(th0(i,j,k)+PT)*PPI                                   
            RTAIR=1./(TAIR-C358)                                  
            QSW=C380*rp0*EXP(C172-C409*RTAIR)                 
            SSW=QV/QSW-1.0                               
            DM=QV-QSW                                      
            RSUB1=C409*AVC*QSW*RTAIR*RTAIR                  
            DD1=MAX(-DM/(1.+RSUB1),0.0)                      
            Y1=.78/ZR**2+RN23A*tem4*SCV    &
                 *exp(alog(ZR)*(-BWH5))
            Y2=RN23B*RHO/(TCA*TAIR**2)+1./(DWV*QSW)   
            ERN=R23T*SSW*Y1/Y2                              
            ERN=MIN(DD1,QR,MAX(ERN,0.))                 
            PT=PT-(lhv*ene1-ene2)*ERN
            QV=QV+ERN                                          
            QR=QR-ERN                                          
          ENDIF                                                         
          dum=ruh(i)*rvh(j)*rmh(i,j,k)
          bud3(k)=bud3(k)+rr(i,j,k)*ERN*dum
          IF(QC.LE.cmin)QC=0.                                     
          IF(QR.LE.cmin)QR=0.                                     
          IF(QI.LE.cmin)QI=0.                                     
          IF(QS.LE.cmin)QS=0.                                     
          IF(QG.LE.cmin)QG=0.                                     
        IF(eqtset.eq.2)THEN
          dum1=abs(qv-qv3d(i,j,k))+abs(qc-qc3d(i,j,k))+abs(qr-qr3d(i,j,k))   &
              +abs(qi-qi3d(i,j,k))+abs(qs-qs3d(i,j,k))+abs(qg-qg3d(i,j,k))
          if(dum1.gt.1.0e-7)then
            pp3d(i,j,k)=((rhod(i,j,k)*(rd+rv*qv)*(th0(i,j,k)+pt)*rp00)**rddcv)-pi0(i,j,k)
            prs(i,j,k)=p00*((pi0(i,j,k)+pp3d(i,j,k))**cpdrd)
          endif
        ELSE
          rhod(i,j,k)=prs(i,j,k)/(rd*(th0(i,j,k)+pt)*(pi0(i,j,k)+pp3d(i,j,k))*(1.0+qv*reps))
        ENDIF

          TH3D(I,J,K)=PT
          QV3D(I,J,K)=QV
          QC3D(I,J,K)=QC
          QR3D(I,J,K)=QR
          QI3D(I,J,K)=QI
          QS3D(I,J,K)=QS
          QG3D(I,J,K)=QG
          vtr(i,j,k)=vfr*0.01
          vts(i,j,k)=vfs*0.01
          vtg(i,j,k)=vfg*0.01
          vti(i,j,k)=0.200

      enddo
      enddo
      enddo

      tem6=dx*dy*dz

      do k=1,nk
        tauto=tauto+bud1(k)*tem6
      enddo

      do k=1,nk
        taccr=taccr+bud2(k)*tem6
      enddo

      do k=1,nk
        tevar=tevar+bud3(k)*tem6
      enddo

      if(timestats.ge.1) time_microphy=time_microphy+mytime()

!C**********************************************************************
!                                                                       
      END subroutine goddard                                                              


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!***********************************************************************
      SUBROUTINE CONSAT
      use input
      use cm1libs , only : gamma
      implicit none
!***********************************************************************
!     (LIN) SPECIFY SOME CONSTANTS IN SATICE ROUTINE                   *
!     LIN ET.AL.  J. CLIM. APPL. METEOR.  22, 1065-1092                *
!     MODIFIED AND CODED BY TAO AND SIMPSON (JAS, 1989; TAO, 1993)     *
!     RECENT CHANGES (6.25.01) BY SCOTT BRAUN INCLUDE
!         * OPTION TO CHOOSE HAIL OR GRAUPEL AS THIRD ICE CATEGORY
!         * OPTION OF SATURATION ADJUSTMENT SCHEMES
!         * INCLUSION OF GAMMA FUNCTION FOR EASIER VARIATION OF FALL
!           SPEED COEFFICIENTS
!***********************************************************************

      real A1(31),A2(31)
      DATA A1/.7939E-7,.7841E-6,.3369E-5,.4336E-5,.5285E-5,.3728E-5,       &
         .1852E-5,.2991E-6,.4248E-6,.7434E-6,.1812E-5,.4394E-5,.9145E-5,   &
         .1725E-4,.3348E-4,.1725E-4,.9175E-5,.4412E-5,.2252E-5,.9115E-6,   &
         .4876E-6,.3473E-6,.4758E-6,.6306E-6,.8573E-6,.7868E-6,.7192E-6,   &
         .6513E-6,.5956E-6,.5333E-6,.4834E-6/
      DATA A2/.4006,.4831,.5320,.5307,.5319,.5249,.4888,.3894,.4047,    &
         .4318,.4771,.5183,.5463,.5651,.5813,.5655,.5478,.5203,.4906,   &
         .4447,.4126,.3960,.4149,.4320,.4506,.4483,.4460,.4433,.4413,   &
         .4382,.4361/

      real cp,cpi,cpi2,grvt,tca,dwv,dva,amw,ars,scv,rw,cw,ci,cd1,cd2,     &
           ga3b,ga4b,ga6b,ga5bh,ga3g,ga4g,ga5gh,ga3d,ga4d,ga5dh,esw,      &
           eri,ami,esr,eiw,ui50,ri50,cmn,y1,rn13,egw,egi,egi2,egr,apri,   &
           bpri,erw,cn0

      integer i,k

!
!
!     USE GRAUPEL OR HAIL (0=GRAUPEL, 1=HAIL)
!
!!!!!! 
!!!!!!      IHAIL=1       ! Now specified in namelist.input
!!!!!! 
!
!*****   IWATER=0 USES A SLIGHT VARIANT OF THE ORIGINAL TAO ET AL METHOD
!*****   THE ONLY DIFFERENCE IS THAT THE WEIGHTING BETWEEN THE LIQUID
!*****   AND ICE SATURATION VALUES IS DONE BY TEMPERATURE RATHER THAN MA
!*****   IWATER=1 USES A SEQUENTIAL METHOD IN WHICH THE ADJUSTMENT IS FI
!*****   DONE FOR LIQUID WATER FOR TEMPERATURES WARMER THAN 253K, THEN F
!*****   ICE ONLY WHERE THE TEMPERATURE IS COLDER THAN 258K. THE MAIN EF
!*****   THIS CHANGE IS TO REDUCE THE AMOUNT OF SUPERCOOLED WATER AT VER
!*****   TEMPERATURES
!!!!!! 
!!!!!!      IWATER=1
!!!!!! 
!*****************************************************************
 
!
      CP=1.0057E7
      CPI=4.*ATAN(1.)
      CPI2=CPI*CPI
      GRVT=981.
      TCA=2.43E3
      DWV=.226
      DVA=1.718E-4
      AMW=18.016
      ARS=8.314E7
      SCV=2.2904487
      T0K=273.15
      T00K=233.15
      ALV=2.501000E10
      ALS=2.834000E10
      ALF=ALS-ALV
      AVC=ALV/CP
      AFC=ALF/CP
      ASC=ALS/CP
      RW=4.615E6
      CW=4.190E7
      CI=2.106000E7
      C76=7.66
      C358=29.65
      C172=17.67
      C409=17.67*(273.15-29.65)
      C218=21.8745584
      C580=21.8745584*(273.15-7.66)
      C380=10.0*611.2*287.04/461.5
      C610=6.1078E3
      C149=1.496286E-5
      C879=8.794142
      C141=1.4144354E7
 
!***   DEFINE THE DENSITY AND SIZE DISTRIBUTION OF PRECIPITATION
!***   DEFINE THE COEFFICIENTS USED IN TERMINAL VELOCITY
!**********   HAIL OR GRAUPEL PARAMETERS   **********
      if (ihail .eq. 1) then
        ROQG=.9
        TNG=.0002
        CD1=6.E-1
        CD2=4.*GRVT/(3.*CD1)
        AGG=SQRT(CD2*ROQG)
        BGG=.5
      else
        ROQG=.4
        TNG=.04
        AGG=351.2
        BGG=.37
      endif
!**********         SNOW PARAMETERS        **********
! Note ... see Potter, 1991, JAM, p. 1040 for more info about these changes
      ROQS=.1
      ROQS_POTTER=1.0
      TNSS=1.
!!!      ASS=152.93
!!!      BSS=.25
      ASS=179.2
      BSS=.42
!**********         RAIN PARAMETERS        **********
      ROQR=1.
      TNW=.08
      AWW=2115.
      BWW=.8
!*****************************************************************
!
!
      BGH=.5*BGG
      BSH=.5*BSS
      BWH=.5*BWW
      BGQ=.25*BGG
      BSQ=.25*BSS
      BWQ=.25*BWW
      GA3B=gamma(3.+BWW)
      GA4B=gamma(4.+BWW)
      GA6B=gamma(6.+BWW)
      GA5BH=gamma((5.+BWW)/2.)
      GA3G=gamma(3.+BGG)
      GA4G=gamma(4.+BGG)
      GA5GH=gamma((5.+BGG)/2.)
      GA3D=gamma(3.+BSS)
      GA4D=gamma(4.+BSS)
      GA5DH=gamma((5.+BSS)/2.)
      ZRC=(CPI*ROQR*TNW)**0.25
      ZSC=(CPI*ROQS_POTTER*TNSS)**0.25
      ZGC=(CPI*ROQG*TNG)**0.25
      VRC=AWW*GA4B/(6.*ZRC**BWW)
      VSC=ASS*GA4D/(6.*ZSC**BSS)
      VGC=AGG*GA4G/(6.*ZGC**BGG)
!-------------------------------
!  fudge check
      if(ihail.eq.1 .and. bgg.ne.0.5)then
        print *,'  BGG must be 0.5 for hail (its fudged into code)!'
        print *,'        (sorry)'
        call stopcm1
      endif
!!!      if(bss.ne.0.25)then
!!!        print *,'  BSS must be 0.25 (its fudged into code)!'
!!!        print *,'        (sorry)'
!!!        call stopcm1
!!!      endif
      if(bww.ne.0.8)then
        print *,'  BWW must be 0.8 (its fudged into code)!'
        print *,'        (sorry)'
        call stopcm1
      endif
!-------------------------------
!     ****************************
      RN1=1.E-3
      RN2=1.E-3
      BND1=5.E-4
      BND2=1.25E-3
      RN3=.25*CPI*TNSS*ASS*GA3D
      ESW=1.
      RN4=.25*CPI*ESW*TNSS*ASS*GA3D
      ERI=1.
      RN5=.25*CPI*ERI*TNW*AWW*GA3B
      AMI=1./(24.*4.19E-10)
      RN6=CPI2*ERI*TNW*AWW*ROQR*GA6B*AMI
      ESR=1.
      RN7=CPI2*ESR*TNW*TNSS*ROQS
      RN8=CPI2*ESR*TNW*TNSS*ROQR
      RN9=CPI2*TNSS*TNG*ROQS
      RN10=2.*CPI*TNSS
      RN101=.31*GA5DH*SQRT(ASS)
      RN10A=ALS*ALS/RW
      RN11=2.*CPI*TNSS/ALF
      RN11A=CW/ALF
      AMI40=2.41e-8
      AMI50=3.76e-8
      AMI100=1.51e-7
      EIW=1.
      UI50=20.
      RI50=10.e-3
      CMN=1.05E-15
      RN12=CPI*EIW*UI50*RI50**2
      DO 10 K=1,31
        Y1=1.-A2(K)
        RN13=A1(K)*Y1/(AMI100**Y1-AMI40**Y1)
        RN12A(K)=RN13/AMI100
        RN12B(K)=A1(K)*AMI100**A2(K)
        RN25A(K)=A1(K)*CMN**A2(K)
   10 CONTINUE
      EGW=1.
      EGI=.1
      EGI2=1.
      RN14=.25*CPI*EGW*TNG*GA3G*agg
      RN15=.25*CPI*EGI*TNG*GA3G*agg
      RN15A=.25*CPI*EGI2*TNG*GA3G*agg
      EGR=1.
      RN16=CPI2*EGR*TNG*TNW*ROQR
      RN17=2.*CPI*TNG
      RN17A=.31*GA5GH*sqrt(agg)
      RN17B=CW-CI
      RN17C=CW
      APRI=.66
      BPRI=1.E-4
      RN18=20.*CPI2*BPRI*TNW*ROQR
      RN18A=APRI
      RN19=2.*CPI*TNG/ALF
      RN19A=.31*GA5GH*sqrt(agg)
      RN19B=CW/ALF
      RN20=2.*CPI*TNG
      RN20A=ALS*ALS/RW
      RN20B=.31*GA5GH*sqrt(agg)
      BND3=2.0E-3
      RN21=1.E3*1.569E-12/0.15
      ERW=1.
      RN22=.25*CPI*ERW*AWW*TNW*GA3B
      RN23=2.*CPI*TNW
      RN23A=.31*GA5BH*SQRT(AWW)
      RN23B=ALV*ALV/RW
      CN0=1.E-6
      RN25=CN0
      RN30A=ALV*ALS*AMW/(TCA*ARS)
      RN30B=ALV/TCA
      RN30C=ARS/(DWV*AMW)
      RN31=1.E-17
      BETA=-.46
      RN32=4.*51.545E-4
!-----------------------------------------------
      RT0=1./(T0K-T00K)
      BW3=BWW+3.
      BS3=BSS+3.
      BG3=BGG+3.
      BWH5=2.5+BWH
      BSH5=2.5+BSH
      BGH5=2.5+BGH
      BW6=BWW+6.
      BS6=BSS+6.
      BETAH=.5*BETA
!-----------------------------------------------
!  lookup tables
!
!  for pwr81, temperature range is from -100 to +50 Celsius in increments
!  of 1 Celsius
      do i=0,151
        pwr81(i)=(173.15+i)**0.81
      enddo
!  rho*q from 0 to 50e-6 in increments of 1e-9
      do i=0,50001
        pwr2(i)=(i*1.0e-9)**bwq
      enddo
!-----------------------------------------------
      END subroutine consat


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!***********************************************************************
      SUBROUTINE CONSAT2(d2t)
      implicit none
      real d2t
      R10T=RN10*D2T
      R11AT=RN11A*D2T
      R19BT=RN19B*D2T
      R20T=-RN20*D2T
      R23T=-RN23*D2T
      END subroutine consat2

  END MODULE goddard_module
