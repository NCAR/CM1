!WRF:MODEL_LAYER:PHYSICS
!
module module_sf_slab

   !---SPECIFY CONSTANTS AND LAYERS FOR SOIL MODEL
   !---SOIL DIFFUSION CONSTANT SET (M^2/S)

   real, parameter :: DIFSL=5.e-7

   !---FACTOR TO MAKE SOIL STEP MORE CONSERVATIVE

   real , parameter :: SOILFAC=1.25

contains

!----------------------------------------------------------------
   subroutine SLAB(T3D,QV3D,P3D,FLHC,FLQC,                      &
                   PSFC,XLAND,TMN,HFX,QFX,LH,TSK,QSFC,CHKLOWQ,  &
                   GSW,GLW,CAPG,THC,SNOWC,EMISS,MAVAIL,         &
                   DELTSM,ROVCP,XLV,lv1,lv2,DTMIN,IFSNOW,       &
                   SVP1,SVP2,SVP3,SVPT0,EP2,                    &
                   KARMAN,EOMEG,STBOLT,                         &
                   TSLB,ZS,DZS,num_soil_layers,radiation,       &
                   P1000mb,                                     &
                   ids,ide, jds,jde, kds,kde,                   &
                   ims,ime, jms,jme, kms,kme,                   &
                   its,ite, jts,jte, kts,kte                    )
!----------------------------------------------------------------
    implicit none
!----------------------------------------------------------------
!                                                                        
!     SUBROUTINE SLAB CALCULATES THE GROUND TEMPERATURE TENDENCY 
!     ACCORDING TO THE RESIDUAL OF THE SURFACE ENERGY BUDGET           
!     (BLACKADAR, 1978B).                                              
!                                                                      
!     CHANGES:                                                         
!          FOR SOIL SUB-TIMESTEPS UPDATE SURFACE HFX AND QFX AS TG     
!          CHANGES TO PREVENT POSSIBLE INSTABILITY FOR LONG MODEL      
!          STEPS (DT > ~200 SEC).                                      
!                                                                      
!          PUT SNOW COVER CHECK ON SOIL SUB-TIMESTEPS                  
!                                                                      
!          MAKE UPPER LIMIT ON SOIL SUB-STEP LENGTH MORE CONSERVATIVE  
!                                                                      
!----------------------------------------------------------------          
!-- T3D         temperature (K)
!-- QV3D        3D water vapor mixing ratio (Kg/Kg)
!-- P3D         3D pressure (Pa)
!-- FLHC        exchange coefficient for heat (m/s)
!-- FLQC        exchange coefficient for moisture (m/s)
!-- PSFC        surface pressure (Pa)
!-- XLAND       land mask (1 for land, 2 for water)
!-- TMN         soil temperature at lower boundary (K)
!-- HFX         upward heat flux at the surface (W/m^2)
!-- QFX         upward moisture flux at the surface (kg/m^2/s)
!-- LH          latent heat flux at the surface (W/m^2)
!-- TSK         surface temperature (K)
!-- GSW         downward short wave flux at ground surface (W/m^2)      
!-- GLW         downward long wave flux at ground surface (W/m^2)
!-- CAPG        heat capacity for soil (J/K/m^3)
!-- THC         thermal inertia (Cal/cm/K/s^0.5)
!-- SNOWC       flag indicating snow coverage (1 for snow cover)
!-- EMISS       surface emissivity (between 0 and 1)
!-- DELTSM      time step (second)
!-- ROVCP       R/CP
!-- XLV         latent heat of melting (J/kg)
!-- DTMIN       time step (minute)
!-- IFSNOW      ifsnow=1 for snow-cover effects
!-- SVP1        constant for saturation vapor pressure (kPa)
!-- SVP2        constant for saturation vapor pressure (dimensionless)
!-- SVP3        constant for saturation vapor pressure (K)
!-- SVPT0       constant for saturation vapor pressure (K)
!-- EP1         constant for virtual temperature (R_v/R_d - 1) (dimensionless)
!-- EP2         constant for specific humidity calculation 
!               (R_d/R_v) (dimensionless)
!-- KARMAN      Von Karman constant
!-- EOMEG       angular velocity of earth's rotation (rad/s)
!-- STBOLT      Stefan-Boltzmann constant (W/m^2/K^4)
!-- TSLB        soil temperature in 5-layer model
!-- ZS          depths of centers of soil layers
!-- DZS         thicknesses of soil layers
!-- num_soil_layers   the number of soil layers
!-- ids         start index for i in domain
!-- ide         end index for i in domain
!-- jds         start index for j in domain
!-- jde         end index for j in domain
!-- kds         start index for k in domain
!-- kde         end index for k in domain
!-- ims         start index for i in memory
!-- ime         end index for i in memory
!-- jms         start index for j in memory
!-- jme         end index for j in memory
!-- kms         start index for k in memory
!-- kme         end index for k in memory
!-- its         start index for i in tile
!-- ite         end index for i in tile
!-- jts         start index for j in tile
!-- jte         end index for j in tile
!-- kts         start index for k in tile
!-- kte         end index for k in tile
!----------------------------------------------------------------
   integer,  intent(in   )   ::     ids,ide, jds,jde, kds,kde,  &
                                    ims,ime, jms,jme, kms,kme,  &
                                    its,ite, jts,jte, kts,kte

   integer, intent(in)       ::     num_soil_layers
   logical, intent(in)       ::     radiation

   integer,  intent(in   )   ::     IFSNOW

!
   real,     intent(in   )   ::     DTMIN,XLV,lv1,lv2,ROVCP,DELTSM

   real,     intent(in )     ::     SVP1,SVP2,SVP3,SVPT0
   real,     intent(in )     ::     EP2,KARMAN,EOMEG,STBOLT
   real,     intent(in )     ::     P1000mb

   real,     dimension( ims:ime , jms:jme , 1:num_soil_layers ), &
             intent(inout)   :: TSLB

   real,     dimension(1:num_soil_layers), intent(in)::ZS,DZS

   real,    dimension( ims:ime, jms:jme, kms:kme )            , &
            intent(in   )    ::                           QV3D, &
                                                           P3D, &
                                                           T3D
!
   real,    dimension( ims:ime, jms:jme )                     , &
            intent(in   )    ::                          SNOWC, &
                                                         XLAND, &
                                                         EMISS, &
                                                        MAVAIL, &
                                                           TMN, &
                                                           GSW, &
                                                           GLW, &
                                                           THC

!CHKLOWQ is declared as memory size
!
   real,    dimension( ims:ime, jms:jme )                     , &
            intent(inout)    ::                            HFX, &
                                                           QFX, &
                                                            LH, &
                                                          CAPG, &
                                                           TSK, &
                                                          QSFC, &
                                                       CHKLOWQ

   real,     dimension( ims:ime, jms:jme )                    , &
             intent(in   )               ::               PSFC
!
   real,    dimension( ims:ime, jms:jme ), intent(inout) ::     &
                                                          FLHC, &
                                                          FLQC

! LOCAL VARS

   real,     dimension( its:ite ) ::                      QV1D, &
                                                           P1D, &
                                                           T1D
   real,     dimension( its:ite , 1:num_soil_layers ) :: TSLB2D
   integer ::  I,J,n

!$omp parallel do default(shared)   &
!$omp private(i,j,n,t1d,qv1d,p1d,tslb2d)
   do J=jts,jte

      do i=its,ite
         T1D(i) =T3D(i,j,1)
         QV1D(i)=QV3D(i,j,1)
         P1D(i) =P3D(i,j,1)
      enddo

      do n=1,num_soil_layers
      do i=its,ite
         TSLB2D(i,n) = TSLB(i,j,n)
      enddo
      enddo

! the indices to the PSFC argument in the following call look
! wrong; however, it is correct to call with its (and not ims)
! because of the way PSFC is defined in SLAB1D. Whether *that*
! is a good idea or not, this commenter cannot comment. JM

      call SLAB1D(J,T1D,QV1D,P1D,FLHC(ims,j),FLQC(ims,j),       &
           PSFC(ims,j),XLAND(ims,j),TMN(ims,j),HFX(ims,j),      &
           QFX(ims,j),TSK(ims,j),QSFC(ims,j),CHKLOWQ(ims,j),    &
           LH(ims,j),GSW(ims,j),GLW(ims,j),                     &
           CAPG(ims,j),THC(ims,j),SNOWC(ims,j),EMISS(ims,j),    &
           MAVAIL(ims,j),DELTSM,ROVCP,XLV,lv1,lv2,DTMIN,IFSNOW, &
           SVP1,SVP2,SVP3,SVPT0,EP2,KARMAN,EOMEG,STBOLT,        &
           TSLB2D       ,ZS,DZS,num_soil_layers,radiation,      &
           P1000mb,                                             &
           ids,ide, jds,jde, kds,kde,                           &
           ims,ime, jms,jme, kms,kme,                           &
           its,ite, jts,jte, kts,kte                            )

      do n=1,num_soil_layers
      do i=its,ite
         TSLB(i,j,n) = TSLB2D(i,n)
      enddo
      enddo

   enddo

   end subroutine SLAB

!----------------------------------------------------------------
   subroutine SLAB1D(J,T1D,QV1D,P1D,FLHC,FLQC,                  &
                   PSFCPA,XLAND,TMN,HFX,QFX,TSK,QSFC,CHKLOWQ,   &
                   LH,GSW,GLW,CAPG,THC,SNOWC,EMISS,MAVAIL,      &
                   DELTSM,ROVCP,XLV,lv1,lv2,DTMIN,IFSNOW,       &
                   SVP1,SVP2,SVP3,SVPT0,EP2,                    &
                   KARMAN,EOMEG,STBOLT,                         &
                   TSLB2D,ZS,DZS,num_soil_layers,radiation,     &
                   P1000mb,                                     &
                   ids,ide, jds,jde, kds,kde,                   &
                   ims,ime, jms,jme, kms,kme,                   &
                   its,ite, jts,jte, kts,kte                    )
!----------------------------------------------------------------
    implicit none
!----------------------------------------------------------------
!                                                                        
!     SUBROUTINE SLAB CALCULATES THE GROUND TEMPERATURE TENDENCY 
!     ACCORDING TO THE RESIDUAL OF THE SURFACE ENERGY BUDGET           
!     (BLACKADAR, 1978B).                                              
!                                                                      
!     CHANGES:                                                         
!          FOR SOIL SUB-TIMESTEPS UPDATE SURFACE HFX AND QFX AS TG     
!          CHANGES TO PREVENT POSSIBLE INSTABILITY FOR LONG MODEL      
!          STEPS (DT > ~200 SEC).                                      
!                                                                      
!          PUT SNOW COVER CHECK ON SOIL SUB-TIMESTEPS                  
!                                                                      
!          MAKE UPPER LIMIT ON SOIL SUB-STEP LENGTH MORE CONSERVATIVE  
!                                                                      
!----------------------------------------------------------------          

   integer,  intent(in   )   ::     ids,ide, jds,jde, kds,kde,  &
                                    ims,ime, jms,jme, kms,kme,  &
                                    its,ite, jts,jte, kts,kte,J 

   integer , intent(in)      ::     num_soil_layers
   logical,  intent(in   )   ::     radiation

   integer,  intent(in   )   ::     IFSNOW
!
   real,     intent(in   )   ::     DTMIN,XLV,lv1,lv2,ROVCP,DELTSM

   real,     intent(in )     ::     SVP1,SVP2,SVP3,SVPT0
   real,     intent(in )     ::     EP2,KARMAN,EOMEG,STBOLT
   real,     intent(in )     ::     P1000mb

   real,     dimension( its:ite , 1:num_soil_layers ),          &
             intent(inout)   :: TSLB2D

   real,     dimension(1:num_soil_layers), intent(in)::ZS,DZS

!
   real,    dimension( ims:ime )                              , &
            intent(inout)    ::                            HFX, &
                                                           QFX, &
                                                            LH, &
                                                          CAPG, &
                                                           TSK, &
                                                          QSFC, &
                                                       CHKLOWQ
!
   real,    dimension( ims:ime )                              , &
            intent(in   )    ::                          SNOWC, &
                                                         XLAND, &
                                                         EMISS, &
                                                        MAVAIL, &
                                                           TMN, &
                                                           GSW, &
                                                           GLW, &
                                                           THC
!
   real,    dimension( its:ite )                              , &
            intent(in   )    ::                           QV1D, &
                                                           P1D, &
                                                           T1D
!
   real,     dimension( ims:ime )                             , &
             intent(in   )               ::             PSFCPA

!
   real,    dimension( ims:ime ), intent(inout) ::              &
                                                          FLHC, &
                                                          FLQC
! LOCAL VARS

   real,    dimension( its:ite )          ::              PSFC

   real,    dimension( its:ite )          ::                    &
                                                           THX, &
                                                            QX, &
                                                          SCR3 

   real,    dimension( its:ite )          ::            DTHGDT, &
                                                           TG0, &
                                                         THTMN, &
                                                          XLD1, &
                                                         TSCVN, &
                                                          OLTG, &
                                                        UPFLUX, &
                                                            HM, &
                                                          RNET, &
                                                         XINET, &
                                                            QS, &
                                                         DTSDT
!
   real, dimension( its:ite, num_soil_layers )        :: FLUX
!
   integer :: I,K,NSOIL,ITSOIL,L,NK,RADSWTCH
   real    :: PS,PS1,XLDCOL,TSKX,RNSOIL,RHOG1,RHOG2,RHOG3,LAMDAG
   real    :: THG,ESG,QSG,HFXT,QFXT,CS,CSW,LAMG(4),THCON,PL
 
!----------------------------------------------------------------------          
!-----DETERMINE IF ANY POINTS IN COLUMN ARE LAND (RATHER THAN OCEAN)             
!       POINTS.  IF NOT, SKIP DOWN TO THE PRINT STATEMENTS SINCE OCEAN           
!       SURFACE TEMPERATURES ARE NOT ALLOWED TO CHANGE.                          
!                                                                                
! from sfcrad   
!----------------------------------------------------------------------
   data CSW/4.183E6/
   data LAMG/1.407E-8, -1.455E-5, 6.290E-3, 0.16857/

   do i=its,ite
! in cmb
      PSFC(I)=PSFCPA(I)/1000.
   enddo


      do I=its,ite
! PL cmb
         PL=P1D(I)/1000.
         SCR3(I)=T1D(I)
!         THCON=(100./PL)**ROVCP
         THCON=(P1000mb*0.001/PL)**ROVCP
         THX(I)=SCR3(I)*THCON
         QX(I)=0.
      enddo

!     IF(IDRY.EQ.1) GOTO 81
      do I=its,ite
         QX(I)=QV1D(I)
      enddo
   81 continue

!
!-----THE SLAB THERMAL CAPACITY CAPG(I) ARE DEPENDENT ON:
!     THC(I) - SOIL THERMAL INERTIAL, ONLY.
!
      do I=its,ite
         CAPG(I)=3.298E6*THC(I)
         if(num_soil_layers .gt. 1)then

! CAPG REPRESENTS SOIL HEAT CAPACITY (J/K/M^3) WHEN DIFSL=5.E-7 (M^2/S)
! TO GIVE A CORRECT THERMAL INERTIA (=CAPG*DIFSL^0.5)

            CAPG(I)=5.9114E7*THC(I)
         endif
      enddo
!        
      XLDCOL=2.0                                                                 
      do 10 I=its,ite
        XLDCOL=AMIN1(XLDCOL,XLAND(I))                                          
   10 continue                                                                   
!                                                                                
      if(XLDCOL.GT.1.5)GOTO 90                                                   
!                                                                                
!                                                                                
!-----CONVERT SLAB TEMPERATURE TO POTENTIAL TEMPERATURE AND                      
!     SET XLD1(I) = 0. FOR OCEAN POINTS:                                         
!                                                                                
!                                                                                
      do 20 I=its,ite
        if((XLAND(I)-1.5).GE.0)then                                            
          XLD1(I)=0.                                                             
        else                                                                     
          XLD1(I)=1.                                                             
        endif                                                                    
   20 continue                                                                   
!                                                                                
!-----CONVERT 'TSK(THETAG)' TO 'TG' FOR 'IUP' CALCULATION ....                   
!       IF WE ARE USING THE BLACKADAR MULTI-LEVEL (HIGH-RESOLUTION)              
!       PBL MODEL                                                                
!                                                                                
      do 50 I=its,ite
        if(XLD1(I).LT.0.5)GOTO 50                                                

! PS cmb
        PS=PSFC(I)

! TSK is Temperature at gound sfc
!       TG0(I)=TSK(I)*(PS*0.01)**ROVCP                                         
        TG0(I)=TSK(I)
   50 continue                                                                   
!                                                                                
!-----COMPUTE THE SURFACE ENERGY BUDGET:                                         
!                                                                                
!     IF(ISOIL.EQ.1)NSOIL=1                                                      
      if(num_soil_layers .gt. 1)NSOIL=1                                                      


      if (radiation) then
        RADSWTCH=1
      else
        RADSWTCH=0
      endif

      do 70 I=its,ite
        if(XLD1(I).LT.0.5)GOTO 70
!        OLTG(I)=TSK(I)*(100./PSFC(I))**ROVCP
        OLTG(I)=TSK(I)*(P1000mb*0.001/PSFC(I))**ROVCP
        UPFLUX(I)=RADSWTCH*STBOLT*TG0(I)**4                            
        XINET(I)=EMISS(I)*(GLW(I)-UPFLUX(I))    
        RNET(I)=GSW(I)+XINET(I)                                                
        HM(I)=1.18*EOMEG*(TG0(I)-TMN(I))                                       
!       MOISTURE FLUX CALCULATED HERE (OVERWRITES SFC LAYER VALUE FOR LAND)
                ESG=SVP1*EXP(SVP2*(TG0(I)-SVPT0)/(TG0(I)-SVP3))
                QSG=EP2*ESG/(PSFC(I)-ESG)
                THG=TSK(I)*(100./PSFC(I))**ROVCP
                HFX(I)=FLHC(I)*(THG-THX(I))
                QFX(I)=FLQC(I)*(QSG-QX(I))
                LH(I)=QFX(I)*( lv1-lv2*tsk(I) )
        QS(I)=HFX(I)+QFX(I)*( lv1-lv2*tsk(I) )
!       IF(ISOIL.EQ.0)THEN                                                       
        if(num_soil_layers .EQ. 1)then                                                       
          DTHGDT(I)=(RNET(I)-QS(I))/CAPG(I)-HM(I)                              
        else
          DTHGDT(I)=0.                                                           
        endif                                                                    
   70 continue                                                                   
!     IF(ISOIL.EQ.1)THEN                                                         
      if(num_soil_layers .gt. 1)then                                                         
        NSOIL=1+IFIX(SOILFAC*4*DIFSL/DZS(1)*DELTSM/DZS(1))   
        RNSOIL=1./FLOAT(NSOIL)                                                   
!                                                                                
!     SOIL SUB-TIMESTEP                                                          
!                                                                                
        do ITSOIL=1,NSOIL                                                        
          do I=its,ite
             do L=1,num_soil_layers-1
              if(XLD1(I).LT.0.5)GOTO 75                                          
              if(L.EQ.1.AND.ITSOIL.GT.1)then                                     
!                PS1=(PSFC(I)*0.01)**ROVCP    
                PS1=(PSFCPA(I)/P1000mb)**ROVCP    

! for rk scheme A and B are the same
                PS=PSFC(I)
                THG=TSLB2D(I,1)/PS1                                              
                ESG=SVP1*EXP(SVP2*(TSLB2D(I,1)-SVPT0)/(TSLB2D(I,1) & 
                    -SVP3))                                                      
                QSG=EP2*ESG/(PS-ESG)                                             
!     UPDATE FLUXES FOR NEW GROUND TEMPERATURE                                   
                HFXT=FLHC(I)*(THG-THX(I))                                     
                QFXT=FLQC(I)*(QSG-QX(I))
                QS(I)=HFXT+QFXT*( lv1-lv2*tsk(I) )
!     SUM HFX AND QFX OVER SOIL TIMESTEPS                                        
                HFX(I)=HFX(I)+HFXT                                           
                QFX(I)=QFX(I)+QFXT                                           
              endif                                                              
              FLUX(I,1)=RNET(I)-QS(I)                                            
              FLUX(I,L+1)=-DIFSL*CAPG(I)*(TSLB2D(I,L+1)-TSLB2D(I,L))/( & 
                          ZS(L+1)-ZS(L))                                         
              DTSDT(I)=-(FLUX(I,L+1)-FLUX(I,L))/(DZS(L)*CAPG(I))               
              TSLB2D(I,L)=TSLB2D(I,L)+DTSDT(I)*DELTSM*RNSOIL                     
              if(IFSNOW.EQ.1.AND.L.EQ.1)then                              
                if((SNOWC(I).GT.0..AND.TSLB2D(I,1).GT.273.16))then             
                  TSLB2D(I,1)=273.16                                             
                endif                                                            
              endif                                                              
              if(L.EQ.1)DTHGDT(I)=DTHGDT(I)+RNSOIL*DTSDT(I)                      
              if(ITSOIL.EQ.NSOIL.AND.L.EQ.1)then                                 
!     AVERAGE HFX AND QFX OVER SOIL TIMESTEPS FOR OUTPUT TO PBL                  
                HFX(I)=HFX(I)*RNSOIL                                         
                QFX(I)=QFX(I)*RNSOIL                                         
                LH(I)=QFX(I)*( lv1-lv2*tsk(I) )
              endif                                                              
   75         continue                                                           
            enddo                                                                
          enddo                                                                  
        enddo                                                                    
      endif                                                                      
!                                                                                
      do 80 I=its,ite
        if(XLD1(I).LT.0.5) GOTO 80                                                
        TSKX=TG0(I)+DELTSM*DTHGDT(I)                                             

! TSK is temperature
!       TSK(I)=TSKX*(100./PS1)**ROVCP                                          
        TSK(I)=TSKX
   80 continue                                                                   

!                                                                                
!-----MODIFY THE THE GROUND TEMPERATURE IF THE SNOW COVER EFFECTS ARE            
!     CONSIDERED: LIMIT THE GROUND TEMPERATURE UNDER 0 C.                        
!                                                                                
      if(IFSNOW.EQ.0)GOTO 90                                              
      do 85 I=its,ite
        if(XLD1(I).LT.0.5)GOTO 85                                                
!       PS1=(PSFC(I)*0.01)**ROVCP             
!       TSCVN(I)=TSK(I)*PS1                                            
        TSCVN(I)=TSK(I)
        if((SNOWC(I).GT.0..AND.TSCVN(I).GT.273.16))then                        
          TSCVN(I)=273.16                                                        
        else                                                                     
          TSCVN(I)=TSCVN(I)                                                      
        endif                                                                    
!       TSK(I)=TSCVN(I)/PS1                                                    
        TSK(I)=TSCVN(I)
   85 continue                                                                   
!                                                                                
   90 continue                                                                   
      do I=its,ite
! QSFC and CHKLOWQ needed by Eta PBL
! WA added check for flqc = 0 to accomodate TEMF (and others?)
        if ( FLQC(I) .ne. 0.) then
           QSFC(I)=QX(I)+QFX(I)/FLQC(I)
        else
           QSFC(I) = QX(I)
        end if
        CHKLOWQ(I)=MAVAIL(I)
      enddo
!                                                                                
  140 continue                                                                   

   end subroutine SLAB1D

!================================================================
   subroutine slabinit(TSK,TMN,                                 &
                       TSLB,ZS,DZS,num_soil_layers,             &
                       allowed_to_read, start_of_simulation,    &
                       ids,ide, jds,jde, kds,kde,               &
                       ims,ime, jms,jme, kms,kme,               &
                       its,ite, jts,jte, kts,kte                )
!----------------------------------------------------------------
   implicit none
!----------------------------------------------------------------
   logical , intent(in)      ::      allowed_to_read
   logical , intent(in)      ::      start_of_simulation
   integer, intent(in   )    ::      ids,ide, jds,jde, kds,kde, &
                                     ims,ime, jms,jme, kms,kme, &
                                     its,ite, jts,jte, kts,kte

   integer, intent(in   )    ::      num_soil_layers
!   
   real,     dimension( ims:ime , 1:num_soil_layers , jms:jme ), intent(inout) :: TSLB

   real,     dimension(1:num_soil_layers), intent(in)  ::  ZS,DZS

   real,    dimension( ims:ime, jms:jme )                     , &
            intent(in)    ::                               TSK, &
                                                           TMN

!  LOCAR VAR

   integer                   ::      L,J,I,itf,jtf
   character*1024 message

!----------------------------------------------------------------
 
   itf=min0(ite,ide-1)
   jtf=min0(jte,jde-1)

   end subroutine slabinit
!-------------------------------------------------------------------          

end module module_sf_slab
