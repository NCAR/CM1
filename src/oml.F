!WRF:MODEL_LAYER:PHYSICS
!
module module_sf_oml

contains

!----------------------------------------------------------------
   subroutine OCEANML(tml,t0ml,hml,h0ml,huml,hvml,ust,u_phy,v_phy, &
                      tmoml,f,g,oml_gamma,                         &
                      OML_RELAXATION_TIME,                         &
                      XLAND,HFX,LH,TSK,GSW,GLW,EMISS,              &
                      DELTSM,STBOLT,                               &
                      ids,ide, jds,jde, kds,kde,                   &
                      ims,ime, jms,jme, kms,kme,                   &
                      its,ite, jts,jte, kts,kte                    )

!----------------------------------------------------------------
   implicit none
!----------------------------------------------------------------
!
!  SUBROUTINE OCEANML CALCULATES THE SEA SURFACE TEMPERATURE (TSK)
!  FROM A SIMPLE OCEAN MIXED LAYER MODEL BASED ON
!  (Pollard, Rhines and Thompson (1973).
!
!-- TML         ocean mixed layer temperature (K)
!-- T0ML        ocean mixed layer temperature (K) at initial time
!-- TMOML       top 200 m ocean mean temperature (K) at initial time
!-- HML         ocean mixed layer depth (m)
!-- H0ML        ocean mixed layer depth (m) at initial time
!-- HUML        ocean mixed layer u component of wind
!-- HVML        ocean mixed layer v component of wind
!-- OML_GAMMA   deep water lapse rate (K m-1)
!-- UAIR,VAIR   lowest model level wind component
!-- UST         frictional velocity
!-- HFX         upward heat flux at the surface (W/m^2)
!-- LH          latent heat flux at the surface (W/m^2)
!-- TSK         surface temperature (K)
!-- GSW         downward short wave flux at ground surface (W/m^2)
!-- GLW         downward long wave flux at ground surface (W/m^2)
!-- EMISS       emissivity of the surface
!-- XLAND       land mask (1 for land, 2 for water)
!-- STBOLT      Stefan-Boltzmann constant (W/m^2/K^4)
!-- F           Coriolis parameter
!-- DT          time step (second)
!-- G           acceleration due to gravity

   integer,  intent(in   )   ::     ids,ide, jds,jde, kds,kde,  &
                                    ims,ime, jms,jme, kms,kme,  &
                                    its,ite, jts,jte, kts,kte

   real,     intent(in   )   ::     DELTSM, STBOLT
   real,    intent(in) ::  OML_RELAXATION_TIME

   real,    dimension( ims:ime, jms:jme )                     , &
            intent(in   )    ::                          EMISS, &
                                                         XLAND, &
                                                           GSW, &
                                                           GLW, &
                                                           HFX, &
                                                            LH

   real,    dimension( ims:ime, jms:jme )                     , &
            intent(inout)    ::                            TSK

   real,    dimension( ims:ime, jms:jme ), intent(inout) ::     &
                                    TML,T0ML,HML,H0ML,HUML,HVML

   real,    dimension( ims:ime, jms:jme, kms:kme ), intent(in   ) ::     &
                                             U_PHY,V_PHY

   real,    dimension( ims:ime, jms:jme ), intent(in   ) ::     &
                                             UST, F, TMOML

   real,    intent(in   )   ::     G
   real,    intent(in   )   ::     OML_GAMMA

! LOCAL VARS

   integer ::  I,J

   do J=jts,jte

!!!      CALL wrf_debug( 100, 'Call OML' )
         do i=its,ite
            if (XLAND(I,J).GT.1.5) then
               call OML1D(I,J,TML(i,j),T0ML(i,j),HML(i,j),H0ML(i,j),           &
                          HUML(i,j),HVML(i,j),TSK(i,j),HFX(i,j),               &
                          LH(i,j),GSW(i,j),GLW(i,j),TMOML(i,j),                &
                          U_PHY(i,j,kts),V_PHY(i,j,kts),UST(i,j),F(i,j),       &
                          EMISS(i,j),STBOLT,G,DELTSM,OML_GAMMA,                &
                          OML_RELAXATION_TIME,                                 &
                          ids,ide, jds,jde, kds,kde,                           &
                          ims,ime, jms,jme, kms,kme,                           &
                          its,ite, jts,jte, kts,kte                            )
            endif
         enddo

   enddo

   end subroutine OCEANML

!----------------------------------------------------------------
   subroutine OML1D(I,J,TML,T0ML,H,H0,HUML,                              &
                    HVML,TSK,HFX,                                        &
                    LH,GSW,GLW,TMOML,                                    &
                    UAIR,VAIR,UST,F,EMISS,STBOLT,G,DT,OML_GAMMA,         &
                    OML_RELAXATION_TIME,                                 &
                    ids,ide, jds,jde, kds,kde,                           &
                    ims,ime, jms,jme, kms,kme,                           &
                    its,ite, jts,jte, kts,kte                            )

!----------------------------------------------------------------
   implicit none
!----------------------------------------------------------------
!
!  SUBROUTINE OCEANML CALCULATES THE SEA SURFACE TEMPERATURE (TSK) 
!  FROM A SIMPLE OCEAN MIXED LAYER MODEL BASED ON 
!  (Pollard, Rhines and Thompson (1973).
!
!-- TML         ocean mixed layer temperature (K)
!-- T0ML        ocean mixed layer temperature (K) at initial time
!-- TMOML       top 200 m ocean mean temperature (K) at initial time
!-- H           ocean mixed layer depth (m)
!-- H0          ocean mixed layer depth (m) at initial time
!-- HUML        ocean mixed layer u component of wind
!-- HVML        ocean mixed layer v component of wind
!-- OML_GAMMA   deep water lapse rate (K m-1)
!-- OMLCALL     whether to call oml model
!-- UAIR,VAIR   lowest model level wind component
!-- UST         frictional velocity
!-- HFX         upward heat flux at the surface (W/m^2)
!-- LH          latent heat flux at the surface (W/m^2)
!-- TSK         surface temperature (K)
!-- GSW         downward short wave flux at ground surface (W/m^2)
!-- GLW         downward long wave flux at ground surface (W/m^2)
!-- EMISS       emissivity of the surface
!-- STBOLT      Stefan-Boltzmann constant (W/m^2/K^4)
!-- F           Coriolis parameter
!-- DT          time step (second)
!-- G           acceleration due to gravity
!-- OML_RELAXATION_TIME  time scale (s) to relax TML to T0ML, H to H0,
!                        HUML and HVML to 0; value <=0 means no relaxation
!
!----------------------------------------------------------------
   integer, intent(in   )    ::      I, J
   integer, intent(in   )    ::      ids,ide, jds,jde, kds,kde, &
                                     ims,ime, jms,jme, kms,kme, &
                                     its,ite, jts,jte, kts,kte

   real,    intent(inout)    :: TML, H, HUML, HVML, TSK

   real,    intent(in   )    :: T0ML, H0, HFX, LH, GSW, GLW,        &
                                UAIR, VAIR, UST, F, EMISS, TMOML

   real,    intent(in) :: STBOLT, G, DT, OML_GAMMA, OML_RELAXATION_TIME

! Local
   real :: rhoair, rhowater, Gam, alp, BV2, A1, A2, B2, u, v, wspd, &
           hu1, hv1, hu2, hv2, taux, tauy, tauxair, tauyair, q, hold, &
           hsqrd, thp, cwater, ust2
   character(LEN=120) :: time_series

      hu1=huml
      hv1=hvml
      rhoair=1.
      rhowater=1000.
      cwater=4200.
! Deep ocean lapse rate (K/m) - from Rich
      Gam=oml_gamma
!     if(i.eq.1 .and. j.eq.1 .or. i.eq.105.and.j.eq.105) print *, 'gamma = ', gam
!     Gam=0.14
!     Gam=5.6/40.
!     Gam=5./100.
! Thermal expansion coeff (/K)
!     alp=.0002
!     temp dependence (/K)
      alp=max((tml-273.15)*1.e-5, 1.e-6)
      BV2=alp*g*Gam
      thp=t0ml-Gam*(h-h0)
      A1=(tml-thp)*h - 0.5*Gam*h*h
      if(h.ne.0.)then
        u=hu1/h
        v=hv1/h
      else
        u=0.
        v=0.
      endif

!  time step

        q=(-hfx-lh+gsw+glw*emiss-stbolt*emiss*tml*tml*tml*tml)/(rhowater*cwater)
!       wspd=max(sqrt(uair*uair+vair*vair),0.1)
        wspd=sqrt(uair*uair+vair*vair)
        if (wspd .lt. 1.e-10 ) then
!          print *, 'i,j,wspd are ', i,j,wspd
           wspd = 1.e-10
        endif
! limit ust to 1.6 to give a value of ust for water of 0.05
!       ust2=min(ust, 1.6)
! new limit for ust: reduce atmospheric ust by half for ocean
        ust2=0.5*ust
        tauxair=ust2*ust2*uair/wspd
        taux=rhoair/rhowater*tauxair
        tauyair=ust2*ust2*vair/wspd
        tauy=rhoair/rhowater*tauyair
! note: forward-backward coriolis force for effective time-centering
        hu2=hu1+dt*( f*hv1 + taux)
        hv2=hv1+dt*(-f*hu2 + tauy)
! consider the flux effect
        A2=A1+q*dt

        huml=hu2
        hvml=hv2

        hold=h
        B2=hu2*hu2+hv2*hv2
        hsqrd=-A2/Gam + sqrt(A2*A2/(Gam*Gam) + 2.*B2/BV2)
        h=sqrt(max(hsqrd,0.0))
! limit to positive h change
        if(h.lt.hold)h=hold

! no change unless tml is warmer than layer mean temp tmol or tsk-5 (see omlinit)
        if(tml.ge.tmoml .and. h.ne.0.)then

! no change unless tml is warmer than layer mean temp tmoml or tsk-5 (see omlinit)
          if(tml.ge.tmoml)then
            tml=max(t0ml - Gam*(h-h0) + 0.5*Gam*h + A2/h, tmoml)
          else 
            tml=tmoml
          endif
          u=hu2/h
          v=hv2/h
        else
          tml=t0ml
          u=0.
          v=0.
        endif

! relax TML T0ML and H to H0, HUML and HVML to 0

        if (oml_relaxation_time .gt. 0.) then
          tml = tml - (tml-t0ml)*dt/oml_relaxation_time
          h = h - (h-h0)*dt/oml_relaxation_time
          huml = huml - huml*dt/oml_relaxation_time
          hvml = hvml - hvml*dt/oml_relaxation_time
        end if

        tsk=tml
!        if(h.gt.100.)print *,i,j,h,tml,' h,tml'

! ww: output point data
!     if( (i.eq.190 .and. j.eq.115) .or. (i.eq.170 .and. j.eq.125) ) then
!        write(jtime,fmt='("TS ",f10.0)') float(itimestep)
!        CALL wrf_message ( TRIM(jtime) )
!        write(time_series,fmt='("OML",2I4,2F9.5,2F8.2,2E15.5,F8.3)') &
!              i,j,u,v,tml,h,taux,tauy,a2
!        CALL wrf_message ( TRIM(time_series) )
!     end if

   end subroutine OML1D

!================================================================
   subroutine omlinit(outfile,oml_hml0, tsk,                   &
                      tml,t0ml,hml,h0ml,huml,hvml,tmoml,       &
                      allowed_to_read, start_of_simulation,    &
                      ids,ide, jds,jde, kds,kde,               &
                      ims,ime, jms,jme, kms,kme,               &
                      its,ite, jts,jte, kts,kte                )
!----------------------------------------------------------------
   implicit none
!----------------------------------------------------------------
   integer, intent(in   )    ::      outfile
   logical , intent(in)      ::      allowed_to_read
   logical , intent(in)      ::      start_of_simulation
   integer, intent(in   )    ::      ids,ide, jds,jde, kds,kde, &
                                     ims,ime, jms,jme, kms,kme, &
                                     its,ite, jts,jte, kts,kte

   real,    dimension( ims:ime, jms:jme )                     , &
            intent(in)    ::                               TSK

   real,    dimension( ims:ime, jms:jme )                     , &
            intent(inout)    ::     TML, T0ML, HML, H0ML, HUML, HVML, TMOML
   real   , intent(in   )    ::     oml_hml0

!  LOCAR VAR

   integer                   ::      L,J,I,itf,jtf
   character*1024 message

!----------------------------------------------------------------
 
   itf=min0(ite,ide-1)
   jtf=min0(jte,jde-1)

   if(start_of_simulation) then
     do J=jts,jtf
     do I=its,itf
       TML(I,J)=TSK(I,J)
       T0ML(I,J)=TSK(I,J)
     enddo
     enddo
     if (oml_hml0 .gt. 0.) then
!!!        WRITE(message,*)'Initializing OML with HML0 = ', oml_hml0
!!!        CALL wrf_debug (0, TRIM(message))
!!!        write(outfile,*) 'Initializing OML with HML0 = ', oml_hml0
        do J=jts,jtf
        do I=its,itf
          HML(I,J)=oml_hml0
          H0ML(I,J)=HML(I,J)
          HUML(I,J)=0.
          HVML(I,J)=0.
          TMOML(I,J)=TSK(I,J)-5.
        enddo
        enddo
     else
!!!        WRITE(message,*)'Initializing OML with real HML0, h(1,1) = ', h0ml(1,1)
!!!        CALL wrf_debug (0, TRIM(message))
!!!        write(outfile,*) 'Initializing OML with real HML0, h(1,1) = ', h0ml(1,1)
        stop 34543
        do J=jts,jtf
        do I=its,itf
          HML(I,J)=H0ML(I,J)
! fill in near coast area with SST: 200 K was set as missing value in ocean pre-processing code
          if(TMOML(I,J).GT.200. .and. TMOML(I,J).LE.201.) TMOML(I,J)=TSK(I,J)
        enddo
        enddo
     endif
   endif

   end subroutine omlinit
!-------------------------------------------------------------------          
end module module_sf_oml
