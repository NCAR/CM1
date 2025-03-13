  MODULE pdef_module

  implicit none

  private
  public :: pdefx1,pdefx2,pdefy1,pdefy2,pdefz

  CONTAINS

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine pdefx1(xh,arh1,arh2,uh,rho0,gz,rgz,rru,advx,dum,mass,s0,s,dt,flag,west,newwest,east,neweast,reqs_s)
      use input
      use bc_module
      use comm_module
      implicit none

      real, intent(in), dimension(ib:ie) :: xh,arh1,arh2,uh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho0,s0,s
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: rru
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: advx,dum,mass
      real, intent(in) :: dt
      logical, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: flag
      real, intent(inout), dimension(cmp,jmp,kmp) :: west,newwest,east,neweast
      integer, intent(inout), dimension(4) :: reqs_s

      integer i,j,k
      real foo1,foo2,foo3,rdt

!----------------------------------------------------------------
! cm1r17:  include divx component

      IF(.not.terrain_flag)THEN
        IF(axisymm.eq.0)THEN
          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            dum(i,j,k)=rho0(1,1,k)*s0(i,j,k)+dt*( advx(i,j,k)  &
                   +s(i,j,k)*(rru(i+1,j,k)-rru(i,j,k))*rdx*uh(i) )
          enddo
          enddo
          enddo
        ELSE
          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            dum(i,j,k)=rho0(1,1,k)*s0(i,j,k)+dt*( advx(i,j,k)  &
                   +s(i,j,k)*(arh2(i)*rru(i+1,j,k)-arh1(i)*rru(i,j,k))*rdx*uh(i) )
          enddo
          enddo
          enddo
        ENDIF
      ELSE
          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            dum(i,j,k)=rho0(i,j,k)*s0(i,j,k)+dt*( advx(i,j,k)  &
                   +s(i,j,k)*(rru(i+1,j,k)-rru(i,j,k))*rdx*uh(i) )*gz(i,j)
          enddo
          enddo
          enddo
      ENDIF

        if(wbc.eq.2 .and. ibw.eq.1)then
          i=1
          !$omp parallel do default(shared)   &
          !$omp private(j,k)
          do k=1,nk
          do j=1,nj
            if(rru(1,j,k).ge.0.0)then
              dum(i,j,k)=rho0(i,j,k)*s0(i,j,k)
            endif
          enddo
          enddo
        endif

        if(ebc.eq.2 .and. ibe.eq.1)then
          i=ni
          !$omp parallel do default(shared)   &
          !$omp private(j,k)
          do k=1,nk
          do j=1,nj
            if(rru(ni+1,j,k).le.0.0)then
              dum(i,j,k)=rho0(i,j,k)*s0(i,j,k)
            endif
          enddo
          enddo
        endif

        if(timestats.ge.1) time_pdef=time_pdef+mytime()

        call bcs(dum)
#ifdef MPI
        call comm_2we_start(dum,west(1,1,1),newwest(1,1,1),east(1,1,1),neweast(1,1,1),reqs_s)
#endif

!----------------------------------------------------------------

      end subroutine pdefx1


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine pdefx2(xh,arh1,arh2,uh,rho0,gz,rgz,rru,advx,dum,mass,s0,s,dt,flag,west,newwest,east,neweast,reqs_s)
      use input
      use comm_module
      implicit none

      real, intent(in), dimension(ib:ie) :: xh,arh1,arh2,uh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho0,s0,s
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: rru
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: advx,dum,mass
      real, intent(in) :: dt
      logical, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: flag
      real, intent(inout), dimension(cmp,jmp,kmp) :: west,newwest,east,neweast
      integer, intent(inout), dimension(4) :: reqs_s

      integer i,j,k
      real foo1,foo2,foo3,rdt

!----------------------------------------------------------------
! cm1r17:  include divx component

#ifdef MPI
        call comm_2we_end(dum,west(1,1,1),newwest(1,1,1),east(1,1,1),neweast(1,1,1),reqs_s)
#endif

        rdt=1.0/dt

    IF(axisymm.eq.0)THEN

      !$omp parallel do default(shared)   &
      !$omp private(i,j,k,foo1,foo2,foo3)
      DO k=1,nk

        do j=1,nj
        do i=-1,ni+2
          mass(i,j,k)=0.0
          flag(i,j,k)=.false.
        enddo
        enddo

        do j=1,nj
        do i=0,ni+1
          if(dum(i,j,k).lt.0.0)then
            foo1=max(0.0,dum(i-1,j,k))
            foo2=max(0.0,dum(i+1,j,k))
            if(foo1+foo2.gt.smeps)then
              foo3=max(dum(i,j,k),-(foo1+foo2))/(foo1+foo2)
              mass(i-1,j,k)=mass(i-1,j,k)+foo1*foo3
              mass(i  ,j,k)=mass(i  ,j,k)-(foo1+foo2)*foo3
              mass(i+1,j,k)=mass(i+1,j,k)+foo2*foo3
              if(dum(i-1,j,k).gt.smeps) flag(i-1,j,k)=.true.
                                        flag(i  ,j,k)=.true.
              if(dum(i+1,j,k).gt.smeps) flag(i+1,j,k)=.true.
            endif
          endif
        enddo
        enddo

      ENDDO

        !-----
        IF(.not.terrain_flag)THEN
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
        if(flag(i,j,k))then
          advx(i,j,k)=(dum(i,j,k)+mass(i,j,k)-rho0(1,1,k)*s0(i,j,k))*rdt  &
                   -s(i,j,k)*(rru(i+1,j,k)-rru(i,j,k))*rdx*uh(i)
        endif
        enddo
        enddo
        enddo
        !-----
        ELSE
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
        if(flag(i,j,k))then
          advx(i,j,k)=(dum(i,j,k)+mass(i,j,k)-rho0(i,j,k)*s0(i,j,k))*rdt*rgz(i,j)  &
                   -s(i,j,k)*(rru(i+1,j,k)-rru(i,j,k))*rdx*uh(i)
        endif
        enddo
        enddo
        enddo
        ENDIF
        !-----

    ELSE

      !$omp parallel do default(shared)   &
      !$omp private(i,j,k,foo1,foo2,foo3)
      DO k=1,nk

        do j=1,nj
        do i=-1,ni+2
          mass(i,j,k)=0.0
          flag(i,j,k)=.false.
        enddo
        enddo

        do j=1,nj
        do i=0,ni+1
          if(dum(i,j,k).lt.0.0)then
            foo1=max(0.0,dum(i-1,j,k))
            foo2=max(0.0,dum(i+1,j,k))
            if(foo1+foo2.gt.smeps)then
              foo3=max(xh(i)*dum(i,j,k),-(xh(i-1)*foo1+xh(i+1)*foo2))   &
                                        /(xh(i-1)*foo1+xh(i+1)*foo2)
              mass(i-1,j,k)=mass(i-1,j,k)+foo1*foo3
              mass(i  ,j,k)=mass(i  ,j,k)-(foo1+foo2)*foo3
              mass(i+1,j,k)=mass(i+1,j,k)+foo2*foo3
              if(dum(i-1,j,k).gt.smeps) flag(i-1,j,k)=.true.
                                        flag(i  ,j,k)=.true.
              if(dum(i+1,j,k).gt.smeps) flag(i+1,j,k)=.true.
            endif
          endif
        enddo
        enddo

      ENDDO

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,foo1,foo2,foo3)
        do k=1,nk
        do j=1,nj
        do i=1,ni
        if(flag(i,j,k))then
          advx(i,j,k)=(dum(i,j,k)+mass(i,j,k)-rho0(1,1,k)*s0(i,j,k))*rdt  &
                   -s(i,j,k)*(arh2(i)*rru(i+1,j,k)-arh1(i)*rru(i,j,k))*rdx*uh(i)
        endif
        enddo
        enddo
        enddo

    ENDIF

!----------------------------------------------------------------

      if(timestats.ge.1) time_pdef=time_pdef+mytime()

      end subroutine pdefx2


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine pdefy1(vh,rho0,gz,rgz,rrv,advy,dum,mass,s0,s,dt,flag,south,newsouth,north,newnorth,reqs_s)
      use input
      use bc_module
      use comm_module
      implicit none

      real, intent(in), dimension(jb:je) :: vh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho0,s0,s
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: rrv
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: advy,dum,mass
      real, intent(in) :: dt
      logical, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: flag
      real, intent(inout), dimension(imp,cmp,kmp) :: south,newsouth,north,newnorth
      integer, intent(inout), dimension(4) :: reqs_s

      integer i,j,k
      real foo1,foo2,foo3,rdt

!----------------------------------------------------------------
! cm1r17:  include divx component

      IF(.not.terrain_flag)THEN
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum(i,j,k)=rho0(1,1,k)*s0(i,j,k)+dt*( advy(i,j,k)  &
                       +s(i,j,k)*(rrv(i,j+1,k)-rrv(i,j,k))*rdy*vh(j) )
        enddo
        enddo
        enddo
      ELSE
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum(i,j,k)=rho0(i,j,k)*s0(i,j,k)+dt*( advy(i,j,k)  &
                       +s(i,j,k)*(rrv(i,j+1,k)-rrv(i,j,k))*rdy*vh(j) )*gz(i,j)
        enddo
        enddo
        enddo
      ENDIF

        if(sbc.eq.2 .and. ibs.eq.1)then
          j=1
          !$omp parallel do default(shared)   &
          !$omp private(i,k)
          do k=1,nk
          do i=1,ni
            if(rrv(i,1,k).ge.0.0)then
              dum(i,j,k)=rho0(i,j,k)*s0(i,j,k)
            endif
          enddo
          enddo
        endif

        if(nbc.eq.2 .and. ibn.eq.1)then
          j=nj
          !$omp parallel do default(shared)   &
          !$omp private(i,k)
          do k=1,nk
          do i=1,ni
            if(rrv(i,nj+1,k).le.0.0)then
              dum(i,j,k)=rho0(i,j,k)*s0(i,j,k)
            endif
          enddo
          enddo
        endif

        if(timestats.ge.1) time_pdef=time_pdef+mytime()

        call bcs(dum)
#ifdef MPI
        call comm_2sn_start(dum,south(1,1,1),newsouth(1,1,1),north(1,1,1),newnorth(1,1,1),reqs_s)
#endif

!----------------------------------------------------------------

      end subroutine pdefy1


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine pdefy2(vh,rho0,gz,rgz,rrv,advy,dum,mass,s0,s,dt,flag,south,newsouth,north,newnorth,reqs_s)
      use input
      use comm_module
      implicit none

      real, intent(in), dimension(jb:je) :: vh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho0,s0,s
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: rrv
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: advy,dum,mass
      real, intent(in) :: dt
      logical, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: flag
      real, intent(inout), dimension(imp,cmp,kmp) :: south,newsouth,north,newnorth
      integer, intent(inout), dimension(4) :: reqs_s

      integer i,j,k
      real foo1,foo2,foo3,rdt

!----------------------------------------------------------------
! cm1r17:  include divx component

#ifdef MPI
        call comm_2sn_end(dum,south(1,1,1),newsouth(1,1,1),north(1,1,1),newnorth(1,1,1),reqs_s)
#endif

      !$omp parallel do default(shared)   &
      !$omp private(i,j,k,foo1,foo2,foo3)
      DO k=1,nk

        do j=-1,nj+2
        do i=1,ni
          mass(i,j,k)=0.0
          flag(i,j,k)=.false.
        enddo
        enddo

        do j=0,nj+1
        do i=1,ni
          if(dum(i,j,k).lt.0.0)then
            foo1=max(0.0,dum(i,j-1,k))
            foo2=max(0.0,dum(i,j+1,k))
            if(foo1+foo2.gt.smeps)then
              foo3=max(dum(i,j,k),-(foo1+foo2))/(foo1+foo2)
              mass(i,j-1,k)=mass(i,j-1,k)+foo1*foo3
              mass(i,j  ,k)=mass(i,j  ,k)-(foo1+foo2)*foo3
              mass(i,j+1,k)=mass(i,j+1,k)+foo2*foo3
              if(dum(i,j-1,k).gt.smeps) flag(i,j-1,k)=.true.
                                        flag(i,j  ,k)=.true.
              if(dum(i,j+1,k).gt.smeps) flag(i,j+1,k)=.true.
            endif
          endif
        enddo
        enddo

      ENDDO

        rdt=1.0/dt

        !-----
        IF(.not.terrain_flag)THEN
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
        if(flag(i,j,k))then
          advy(i,j,k)=(dum(i,j,k)+mass(i,j,k)-rho0(1,1,k)*s0(i,j,k))*rdt  &
                       -s(i,j,k)*(rrv(i,j+1,k)-rrv(i,j,k))*rdy*vh(j)
        endif
        enddo
        enddo
        enddo
        !-----
        ELSE
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
        if(flag(i,j,k))then
          advy(i,j,k)=(dum(i,j,k)+mass(i,j,k)-rho0(i,j,k)*s0(i,j,k))*rdt*rgz(i,j)  &
                       -s(i,j,k)*(rrv(i,j+1,k)-rrv(i,j,k))*rdy*vh(j)
        endif
        enddo
        enddo
        enddo
        ENDIF
        !-----

!----------------------------------------------------------------

      if(timestats.ge.1) time_pdef=time_pdef+mytime()

      end subroutine pdefy2


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine pdefz(mh,rho0,gz,rgz,rdsf,rrw,advz,dum,mass,s0,s,dt,flag)
      use input
      implicit none

      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh,rho0,s0,s
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz
      real, intent(in), dimension(kb:ke+1) :: rdsf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: rrw
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: advz,dum,mass
      real, intent(in) :: dt
      logical, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: flag

      integer i,j,k
      real foo1,foo2,foo3,rdt

!----------------------------------------------------------------

      IF(.not.terrain_flag)THEN
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum(i,j,k)=rho0(1,1,k)*s0(i,j,k)+dt*( advz(i,j,k)  &
                       +s(i,j,k)*(rrw(i,j,k+1)-rrw(i,j,k))*rdz*mh(1,1,k) )
        enddo
        enddo
        enddo
      ELSE
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum(i,j,k)=rho0(i,j,k)*s0(i,j,k)+dt*( advz(i,j,k)  &
                       +s(i,j,k)*(rrw(i,j,k+1)-rrw(i,j,k))*rdsf(k) )*gz(i,j)
        enddo
        enddo
        enddo
      ENDIF

        !$omp parallel do default(shared)   &
        !$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum(i,j, 0)=0.0
          dum(i,j,nk+1)=0.0
        enddo
        enddo

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=0,nk+1
        do j=1,nj
        do i=1,ni
          mass(i,j,k)=0.0
          flag(i,j,k)=.false.
        enddo
        enddo
        enddo

        do k=1,nk

        !$omp parallel do default(shared)   &
        !$omp private(i,j,foo1,foo2,foo3)
        do j=1,nj
        do i=1,ni
          if(dum(i,j,k).lt.0.0)then
            foo1=max(0.0,dum(i,j,k-1))
            foo2=max(0.0,dum(i,j,k+1))
            if(foo1+foo2.gt.smeps)then
              foo3=max(dum(i,j,k),-(foo1+foo2))/(foo1+foo2)
              mass(i,j,k-1)=mass(i,j,k-1)+foo1*foo3
              mass(i,j,k  )=mass(i,j,k  )-(foo1+foo2)*foo3
              mass(i,j,k+1)=mass(i,j,k+1)+foo2*foo3
              if(dum(i,j,k-1).gt.smeps) flag(i,j,k-1)=.true.
                                        flag(i,j,k  )=.true.
              if(dum(i,j,k+1).gt.smeps) flag(i,j,k+1)=.true.
            endif
          endif
        enddo
        enddo

        enddo

        rdt=1.0/dt

        !-----
        IF(.not.terrain_flag)THEN
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
        if(flag(i,j,k))then
          advz(i,j,k)=(dum(i,j,k)+mass(i,j,k)-rho0(1,1,k)*s0(i,j,k))*rdt  &
                       -s(i,j,k)*(rrw(i,j,k+1)-rrw(i,j,k))*rdz*mh(1,1,k)
        endif
        enddo
        enddo
        enddo
        !-----
        ELSE
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
        if(flag(i,j,k))then
          advz(i,j,k)=(dum(i,j,k)+mass(i,j,k)-rho0(i,j,k)*s0(i,j,k))*rdt*rgz(i,j)  &
                       -s(i,j,k)*(rrw(i,j,k+1)-rrw(i,j,k))*rdsf(k)
        endif
        enddo
        enddo
        enddo
        ENDIF
        !-----

!----------------------------------------------------------------

      if(timestats.ge.1) time_pdef=time_pdef+mytime()

      end subroutine pdefz


  END MODULE pdef_module
