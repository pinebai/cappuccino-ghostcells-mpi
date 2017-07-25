!***********************************************************************
!
      subroutine calcstress
!
!***********************************************************************
!
      use types
      use parameters
      use indexes
      use geometry
      use coef
      use coefb
      use variables
      use buoy
      use time_mod
      use gradients
      use omega_turb_models

      implicit none 
!
!***********************************************************************
!
      integer :: i, j, k, inp
      real(prec) :: vist, &                               ! turbulent viscosity
                    dudx, dudy, dudz, &                   ! gradu - the velocity gradient
                    dvdx, dvdy, dvdz, &                   ! gradv - the velocity gradient
                    dwdx, dwdy, dwdz, &                   ! gradw - the velocity gradient
                    csasa, &                              !
                    ed11, ed22, ed33, ed12, ed13, ed23, & ! turb. dissipation tensor components
                    p11, p22, p33, p12, p13, p23, p123, & !
                    g11, g22, g33, g12, g13, g23, g123, & ! bouyancy contribution terms
                    uuold,vvold, wwold, &                 ! reynolds stress tensor 
                    uvold, uwold,vwold                    ! components 
!
!=========================================
!     [turbulent stresses calculations: ]
!=========================================
!
      do k=3,nkmm
      do i=3,nimm
      do j=3,njmm

      inp=lk(k)+li(i)+j

      vist=(vis(inp)-viscos)/densit

      uuold=uu(inp)
      vvold=vv(inp)
      wwold=ww(inp)
      uvold=uv(inp)
      uwold=uw(inp)
      vwold=vw(inp)
!
!--------------------------------
!      [velocity gradients: ]
!--------------------------------
      dudx = gradu(1,inp)
      dudy = gradu(2,inp)
      dudz = gradu(3,inp)
      
      dvdx = gradv(1,inp)
      dvdy = gradv(2,inp)
      dvdz = gradv(3,inp)

      dwdx = gradw(1,inp)
      dwdy = gradw(2,inp)
      dwdz = gradw(3,inp)
!
!==============================================
!---[evm approach: ]
!==============================================
      if(levm) then

      uu(inp)=2.0d0/3.0d0*te(inp)-vist*(dudx+dudx)
      vv(inp)=2.0d0/3.0d0*te(inp)-vist*(dvdy+dvdy)
      ww(inp)=2.0d0/3.0d0*te(inp)-vist*(dwdz+dwdz)

      uv(inp)=-vist*(dudy+dvdx)
      uw(inp)=-vist*(dudz+dwdx)
      vw(inp)=-vist*(dvdz+dwdy)
!
!==============================================

!==============================================
!-----------------------------------------------------
!     don't forget to include exact production gen !!
!-----------------------------------------------------
      else if(lasm) then

!.....wallin-johansson earsm
      uu(inp) = 2.0d0/3.0d0*te(inp)-vist*(dudx+dudx) + 2.*bij(1,inp)*te(inp)
      vv(inp) = 2.0d0/3.0d0*te(inp)-vist*(dvdy+dvdy) + 2.*bij(4,inp)*te(inp)
!.....it seems that b(3,3)=-b(1,1)-b(2,2):
      ww(inp) = 2.0d0/3.0d0*te(inp)-vist*(dwdz+dwdz)  &
              - 2.*(bij(1,inp)+bij(4,inp))*te(inp) 

      uv(inp) = -vist*(dudy+dvdx) + 2.*bij(2,inp)*te(inp)
      uw(inp) = -vist*(dudz+dwdx) + 2.*bij(3,inp)*te(inp)
      vw(inp) = -vist*(dvdz+dwdy) + 2.*bij(5,inp)*te(inp)

!---[asm approach: ip-model]
      if (1.eq.0) then ! switched-off for now

      csasa=te(inp)/(c1asm*ed(inp)+small)

      ed11=2.0d0/3.0d0*ed(inp)
      ed22=2.0d0/3.0d0*ed(inp)
      ed33=2.0d0/3.0d0*ed(inp)
      ed12=0.
      ed13=0.
      ed23=0.

      p11=-2.*(uu(inp)*dudx+uv(inp)*dudy+uw(inp)*dudz)
      p22=-2.*(uv(inp)*dvdx+vv(inp)*dvdy+vw(inp)*dvdz)
      p33=-2.*(uw(inp)*dwdx+vw(inp)*dwdy+ww(inp)*dwdz)

      p12=-(uu(inp)*dvdx+uv(inp)*(dvdy+dudx)+uw(inp)*dvdz+ &
            vv(inp)*dudy+vw(inp)*dudz)
      p13=-(uu(inp)*dwdx+uw(inp)*(dwdz+dudx)+uv(inp)*dwdy+ &
            vw(inp)*dudy+ww(inp)*dudz)
      p23=-(uv(inp)*dwdx+vw(inp)*(dvdy+dwdz)+vv(inp)*dwdy+ &
            uw(inp)*dvdx+ww(inp)*dvdz)

      p123=0.5*(p11+p22+p33)

!-----
      uu(inp)=2*te(inp)/3.+ &
              csasa*((1-c2asm)*p11+2/3.*c2asm*p123-ed11)
      vv(inp)=2*te(inp)/3.+ &
              csasa*((1-c2asm)*p22+2/3.*c2asm*p123-ed22)
      ww(inp)=2*te(inp)/3.+ &
              csasa*((1-c2asm)*p33+2/3.*c2asm*p123-ed33)

      uv(inp)=csasa*((1-c2asm)*p12-ed12)
      uw(inp)=csasa*((1-c2asm)*p13-ed13)
      vw(inp)=csasa*((1-c2asm)*p23-ed23)

!
!--------
      if(lbuoy) then

      if(boussinesq) then
            g11=-2.*beta*gravx*utt(inp)
            g22=-2.*beta*gravy*vtt(inp)
            g33=-2.*beta*gravz*wtt(inp)
            g12=-beta*(gravx*vtt(inp)+gravy*utt(inp))
            g13=-beta*(gravx*wtt(inp)+gravz*utt(inp))
            g23=-beta*(gravy*wtt(inp)+gravz*vtt(inp))
      else !if(boussinesq.eq.0) then
            g11=-2.*gravx*utt(inp)/(t(inp)+273.)
            g22=-2.*gravy*vtt(inp)/(t(inp)+273.)
            g33=-2.*gravz*wtt(inp)/(t(inp)+273.)
            g12=-(gravx*vtt(inp)+gravy*utt(inp))/(t(inp)+273.)
            g13=-(gravx*wtt(inp)+gravz*utt(inp))/(t(inp)+273.)
            g23=-(gravy*wtt(inp)+gravz*vtt(inp))/(t(inp)+273.)
      end if

      g123=0.5*(g11+g22+g33)

      uu(inp)=uu(inp)+csasa*((1-c3asm)*g11+2/3.*c3asm*g123)
      vv(inp)=vv(inp)+csasa*((1-c3asm)*g22+2/3.*c3asm*g123)
      ww(inp)=ww(inp)+csasa*((1-c3asm)*g33+2/3.*c3asm*g123)

      uv(inp)=uv(inp)+csasa*(1-c3asm)*g12
      uw(inp)=uw(inp)+csasa*(1-c3asm)*g13
      vw(inp)=vw(inp)+csasa*(1-c3asm)*g23

      end if !![for buoyancy asm contribution]
!--------
 
!==========================================================
!     reynolds stresses using the anisotropy tensor:
!     ______
!     u_iu_j = k(a_ij + 2/3 \delta_ij)
!     a_ij = 2* b_ij   
!==========================================================
      uu(inp)=2*bij(1,inp)*te(inp) + 2*te(inp)/3.
      vv(inp)=2*bij(4,inp)*te(inp) + 2*te(inp)/3.
      ww(inp)=0.                   + 2*te(inp)/3.

      uv(inp)=2*bij(2,inp)*te(inp)
      uw(inp)=2*bij(3,inp)*te(inp)
      vw(inp)=2*bij(5,inp)*te(inp)

      endif ! 1.eq.0

      end if !![for evm or asm approach]
!----------------------------------------------

      uu(inp)=max(uu(inp),small)
      vv(inp)=max(vv(inp),small)
      ww(inp)=max(ww(inp),small)
!---
      uu(inp)=facnap*uu(inp)+(1.-facnap)*uuold
      vv(inp)=facnap*vv(inp)+(1.-facnap)*vvold
      ww(inp)=facnap*ww(inp)+(1.-facnap)*wwold

      uv(inp)=facnap*uv(inp)+(1.-facnap)*uvold
      uw(inp)=facnap*uw(inp)+(1.-facnap)*uwold
      vw(inp)=facnap*vw(inp)+(1.-facnap)*vwold
!---
      end do ! i-loop
      end do ! j-loop
      end do ! k-loop
!-------------------------------
      return
      end
