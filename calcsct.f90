!***********************************************************************
!
      subroutine calcsct(fi,gradfi,ifi)
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
      use obstacle
      use gradients

      implicit none
!
!***********************************************************************
!
      integer, intent(in) :: ifi
      real(prec), dimension(nxyza) :: fi
      real(prec), dimension(3,nxyza) :: gradfi

!
!     local variables
!
      integer :: i, j, k, inp
      real(prec) :: gam, prtr, apotime, sut, urfrs, urfms, &
                    dtdx,dtdy, dtdz, &
                    term1e, term1w, term1n, term1s, term1t, term1b, dterm1dx, &
                    term2e, term2w, term2n, term2s, term2t, term2b, dterm2dy, &
                    term3e, term3w, term3n, term3s, term3t, term3b, dterm3dz, & 
                    edte, stani, dva, att, btt, ctt, &
                    attp, attn, bttp, bttn, cttp, cttn, svip, svin

      gam = gds(ifi)
      prtr = prtinv(ifi)

!.....Calculate gradient: 
      if (lstsq) then
        call grad_lsq_qr(fi,gradfi,2)
      elseif (gauss) then
        call grad_gauss(fi,gradfi(1,:),gradfi(2,:),gradfi(3,:))
      endif

!
!****************************************
      if(ifi.eq.ien) then
!****************************************
      do k = 3,nkmm
      do i = 3,nimm
      do j = 3,njmm

      inp = lk(k)+li(i)+j

!.....unsteady term
      apotime = den(inp)*vol(inp)/timestep
      sut = apotime*((1+btime)*to(inp)-0.5*btime*too(inp))

      sv(inp) = 0.0d0
      su(inp) = 0.0d0
      sp(inp) = 0.0d0

      su(inp) = su(inp)+sut
      sp(inp) = sp(inp)+apotime*(1+0.5*btime)

      end do
      end do
      end do

!.....[Additional terms if turbulent and bouyant flow: ]
      if(lturb.and.lbuoy) then

      call calcheatflux

      do k = 3,nkmm
      do i = 3,nimm
      do j = 3,njmm

      inp = lk(k)+li(i)+j



      term1e = (den(inp+nj)*utt(inp+nj)+gradt(1,inp+nj)*prt1*(vis(inp+nj)-visob(inp+nj)))*fx(inp)+ &
             (den(inp)*utt(inp)+gradt(1,inp)*prt1*(vis(inp)-visob(inp)))*(1.-fx(inp))

      term1w = (den(inp)*utt(inp)+gradt(1,inp)*prt1*(vis(inp)-visob(inp)))*fx(inp-nj)+ &
             (den(inp-nj)*utt(inp-nj)+gradt(1,inp-nj)*prt1*(vis(inp-nj)-visob(inp-nj)))*(1.-fx(inp-nj))

      term1n = (den(inp+1)*utt(inp+1)+gradt(1,inp+1)*prt1*(vis(inp+1)-visob(inp+1)))*fy(inp)+ &
             (den(inp)*utt(inp)+gradt(1,inp)*prt1*(vis(inp)-visob(inp)))*(1.-fy(inp))

      term1s = (den(inp)*utt(inp)+gradt(1,inp)*prt1*(vis(inp)-visob(inp)))*fy(inp-1)+ &
             (den(inp-1)*utt(inp-1)+gradt(1,inp-1)*prt1*(vis(inp-1)-visob(inp-1)))*(1.-fy(inp-1))

      term1t = (den(inp+nij)*utt(inp+nij)+gradt(1,inp+nij)*prt1*(vis(inp+nij)-visob(inp+nij)))*fz(inp)+ &
             (den(inp)*utt(inp)+gradt(1,inp)*prt1*(vis(inp)-visob(inp)))*(1.-fz(inp))

      term1b = (den(inp)*utt(inp)+gradt(1,inp)*prt1*(vis(inp)-visob(inp)))*fz(inp-nij)+ &
             (den(inp-nij)*utt(inp-nij)+gradt(1,inp-nij)*prt1*(vis(inp-nij)-visob(inp-nij)))*(1.-fz(inp-nij))

      dterm1dx = ((term1e-term1w)*ar1x(inp)+ &
                (term1n-term1s)*ar2x(inp)+ &
                (term1t-term1b)*ar3x(inp))
!------------------------------------------------------------
!
      term2e = (den(inp+nj)*vtt(inp+nj)+an(inp+nj)*prt1*(vis(inp+nj)-visob(inp+nj)))*fx(inp)+ &
             (den(inp)*vtt(inp)+gradt(2,inp)*prt1*(vis(inp)-visob(inp)))*(1.-fx(inp))

      term2w = (den(inp)*vtt(inp)+gradt(2,inp)*prt1*(vis(inp)-visob(inp)))*fx(inp-nj)+ &
             (den(inp-nj)*vtt(inp-nj)+gradt(2,inp-nj)*prt1*(vis(inp-nj)-visob(inp-nj)))*(1.-fx(inp-nj))

      term2n = (den(inp+1)*vtt(inp+1)+gradt(2,inp+1)*prt1*(vis(inp+1)-visob(inp+1)))*fy(inp)+ &
             (den(inp)*vtt(inp)+gradt(2,inp)*prt1*(vis(inp)-visob(inp)))*(1.-fy(inp))

      term2s = (den(inp)*vtt(inp)+gradt(2,inp)*prt1*(vis(inp)-visob(inp)))*fy(inp-1)+ &
             (den(inp-1)*vtt(inp-1)+gradt(2,inp-1)*prt1*(vis(inp-1)-visob(inp-1)))*(1.-fy(inp-1))

      term2t = (den(inp+nij)*vtt(inp+nij)+gradt(2,inp+nij)*prt1*(vis(inp+nij)-visob(inp+nij)))*fz(inp)+ &
             (den(inp)*vtt(inp)+gradt(2,inp)*prt1*(vis(inp)-visob(inp)))*(1.-fz(inp))

      term2b = (den(inp)*vtt(inp)+gradt(2,inp)*prt1*(vis(inp)-visob(inp)))*fz(inp-nij)+ &
             (den(inp-nij)*vtt(inp-nij)+gradt(2,inp-nij)*prt1*(vis(inp-nij)-visob(inp-nij)))*(1.-fz(inp-nij))

      dterm2dy = ((term2e-term2w)*ar1y(inp)+ &
                (term2n-term2s)*ar2y(inp)+ &
                (term2t-term2b)*ar3y(inp))
!-----------------------------------------------------------------
!
      term3e = (den(inp+nj)*wtt(inp+nj)+gradt(3,inp+nj)*prt1*(vis(inp+nj)-visob(inp+nj)))*fx(inp)+ &
             (den(inp)*wtt(inp)+gradt(3,inp)*prt1*(vis(inp)-visob(inp)))*(1.-fx(inp))

      term3w = (den(inp)*wtt(inp)+gradt(3,inp)*prt1*(vis(inp)-visob(inp)))*fx(inp-nj)+ &
             (den(inp-nj)*wtt(inp-nj)+gradt(3,inp-nj)*prt1*(vis(inp-nj)-visob(inp-nj)))*(1.-fx(inp-nj))

      term3n = (den(inp+1)*wtt(inp+1)+gradt(3,inp+1)*prt1*(vis(inp+1)-visob(inp+1)))*fy(inp)+ &
             (den(inp)*wtt(inp)+gradt(3,inp)*prt1*(vis(inp)-visob(inp)))*(1.-fy(inp))

      term3s = (den(inp)*wtt(inp)+gradt(3,inp)*prt1*(vis(inp)-visob(inp)))*fy(inp-1)+ &
             (den(inp-1)*wtt(inp-1)+gradt(3,inp-1)*prt1*(vis(inp-1)-visob(inp-1)))*(1.-fy(inp-1))

      term3t = (den(inp+nij)*wtt(inp+nij)+gradt(3,inp+nij)*prt1*(vis(inp+nij)-visob(inp+nij)))*fz(inp)+ &
             (den(inp)*wtt(inp)+gradt(3,inp)*prt1*(vis(inp)-visob(inp)))*(1.-fz(inp))

      term3b = (den(inp)*wtt(inp)+gradt(3,inp)*prt1*(vis(inp)-visob(inp)))*fz(inp-nij)+ &
             (den(inp-nij)*wtt(inp-nij)+gradt(3,inp-nij)*prt1*(vis(inp-nij)-visob(inp-nij)))*(1.-fz(inp-nij))

      dterm3dz = ((term3e-term3w)*ar1z(inp)+ &
                (term3n-term3s)*ar2z(inp)+ &
                (term3t-term3b)*ar3z(inp))
!------------------------------------------------------------------
!
      su(inp) = su(inp)-dterm1dx-dterm2dy-dterm3dz

      end do !i-loop
      end do !j-loop
      end do !k-loop

      end if !(condition loop for lturb & lbuoy)
!
!****************************************
      elseif(ifi.eq.ivart) then
!****************************************
      do k = 3,nkmm
      do i = 3,nimm
      do j = 3,njmm

      inp = lk(k)+li(i)+j

      dtdx = gradfi(1,inp)
      dtdy = gradfi(2,inp)
      dtdz = gradfi(3,inp)

      edte = dabs(ed(inp)/(te(inp)+small))
      stani = den(inp)*edte*vol(inp)/rcost

      dva = -2.*den(inp)*vol(inp)

      att = dva*utt(inp)*dtdx
      btt = dva*vtt(inp)*dtdy
      ctt = dva*wtt(inp)*dtdz

      attp = dmax1(att,zero)
      attn = dmin1(att,zero)
      bttp = dmax1(btt,zero)
      bttn = dmin1(btt,zero)
      cttp = dmax1(ctt,zero)
      cttn = dmin1(ctt,zero)

      svip = attp+bttp+cttp
      svin = attn+bttn+cttn

!.....unsteady term
      apotime = den(inp)*vol(inp)/timestep
      sut = apotime*((1+btime)*varto(inp)-0.5*btime*vartoo(inp))
      sv(inp) = 0.0d0
      su(inp) = 0.0d0
      sp(inp) = 0.0d0

      su(inp) = su(inp)+sut
      sp(inp) = sp(inp)+apotime*(1+0.5*btime)

      sp(inp) = sp(inp)+stani

      su(inp) = su(inp)+svip
      sp(inp) = sp(inp)-svin/dabs(vart(inp)+small)

      enddo
      enddo
      enddo
!
!****************************************
      elseif(ifi.eq.icon) then
!****************************************
      do k = 3,nkmm
      do i = 3,nimm
      do j = 3,njmm

      inp = lk(k)+li(i)+j

!.....unsteady term 
      apotime = den(inp)*vol(inp)/timestep
      sut = apotime*((1+btime)*cono(inp)-0.5*btime*conoo(inp))

      sv(inp) = 0.0d0
      su(inp) = 0.0d0
      sp(inp) = 0.0d0

      su(inp) = su(inp)+sut
      sp(inp) = sp(inp)+apotime*(1+0.5*btime)

      enddo
      enddo
      enddo

      end if

!
!.....calculate terms integrated over faces

!.....east cell - face
      call fluxsct(nj,1,nij,fi,gradfi,ifi, &
                   ar1x,ar1y,ar1z, &
                   fx,ae,aw,f1)
!.....north cell - face
      call fluxsct(1,nij,nj,fi,gradfi,ifi, &
                   ar2x,ar2y,ar2z, &
                   fy,an,as,f2)
!.....top   cell - face
      call fluxsct(nij,nj,1,fi,gradfi,ifi, &
                   ar3x,ar3y,ar3z, &
                   fz,at,ab,f3)

!.....main diagonal coefficient and element of rhs vector; under-relaxation
      urfrs = urfr(ifi)
      urfms = urfm(ifi)

      do k = 3,nkmm
      do i = 3,nimm
      do j = 3,njmm

      inp = lk(k)+li(i)+j

      ap(inp) = ae(inp)+aw(inp)+an(inp)+as(inp)+at(inp)+ab(inp)+sp(inp)
      ap(inp) = ap(inp)*urfrs

      su(inp) = su(inp)+urfms*ap(inp)*fi(inp)

      end do
      end do
      end do
!
!.....solving linear system:
      call sipsol(fi,ifi)
      !call cgstab(fi,ifi)
      !call iccg(fi,ifi)
      !call pcgsip(fi,ifi)
      !call cgstab_sip(fi,ifi)  


      if(ifi.eq.ivart) then
        do inp = icst,icen
          fi(inp) = max(fi(inp),small)
        enddo
      endif

      return
      end
