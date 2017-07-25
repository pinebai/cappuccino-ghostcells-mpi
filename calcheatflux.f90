      subroutine calcheatflux
!***********************************************************************
!
!
!     Turbulent heat flux calculations
!
!***********************************************************************
!
      use types
      use parameters
      use indexes
      use geometry
      use variables
      use buoy
      use gradients

      implicit none 
!
!***********************************************************************
!
      integer :: i, j, k, inp
      real(prec) :: vist, &
                    tedi, &                                     ! k/eps
                    dtdx, dtdy, dtdz, &                         ! temperature gradient vector components
                    ut1, ut2, ut3, ut4, ut5, ut6, ut7, &        ! needed for afm approach
                    vt1, vt2, vt3, vt4, vt5, vt6, vt7, &
                    wt1, wt2, wt3, wt4, wt5, wt6, wt7, &
                    dudx, dudy, dudz, &                         ! gradu - the velocity gradient
                    dvdx, dvdy, dvdz, &                         !gradv - the velocity gradient
                    dwdx, dwdy, dwdz                            !gradw - the velocity gradient
      real(prec) :: uttold, vttold, wttold                !heat flux vector compnnts


      do k = 3,nkmm
      do i = 3,nimm
      do j = 3,njmm

      inp = lk(k)+li(i)+j

      vist = (vis(inp)-viscos)/densit
      tedi = te(inp)/(ed(inp)+small)

      uttold = utt(inp)
      vttold = vtt(inp)
      wttold = wtt(inp)
!
!--------------------------------
!      [temperature gradients: ]
!--------------------------------
      dtdx = gradt(1,inp)
      dtdy = gradt(2,inp)
      dtdz = gradt(3,inp)

!.....[sgdh approach: ] 
      if(lsgdh) then
      utt(inp) = -(vist*prt1)*dtdx
      vtt(inp) = -(vist*prt1)*dtdy
      wtt(inp) = -(vist*prt1)*dtdz

!.....[ggdh approach: ]
      else if(lggdh) then
      ut1 = uu(inp)*dtdx
      ut2 = uv(inp)*dtdy
      ut3 = uw(inp)*dtdz

      vt1 = uv(inp)*dtdx
      vt2 = vv(inp)*dtdy
      vt3 = vw(inp)*dtdz

      wt1 = uw(inp)*dtdx
      wt2 = vw(inp)*dtdy
      wt3 = ww(inp)*dtdz

      utt(inp) = -phit*tedi*(ut1+ut2+ut3)
      vtt(inp) = -phit*tedi*(vt1+vt2+vt3)
      wtt(inp) = -phit*tedi*(wt1+wt2+wt3)

!.....[afm: include all]
      else if(lafm) then

!      [velocity gradients: ]
      dudx = gradu(1,inp)
      dudy = gradu(2,inp)
      dudz = gradu(3,inp)

      dvdx = gradv(1,inp)
      dvdy = gradv(2,inp)
      dvdz = gradv(3,inp)

      dwdx = gradw(1,inp)
      dwdy = gradw(2,inp)
      dwdz = gradw(3,inp)

!----------------------------------------------------------------
!
      ut1 = uu(inp)*dtdx
      ut2 = uv(inp)*dtdy
      ut3 = uw(inp)*dtdz
      ut4 = sksi*utt(inp)*dudx
      ut5 = sksi*vtt(inp)*dudy
      ut6 = sksi*wtt(inp)*dudz
      ut7 = gravx*eta*vart(inp)*beta
      if(.not.boussinesq) ut7 = gravx*eta*vart(inp)/(t(inp)+273.)

      vt1 = uv(inp)*dtdx
      vt2 = vv(inp)*dtdy
      vt3 = vw(inp)*dtdz
      vt4 = sksi*utt(inp)*dvdx
      vt5 = sksi*vtt(inp)*dvdy
      vt6 = sksi*wtt(inp)*dvdz
      vt7 = gravy*eta*vart(inp)*beta
      if(.not.boussinesq) vt7 = gravy*eta*vart(inp)/(t(inp)+273.)

      wt1 = uw(inp)*dtdx
      wt2 = vw(inp)*dtdy
      wt3 = ww(inp)*dtdz
      wt4 = sksi*utt(inp)*dwdx
      wt5 = sksi*vtt(inp)*dwdy
      wt6 = sksi*wtt(inp)*dwdz
      wt7 = gravz*eta*vart(inp)*beta
      if(.not.boussinesq) wt7 = gravz*eta*vart(inp)/(t(inp)+273.)

      utt(inp) = -phit*tedi*(ut1+ut2+ut3+ut4+ut5+ut6+ut7)
      vtt(inp) = -phit*tedi*(vt1+vt2+vt3+vt4+vt5+vt6+vt7)
      wtt(inp) = -phit*tedi*(wt1+wt2+wt3+wt4+wt5+wt6+wt7)

!      if(i.eq.10.and.j.eq.10.and.k.eq.2) write(6,*)'wtt: ',wtt(inp)
!      if(i.eq.10.and.j.eq.10.and.k.eq.2) write(6,*)'vart: ',vart(inp)
!      if(i.eq.10.and.j.eq.10.and.k.eq.2) write(6,*)'t: ',t(inp)
!      if(i.eq.10.and.j.eq.10.and.k.eq.2) write(6,*)'w: ',w(inp)

      end if !! [conditional for sgdh,ggdh,afm]
!-------------------------------------------------------
!
      utt(inp) = facflx*utt(inp)+(1.-facflx)*uttold
      vtt(inp) = facflx*vtt(inp)+(1.-facflx)*vttold
      wtt(inp) = facflx*wtt(inp)+(1.-facflx)*wttold

      end do ! j-loop
      end do ! i-loop
      end do ! k-loop

      return
      end
