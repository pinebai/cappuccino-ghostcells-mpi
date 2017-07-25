!.....Velocity gradients: 
      if (lstsq) then
        call grad_lsq_qr(u,gradu,2)
        call grad_lsq_qr(v,gradv,2)
        call grad_lsq_qr(w,gradw,2)
      elseif (gauss) then
        call grad_gauss(u,gradu(1,:),gradu(2,:),gradu(3,:))
        call grad_gauss(v,gradv(1,:),gradv(2,:),gradv(3,:))
        call grad_gauss(w,gradw(1,:),gradw(2,:),gradw(3,:))
      endif

!.....only inner faces
!.....east cell - face
      call fluxmass_plain(nj,1,nij, &
                          ar1x,ar1y,ar1z, &
                          fx,f1)
!.....north cell - face
      call fluxmass_plain(1,nij,nj, &
                          ar2x,ar2y,ar2z, &
                          fy,f2)
!.....top   cell - face
      call fluxmass_plain(nij,nj,1, &
                          ar3x,ar3y,ar3z, &
                          fz,f3)

!// adjusts the inlet and outlet fluxes to obey continuity, which is necessary for creating a well-posed
!// problem where a solution for pressure exists.
!     adjustPhi(phi, U, p);
      if(.not.const_mflux) call outbc
