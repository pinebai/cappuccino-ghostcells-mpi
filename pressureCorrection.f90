!
!.....Divergence of the Hi operator (i=1,2,3) which is a source term of explicit Convection + explicit diffusion
!     contained in su, sv, sw arrays.



! !.....Velocity gradients (will be used for interpolation in peqn_flux function: 
!       if (lstsq) then
! !=====lstsq=============================================================  
!       call grad_lsq_qr(u,gradu,2)
!       call grad_lsq_qr(v,gradv,2)
!       call grad_lsq_qr(w,gradw,2)

!       elseif (gauss) then
! !=====gauss=============================================================
!       call grad_gauss(u,gradu(1,:),gradu(2,:),gradu(3,:))
!       call grad_gauss(v,gradv(1,:),gradv(2,:),gradv(3,:))
!       call grad_gauss(w,gradw(1,:),gradw(2,:),gradw(3,:))
!       endif

!.....ASSEMBLE OFF DIAGONAL ENTRIES OF Pressure Eq. SYSTEM MATRIX,
!// calculate the fluxes F1,F2,F3 by dotting the interpolated velocity (to cell faces) with face normals

!.....only inner surfaces
!.....East cell - face
      call assemble_pressure_eq(nj,1,nij, &
                                ar1x,ar1y,ar1z, &
                                fx,ae,aw,f1)
!.....North cell - face
      call assemble_pressure_eq(1,nij,nj, &
                                ar2x,ar2y,ar2z, &
                                fy,an,as,f2)
!.....Top cell - face
      call assemble_pressure_eq(nij,nj,1, &
                                ar3x,ar3y,ar3z, &
                                fz,at,ab,f3)




!// adjusts the inlet and outlet fluxes to obey continuity, which is necessary for creating a well-posed
!// problem where a solution for pressure exists.
!     adjustPhi(phi, U, p);
      ! if(.not.const_mflux) call outbc





!.....Assemble The System's Main Diagonal Entry and rhs Vector Stored in Su Array 
      idew=nj
      idns=1
      idtb=nij

      do k=3,nkmm
        do i=3,nimm
          do j=3,njmm

            inp=lk(k)+li(i)+j

            ap(inp)=ae(inp)+aw(inp)+an(inp)+as(inp)+at(inp)+ab(inp)

            su(inp)=-f1(inp)+f1(inp-idew) &
                    -f2(inp)+f2(inp-idns) &
                    -f3(inp)+f3(inp-idtb)
          enddo
        enddo
      enddo


!!  "If you have a pressure equations with boundaries that do not fix pressure level, you have to fix a reference pressure." H.Jasak cfd-online forum
!// In incompressible flow, only relative pressure matters.  Unless there is a pressure BC present,
!// one cell's pressure has to be set to produce a unique pressure solution
!     pEqn.setReference(pRefCell, pRefValue);
!//
      inpr=lk(kpr)+li(ipr)+jpr

      ap(inpr)=1.0d0
      ae(inpr)=0.0d0
      aw(inpr)=0.0d0
      an(inpr)=0.0d0
      as(inpr)=0.0d0
      at(inpr)=0.0d0
      ab(inpr)=0.0d0

!.....Reference pressure
      ppref = p(inpr)
      su(inpr) = ppref
                                                                                                                                                 
                                                                                                                                                                                 
      ! Initialize
      ! p = 0.0d0

!.....Solving pressure correction equation                                               
      ! call sipsol(p,ip)                                                   
      ! call cgstab(p,ip)                                                   
      call iccg(p,ip)                                                     
      ! call pcgsip(p,ip)                                                   
      ! call cgstab_sip(p,ip)                                               
                                                                                      
!.....Calculate pressure gradients                                                       
      call bpres(p)                                                                      
      if (lstsq) then   
        call grad_lsq_qr(p,gradp,2)
      elseif (gauss) then
        call grad_gauss(p,gradp(1,:),gradp(2,:),gradp(3,:))
      endif                                                                      


!.....Make divergence free velocity field
      do k=3,nkmm
        do i=3,nimm
          do j=3,njmm

            inp=lk(k)+li(i)+j     

            u(inp) = u(inp) - dt/den(inp) * (gradP(1,inp)*vol(inp))
            v(inp) = v(inp) - dt/den(inp) * (gradP(2,inp)*vol(inp))
            w(inp) = w(inp) - dt/den(inp) * (gradP(3,inp)*vol(inp))

          enddo
        enddo
      enddo

      ! ! For calcuvw subroutine - copy pressure field to pp
      ! pp = p

!.....Explicit correction of boundary conditions 
      !call correctBoundaryConditions
