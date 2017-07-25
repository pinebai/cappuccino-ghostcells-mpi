!***********************************************************************
!
      SUBROUTINE PIMPLE_multiple_correction
!
!***********************************************************************
!
      use types
      use parameters
      use indexes 
      use geometry
      use coef
      use variables
      use buoy
      use time_mod
      use title_mod
      use gradients
      use hcoef
      use fieldManipulation

      implicit none
!
!***********************************************************************
!
      integer :: i, j, k, inp, &
                 ine, inn, int, &
                 idew, idns, idtb,inpr
      real(prec) :: ppref,fmcor

!.....Before entering the corection loop backup a_nb coefficient arrays:
      he=ae
      hw=aw
      hn=an
      hs=as
      ht=at
      hb=ab     

!+++++PISO Corrector loop++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO icorr=1,ncorr

!// From the last solution of velocity, extract the diag. term from the matrix and store the reciprocal
!// note that the matrix coefficients are functions of U due to the non-linearity of convection.
!            volScalarField rUA = 1.0/UEqn.A();
!// take a Jacobi pass and update U.  See Hrv Jasak's thesis eqn. 3.137 and Henrik Rusche's thesis, eqn. 2.43
!// UEqn.H is the right-hand side of the UEqn minus the product of (the off-diagonal terms and U).
!// Note that since the pressure gradient is not included in the UEqn. above, 
!// this gives us U without the pressure gradient.  Also note that UEqn.H() is a function of U.
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Posle ovoga imamo novo H(u)/ap, H(v)/ap ,i H(w)/ap A.K.A. "HbyA" smesteno u U,V, i W. To je polje brzina 
! bez uticaja gradijenta pritiska!
!
      call PISO_getHbyA 
 

!.....ASSEMBLE OFF DIAGONAL ENTRIES OF Pressure Eq. SYSTEM MATRIX,
!// calculate the fluxes F1,F2,F3 by dotting the interpolated velocity (to cell faces) with face normals
!// {The ddtPhiCorr term accounts for the divergence of the face velocity field by taking out the 
!//  difference between the interpolated velocity and the flux.} NOT ACCOUNTED FOR HERE.
!.....ONLY INNER SURFACES
!.....EAST CELL - FACE
      CALL PISO_assemble_pressure_eq(nj,1,nij, &
                                     ar1x,ar1y,ar1z, &
                                     fx,ae,aw,f1)
!.....NORTH CELL - FACE
      CALL PISO_assemble_pressure_eq(1,nij,nj, &
                                     ar2x,ar2y,ar2z, &
                                     fy,an,as,f2)
!.....TOP   CELL - FACE
      CALL PISO_assemble_pressure_eq(nij,nj,1, &
                                     ar3x,ar3y,ar3z, &
                                     fz,at,ab,f3)

!// adjusts the inlet and outlet fluxes to obey continuity, which is necessary for creating a well-posed
!// problem where a solution for pressure exists.
!     adjustPhi(phi, U, p);
      if(.not.const_mflux) call outbc

!.....ASSEMBLE THE SYSTEM's MAIN DIAGONAL ENTRY AND RHS VECTOR STORED IN SU ARRAY
      idew=nj
      idns=1
      idtb=nij

      do k=3,nkmm
        do i=3,nimm
          do j=3,njmm

            inp=lk(k)+li(i)+j

            su(inp) = -f1(inp)+f1(inp-idew) &
                      -f2(inp)+f2(inp-idns) &
                      -f3(inp)+f3(inp-idtb)

            ap(inp) = ae(inp)+aw(inp)+an(inp)+as(inp)+at(inp)+ab(inp)

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
      ppref=pp(inpr)
      su(inpr)=ppref


!=====Multiple pressure corrections======================================================.
      DO ipcorr=1,npcor                                                                  
!                                                                                        
!.....Initialize pressure                                                                
      pp = 0.0d0                                                                         
                                                                                         
!                                                                                        
!.....SOLVING PRESSURE CORRECTION EQUATION                                               
      ! call sipsol(pp,ip)                                                              
      ! call cgstab(pp,ip)                                                               
      ! call iccg(pp,ip)                                                                 
      call pcgsip(pp,ip)                                                                 
      ! call cgstab_sip(pp,ip)                                                           
!                                                                                        
!.....Calculate pressure gradients                                                       
      call bpres(pp)
      if (lstsq) then   
        call grad_lsq_qr(pp,gradp,2)
      elseif (gauss) then
        call grad_gauss(pp,gradp(1,:),gradp(2,:),gradp(3,:))
      endif
                                                                                         
                                                                                         
!// Correct the flux using the most up-to-date pressure.                                 
!// The .flux method includes contributions from all implicit terms of the               
!// pEqn (the Laplacian):                                                                
!                    phi -= pEqn.flux();                                                 
!//                                                                                      
                                                                                         
!.....Correct mass fluxes to get conservative ones:       
      idew=nj                                                                          
      idns=1                                                                            
      idtb=nij                                                                          
                                                                                        
      do k=2,nkmm                                                                       
        do i=2,nimm                                                                       
          do j=2,njmm                                                                       
                                                                                           
            inp=lk(k)+li(i)+j                                                                 
                                                                                              
            ine=inp+idew                                                                      
            inn=inp+idns                                                                     
            int=inp+idtb     

    !.....Correct mass fluxes at inner cv-faces only (only inner flux)                      
            f1(inp)=f1(inp)-ae(inp)*(pp(ine)-pp(inp))                                         
            f2(inp)=f2(inp)-an(inp)*(pp(inn)-pp(inp))                                         
            f3(inp)=f3(inp)-at(inp)*(pp(int)-pp(inp))                                         
                                                                                              
          enddo                                                                             
        enddo                                                                             
      enddo                                                                             

!.....Correct velocities and pressures                                                                                                                                   
      do k=3,nkmm                                                                        
        do i=3,nimm                                                                        
          do j=3,njmm                                                                                                                            
                                                                                        
!// Add pressure gradient to interior velocity and BC's.                                
!// Note that this pressure is not just a small correction                              
!// to a previous pressure, but is the entire pressure field.                           
!// Contrast this to the use of p' in Ferziger & Peric, Eqn. 7.37.                      
!//                                                                                          
            u(inp)=u(inp)-apu(inp)*gradP(1,inp)*vol(inp)                                      
            v(inp)=v(inp)-apv(inp)*gradP(2,inp)*vol(inp)                                     
            w(inp)=w(inp)-apw(inp)*gradP(3,inp)*vol(inp)                                      
            p(inp)=p(inp)+urf(ip)*(pp(inp)-p(inp))                                            
          enddo                                                                             
        enddo                                                                             
      enddo                                                                             

!.....Update pressure at ghost cells
      call bpres(p)  
                                                                                         
      ! Clean rhs vector                                                                 
      su = 0.0_dp                                                                        
                                                                                                   
!.....Source term modification for the ipcorr-th corrector                              
      if(ipcorr.ne.npcor) then !.........................................               
                                                                        !
      idew=nj                                                           !                
      idns=1                                                            !                
      idtb=nij                                                          !                
                                                                        !
      do k=2,nkmm                                                       !                
      do i=2,nimm                                                       !                
      do j=2,njmm                                                       !                
                                                                        !
      inp=lk(k)+li(i)+j                                                 !                
                                                                        !
      ine=inp+idew                                                      !                
      inn=inp+idns                                                      !                
      int=inp+idtb                                                      !                
!                                                                       !                
!.....correct mass fluxes at inner cv-faces for second corr.            !                
!.....east cell - face                                                  !                
      call fluxmc(inp,nj,1,nij,fx,fmcor)                                !                
        f1(inp) = f1(inp)-fmcor                                         !
        su(inp) = su(inp)+fmcor                                         !                
        su(ine) = su(ine)-fmcor                                         !                
!.....north cell - face                                                 !                
      call fluxmc(inp,1,nij,nj,fy,fmcor)                                !                
        f2(inp) = f2(inp)-fmcor                                         !                
        su(inp) = su(inp)+fmcor                                         !                
        su(inn) = su(inn)-fmcor                                         !                
!.....top   cell - face                                                 !                
      call fluxmc(inp,nij,nj,1,fz,fmcor)                                !                
        f3(inp) = f3(inp)-fmcor                                         !                
        su(inp) = su(inp)+fmcor                                         !                
        su(int) = su(int)-fmcor                                         !                
!.....end: source term modification                                     !                
      enddo                                                             !                
      enddo                                                             !                
      enddo                                                             !                
                                                                        !                
      endif!............................................................!                
                                                                                         
                                                                                         
! We have hit the last iteration of nonorthogonality correction:                         
      if(ipcorr.eq.npcor.and.npcor.gt.1) then !..........................                
      idew=nj                                                           !                
      idns=1                                                            !                
      idtb=nij                                                          !                
      do k=2,nkmm                                                       !                
      do i=2,nimm                                                       !                
      do j=2,njmm                                                       !                
      inp=lk(k)+li(i)+j                                                 !                
      ine=inp+idew                                                      !                
      inn=inp+idns                                                      !                
      int=inp+idtb                                                      !                
!                                                                       !                
!.....correct mass fluxes at inner cv-faces for second corr.            !                
!.....east cell - face                                                  !                
      call fluxmc(inp,nj,1,nij,fx,fmcor)                                !                
        f1(inp) = f1(inp)-fmcor                                         !                        
!.....north cell - face                                                 !                
      call fluxmc(inp,1,nij,nj,fy,fmcor)                                !                
        f2(inp) = f2(inp)-fmcor                                         !                
!.....top   cell - face                                                 !                
      call fluxmc(inp,nij,nj,1,fmcor)                                   !                
        f3(inp) = f3(inp)-fmcor                                         !                
      enddo                                                             !                
      enddo                                                             !                
      enddo                                                             !                
                                                                        !                
      endif!............................................................!                 
                                                                                         
!.....Explicit correction of boundary conditions                                         
      !call correctBoundaryConditions                                                    
                                                                                         
!=====END:Multiple pressure corrections==================================================
      enddo

!.....Write continuity error report:
      include 'continuityErrors.h'

!+++++PISO Corrector loop++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      enddo

      ! For calcuvw subroutine - copy pressure field to pp
      pp = p

      return
      end
