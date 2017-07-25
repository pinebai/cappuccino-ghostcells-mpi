!     Contains:
!     fvm_ddt(phi)
!     fvm_ddt(rho,U)

      SUBROUTINE fvm_ddt(rho,U)
!  
!******************************************************************************
!
!     Fills matrix with coefficients representing implicit FVM discretization 
!     of Nonstationary term: \frac{\partial \rho U}{\partial t}.
!
!     Inteded for structured 3D grids.
!
!     Matrix is stored in diagonal format with seven diagonals.
!     Each diagonal is stored in an array,
!     Off-main diagonnal: ae,aw,an,as,at,ab,
!     Main diagonal: ap.
!     RHS vector is SU.
!     Implicitely treated part of the source is in SP array.
!     System of linear equations is written as:
!     $ a_p^{(i)}*\phi_p^(i)-\sum_{j=1}^{nb} a_j^{(i)}*\phi_j^(i) = b_p{(i)}, i=1,ncells $
!
!******************************************************************************
!
      Use Types
      Use Parameters
      Use Indexes
      Use Geometry
      Use Coef
      Use Variables
      Use Buoy
      Use Time_Mod
      Use Gradients

      Implicit None

      Real(prec), Dimension(nxyza) :: rho, U
!
!     Local Variables
!
      Integer :: I, J, K, Inp

      Do K=2,Nkm
      Do I=2,Nim
      Do J=2,Njm
      Inp=Lk(K)+Li(I)+J

      IF(BDF) THEN
!-----------------------------------------------------------------------
!    Three Level Implicit Time Integration Method:
!    in case that BTIME=0. --> Implicit Euler
!-----------------------------------------------------------------------
      apotime=rho(inp)*vol(inp)/timestep
      sut=apotime*((1+btime)*phio(inp)-0.5*btime*phioo(inp))
      su(inp)=su(inp)+sut
      sp(inp) = sp(inp)+apotime*(1+0.5*btime)
!-----------------------------------------------------------------------
      ENDIF

      IF(CN) THEN
!.....Crank-Nicolson stuff:
      ae(inp)=0.5d0*ae(inp)
      an(inp)=0.5d0*an(inp)
      at(inp)=0.5d0*at(inp)
      aw(inp)=0.5d0*aw(inp)
      as(inp)=0.5d0*as(inp)
      ab(inp)=0.5d0*ab(inp)
!.....Crank-Nicolson time stepping source terms
      apotime=rho(inp)*vol(inp)/timestep
      su(inp)=su(inp)+(ae(inp)*phio(inp+nj)  + aw(inp)*phio(inp-nj)+    &
                       an(inp)*phio(inp+1)   + as(inp)*phio(inp-1)+     &
                       at(inp)*phio(inp+nij) + ab(inp)*phio(inp-nij))+  &
              (apotime-ae(inp)-aw(inp)                              &
                      -an(inp)-as(inp)                              &
                      -at(inp)-ab(inp))*phio(inp)      
      sp(inp)=sp(inp)+apotime
!.....End of Crank-Nicolson time stepping source terms.
!.....End of Crank-Nicolson stuff.
      ENDIF

      enddo
      enddo
      enddo


      Return
      End

      SUBROUTINE fvm_ddt(phi)
!  
!******************************************************************************
!
!     Fills matrix with coefficients representing implicit FVM discretization 
!     of Nonstationary term: \frac{\partial \rho U}{\partial t}.
!
!     Inteded for structured 3D grids.
!
!     Matrix is stored in diagonal format with seven diagonals.
!     Each diagonal is stored in an array,
!     Off-main diagonnal: ae,aw,an,as,at,ab,
!     Main diagonal: ap.
!     RHS vector is SU.
!     Implicitely treated part of the source is in SP array.
!     System of linear equations is written as:
!     $ a_p^{(i)}*\phi_p^(i)-\sum_{j=1}^{nb} a_j^{(i)}*\phi_j^(i) = b_p{(i)}, i=1,ncells $
!
!******************************************************************************
!
      Use Types
      Use Parameters
      Use Indexes
      Use Geometry
      Use Coef
      Use Variables
      Use Buoy
      Use Time_Mod
      Use Gradients

      Implicit None

      Real(prec), Dimension(nxyza) :: phi
!
!     Local Variables
!
      Integer :: I, J, K, Inp

      Do K=2,Nkm
      Do I=2,Nim
      Do J=2,Njm
      Inp=Lk(K)+Li(I)+J

      IF(BDF) THEN
!-----------------------------------------------------------------------
!    Three Level Implicit Time Integration Method:
!    in case that BTIME=0. --> Implicit Euler
!-----------------------------------------------------------------------
      apotime = vol(inp)/timestep
      sut = apotime*((1+btime)*phio(inp)-0.5*btime*phioo(inp))
      su(inp) = su(inp)+sut
      sp(inp) = sp(inp)+apotime*(1+0.5*btime)
!-----------------------------------------------------------------------
      ENDIF

      IF(CN) THEN
!.....Crank-Nicolson stuff:
      ae(inp)=0.5d0*ae(inp)
      an(inp)=0.5d0*an(inp)
      at(inp)=0.5d0*at(inp)
      aw(inp)=0.5d0*aw(inp)
      as(inp)=0.5d0*as(inp)
      ab(inp)=0.5d0*ab(inp)
!.....Crank-Nicolson time stepping source terms
      apotime = vol(inp)/timestep
      su(inp) = su(inp)+(ae(inp)*phio(inp+nj)  + aw(inp)*phio(inp-nj)+    &
                         an(inp)*phio(inp+1)   + as(inp)*phio(inp-1)+     &
                         at(inp)*phio(inp+nij) + ab(inp)*phio(inp-nij))+  &
                (apotime-ae(inp)-aw(inp)                                  &
                        -an(inp)-as(inp)                                  &
                        -at(inp)-ab(inp))*phio(inp)      
      sp(inp) = sp(inp)+apotime
!.....End of Crank-Nicolson time stepping source terms.
!.....End of Crank-Nicolson stuff.
      ENDIF

      enddo
      enddo
      enddo


      Return
      End
