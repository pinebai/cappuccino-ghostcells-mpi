!***********************************************************************
!
      SUBROUTINE adaptive_wallbc_omega(inhlp,idew)
!
!***********************************************************************
!
!     Wall BC for Omega [1/s] - Automatic wall treatment
!
!***********************************************************************
!
      USE TYPES
      USE PARAMETERS
      USE INDEXES
      USE GEOMETRY
      USE COEF
      USE VARIABLES
      USE BC
      USE BOUNDC

      IMPLICIT NONE
!
!***********************************************************************
!

      INTEGER :: INHLP,IDEW,INP
      REAL(PREC) :: FAC,DXF,DYF,DZF,ARE,ARER,ALF,BET,GAM,DELNR
      REAL(PREC) :: Utot2,Un,UN2,Utan,UtauVis,Ustar,Wvis,Wlog

!---------------------------------------------------------
!.....CALCULATE INDEX
!---------------------------------------------------------
      INP=INHLP-IDEW
!---------------------------------------------------------
!.....FACTOR
!---------------------------------------------------------
      FAC=SIGN(1.,FLOAT(-IDEW))
                                       
!.....FIRST VECTOR-DISTANCE FROM GHOST CELL CENTER TO THIS CELL'S CENTER                                                   
      DXF=FAC*(XC(INP)-XC(INHLP))                                
      DYF=FAC*(YC(INP)-YC(INHLP))    
      DZF=FAC*(ZC(INP)-ZC(INHLP))
!.....SECOND X THIRD DEFINE  ALWAYS CF AREA
      ARE=DSQRT(AR1X(INP)**2+AR1Y(INP)**2+AR1Z(INP)**2)
!.....COMPONENTS OF THE UNIT NORMAL VECTOR
      ARER=1./ARE
      ALF=AR1X(INP)*ARER
      BET=AR1Y(INP)*ARER
      GAM=AR1Z(INP)*ARER
!.....NORMAL DISTANCE FROM SCALAR PRODUCT
      DELN=(DXF*ALF+DYF*BET+DZF*GAM)
      DELNR=1./DELN
!.....1.============================================                  
!      TE(INHLP)=DABS(TE(INHLP))                     !
!      TEPR=DSQRT(TE(INHLP))                         !
!      CONST=DENSIT*CMU25*TEPR                       !const=rho*c_mu^0.25*sqrt(k)=rho*u_tau
!      Ystar=CONST*DELN/VISCOS
!      Ystar=max(Ystar,11.63)                        ! Scalable wall function Grotjans&Menter
!.....2.===========================================   
      Utot2=U(inhlp)*U(inhlp) + V(inhlp)*V(inhlp) + W(inhlp)*W(inhlp)  
      Un = ( U(inhlp)*ALF + V(inhlp)*BET + W(inhlp)*GAM )              
      Un2 = Un*Un                                                      
      Utan = sqrt(Utot2 - Un2)                                          

!.....3.===========================================
!-----Adaptive Wall Treatment---------------------------------------- 
      UtauVis=sqrt(Utan*DELNR*VISCOS/DENSIT)
!      UtauLog=Utan*CAPPA/DLOG(ELOG*Ystar) 
     
!      UTAU=sqrt(sqrt(UtauVis**4+UtauLog**4)) 
      Ustar=sqrt(sqrt(UtauVis**4 + sqrt(0.31*TE(inhlp))**4)) 
 
!      YPL=sqrt(UTAU*USTAR)*DELN*DENSIT/VISCOS

!.....Menter&Esch u_tau blending.....................................    
!      UTAU=Utan*sqrt(sqrt( (1./YPL)**4 + (CAPPA/DLOG(ELOG*YPL))**4 ))     


!.....Viscous sub-layer..........................................
!!      IF (YPL.lt.2.5) THEN                                    !
!.....Ovako je originalno kod Wilcox-a, Od Mentera potice faktor 10  
!      SU(INHLP)=GREAT*6.*(VISCOS/DEN(INHLP))                   !
!     &                /(0.075*DELN**2)                         !
!.....Preko law-of-the-wall.Eliminating u_tau from equation:    !
!      SU(INHLP)=GREAT*6.*U(INHLP)**2                           !
!     &               /(0.075*(VISCOS/DEN(INHLP))*YPL**4) ! 
!.....Using y+ in equation:                                     !<--Pravilno za viskozni podsloj
!!!      SU(INHLP)=GREAT*6.*UTAU**2                               !
!!!     &     /(0.075*(VISCOS/DEN(INHLP))*YPL**2)               !
!      SP(INHLP)=GREAT                                          !                   
!...............................................................!

!.....Log-layer.........................................................
!      ELSEIF (YPL.gt.11.63) THEN                                      !
!.....Preko kin.ener.turb. - jer u log sloju znamo vezu u_tau=f(k)     !
!      SU(INHLP)=GREAT*TE(INHLP)*SQRT(BETTAST)                         !
!     &     /(0.3*CAPPA*(VISCOS/DEN(INHLP))*YPL)                       !
!      SP(INHLP)=GREAT                                                 !
!.....Preko law-of-the-wall.                                           !
!      SU(INHLP)=GREAT*UTAU**2                                         ! <--Pravilno za log-sloj umesto 0.3 moze da stoji sqrt(BETTAST)
!     &     /(0.3*CAPPA*(VISCOS/DEN(INHLP))*YPL)                       !
!.....Standard wall function for log-layer:                            ! <--- Odlicno funckionise                                        
!!      SU(INHLP)=GREAT*DSQRT(TE(INP))                                !
!!     &               /(CMU25*CAPPA*DELN)                            !
!!      SP(INHLP)=GREAT                                               !
!......................................................................!

!.....Buffer layer..........................................................
!      ELSE                                                                 !
!.....FIND OMEGA IN FIRST NEAR-WALL CELL BY BLENDING                       !
!.....Preko law-of-the-wall                                                !1 <--Ne koristis UTAU vec U(inhlp)
!      WVIS=6.*U(INHLP)**2                                                 !1
!     &     /(0.075*(VISCOS/DEN(INHLP))*YPL**4)                            !1
!.....Eliminating y+ from the equation:                                    !2 <--Ovako je originalno kod Wilcox-a
      WVIS=6.*(VISCOS/DEN(INHLP))  &                                       !2
           /(0.075*DELN**2)                                                !2
!.....Preko Ustar koje nalazimo blendingom                                 !3
!      WVIS=6.*Ustar**2*DEN(INHLP)                                         !3 <--Ovo je najbolja verzija
!     &     /(0.075*VISCOS*YPL**2)                                         !3
!.....Preko pretpostavke da i ovde vazi relacija k=u_tau**2/betta_star**0.5!4
!c!      WVIS=6.*TE(INHLP)*SQRT(BETTAST)                                   !4
!c!     &     /(0.075*(VISCOS/DEN(INHLP))*YPL**2)                          !4
!--------------------------------------------------------------------------!
!.....Eliminating y+ from the equation:                                    !2
!      WLOG=DSQRT(TE(INP))                                                 !2
!     &     /(CMU25*CAPPA*DELN)                                            !2
!.....Preko UTAU koje nalazimo blendingom                                  !3  <--Pravilno za log-sloj umesto 0.3 moze da stoji sqrt(BETTAST)
      WLOG=Ustar/(0.3*CAPPA*DELN)
!      WLOG=UTAU**2                                                       !3 <--Ovo je najbolja verzija
!     &     /(0.3*CAPPA*(VISCOS/DEN(INHLP))*YPL)                           !3
!.....Preko kin.ener.turb. - jer je u_tau=f(k): k=u_tau**2/betta_star**0.5 !4
!c!      WLOG=TE(INHLP)*SQRT(BETTAST)                                      !4
!c!     &     /(0.3*CAPPA*(VISCOS/DEN(INHLP))*YPL)                         !4
                                                                           !
!.....Blending:                                                            !
      SU(INHLP)=DSQRT(WVIS**2+WLOG**2)                                     !
      sp(inhlp)=1.0d0                                                      !
      ae(inhlp)=0.0d0                                                      !
      aw(inhlp)=0.0d0                                                      !
      an(inhlp)=0.0d0                                                      !
      as(inhlp)=0.0d0                                                      !
      at(inhlp)=0.0d0                                                      !
      ab(inhlp)=0.0d0                                                      !
!      ENDIF                                                               !
!..........................................................................!
!=====END OF AUTOMATIC WALL TREATMENT==================================
      RETURN
      END
