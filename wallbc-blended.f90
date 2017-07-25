!***********************************************************************
!
      SUBROUTINE WALLBC(INHLP,IDEW,IDNS,IDTB)
!
!***********************************************************************
!
!     Law of the wall:
!     U* = 1/cappa * ln(Elog*y*)  in log-layer                     (1)
!     U* = y*                     in viscous sublayer              (2)
!
!     By definition:
!     U* = Up * Cmu**0.25 * kp**0.5 / (Tau_wall / Densit)          (3)
!     y* = Densit * Cmu**0.25 * kp**0.5 * yp / mu                  (4)
!
!     Inserting (3) and (4) into (1) and (2) we get
!     for viscus sublayer:
!     Tau_wall = Up*mu/yp                                          (5)
!     for logarithmic layer:
!     Tau_wall = cappa*Densit*Cmu**0.25*kp**0.5*Up/ln(Elog*y*)     (6)
!
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
      INTEGER :: INHLP,IDEW,IDNS,IDTB,INP
      REAL(PREC) :: FAC,DXF,DYF,DZF,ARE,ARER,ALF,BET,GAM,DELNR
      REAL(PREC) :: TEPR, CONST, Ystar
      REAL(PREC) :: Utot2,Un,UN2,Utan,UtauVis,UtauLog,Ustar
      REAL(PREC) :: ESTAR,TCOEF,TAR
!---------------------------------------------------------
!.....CALCULATE INDEX
!---------------------------------------------------------
      INP=INHLP-IDEW
!---------------------------------------------------------
!.....FACTOR
!---------------------------------------------------------
      FAC=SIGN(1.,FLOAT(-IDEW))
!---------------------------------------------------------
!.....FIRST VECTOR
!---------------------------------------------------------
      DXF=FAC*(XC(INP)-XC(INHLP))   
      DYF=FAC*(YC(INP)-YC(INHLP))    
      DZF=FAC*(ZC(INP)-ZC(INHLP))
!
!---------------------------------------------------------
!.....SECOND X THIRD DEFINE  ALWAYS CF AREA
!.....AREA
!---------------------------------------------------------
      ARE=DSQRT(AR1X(INP)**2+AR1Y(INP)**2+AR1Z(INP)**2)
!---------------------------------------------------------
!.....COMPONENTS OF THE UNIT NORMAL VECTOR
!---------------------------------------------------------
      ARER=1./ARE
      ALF=AR1X(INP)*ARER
      BET=AR1Y(INP)*ARER
      GAM=AR1Z(INP)*ARER
!---------------------------------------------------------
!.....NORMAL DISTANCE FROM SCALAR PRODUCT
!---------------------------------------------------------
      DELN=(DXF*ALF+DYF*BET+DZF*GAM)
      DELNR=1./DELN
!---------------------------------------------------------

!-----Find Ystar--------------------
      TE(INHLP)=DABS(TE(INHLP))                      
      TEPR=DSQRT(TE(INHLP))                          
      CONST=DENSIT*CMU25*TEPR  
      Ystar=CONST*DELN/VISCOS
!      Ystar=max(Ystar,11.63)   ! Scalable wall function Grotjans&Menter
!-----When the grid is fine enough:----------------------------------- 
      Utot2=U(inhlp)*U(inhlp) + V(inhlp)*V(inhlp) + W(inhlp)*W(inhlp)  
      Un = ( U(inhlp)*ALF + V(inhlp)*BET + W(inhlp)*GAM )              
      Un2 = Un*Un                                                      
      Utan = sqrt(Utot2 - Un2)                                       
!-----
      UtauVis=sqrt(Utan*DELNR*VISCOS/DENSIT)
      UtauLog=Utan*CAPPA/DLOG(ELOG*Ystar) 
     
      UTAU=sqrt(sqrt(UtauVis**4+UtauLog**4)) 
      Ustar=sqrt(sqrt( UtauVis**4 + sqrt(0.31*TE(INP))**4))  
      CONST=DENSIT*Ustar
 
      YPL=sqrt(UTAU*USTAR)*DELN*DENSIT/VISCOS
!      Write (*,*) YPL
!---------------------------------------------------------
!     [VISCOUS SUBLAYER: ]
!---------------------------------------------------------
      TCOEF=VISCOS*DELNR                                
!---------------------------------------------------------
!     [TURBULENT: ] SMOOTH WALL
!---------------------------------------------------------
      IF(Ystar.GT.5) THEN 
!.....Blending - Automatic wall treatment:                             
        TCOEF=((VISCOS/DELNR)**4+(CONST*CAPPA/DLOG(ELOG*YPL))**4)**0.25 
      ELSEIF(Ystar.GT.11.63) THEN 
!.....Log-layer-Wall Function:                                                                         
        TCOEF=CONST*CAPPA/DLOG(ELOG*Ystar)                                   
      ENDIF

!==========================================================================
!      [TURBULENT: ] ROUGH WALL        [ADDED BY REMCO]	
!---------------------------------------------------------
      USTAR=CMU25*TEPR
!---------------------------------------------------------
!      ROUGH WALL        [ADDED BY Nikola/scaled top hill]	
!---------------------------------------------------------
!     [scaling factor = 1/120]
!      XC--> 327 m --> 2.725 m
!      ZC--> 0.75 m --> 0.00625  m
!-----Uvesti da se ZZERO cita iz fajla--------------------
!      IF(XC(INP).LT.2.725) THEN
!      ZZERO=3.E-4/11.
!         IF(ZC(INP).GT.0.00625) ZZERO=0.015/11.
!      END IF
!      IF(XC(INP).GE.2.725) THEN
!      ZZERO=0.015/11.
!      END IF
!---------------------------------------------------------
      if (roughWall) then
        ESTAR=USTAR*ZZERO/VISCOS
        IF(YPL.GT.11.63) TCOEF=CONST*CAPPA/DLOG((EROUGH/ESTAR)*YPL)
        TCOEF=CONST*CAPPA/DLOG((EROUGH/ESTAR)*YPL)
      endif
!==========================================================================
      TAR=TCOEF*ARE 
!---------------------------------------------------------
!.....SOURCE TERMS
!---------------------------------------------------------
      SUU=TAR*(ALF*(BET*(V(INHLP)-V(INP))+GAM*(W(INHLP)-W(INP))) &
         +(1.-ALF**2)*U(INP))
      SVU=TAR*(BET*(ALF*(U(INHLP)-U(INP))+GAM*(W(INHLP)-W(INP))) &
         +(1.-BET**2)*V(INP))
      SWU=TAR*(GAM*(ALF*(U(INHLP)-U(INP))+BET*(V(INHLP)-V(INP))) &
         +(1.-BET**2)*W(INP))
      SUP=TAR*(1.-ALF**2)
      SVP=TAR*(1.-BET**2)
      SWP=TAR*(1.-GAM**2)
!---------------------------------------------------------
!.....WALL SHEAR STRESS AND GENER. TERM
!.....PARALEL COMPONENTS
!---------------------------------------------------------
!      UVWP=DSQRT((SUP*U(INHLP)-SUU)**2+    &
!                 (SVP*V(INHLP)-SVU)**2+    &
!                (SWP*W(INHLP)-SWU)**2)/(TAR+SMALL)
!      TAU=UVWP*TCOEF
      TAU = Utan*TCOEF

!.....When first cell lies in viscous sublayer
      GENT = DABS(TAU)*CONST/DENSIT
!.....When first cell lies in log-layer
      IF(Ystar.GT.11.63) THEN 
      GENT=DABS(TAU)*CONST*DELNR/(CAPPA*DENSIT) ! Original - Sasa
      ENDIF

!      SUED=CMU75*DELNR/CAPPA
      SUED = TAU !Trebace mi u writefiles.f

!      write(6,*) 'UVWP: ',UVWP, ' TAR:  ',TAR
!      write(6,*) 'GENT: ',GENT, ' SUED: ',SUED

      RETURN
      END

