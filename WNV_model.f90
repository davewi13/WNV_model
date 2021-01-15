MODULE define_DDEs

  IMPLICIT NONE

  ! Set number of equations, number of delays and event functions
  INTEGER, PARAMETER :: NEQN=21,NLAGS=7,NEF=1
  ! Set length of the temperature datasets
  INTEGER, PARAMETER :: TEMPNair = 2161
  INTEGER :: K
  DOUBLE PRECISION :: PI=3.1415927D0
  ! Set predation and competition parameters
  DOUBLE PRECISION :: densa=1.03D0,densr=0.021298D0,densh=0.043D0,UPS=19.841D0,SHP=0D0,VOL=500D0,B0=0.0031946D0,B1=0.0046884D0
  ! Set egg raft size
  DOUBLE PRECISION :: MAXEGG=200D0
  ! Set WNV transmission parameters
  DOUBLE PRECISION :: DEATHbird=0.000685D0,DEATHbirdWNV=0.167D0,INF_ARR=120D0
  DOUBLE PRECISION :: RECOVERY=0.25D0,PHH=0.33D0,VERTICAL=0.004D0
  DOUBLE PRECISION :: tranMB=0.88D0,tranBM=0.4D0,CAR=200D0,INOCB=0.5D0
  ! Set time points for temperature data
  double precision, DIMENSION(TEMPNair) :: Xtempair= (/ ((K-1)/2D0*365D0/360D0, K= 1,TEMPNair) /)
  character*12 :: filename1,filename2
  INTEGER :: COUNTER
  DOUBLE PRECISION :: LAT


  
CONTAINS

  SUBROUTINE DDES(T,Y,Z,DY)

	! Set variables for DDE solver
    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(NEQN) :: Y,DY
    DOUBLE PRECISION, DIMENSION(NEQN,NLAGS) :: Z
    INTENT(IN)  :: T,Y,Z
    INTENT(OUT) :: DY
    INTEGER :: K,L
    ! Set values used in extracting temperature data and interpolating between values
    double precision :: ytempp1 = 1D0, ytemppn = 1D0
    double precision, DIMENSION(TEMPNair) :: airtemp, airtemp2

    ! Define variable names
    DOUBLE PRECISION :: TEMPnowair,TEMPGCair,TEMPEIPair,TEMPEair,TEMPLair,TEMPPair,TEMPELair,TEMPLPair,TEMPELPair
    DOUBLE PRECISION :: TEMPnowwater,TEMPEwater,TEMPLwater,TEMPPwater,TEMPELwater,TEMPLPwater,TEMPELPwater
    DOUBLE PRECISION :: BIRTHnow,BIRTHE,BIRTHEL,BIRTHELP,BIRTHbird
    DOUBLE PRECISION :: EGGMATnow,EGGMATE,EGGMATL,EGGMATEL,EGGMATLP,EGGMATELP
    DOUBLE PRECISION :: LARMATnow,LARMATL,LARMATP,LARMATLP
    DOUBLE PRECISION :: REt,RLt,RPt,RSAt,REAt,RIAt,MEt,MLt,MPt,MSAt,MEAt,DEt,DLt,DPt,DSAt,DEAt,DIAt
    DOUBLE PRECISION :: dEdt,dLdt,dPdt,dSAdt,dEAdt,dIAdt,dSEdt,dSLdt,dSPdt,dDEdt,dDLdt,dDELdt,dDPdt,dDLPdt,dDELPdt,dGCdt
    DOUBLE PRECISION :: dSBdt,dIBdt,dRBdt,dEIPdt,dSEAdt
    DOUBLE PRECISION :: E,LAR,PUP,SA,EA,IA,SE,SL,SP,DE,DL,DP,GC,EIP,SEA,NB,SB,IB,RB,NBEIP
    DOUBLE PRECISION :: INOCM,INOCI,BITINGnow,BITINGEIP
    DOUBLE PRECISION :: DEATHeggnow,DEATHeggE,DEATHlarnow,DEATHlarL,DEATHpupnow,DEATHpupP,DEATHadunow,DEATHaduEIP
    DOUBLE PRECISION :: PUPMATnow,PUPMATP
    DOUBLE PRECISION :: DIAPAUSEE,DIAPAUSEEL,DIAPAUSEELP,DIAPAUSEnow,DIAPAUSEEIP
    DOUBLE PRECISION :: PPnow,PPE,PPEL,PPELP,PPEIP
    DOUBLE PRECISION :: GCnow,GCC,GCEIP,EIPnow,EIPdelay

    OPEN (UNIT=COUNTER, FILE='tasminmax_rcp85_land_rcm_uk_12km_15_day_2077_2079_lon'//trim(adjustl(filename1))//'lat'//trim(adjustl(filename2))//'.csv', STATUS="OLD", ACTION="READ")
    READ(COUNTER,*)
    DO L = 1, TEMPNair
      IF(L==1) THEN
        READ(COUNTER,*) LAT
      ELSE
      	READ (COUNTER,*) airtemp(L-1)
      END IF
    END DO
    CLOSE(COUNTER)

	! Interpolate between temperature values using spline
    CALL cubic_spline_air(xtempair, airtemp, TEMPnair, ytempp1, ytemppn, airtemp2)

    ! Give names to responses of system of equations    
    E = Y(1)
    LAR = Y(2)
    PUP = Y(3)
    SA = Y(4)
    EA = Y(5)
    IA = Y(6)
    SE = Y(7)
    SL = Y(8)
    SP = Y(9)
    SEA = Y(10)
    DE = Y(11)
    DL = Y(12)
    DP = Y(13)
    GC = Y(17)
    EIP = Y(18)
    SB = Y(19)
    IB = Y(20)
    RB = Y(21)
    NB = (SB+IB+RB)
    NBEIP = Z(19,7)+Z(20,7)+Z(21,7)

    ! Extract temperature values at various time points
    TEMPnowair = splint(xtempair,airtemp,airtemp2,TEMPnair,T)
    TEMPnowwater = TEMPWATER(TEMPnowair)

	! Define temperature values for time points before t=0.  Air temperature
    ! is the mean daily temperature on the first day of observation.
    IF ((T-DE) .LE. 0) THEN
      TEMPEair = splint(xtempair,airtemp,airtemp2,TEMPnair,0D0)
      TEMPEwater = TEMPWATER(TEMPEair)
    ELSE
      TEMPEair = splint(xtempair,airtemp,airtemp2,TEMPnair,T-DE)
      TEMPEwater = TEMPWATER(TEMPEair)
    END IF
    IF ((T-DL) .LE. 0) THEN
      TEMPLair = splint(xtempair,airtemp,airtemp2,TEMPnair,0D0)
      TEMPLwater = TEMPWATER(TEMPLair)
    ELSE
      TEMPLair = splint(xtempair,airtemp,airtemp2,TEMPnair,T-DL)
      TEMPLwater = TEMPWATER(TEMPLair)
    END IF
    IF ((T-DP) .LE. 0) THEN
      TEMPPair = splint(xtempair,airtemp,airtemp2,TEMPnair,0D0)
      TEMPPwater = TEMPWATER(TEMPPair)
    ELSE
      TEMPPair = splint(xtempair,airtemp,airtemp2,TEMPnair,T-DP)
      TEMPPwater = TEMPWATER(TEMPPair)
    END IF
    IF ((T-DL-Z(11,4)) .LE. 0) THEN
      TEMPELair = splint(xtempair,airtemp,airtemp2,TEMPnair,0D0)
      TEMPELwater = TEMPWATER(TEMPELair)
    ELSE
      TEMPELair = splint(xtempair,airtemp,airtemp2,TEMPnair,T-DL-Z(11,4))
      TEMPELwater = TEMPWATER(TEMPELair)    
    END IF
    IF ((T-DP-Z(12,6)-Z(11,5)) .LE. 0) THEN
      TEMPELPair = splint(xtempair,airtemp,airtemp2,TEMPnair,0D0)
      TEMPELPwater = TEMPWATER(TEMPELPair)
    ELSE
      TEMPELPair = splint(xtempair,airtemp,airtemp2,TEMPnair,T-DP-Z(12,6)-Z(11,5))
      TEMPELPwater = TEMPWATER(TEMPELPair)
    END IF
    IF ((T-DP-Z(12,6)) .LE. 0) THEN
      TEMPLPair = splint(xtempair,airtemp,airtemp2,TEMPnair,0D0)
      TEMPLPwater = TEMPWATER(TEMPLPair)
    ELSE
      TEMPLPair = splint(xtempair,airtemp,airtemp2,TEMPnair,T-DP-Z(12,6))
      TEMPLPwater = TEMPWATER(TEMPLPair)
    END IF
    IF ((T-GC) .LE. 0) THEN
      TEMPGCair = splint(xtempair,airtemp,airtemp2,TEMPnair,0D0)
    ELSE
      TEMPGCair = splint(xtempair,airtemp,airtemp2,TEMPnair,T-GC)
    END IF
    IF ((T-EIP) .LE. 0) THEN
      TEMPEIPair = splint(xtempair,airtemp,airtemp2,TEMPnair,0D0)
    ELSE
      TEMPEIPair = splint(xtempair,airtemp,airtemp2,TEMPnair,T-EIP)
    END IF  
    
	! Calculate photoperiod values for a range of time points
    PPnow = DAYLIGHT(T,LAT)
    PPE = DAYLIGHT(T-DE,LAT)
    PPEL = DAYLIGHT(T-DL-Z(11,4),LAT)
    PPELP = DAYLIGHT(T-DP-Z(12,6)-Z(11,5),LAT)
    PPEIP = DAYLIGHT(T-EIP,LAT)

	! Calculate gonotrophic cycle length for various time points
    GCnow = GONOTROPHIC(TEMPnowair)
    GCC = GONOTROPHIC(TEMPGCair)
    GCEIP = GONOTROPHIC(TEMPEIPair)

    EIPnow = EXTRINSIC_INCUBATION(TEMPnowair)
    EIPdelay = EXTRINSIC_INCUBATION(TEMPEIPair)

	! Calculate diapause percentage for given photoperiod values.
    ! Use these diapause percentages to caluclate birth rates.
    ! Do this for both spring and autumn diapause thresholds.
    IF (DAYLIGHT(T,LAT) > DAYLIGHT(T-1,LAT)) THEN
      DIAPAUSEnow = DIAPAUSE_SPRING(PPnow)
      DIAPAUSEE = DIAPAUSE_SPRING(PPE)
      DIAPAUSEEL = DIAPAUSE_SPRING(PPEL)
      DIAPAUSEELP = DIAPAUSE_SPRING(PPELP)
      DIAPAUSEEIP = DIAPAUSE_SPRING(PPEIP)
      BIRTHnow = BIRTH(DIAPAUSEnow,GC)
      BIRTHE = BIRTH(DIAPAUSEE,Z(17,1))
      BIRTHEL = BIRTH(DIAPAUSEEL,Z(17,2))
      BIRTHELP = BIRTH(DIAPAUSEELP,Z(17,3))
    ELSE
      DIAPAUSEnow = DIAPAUSE_AUTUMN(PPnow)
      DIAPAUSEE = DIAPAUSE_AUTUMN(PPE)
      DIAPAUSEEL = DIAPAUSE_AUTUMN(PPEL)
      DIAPAUSEELP = DIAPAUSE_AUTUMN(PPELP)
      DIAPAUSEEIP = DIAPAUSE_AUTUMN(PPEIP)
      BIRTHnow = BIRTH(DIAPAUSEnow,GC)
      BIRTHE = BIRTH(DIAPAUSEE,Z(17,1))
      BIRTHEL = BIRTH(DIAPAUSEEL,Z(17,2))
      BIRTHELP = BIRTH(DIAPAUSEELP,Z(17,3))
    END IF

	! Calculate death rates for each life stage
	DEATHeggnow = DEATHegg(TEMPnowwater)
    DEATHeggE = DEATHegg(TEMPEwater)

    DEATHlarnow = DEATHlar(TEMPnowwater)
    DEATHlarL = DEATHlar(TEMPLwater)

    DEATHpupnow = DEATHpup(TEMPnowwater)
    DEATHpupP = DEATHpup(TEMPPwater)
    
    DEATHadunow = DEATHadu(TEMPnowair,GC,T)
    DEATHaduEIP = DEATHadu(TEMPEIPair,Z(17,7),T-EIP)

    BIRTHbird = BIRD_BIRTH_FUNC(T)
    
	! Calculate development rates for the immature stages
    LARMATnow = LARMATURATION(TEMPnowwater)
    LARMATL = LARMATURATION(TEMPLwater)
    LARMATP = LARMATURATION(TEMPPwater)
    LARMATLP = LARMATURATION(TEMPLPwater)

    EGGMATnow = EGGMATURATION(TEMPnowwater)
    EGGMATE = EGGMATURATION(TEMPEwater)
    EGGMATL = EGGMATURATION(TEMPLwater)
    EGGMATEL = EGGMATURATION(TEMPELwater)
	EGGMATLP = EGGMATURATION(TEMPLPwater)
	EGGMATELP = EGGMATURATION(TEMPELPwater)

    PUPMATnow = PUPMATURATION(TEMPnowwater)
    PUPMATP = PUPMATURATION(TEMPPwater)

    BITINGnow = 0.5D0*DIAPAUSEnow*GONOTROPHIC(TEMPnowair)
    BITINGEIP = 0.5D0*DIAPAUSEEIP*GONOTROPHIC(TEMPEIPair)

	! Innoculate the system with a given number of adults.
    INOCM = INOCCULATEM(T)

	! Delay differential equations for immature stage durations
    dDEdt = 1D0 - EGGMATnow/EGGMATE
    dDLdt = 1D0 - LARMATnow/LARMATL
    dDPdt = 1D0 - PUPMATnow/PUPMATP

    ! Delay differential equation for gonotrophic cycle duration 
    dGCdt = 1D0 - GCnow/GCC
    dEIPdt = 1D0 - EIPnow/EIPdelay

    ! Delay differential equations for stage durations referenced back through previous stages
    dDELdt = (1D0 - dDLdt) * (1D0 - EGGMATL/EGGMATEL)
    dDLPdt = (1D0 - dDPdt) * (1D0 - LARMATP/LARMATLP)
    dDELPdt = (1D0 - dDPdt - dDLPdt) * (1D0 - EGGMATLP/EGGMATELP)
          
	! Recruitment equations for each stage
	REt = BIRTHnow * (SA + EA + IA)
	RLt = BIRTHE * (Z(4,1) +Z(5,1) +Z(6,1)) * SE * EGGMATnow/EGGMATE
    RPt = BIRTHEL * (Z(4,2) +Z(5,2) +Z(6,2)) * Z(7,4) * SL * LARMATnow/LARMATL * EGGMATL/EGGMATEL
    
    RSAt = BIRTHELP * (Z(4,3) +Z(5,3) +Z(6,3) * (1D0-VERTICAL)) * Z(7,5) * Z(8,6) * SP * PUPMATnow/PUPMATP *&
    LARMATP/LARMATLP * EGGMATLP/EGGMATELP + INOCM
    
    REAt = BITINGnow * SA * (tranBM * IB) / NB
    
    RIAt = BITINGEIP * Z(4,7) * SEA * EIPnow/EIPdelay * (tranBM * Z(20,7)) /&
    NBEIP + BIRTHELP * Z(6,3) * VERTICAL * Z(7,5) * Z(8,6) * SP * PUPMATnow/PUPMATP *&
    LARMATP/LARMATLP * EGGMATLP/EGGMATELP

    
	! Maturation equations for each stage
    MEt = RLt
    MLt = RPt
    MPt = BIRTHELP * Z(4,3) * Z(7,5) * Z(8,6) * SP * PUPMATnow/PUPMATP * LARMATP/LARMATLP * EGGMATLP/EGGMATELP
    MSAt = REAt
    MEAt = BITINGEIP * Z(4,7) * SEA * EIPnow/EIPdelay * (tranBM * Z(20,7)) / NBEIP

	! Death equations for each stage
	DEt = DEATHeggnow * E
	DLt = (densa*densr*(((1D0+COS(2D0*PI*(T-182.5D0-UPS)/365D0))/2D0)**SHP)*(LAR/(VOL/5D0))/(1D0+densh*densa*LAR/(VOL/5D0))&
    + DEATHlarnow + B0*EXP(B1*LAR/VOL))*LAR
    DPt = DEATHpupnow * PUP
    DSAt = DEATHadunow * SA
    DEAt = DEATHadunow * EA
    DIAt = DEATHadunow * IA

    ! Balance equations for each life stage
	dEdt = REt - MEt - DEt
	dLdt = RLt - MLt - DLt
    dPdt = RPt - MPt - DPt
    
	dSAdt = RSAt - MSAt - DSAt
    dEAdt = REAt - MEAt - DEAt
    dIAdt = RIAt - DIAt
    
    dSBdt = BIRTHbird * NB - BITINGnow * tranMB * IA * SB / NB - PHH * IB * SB / NB - (DEATHbird + BIRTHbird*NB/CAR) * SB
    
    dIBdt = BITINGnow*tranMB*IA*SB/NB + PHH*IB*SB/NB - (DEATHbird+DEATHbirdWNV+BIRTHbird*NB/CAR)*IB - RECOVERY*IB
    
    dRBdt = RECOVERY * IB - (DEATHbird + BIRTHbird*NB/CAR) * RB

	! Survival equations for the immature stages
    dSEdt = SE * ((EGGMATnow * DEATHeggE / EGGMATE) - DEATHeggnow)
	dSLdt =SL*(((densa*densr*(((1D0+COS(2D0*PI*(T-DL-182.5D0-UPS)/365D0))/2D0)**SHP)*Z(2,4)/(VOL/5D0)/(1D0+densa*densh*Z(2,4)&
    /(VOL/5D0)))+DEATHlarL + B0*EXP(B1*Z(2,4)/VOL)) * (1-dDLdt) - (densa*densr*(((1D0+COS(2D0*PI*(T-182.5D0-UPS)/365D0))/2D0)**SHP)&
    *LAR/(VOL/5D0) / (1D0+densa*densh*LAR/(VOL/5D0))) - DEATHlarnow - B0*EXP(B1*LAR/VOL))
    dSPdt = SP * ((PUPMATnow * DEATHpupP / PUPMATP) - DEATHpupnow)
    dSEAdt = SEA * ((EIPnow * DEATHaduEIP / EIPdelay) - DEATHadunow)

    ! Derivatives for the integrator:
    DY = (/ dEdt, dLdt, dPdt, dSAdt, dEAdt, dIAdt, dSEdt, dSLdt, dSPdt, dSEAdt, dDEdt, dDLdt, dDPdt, dDELdt, dDLPdt, dDELPdt,&
    dGCdt, dEIPdt, dSBdt, dIBdt, dRBdt /)

    RETURN
    END SUBROUTINE DDES

  SUBROUTINE BETA(T,Y,BVAL)
    
    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(NEQN) :: Y
    DOUBLE PRECISION, DIMENSION(NLAGS) :: BVAL
    INTENT(IN)  :: T,Y
    INTENT(OUT) :: BVAL

    ! Set the delay values
	! T - Eggdelay(T)
	BVAL(1) = T-Y(11)
    ! T - Lardelay(T) - Eggdelay(T-Lardelay(T))				
    BVAL(2) = T-Y(12)-Y(14)
    !T - Pupdelay(T) - Lardelay(T-Pupdelay(T)) - Eggdelay(T-Pupdelay(T)-Lardelay(T-Pupdelay(T)))
    BVAL(3) = T-Y(13)-Y(15)-Y(16)	
    !T - Lardelay(T)
    BVAL(4) = T-Y(12)				
    !T - Pupdelay(T) - Lardelay(T-Pupdelay(T))
    BVAL(5) = T-Y(13)-Y(15)			
    !T - Pupdelay(T)
    BVAL(6) = T-Y(13)				
    !T - EIPdelay(T)
    BVAL(7) = T-Y(18)
   	  
    RETURN
  END SUBROUTINE BETA

  SUBROUTINE HISTORY(T,Y)
    DOUBLE PRECISION :: T,TEMPhist,TEMPhistair
    DOUBLE PRECISION, DIMENSION(2*NEQN) :: Y
    INTENT(IN)  :: T
    INTENT(OUT) :: Y
    INTEGER :: L
    ! Set values used in extracting temperature data and interpolating between values
    double precision :: ytempp1 = 1D0, ytemppn = 1D0
    double precision, DIMENSION(TEMPNair) :: airtemp, airtemp2

    OPEN (UNIT=COUNTER, FILE='tasminmax_rcp85_land_rcm_uk_12km_15_day_2077_2079_lon'//trim(adjustl(filename1))//'lat'//trim(adjustl(filename2))//'.csv', STATUS="OLD", ACTION="READ")
    READ(COUNTER,*)
    DO L = 1, TEMPNair
      IF(L==1) THEN
        READ(COUNTER,*) LAT
      ELSE
      	READ (COUNTER,*) airtemp(L-1)
      END IF
    END DO
    CLOSE(COUNTER)

	! Interpolate between temperature values using spline
    CALL cubic_spline_air(xtempair, airtemp, TEMPnair, ytempp1, ytemppn, airtemp2)

	! Set the temperatures for T < 0
    TEMPhist = TEMPWATER(splint(xtempair,airtemp,airtemp2,TEMPnair,0D0))
    TEMPhistair = splint(xtempair,airtemp,airtemp2,TEMPnair,0D0)

	! Set historical values for all stages to be zero
    Y(1) = 0D0
    Y(2) = 0D0
    Y(3) = 0D0
    Y(4) = 0D0
    Y(5) = 0D0
    Y(6) = 0D0

    ! Calculate historical survival rates based on temperature
    Y(7) = EXP(-DEATHegg(TEMPhist)*(1D0/EGGMATURATION(TEMPhist)))
    Y(8) = EXP(-DEATHlar(TEMPhist)*(1D0/LARMATURATION(TEMPhist)))
    Y(9) = EXP(-DEATHpup(TEMPhist)*(1D0/PUPMATURATION(TEMPhist)))

	! Calculate historical development rates based on temperature
    Y(11) = 1D0/EGGMATURATION(TEMPhist)
    Y(12) = 1D0/LARMATURATION(TEMPhist)
    Y(13) = 1D0/PUPMATURATION(TEMPhist)
    Y(14) = 1D0/EGGMATURATION(TEMPhist)
    Y(15) = 1D0/LARMATURATION(TEMPhist)
    Y(16) = 1D0/EGGMATURATION(TEMPhist)
    Y(17) = 1D0/GONOTROPHIC(TEMPhistair)

    Y(10) = EXP(-DEATHadu(TEMPhistair,Y(17),T)*(1D0/EXTRINSIC_INCUBATION(TEMPhistair)))

    ! Set historical EIP
    Y(18) = 1D0/EXTRINSIC_INCUBATION(TEMPhistair)
    
	! Set historical values for birds to zero
	Y(19) = 0.875D0*CAR      
    Y(20) = 0D0
    Y(21) = 0D0
    
    RETURN
  END SUBROUTINE HISTORY


  SUBROUTINE EF(T,Y,DY,Z,G)

    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(NEQN) :: Y,DY
    DOUBLE PRECISION, DIMENSION(NEQN,NLAGS) :: Z
    DOUBLE PRECISION, DIMENSION(NEF) :: G
    INTENT(IN) :: T,Y,DY,Z
    INTENT(OUT) :: G

    ! Set events as turning points in adult time series and diapause entry/exit
    G(1) = T-730D0-INF_ARR

     RETURN
  END SUBROUTINE EF

  SUBROUTINE CHNG(NEVENT,TEVENT,YEVENT,DYEVENT,HINIT, &
                  DIRECTION,ISTERMINAL,QUIT)
  ! Function to change a flag so that the DDE model will
  ! be evaluated in subroutine DDES instead of the ODE model.

     INTEGER :: NEVENT
     INTEGER, DIMENSION(NEF) :: DIRECTION
     DOUBLE PRECISION :: TEVENT,HINIT
     DOUBLE PRECISION, DIMENSION(NEQN) :: YEVENT,DYEVENT
     LOGICAL :: QUIT
     LOGICAL, DIMENSION(NEF) :: ISTERMINAL
     INTENT(IN) :: NEVENT,TEVENT
     INTENT(INOUT) :: YEVENT,DYEVENT,HINIT,DIRECTION,&
                      ISTERMINAL,QUIT

     YEVENT(20) = 0.0025D0*CAR

    RETURN
  END SUBROUTINE CHNG

   SUBROUTINE cubic_spline_air(x,airtemp,n,airtempp1,airtemppn,airtemp2)
    INTEGER n,NMAX
    double precision airtempp1,airtemppn,x(n),airtemp(n),airtemp2(n)
    PARAMETER (NMAX=26304)
    
    !Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., yi = f(xi), with
    !x1 < x2 < ::: < xN, and given values yp1 and ypn for the rst derivative of the interpolating
    !function at points 1 and n, respectively, this routine returns an array y2(1:n) of
    !length n which contains the second derivatives of the interpolating function at the tabulated
    !points xi. If yp1 and/or ypn are equal to 1  1030 or larger, the routine is signaled to set
    !the corresponding boundary condition for a natural spline, with zero second derivative on
    !that boundary.
    
    !Parameter: NMAX is the largest anticipated value of n.
    INTEGER i,k
    double precision p,qn,sig,un,u(NMAX)
    if (airtempp1.gt..99e30) then !The lower boundary condition is set either to be
    airtemp2(1)=0. !\natural"
    u(1)=0.
    else !or else to have a specified first derivative.
    airtemp2(1)=-0.5
    u(1)=(3./(x(2)-x(1)))*((airtemp(2)-airtemp(1))/(x(2)-x(1))-airtempp1)
    endif
    do i=2,n-1 
    !This is the decomposition loop of the tridiagonal
    !algorithm. y2 and u are used for temporary
    !storage of the decomposed factors.
    sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
    p=sig*airtemp2(i-1)+2.
    airtemp2(i)=(sig-1.)/p
    u(i)=(6.*((airtemp(i+1)-airtemp(i))/(x(i+1)-x(i))-(airtemp(i)-airtemp(i-1)) &
    /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    enddo
    if (airtemppn.gt..99e30) then !The upper boundary condition is set either to be
    qn=0. !\natural"
    un=0.
    else !or else to have a specified first derivative.
    qn=0.5
    un=(3./(x(n)-x(n-1)))*(airtemppn-(airtemp(n)-airtemp(n-1))/(x(n)-x(n-1)))
    endif
    airtemp2(n)=(un-qn*u(n-1))/(qn*airtemp2(n-1)+1.)
    do k=n-1,1,-1 !This is the backsubstitution loop of the tridiagonal algorithm.
    airtemp2(k)=airtemp2(k)*airtemp2(k+1)+u(k) 
    enddo
    return

   END SUBROUTINE

   double precision FUNCTION SPLINT(xa,ytempa,ytemp2a,tempn,x)
    INTEGER tempn
    double precision xa(tempn),ytemp2a(tempn),ytempa(tempn)
    double precision x
    !Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the
    !xai 's in order), and given the array y2a(1:n), which is the output from spline above,
    !and given a value of x, this routine returns a cubic-spline interpolated value y.
    INTEGER k,khi,klo
    double precision a1,b,h1

    !We will find the right place in the table by means of bisection.
    !This is optimal if sequential calls to this routine are at random
    !values of x. If sequential calls are in order, and closely
    !spaced, one would do better to store previous values of
    !klo and khi and test if they remain appropriate on the
    !next call.
    
    klo=1 
    khi=tempn
    do
    if (khi-klo.le.1) exit
    k=(khi+klo)/2
    if(xa(k).gt.x) then
        khi=k
    else
        klo=k
    endif
    enddo
    !klo and khi now bracket the input value of x.
    
    h1=xa(khi)-xa(klo)
    if (h1.eq.0) pause 'bad xa input in splint' !The xa's must be distinct.
    a1=(xa(khi)-x)/h1 !Cubic spline polynomial is now evaluated.
    b=(x-xa(klo))/h1
    if (x .le. 0) then
      splint=0.1d0
      else
    splint=a1*ytempa(klo)+b*ytempa(khi)+ ((a1**3-a1)*ytemp2a(klo)+(b**3-b)*ytemp2a(khi))*(h1**2)/6.
    end if
    return
  END FUNCTION


  DOUBLE PRECISION FUNCTION TEMPWATER(AIRTEMP)
     
    DOUBLE PRECISION :: AIRTEMP
    
    TEMPWATER = 0.9505D0*AIRTEMP + 3.8887D0
    
  RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION INOCCULATEM(T)

    DOUBLE PRECISION :: T

    IF (T < 1D0 .AND. T > 0D0) THEN
	   INOCCULATEM = 20000D0
    ELSE
       INOCCULATEM = 0D0
    END IF

	RETURN
  END FUNCTION
  
  DOUBLE PRECISION FUNCTION DAYLIGHT(T,LAT)
     
    DOUBLE PRECISION :: T,EPS,NUM,DEN,TIME360,LAT
    REAL, PARAMETER :: Pi = 3.1415927D0
    TIME360 = T
    EPS = ASIN(0.39795D0 * COS(0.2163108D0 + 2 * ATAN(0.9671396D0 * TAN(0.00860D0 * (TIME360-185.5D0)))))
    NUM = SIN(0.8333D0*Pi/180D0) + (SIN(LAT*Pi/180D0) * SIN(EPS))
    DEN = COS(LAT*Pi/180D0) * COS(EPS)
    DAYLIGHT = 24D0 - (24D0/Pi) * ACOS(NUM / DEN)
    RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION DIAPAUSE_SPRING(PP)
     
    DOUBLE PRECISION :: PP
    DIAPAUSE_SPRING = 1D0 / (1D0 + EXP(1.5D0*(13.7D0-PP)))
    RETURN
    
  END FUNCTION
  
  DOUBLE PRECISION FUNCTION DIAPAUSE_AUTUMN(PP)
     
    DOUBLE PRECISION :: PP
    DIAPAUSE_AUTUMN = 1D0 / (1D0 + EXP(1D0*(15D0-PP)))
    RETURN
    
  END FUNCTION

  DOUBLE PRECISION FUNCTION BIRTH(DIAPAUSE,GONOTROPHICtime)
     
    DOUBLE PRECISION :: GONOTROPHICtime,EGGRAFT,DIAPAUSE

    EGGRAFT = DIAPAUSE*MAXEGG
    BIRTH = 0.5D0*EGGRAFT/GONOTROPHICtime
  
    RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION BIRD_BIRTH_FUNC(T)
     
    DOUBLE PRECISION :: T
    DOUBLE PRECISION :: K_BIRD=0.15D0,PHI_BIRD=50D0,S_BIRD=10D0

    BIRD_BIRTH_FUNC = 0.5D0 * K_BIRD * EXP(-S_BIRD * (COS((PI*T+PHI_BIRD)/365D0))**2)
  
    RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION GONOTROPHIC(TEMP)
     
	DOUBLE PRECISION :: TEMP
	DOUBLE PRECISION :: KG=0.2024D0,QG=74.48D0,BG=0.2456D0

	IF (TEMP <= 0D0) THEN
      GONOTROPHIC = 0.0333D0
    ELSE
	  GONOTROPHIC = KG / (1+QG*EXP(-BG*TEMP))
    END IF
	IF(GONOTROPHIC < 0.0333D0) THEN
      GONOTROPHIC = 0.0333D0
    END IF
    
    RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION EXTRINSIC_INCUBATION(TEMP)
     
	DOUBLE PRECISION :: TEMP
    DOUBLE PRECISION :: Tmax = 45.2D0, Tmin = 11.4D0, q=7.38D-5
    
	IF(TEMP > Tmax) THEN
	  EXTRINSIC_INCUBATION = 0.005D0
	ELSE IF(TEMP < Tmin) THEN
	  EXTRINSIC_INCUBATION = 0.005D0
	ELSE
	  EXTRINSIC_INCUBATION = q*TEMP*(TEMP-Tmin)*SQRT(Tmax-TEMP)
    END IF
    IF (EXTRINSIC_INCUBATION .LE. 0.005D0) THEN
      EXTRINSIC_INCUBATION = 0.005D0
    END IF

    RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION DEATHegg(TEMP)
    
	DOUBLE PRECISION :: TEMP
    DOUBLE PRECISION :: U3=0.0157D0,U4=20.5D0,U5=7D0
    
	DEATHegg = U3 * EXP(((TEMP-U4)/U5)**2)
	IF (DEATHegg > 1D0) THEN
      DEATHegg = 1D0
    END IF
        
    RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION DEATHlar(TEMP)
    
	DOUBLE PRECISION :: TEMP
    DOUBLE PRECISION :: U3=0.0157D0,U4=20.5D0,U5=7D0
    
	DEATHlar = U3 * EXP(((TEMP-U4)/U5)**2)
    IF (DEATHlar > 1D0) THEN
      DEATHlar = 1D0
    END IF
        
    RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION DEATHpup(TEMP)
    
	DOUBLE PRECISION :: TEMP
    DOUBLE PRECISION :: U3=0.0157D0,U4=20.5D0,U5=7D0
    
	DEATHpup = U3 * EXP(((TEMP-U4)/U5)**2)
    IF (DEATHpup > 1D0) THEN
      DEATHpup = 1D0
    END IF
        
    RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION DEATHadu(TEMP,GONOTROPHICtime,T)
    
	DOUBLE PRECISION :: TEMP,GONOTROPHICtime,T
    DOUBLE PRECISION :: ALPHA=2.166D-8,BETA=4.483D0,PI=3.1415927D0,MULTIPLIER=8D0,SIGMASQ=4D0
	
	!Threshold set to 80% exit from diapause
	IF (TEMP <= 0D0) THEN
    	DEATHadu = 0.003D0
    ELSE
		DEATHadu = ALPHA*(TEMP**BETA)
    END IF
	IF (DEATHadu < 0.003D0) THEN
      DEATHadu = 0.003D0
    END IF
   	DEATHadu = DEATHadu + (MULTIPLIER/SQRT(SIGMASQ*2D0*PI))*EXP((-1D0/(SIGMASQ*2D0))*&
    (MOD(T,365D0)-GONOTROPHICtime-109D0)**2)
   
    RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION EGGMATURATION(TEMP)
     
	DOUBLE PRECISION :: TEMP
    DOUBLE PRECISION :: ALPHA=0.0022D0,BETA=1.77D0
	
	IF (TEMP <= 0D0) THEN
    	EGGMATURATION = 0.016667D0
    ELSE
		EGGMATURATION = ALPHA*(TEMP**BETA)
	END IF
	IF (EGGMATURATION < 0.016667D0) THEN
    	EGGMATURATION = 0.016667D0
    END IF   
	IF (EGGMATURATION > 1D0) THEN
    	EGGMATURATION = 1D0
    END IF 	

    RETURN

  END FUNCTION

  DOUBLE PRECISION FUNCTION LARMATURATION(TEMP)
     
	DOUBLE PRECISION :: TEMP
    DOUBLE PRECISION :: ALPHA=0.00315D0,BETA=1.12D0
	
	IF (TEMP <= 0D0) THEN
    	LARMATURATION = 0.016667D0
    ELSE
		LARMATURATION = ALPHA*(TEMP**BETA)
	END IF
	IF (LARMATURATION < 0.016667D0) THEN
    	LARMATURATION = 0.016667D0
    END IF 
	IF (LARMATURATION > 1D0) THEN
    	LARMATURATION = 1D0
    END IF 	

    RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION PUPMATURATION(TEMP)
     
	DOUBLE PRECISION :: TEMP
    DOUBLE PRECISION :: ALPHA=0.0007109D0,BETA=1.8865648D0
	
	IF (TEMP <= 0D0) THEN
    	PUPMATURATION = 0.016667D0
    ELSE
		PUPMATURATION = ALPHA*(TEMP**BETA)
	END IF
    IF (PUPMATURATION < 0.016667D0) THEN
    	PUPMATURATION = 0.016667D0
    END IF 
	IF (PUPMATURATION > 1D0) THEN
    	PUPMATURATION = 1D0
    END IF 
    
    RETURN
    
  END FUNCTION

END MODULE define_DDEs

!******************************************************************

PROGRAM runmodel

  USE define_DDEs
  USE DDE_SOLVER_M

  IMPLICIT NONE

  INTEGER :: I,J ! Local variables

  INTEGER, DIMENSION(3) :: NVAR = (/NEQN,NLAGS,NEF/)

  
  INTEGER, PARAMETER :: NOUT=1096D0
  DOUBLE PRECISION, PARAMETER :: T0=0D0,TFINAL=1095D0
  DOUBLE PRECISION, DIMENSION(NOUT) :: TSPAN= &
  (/ (T0+(I-1)*((TFINAL - T0)/(NOUT-1)), I=1,NOUT) /)

  INTEGER :: LONSPAN(112) = (/ (I, I=1,112) /)
  INTEGER :: LATSPAN(112) = (/ (I, I=1,112) /)

  INTEGER :: VEC(12544) = (/ (I, I = 0, 12543) /)
  
  INTEGER, DIMENSION(12544) :: VEC1,VEC2

  TYPE(DDE_SOL) :: SOL
  TYPE(DDE_OPTS) :: OPTS
  
  DOUBLE PRECISION :: MAXDELAY = 200D0
  INTEGER :: ITERSTART,ITEREND
  CHARACTER*10 :: BUFFER,BUFFER2
  INTEGER :: POINTER1, POINTER2

  CALL GETARG(1,BUFFER)
  CALL GETARG(2,BUFFER2)
  READ(BUFFER,*) ITERSTART
  READ(BUFFER2,*) ITEREND

  VEC1 = INT(VEC/112)+1
  VEC2 = MOD(VEC,112)+1

  DO COUNTER=ITERSTART,ITEREND

      POINTER1 = VEC1(COUNTER)
      POINTER2 = VEC2(COUNTER)
      WRITE(FILENAME1, '(I3)') LONSPAN(POINTER1)
      WRITE(FILENAME2, '(I3)') LATSPAN(POINTER2)

      OPTS = DDE_SET(RE=1D-5,AE=1D-5,MAX_STEPS=100000000,MAX_DELAY=MAXDELAY)  

      SOL = DDE_SOLVER(NVAR,DDES,BETA,HISTORY,TSPAN,EVENT_FCN=EF,CHANGE_FCN=CHNG,OPTIONS=OPTS)
  
      IF (SOL%FLAG == 0) THEN

        OPEN(unit=10,FILE='tasminmax_rcp85_land_rcm_uk_12km_15_day_2077_2079_lon'//trim(adjustl(filename1))//'lat'//trim(adjustl(filename2))//'.dat')
        DO I = 1,SOL%NPTS
          WRITE(UNIT=10,FMT='(28E14.5E3)') SOL%T(I),(SOL%Y(I,J),J=1,NEQN)
        END DO
        CLOSE(unit=10)

      ELSE

        PRINT *,' Abnormal return from DDE_SOLVER with FLAG = ',&
        SOL%FLAG

      END IF

  END DO

  STOP
END PROGRAM runmodel