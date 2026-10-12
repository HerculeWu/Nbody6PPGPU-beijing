      SUBROUTINE UNITS
*
*
*       Initialization of units & scaling factors.
*       ------------------------------------------
*
      INCLUDE 'common6.h'
      INCLUDE 'tidaltiles_units.h'
*
*
*       Derive cgs constants from the common pc-Msun-Myr contract.
      PC = TT_PC_CM
      GM = TT_G_KMS*PC*1.0D10
      AU = PC/1.4959787D+13  ! bug fix from D+12 4/7/11 (was OK in 2009).
*
*       Form scaling factors for binary periods A*SQRT(A/M) to yrs and days.
      YRS = (RBAR*AU)**1.5/SQRT(ZMBAR)
      DAYS = 365.24*YRS
*
*       Specify conversion factors for lengths to solar radii & AU.
      SU = PC/(AU*6.955D+10)*RBAR*AU
      RAU = RBAR*AU
*
*       Copy solar mass scaling to new variable (M = BODY*<M>).
      SMU = ZMBAR
*
*     Time in Julian Myr and velocity in km/s; R/T = V*c exactly.
      TSTAR = 1.0D0/SQRT(TT_G_PCMYR)
      VSTAR = SQRT(TT_G_KMS)

      IF (KZ(22).LT.10.AND.KZ(19).LT.4) THEN
*
*       Ensure ZMBAR & RBAR > 0 (=0: assume <M>/Sun = 1, RBAR = 1 pc).
         IF (ZMBAR.LE.0.0D0) ZMBAR = FLOAT(N)/ZMASS
         IF (RBAR.LE.0.0D0) RBAR = 1.0
*
*       Scale to working units of RBAR in pc & ZMBAR in solar masses.
        TSTAR = TSTAR*SQRT(RBAR**3/(ZMASS*ZMBAR))
        VSTAR = VSTAR*SQRT(ZMASS*ZMBAR/RBAR)
*
      ELSE
*     USE scale factor from original input data instead
         TSTAR = TSTAR*SQRT(RBAR**3/ZMBAR)
         VSTAR = VSTAR*SQRT(ZMBAR/RBAR)
      END IF
*
*       Copy TSTAR to secondary time-scale factor.
      TSCALE = TSTAR
*     Runtime provenance: never infer this version from the runner.
      if(rank.eq.0) WRITE (6,*) 'TIDALTILES_UNIT_CONTRACT ',
     &     TT_UNIT_VERSION,' G_PCMYR=',TT_G_PCMYR,
     &     ' KMS_TO_PCMYR=',TT_KMS_PCMYR
*
*       Physical scaling: X, M, V, T from RBAR*X, ZMBAR*M, VSTAR*V, TSTAR*T.
      if(rank.eq.0)
     &WRITE (6,10)  RBAR, ZMBAR, VSTAR, TSTAR, BODYM*ZMBAR, SU
 10   FORMAT (/,12X,1P,'PHYSICAL SCALING:    R* =',E21.12,'  M* =',
     &     E21.12,'  V* =',E21.12,'  T* =',E21.12,'  <M> =',E21.12,
     &     '  SU =',E21.12)
*
*       Define relevant parameter for the GR case (RZ = 6*<m>/c^2).
          CLIGHT = 3.0D+05/VSTAR
          RZ = 6.0*ZMASS/(FLOAT(N)*CLIGHT**2)
          if(rank.eq.0)
     &    WRITE (6,20)  VSTAR, CLIGHT, RZ
   20     FORMAT (/,12X,'GR SCALING:    V* =',1P,E10.2,'  C =',E10.2,
     &                                  '  RZ =',E10.2)
*
      RETURN
*
      END
