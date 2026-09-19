      SUBROUTINE TAIL0(I,EI,RI2,VI2,POTI)
*
*
*      Initialization of tidal tail member.
*      ------------------------------------
*
      INCLUDE 'common6.h'
      INCLUDE 'galaxy.h'
      REAL*8  XI(3),XIDOT(3),FIRR(3),FD(3)
      SAVE IWARN
      DATA IWARN /0/
*
*
*       Transform to galactic coordinates.
      DO 5 K = 1,3
          XI(K) = X(K,I) + RG(K)
          XIDOT(K) = XDOT(K,I) + VG(K)
          FIRR(K) = 0.0
          FD(K) = 0.0
    5 CONTINUE

*       Obtain relevant force components.
      CALL XTRNLT(XI,XIDOT,FIRR,FD)
*
*       Update membership (first member in NTTOT = ITAIL0).
      IF (NTAIL.EQ.0) THEN
          NSTAIL = 0
*       Reserve the full KS c.m. range: new pairs are appended at
*       NTOT = N + NPAIRS (ksreg.F) and NPAIRS may grow up to KMAX,
*       so MIN(KMAX,NBIN0+10) would let late dynamical pairs overwrite
*       the tail block.
          ITAIL0 = NZERO + KMAX + 1
          NTTOT = ITAIL0 - 1
      END IF
*
*       Include safety test and 10 warnings for maximum membership.
      IF (NTTOT.GE.NMAX) THEN
          IWARN = IWARN + 1
          IF (IWARN.LE.10) THEN
              if(rank.eq.0)WRITE (6,10)  NTTOT
   10         FORMAT (' WARNING!    MAXIMUM TIDAL TAIL ',I6)
          END IF
          RETURN
      END IF

*       Record the exact one-way hand-off state before array compaction.
      IF (RANK.EQ.0) THEN
          OPEN (UNIT=53,STATUS='UNKNOWN',FORM='FORMATTED',
     &          FILE='handoff.53',POSITION='APPEND')
          WRITE (53,15) TTOT,TTOT*TSTAR,NAME(I),BODY(I)*ZMBAR,
     &          EI,SQRT(RI2),SQRT(VI2),POTI,RTIDE,ZMASS*ZMBAR,
     &          RBAR,TSTAR,
     &          (X(K,I),K=1,3),(XDOT(K,I),K=1,3),
     &          (XI(K),K=1,3),(XIDOT(K),K=1,3)
   15     FORMAT(1X,1P,2E24.16,I12,21E24.16)
          CALL FLUSH(53)
          CLOSE (53)
      END IF
*
*       Increase the memberships (note initial value ITAIL0 - 1).
      NTTOT = NTTOT + 1
      NTAIL = NTAIL + 1
      J = NTTOT
*
*       Copy escaper data and initialize integration variables.
      NAME(J) = NAME(I)
      BODY(J) = BODY(I)
      FF = 0.0
      FFD = 0.0
      DO 20 K = 1,3
          X(K,J) = XI(K)
          X0(K,J) = XI(K)
          XDOT(K,J) = XIDOT(K)
          X0DOT(K,J) = XIDOT(K)
          F(K,J) = 0.5*FIRR(K)
          FDOT(K,J) = ONE6*FD(K)
          FI(K,J) = FIRR(K)
          FIDOT(K,J) = FD(K)
          FF = FF + FIRR(K)**2
          FFD = FFD + FD(K)**2
   20 CONTINUE
*
*       Form a power-of-two step capped independently from the main system.
      DT = 0.5*ETAI*SQRT(FF/FFD)
      DTN = SMAXT
      IF (DT.GE.SMAXT) GO TO 30
      DO 25 K = 1,63
          DTN = 0.5D0*DTN
          IF (DTN.LT.DT) GO TO 30
   25 CONTINUE
   30 CONTINUE
      STEP(J) = DTN
      T0(J) = TIME
      TIMENW(J) = T0(J) + STEP(J)
*
      RETURN
*
      END
