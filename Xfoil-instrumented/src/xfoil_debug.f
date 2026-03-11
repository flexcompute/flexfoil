C***********************************************************************
C    Module:  xfoil_debug.f
C 
C    Debug output module for XFOIL instrumentation
C    Outputs JSON-formatted data for comparison with RustFoil
C
C    This module provides subroutines to dump XFOIL internal state
C    in JSON format for side-by-side comparison with RustFoil.
C***********************************************************************

C---- Block data to initialize debug common block
      BLOCK DATA XDBGINIT
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      COMMON /XDBGCTX/ ISDBG, IBLDBG
      COMMON /XDBGLOCK/ LDBGLOCK
      LOGICAL LDBG, LDBGLOCK
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER ISDBG, IBLDBG
      DATA LDBG /.TRUE./
      DATA LUDBG /77/
      DATA IDBGCALL /0/
      DATA IDBGITER /0/
      DATA ISDBG /0/
      DATA IBLDBG /0/
      DATA LDBGLOCK /.FALSE./
      END

C---- Dump TRDIF chain-rule derivatives (XT/WF/TT/DT/UT)
      SUBROUTINE DBGTRDIF_DERIVS(ISIDE, IBL, ITBL,
     &  WF1, WF2, XT,
     &  XT_A1, XT_X1, XT_X2, XT_T1, XT_T2, XT_D1, XT_D2, XT_U1, XT_U2,
     &  XT_MS, XT_RE,
     &  TT_A1, TT_X1, TT_X2, TT_T1, TT_T2, TT_D1, TT_D2, TT_U1, TT_U2,
     &  TT_MS, TT_RE,
     &  DT_A1, DT_X1, DT_X2, DT_T1, DT_T2, DT_D1, DT_D2, DT_U1, DT_U2,
     &  DT_MS, DT_RE,
     &  UT_A1, UT_X1, UT_X2, UT_T1, UT_T2, UT_D1, UT_D2, UT_U1, UT_U2,
     &  UT_MS, UT_RE)
      INTEGER ISIDE, IBL, ITBL
      REAL WF1, WF2, XT
      REAL XT_A1, XT_X1, XT_X2, XT_T1, XT_T2, XT_D1, XT_D2, XT_U1, XT_U2
      REAL XT_MS, XT_RE
      REAL TT_A1, TT_X1, TT_X2, TT_T1, TT_T2, TT_D1, TT_D2, TT_U1, TT_U2
      REAL TT_MS, TT_RE
      REAL DT_A1, DT_X1, DT_X2, DT_T1, DT_T2, DT_D1, DT_D2, DT_U1, DT_U2
      REAL DT_MS, DT_RE
      REAL UT_A1, UT_X1, UT_X2, UT_T1, UT_T2, UT_D1, UT_D2, UT_U1, UT_U2
      REAL UT_MS, UT_RE
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "TRDIF_DERIVS",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,I2,A)') '  "side": ', ISIDE, ','
      WRITE(LUDBG,'(A,I4,A)') '  "ibl": ', IBL, ','
      WRITE(LUDBG,'(A,I3,A)') '  "newton_iter": ', ITBL, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "WF1": ', WF1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "WF2": ', WF2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "XT": ', XT, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "XT_A1": ', XT_A1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "XT_X1": ', XT_X1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "XT_X2": ', XT_X2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "XT_T1": ', XT_T1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "XT_T2": ', XT_T2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "XT_D1": ', XT_D1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "XT_D2": ', XT_D2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "XT_U1": ', XT_U1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "XT_U2": ', XT_U2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "XT_MS": ', XT_MS, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "XT_RE": ', XT_RE, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "TT_A1": ', TT_A1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "TT_X1": ', TT_X1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "TT_X2": ', TT_X2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "TT_T1": ', TT_T1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "TT_T2": ', TT_T2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "TT_D1": ', TT_D1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "TT_D2": ', TT_D2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "TT_U1": ', TT_U1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "TT_U2": ', TT_U2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "TT_MS": ', TT_MS, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "TT_RE": ', TT_RE, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "DT_A1": ', DT_A1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "DT_X1": ', DT_X1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "DT_X2": ', DT_X2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "DT_T1": ', DT_T1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "DT_T2": ', DT_T2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "DT_D1": ', DT_D1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "DT_D2": ', DT_D2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "DT_U1": ', DT_U1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "DT_U2": ', DT_U2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "DT_MS": ', DT_MS, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "DT_RE": ', DT_RE, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "UT_A1": ', UT_A1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "UT_X1": ', UT_X1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "UT_X2": ', UT_X2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "UT_T1": ', UT_T1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "UT_T2": ', UT_T2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "UT_D1": ', UT_D1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "UT_D2": ', UT_D2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "UT_U1": ', UT_U1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "UT_U2": ', UT_U2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "UT_MS": ', UT_MS, ','
      WRITE(LUDBG,'(A,E15.8)') '  "UT_RE": ', UT_RE
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END

C---- Dump TRDIF Jacobian and residuals (combined laminar+turbulent)
      SUBROUTINE DBGTRDIF(ISIDE, IBL, ITBL, VS1, VS2, VSREZ)
      INTEGER ISIDE, IBL, ITBL
      REAL VS1(4,5), VS2(4,5), VSREZ(4)
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER I, J
C
      IF(.NOT.LDBG) RETURN
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "TRDIF",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,I2,A)') '  "side": ', ISIDE, ','
      WRITE(LUDBG,'(A,I4,A)') '  "ibl": ', IBL, ','
      WRITE(LUDBG,'(A,I3,A)') '  "newton_iter": ', ITBL, ','
C---- VS1 matrix
      WRITE(LUDBG,'(A)') '  "VS1": ['
      DO 10 I=1,4
        IF(I.LT.4) THEN
          WRITE(LUDBG,'(A,5(E14.7,A),A)')
     &      '    [', (VS1(I,J), ',', J=1,4), VS1(I,5), '],'
        ELSE
          WRITE(LUDBG,'(A,5(E14.7,A),A)')
     &      '    [', (VS1(I,J), ',', J=1,4), VS1(I,5), ']'
        ENDIF
   10 CONTINUE
      WRITE(LUDBG,'(A)') '  ],'
C---- VS2 matrix
      WRITE(LUDBG,'(A)') '  "VS2": ['
      DO 20 I=1,4
        IF(I.LT.4) THEN
          WRITE(LUDBG,'(A,5(E14.7,A),A)')
     &      '    [', (VS2(I,J), ',', J=1,4), VS2(I,5), '],'
        ELSE
          WRITE(LUDBG,'(A,5(E14.7,A),A)')
     &      '    [', (VS2(I,J), ',', J=1,4), VS2(I,5), ']'
        ENDIF
   20 CONTINUE
      WRITE(LUDBG,'(A)') '  ],'
C---- Residuals
      WRITE(LUDBG,'(A,4(E14.7,A),A)')
     &  '  "VSREZ": [', (VSREZ(I), ',', I=1,3), VSREZ(4), ']'
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump BLPRV state (T2/D2/DW2 after DSWAKI subtraction)
      SUBROUTINE DBGBLPRV(IS, IBL, ITBL)
      INTEGER IS, IBL, ITBL
      INCLUDE 'XBL.INC'
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      COMMON /XDBGIT/ ITBLDBG
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER ITBLDBG
C
      IF(.NOT.LDBG) RETURN
C---- Only dump upper surface IBL 30-35 to keep output manageable
      IF(IS.EQ.1) THEN
        IF(IBL.LT.24 .OR. IBL.GT.35) RETURN
      ELSE IF(IS.EQ.2) THEN
        IF(IBL.LT.70 .OR. IBL.GT.75) RETURN
      ELSE
        RETURN
      ENDIF
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "BLPRV_STATE",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,I2,A)') '  "side": ', IS, ','
      WRITE(LUDBG,'(A,I4,A)') '  "ibl": ', IBL, ','
      WRITE(LUDBG,'(A,I3,A)') '  "newton_iter": ', ITBL, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "T2": ', T2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "D2": ', D2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "DW2": ', DW2, ','
      WRITE(LUDBG,'(A,E15.8)') '  "U2": ', U2
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Set the current station context for debug output inside TRCHEK2
      SUBROUTINE DBGSETCTX(IS, IBL)
      INTEGER IS, IBL
      COMMON /XDBGCTX/ ISDBG, IBLDBG
      INTEGER ISDBG, IBLDBG
C
      ISDBG = IS
      IBLDBG = IBL
      RETURN
      END


C---- Open debug output file
      SUBROUTINE DBGOPEN(FNAME)
      CHARACTER*(*) FNAME
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
      OPEN(LUDBG, FILE=FNAME, STATUS='REPLACE')
      WRITE(LUDBG,'(A)') '{"events": ['
      IDBGCALL = 0
      IDBGITER = 0
      RETURN
      END


C---- Close debug output file
      SUBROUTINE DBGCLOSE()
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
      WRITE(LUDBG,'(A)') ']}'
      CLOSE(LUDBG)
      RETURN
      END


C---- Write comma separator if not first call
      SUBROUTINE DBGCOMMA()
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
      IF(IDBGCALL .GT. 0) THEN
        WRITE(LUDBG,'(A)') ','
      ENDIF
      IDBGCALL = IDBGCALL + 1
      RETURN
      END


C---- Set current iteration number
      SUBROUTINE DBGSETITER(ITER)
      INTEGER ITER
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IDBGITER = ITER
      RETURN
      END


C---- Dump VISCAL iteration start
      SUBROUTINE DBGVISCAL(ITER, ALFA, REINF, MINF, ACRIT)
      INTEGER ITER
      REAL ALFA, REINF, MINF, ACRIT
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "VISCAL",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', ITER, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "alpha_rad": ', ALFA, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "reynolds": ', REINF, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "mach": ', MINF, ','
      WRITE(LUDBG,'(A,E15.8)') '  "ncrit": ', ACRIT
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump VISCAL iteration result
      SUBROUTINE DBGVISCRES(ITER, RMSBL, RMXBL, CL, CD, CM, CDV, CDF)
      INTEGER ITER
      REAL RMSBL, RMXBL, CL, CD, CM, CDV, CDF
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "VISCAL_RESULT",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', ITER, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "rms_residual": ', RMSBL, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "max_residual": ', RMXBL, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "CL": ', CL, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "CD": ', CD, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "CM": ', CM, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "CD_viscous": ', CDV, ','
      WRITE(LUDBG,'(A,E15.8)') '  "CD_friction": ', CDF
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump BLVAR inputs and outputs
      SUBROUTINE DBGBLVAR(ISIDE, IBL, ITYP,
     &                    X, U, T, D, S, AMPL,
     &                    H, HK, HS, HC, RT, CF, DI, US, CQ, DE)
      INTEGER ISIDE, IBL, ITYP
      REAL X, U, T, D, S, AMPL
      REAL H, HK, HS, HC, RT, CF, DI, US, CQ, DE
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "BLVAR",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,I2,A)') '  "side": ', ISIDE, ','
      WRITE(LUDBG,'(A,I4,A)') '  "ibl": ', IBL, ','
      WRITE(LUDBG,'(A,I2,A)') '  "flow_type": ', ITYP, ','
      WRITE(LUDBG,'(A)') '  "input": {'
      WRITE(LUDBG,'(A,E15.8,A)') '    "x": ', X, ','
      WRITE(LUDBG,'(A,E15.8,A)') '    "u": ', U, ','
      WRITE(LUDBG,'(A,E15.8,A)') '    "theta": ', T, ','
      WRITE(LUDBG,'(A,E15.8,A)') '    "delta_star": ', D, ','
      WRITE(LUDBG,'(A,E15.8,A)') '    "ctau": ', S, ','
      WRITE(LUDBG,'(A,E15.8)') '    "ampl": ', AMPL
      WRITE(LUDBG,'(A)') '  },'
      WRITE(LUDBG,'(A)') '  "output": {'
      WRITE(LUDBG,'(A,E15.8,A)') '    "H": ', H, ','
      WRITE(LUDBG,'(A,E15.8,A)') '    "Hk": ', HK, ','
      WRITE(LUDBG,'(A,E15.8,A)') '    "Hs": ', HS, ','
      WRITE(LUDBG,'(A,E15.8,A)') '    "Hc": ', HC, ','
      WRITE(LUDBG,'(A,E15.8,A)') '    "Rtheta": ', RT, ','
      WRITE(LUDBG,'(A,E15.8,A)') '    "Cf": ', CF, ','
      WRITE(LUDBG,'(A,E15.8,A)') '    "Cd": ', DI, ','
      WRITE(LUDBG,'(A,E15.8,A)') '    "Us": ', US, ','
      WRITE(LUDBG,'(A,E15.8,A)') '    "Cq": ', CQ, ','
      WRITE(LUDBG,'(A,E15.8)') '    "De": ', DE
      WRITE(LUDBG,'(A)') '  }'
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump BLDIF Jacobian and residuals
      SUBROUTINE DBGBLDIF(ISIDE, IBL, ITYP,
     &                    X1, U1, T1, D1, S1, A1,
     &                    X2, U2, T2, D2, S2, A2,
     &                    VS1, VS2, VSREZ)
      INTEGER ISIDE, IBL, ITYP
      REAL X1, U1, T1, D1, S1, A1, X2, U2, T2, D2, S2, A2
      REAL VS1(4,5), VS2(4,5), VSREZ(4)
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      COMMON /XDBGIT/ ITBLDBG
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER ITBLDBG
      INTEGER I, J
C
      IF(.NOT.LDBG) RETURN
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "BLDIF",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,I2,A)') '  "side": ', ISIDE, ','
      WRITE(LUDBG,'(A,I4,A)') '  "ibl": ', IBL, ','
      WRITE(LUDBG,'(A,I3,A)') '  "newton_iter": ', ITBLDBG, ','
      WRITE(LUDBG,'(A,I2,A)') '  "flow_type": ', ITYP, ','
C---- primary state used in BLDIF
      WRITE(LUDBG,'(A,E15.8,A)') '  "X1": ', X1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "U1": ', U1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "T1": ', T1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "D1": ', D1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "S1": ', S1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "A1": ', A1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "X2": ', X2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "U2": ', U2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "T2": ', T2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "D2": ', D2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "S2": ', S2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "A2": ', A2, ','
C---- VS1 matrix
      WRITE(LUDBG,'(A)') '  "VS1": ['
      DO 10 I=1,4
        IF(I.LT.4) THEN
          WRITE(LUDBG,'(A,5(E14.7,A),A)')
     &      '    [', (VS1(I,J), ',', J=1,4), VS1(I,5), '],'
        ELSE
          WRITE(LUDBG,'(A,5(E14.7,A),A)')
     &      '    [', (VS1(I,J), ',', J=1,4), VS1(I,5), ']'
        ENDIF
   10 CONTINUE
      WRITE(LUDBG,'(A)') '  ],'
C---- VS2 matrix
      WRITE(LUDBG,'(A)') '  "VS2": ['
      DO 20 I=1,4
        IF(I.LT.4) THEN
          WRITE(LUDBG,'(A,5(E14.7,A),A)')
     &      '    [', (VS2(I,J), ',', J=1,4), VS2(I,5), '],'
        ELSE
          WRITE(LUDBG,'(A,5(E14.7,A),A)')
     &      '    [', (VS2(I,J), ',', J=1,4), VS2(I,5), ']'
        ENDIF
   20 CONTINUE
      WRITE(LUDBG,'(A)') '  ],'
C---- Residuals
      WRITE(LUDBG,'(A,4(E14.7,A),A)')
     &  '  "VSREZ": [', (VSREZ(I), ',', I=1,3), VSREZ(4), ']'
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump MRCHUE station state
      SUBROUTINE DBGMRCHUE(ISIDE, IBL, X, UE, THETA, DSTAR, CTAU,
     &                     HK, CF, TRAN)
      INTEGER ISIDE, IBL
      REAL X, UE, THETA, DSTAR, CTAU, HK, CF
      LOGICAL TRAN
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "MRCHUE",'
      WRITE(LUDBG,'(A,I2,A)') '  "side": ', ISIDE, ','
      WRITE(LUDBG,'(A,I4,A)') '  "ibl": ', IBL, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "x": ', X, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Ue": ', UE, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "theta": ', THETA, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "delta_star": ', DSTAR, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "ctau": ', CTAU, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Hk": ', HK, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Cf": ', CF, ','
      IF(TRAN) THEN
        WRITE(LUDBG,'(A)') '  "transitional": true'
      ELSE
        WRITE(LUDBG,'(A)') '  "transitional": false'
      ENDIF
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump MRCHUE Newton iteration state (per-iteration debugging)
C     Computes Hk and Rtheta from output state to ensure consistency
      SUBROUTINE DBGMRCHUE_ITER(IS, IBL, ITBL, XSI, UEI, DSWAKI,
     &                          THI_IN, DSI_IN, CTI_IN, AMI_IN,
     &                          THI_OUT, DSI_OUT, CTI_OUT, AMI_OUT,
     &                          VS2, VSREZ, DMAX, RLX, CONV,
     &                          HSTINV_V, GM1BL_V, REYBL_V)
      INTEGER IS, IBL, ITBL
      REAL XSI, UEI, DSWAKI
      REAL THI_IN, DSI_IN, CTI_IN, AMI_IN
      REAL THI_OUT, DSI_OUT, CTI_OUT, AMI_OUT
      REAL VS2(4,5), VSREZ(4), DMAX, RLX
      REAL HSTINV_V, GM1BL_V, REYBL_V
      LOGICAL CONV
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER I, J
C---- Local variables for computed Hk and Rtheta
      REAL H_OUT, HK_OUT, RT_OUT, MSQ
C
      IF(.NOT.LDBG) RETURN
C
C---- Compute H = delta_star / theta from output state
      H_OUT = DSI_OUT / THI_OUT
C
C---- Compute MSQ (Mach^2) from edge velocity
C     Using formula from BLKIN: MSQ = UEI^2 * HSTINV / (GM1BL * (1 - 0.5*UEI^2*HSTINV))
      MSQ = UEI*UEI*HSTINV_V / (GM1BL_V*(1.0 - 0.5*UEI*UEI*HSTINV_V))
C
C---- Compute Hk from H using HKIN correlation (inline for simplicity)
C     At M=0: Hk = H
C     With compressibility: Hk = (H - 0.29*MSQ) / (1 + 0.113*MSQ)
      HK_OUT = (H_OUT - 0.29*MSQ) / (1.0 + 0.113*MSQ)
C
C---- Compute Rtheta = Re * Ue * theta
      RT_OUT = REYBL_V * UEI * THI_OUT
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "MRCHUE_ITER",'
      WRITE(LUDBG,'(A,I2,A)') '  "side": ', IS, ','
      WRITE(LUDBG,'(A,I4,A)') '  "ibl": ', IBL, ','
      WRITE(LUDBG,'(A,I3,A)') '  "newton_iter": ', ITBL, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "x": ', XSI, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Ue": ', UEI, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "d_s_wake": ', DSWAKI, ','
C---- Input state (before update)
      WRITE(LUDBG,'(A)') '  "input": {'
      WRITE(LUDBG,'(A,E15.8,A)') '    "theta": ', THI_IN, ','
      WRITE(LUDBG,'(A,E15.8,A)') '    "delta_star": ', DSI_IN, ','
      WRITE(LUDBG,'(A,E15.8,A)') '    "ctau": ', CTI_IN, ','
      WRITE(LUDBG,'(A,E15.8)') '    "ampl": ', AMI_IN
      WRITE(LUDBG,'(A)') '  },'
C---- Output state (after update) - Hk and Rtheta computed from output theta, delta_star
      WRITE(LUDBG,'(A)') '  "output": {'
      WRITE(LUDBG,'(A,E15.8,A)') '    "theta": ', THI_OUT, ','
      WRITE(LUDBG,'(A,E15.8,A)') '    "delta_star": ', DSI_OUT, ','
      WRITE(LUDBG,'(A,E15.8,A)') '    "ctau": ', CTI_OUT, ','
      WRITE(LUDBG,'(A,E15.8,A)') '    "ampl": ', AMI_OUT, ','
      WRITE(LUDBG,'(A,E15.8,A)') '    "Hk": ', HK_OUT, ','
      WRITE(LUDBG,'(A,E15.8)') '    "Rtheta": ', RT_OUT
      WRITE(LUDBG,'(A)') '  },'
C---- Newton update vector
      WRITE(LUDBG,'(A,4(E14.7,A),A)')
     &  '  "vsrez": [', (VSREZ(I), ',', I=1,3), VSREZ(4), '],'
C---- VS2 matrix (Jacobian for current station)
      WRITE(LUDBG,'(A)') '  "VS2": ['
      DO 10 I=1,4
        IF(I.LT.4) THEN
          WRITE(LUDBG,'(A,5(E14.7,A),A)')
     &      '    [', (VS2(I,J), ',', J=1,4), VS2(I,5), '],'
        ELSE
          WRITE(LUDBG,'(A,5(E14.7,A),A)')
     &      '    [', (VS2(I,J), ',', J=1,4), VS2(I,5), ']'
        ENDIF
   10 CONTINUE
      WRITE(LUDBG,'(A)') '  ],'
C---- Convergence info
      WRITE(LUDBG,'(A,E15.8,A)') '  "dmax": ', DMAX, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "relaxation": ', RLX, ','
      IF(CONV) THEN
        WRITE(LUDBG,'(A)') '  "converged": true'
      ELSE
        WRITE(LUDBG,'(A)') '  "converged": false'
      ENDIF
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump UPDATE delta values
      SUBROUTINE DBGUPDATE(ISIDE, IBL, DCTAU, DTHET, DMASS, DUEDG, RLX)
      INTEGER ISIDE, IBL
      REAL DCTAU, DTHET, DMASS, DUEDG, RLX
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "UPDATE",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,I2,A)') '  "side": ', ISIDE, ','
      WRITE(LUDBG,'(A,I4,A)') '  "ibl": ', IBL, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "delta_ctau": ', DCTAU, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "delta_theta": ', DTHET, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "delta_mass": ', DMASS, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "delta_Ue": ', DUEDG, ','
      WRITE(LUDBG,'(A,E15.8)') '  "relaxation": ', RLX
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump UPDATE detailed before/after state
C     Captures complete state of BL variables before and after update
C     for tracking where divergence occurs
C     Enhanced version includes mass_defect, H, and Hk
      SUBROUTINE DBGUPDATEDETAIL(ISIDE, IBL,
     &   CTAU_B, THETA_B, DSTAR_B, UE_B,
     &   MASS_B, H_B, HK_B,
     &   DCTAU, DTHETA, DMASS, DUE,
     &   RLX,
     &   CTAU_A, THETA_A, DSTAR_A, UE_A,
     &   MASS_A, H_A, HK_A)
      INTEGER ISIDE, IBL
      REAL CTAU_B, THETA_B, DSTAR_B, UE_B, MASS_B, H_B, HK_B
      REAL DCTAU, DTHETA, DMASS, DUE
      REAL RLX
      REAL CTAU_A, THETA_A, DSTAR_A, UE_A, MASS_A, H_A, HK_A
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
C---- Limit output to first 30 stations per side (extended from 15)
      IF(IBL.GT.30) RETURN
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "UPDATE_DETAILED",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,I2,A)') '  "side": ', ISIDE, ','
      WRITE(LUDBG,'(A,I4,A)') '  "ibl": ', IBL, ','
C---- Before state
      WRITE(LUDBG,'(A,E15.8,A)') '  "ctau_before": ', CTAU_B, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "theta_before": ', THETA_B, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "delta_star_before": ', DSTAR_B, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "ue_before": ', UE_B, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "mass_defect_before": ', MASS_B, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "h_before": ', H_B, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "hk_before": ', HK_B, ','
C---- Delta values (before relaxation)
      WRITE(LUDBG,'(A,E15.8,A)') '  "delta_ctau": ', DCTAU, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "delta_theta": ', DTHETA, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "delta_mass": ', DMASS, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "delta_ue": ', DUE, ','
C---- Relaxation factor
      WRITE(LUDBG,'(A,E15.8,A)') '  "relaxation": ', RLX, ','
C---- After state
      WRITE(LUDBG,'(A,E15.8,A)') '  "ctau_after": ', CTAU_A, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "theta_after": ', THETA_A, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "delta_star_after": ', DSTAR_A, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "ue_after": ', UE_A, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "mass_defect_after": ', MASS_A, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "h_after": ', H_A, ','
      WRITE(LUDBG,'(A,E15.8)') '  "hk_after": ', HK_A
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump UESET edge velocity update via DIJ matrix
C     Called before and after UESET to capture the edge velocity changes
C     resulting from mass defect influence
      SUBROUTINE DBGUESET(NBL1, NBL2, UINV, MASS, UEDG_B, UEDG_A, IVX)
      INTEGER NBL1, NBL2, IVX
      REAL UINV(IVX,2), MASS(IVX,2), UEDG_B(IVX,2), UEDG_A(IVX,2)
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER I, IS, NOUT
C
      IF(.NOT.LDBG) RETURN
C
C---- Limit output to first 20 stations per side
      NOUT = MIN(20, MAX(NBL1, NBL2))
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "UESET",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,I4,A)') '  "nbl_upper": ', NBL1, ','
      WRITE(LUDBG,'(A,I4,A)') '  "nbl_lower": ', NBL2, ','
C
C---- Upper surface (IS=1)
      WRITE(LUDBG,'(A)') '  "upper_surface": {'
C---- Ue inviscid
      WRITE(LUDBG,'(A)') '    "ue_inviscid": ['
      DO 10 I=2, MIN(NOUT, NBL1)
        IF(I.LT.MIN(NOUT, NBL1)) THEN
          WRITE(LUDBG,'(A,E15.8,A)') '      ', UINV(I,1), ','
        ELSE
          WRITE(LUDBG,'(A,E15.8)') '      ', UINV(I,1)
        ENDIF
   10 CONTINUE
      WRITE(LUDBG,'(A)') '    ],'
C---- Mass defect
      WRITE(LUDBG,'(A)') '    "mass_defect": ['
      DO 20 I=2, MIN(NOUT, NBL1)
        IF(I.LT.MIN(NOUT, NBL1)) THEN
          WRITE(LUDBG,'(A,E15.8,A)') '      ', MASS(I,1), ','
        ELSE
          WRITE(LUDBG,'(A,E15.8)') '      ', MASS(I,1)
        ENDIF
   20 CONTINUE
      WRITE(LUDBG,'(A)') '    ],'
C---- Ue before
      WRITE(LUDBG,'(A)') '    "ue_before": ['
      DO 30 I=2, MIN(NOUT, NBL1)
        IF(I.LT.MIN(NOUT, NBL1)) THEN
          WRITE(LUDBG,'(A,E15.8,A)') '      ', UEDG_B(I,1), ','
        ELSE
          WRITE(LUDBG,'(A,E15.8)') '      ', UEDG_B(I,1)
        ENDIF
   30 CONTINUE
      WRITE(LUDBG,'(A)') '    ],'
C---- Ue after
      WRITE(LUDBG,'(A)') '    "ue_after": ['
      DO 40 I=2, MIN(NOUT, NBL1)
        IF(I.LT.MIN(NOUT, NBL1)) THEN
          WRITE(LUDBG,'(A,E15.8,A)') '      ', UEDG_A(I,1), ','
        ELSE
          WRITE(LUDBG,'(A,E15.8)') '      ', UEDG_A(I,1)
        ENDIF
   40 CONTINUE
      WRITE(LUDBG,'(A)') '    ]'
      WRITE(LUDBG,'(A)') '  },'
C
C---- Lower surface (IS=2)
      WRITE(LUDBG,'(A)') '  "lower_surface": {'
C---- Ue inviscid
      WRITE(LUDBG,'(A)') '    "ue_inviscid": ['
      DO 50 I=2, MIN(NOUT, NBL2)
        IF(I.LT.MIN(NOUT, NBL2)) THEN
          WRITE(LUDBG,'(A,E15.8,A)') '      ', UINV(I,2), ','
        ELSE
          WRITE(LUDBG,'(A,E15.8)') '      ', UINV(I,2)
        ENDIF
   50 CONTINUE
      WRITE(LUDBG,'(A)') '    ],'
C---- Mass defect
      WRITE(LUDBG,'(A)') '    "mass_defect": ['
      DO 60 I=2, MIN(NOUT, NBL2)
        IF(I.LT.MIN(NOUT, NBL2)) THEN
          WRITE(LUDBG,'(A,E15.8,A)') '      ', MASS(I,2), ','
        ELSE
          WRITE(LUDBG,'(A,E15.8)') '      ', MASS(I,2)
        ENDIF
   60 CONTINUE
      WRITE(LUDBG,'(A)') '    ],'
C---- Ue before
      WRITE(LUDBG,'(A)') '    "ue_before": ['
      DO 70 I=2, MIN(NOUT, NBL2)
        IF(I.LT.MIN(NOUT, NBL2)) THEN
          WRITE(LUDBG,'(A,E15.8,A)') '      ', UEDG_B(I,2), ','
        ELSE
          WRITE(LUDBG,'(A,E15.8)') '      ', UEDG_B(I,2)
        ENDIF
   70 CONTINUE
      WRITE(LUDBG,'(A)') '    ],'
C---- Ue after
      WRITE(LUDBG,'(A)') '    "ue_after": ['
      DO 80 I=2, MIN(NOUT, NBL2)
        IF(I.LT.MIN(NOUT, NBL2)) THEN
          WRITE(LUDBG,'(A,E15.8,A)') '      ', UEDG_A(I,2), ','
        ELSE
          WRITE(LUDBG,'(A,E15.8)') '      ', UEDG_A(I,2)
        ENDIF
   80 CONTINUE
      WRITE(LUDBG,'(A)') '    ]'
      WRITE(LUDBG,'(A)') '  }'
C
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump QDCALC DIJ matrix (condensed - just dimensions and sample)
      SUBROUTINE DBGQDCALC(N, NW, DIJ, IZX)
      INTEGER N, NW, IZX
      REAL DIJ(IZX,IZX)
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER I, J, NTOT
C
      IF(.NOT.LDBG) RETURN
      CALL DBGCOMMA()
      NTOT = N + NW
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "QDCALC",'
      WRITE(LUDBG,'(A,I4,A)') '  "n_airfoil": ', N, ','
      WRITE(LUDBG,'(A,I4,A)') '  "n_wake": ', NW, ','
      WRITE(LUDBG,'(A,I4,A)') '  "n_total": ', NTOT, ','
C---- Output first 10 diagonal elements as sample
      WRITE(LUDBG,'(A)') '  "DIJ_diagonal_sample": ['
      DO 10 I=1, MIN(10, NTOT)
        IF(I.LT.MIN(10,NTOT)) THEN
          WRITE(LUDBG,'(A,E14.7,A)') '    ', DIJ(I,I), ','
        ELSE
          WRITE(LUDBG,'(A,E14.7)') '    ', DIJ(I,I)
        ENDIF
   10 CONTINUE
      WRITE(LUDBG,'(A)') '  ],'
C---- Output first row as sample
      WRITE(LUDBG,'(A)') '  "DIJ_row1_sample": ['
      DO 20 J=1, MIN(10, NTOT)
        IF(J.LT.MIN(10,NTOT)) THEN
          WRITE(LUDBG,'(A,E14.7,A)') '    ', DIJ(1,J), ','
        ELSE
          WRITE(LUDBG,'(A,E14.7)') '    ', DIJ(1,J)
        ENDIF
   20 CONTINUE
      WRITE(LUDBG,'(A)') '  ]'
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump DQDM values from PSILIN call in QDCALC (wake-airfoil block)
C     This captures the source influence on tangential velocity at wake points
C     Note: QDCALC is called before debug file is opened, so write to separate file
      SUBROUTINE DBGDQDM_PSILIN(IWAKE, IW, N, XI, YI, NXI, NYI, DQDM,
     &                          XPANEL, YPANEL)
      INTEGER IWAKE, IW, N
      REAL XI, YI, NXI, NYI, DQDM(*)
      REAL XPANEL(*), YPANEL(*)
      INTEGER J
      LOGICAL FIRST
      SAVE FIRST
      DATA FIRST /.TRUE./
C
C---- Open file on first call
      IF(FIRST) THEN
        OPEN(UNIT=89, FILE='xfoil_dqdm.json', STATUS='REPLACE')
        WRITE(89,'(A)') '{"events": ['
        FIRST = .FALSE.
C
C------ Also dump panel geometry on first call
        OPEN(UNIT=90, FILE='xfoil_panel_geom.json', STATUS='REPLACE')
        WRITE(90,'(A)') '{'
        WRITE(90,'(A,I4,A)') '  "n": ', N, ','
        WRITE(90,'(A)') '  "x": ['
        DO 5 J=1, MIN(20, N)
          IF(J.LT.MIN(20, N)) THEN
            WRITE(90,'(A,E17.10,A)') '    ', XPANEL(J), ','
          ELSE
            WRITE(90,'(A,E17.10)') '    ', XPANEL(J)
          ENDIF
    5   CONTINUE
        WRITE(90,'(A)') '  ],'
        WRITE(90,'(A)') '  "y": ['
        DO 6 J=1, MIN(20, N)
          IF(J.LT.MIN(20, N)) THEN
            WRITE(90,'(A,E17.10,A)') '    ', YPANEL(J), ','
          ELSE
            WRITE(90,'(A,E17.10)') '    ', YPANEL(J)
          ENDIF
    6   CONTINUE
        WRITE(90,'(A)') '  ]'
        WRITE(90,'(A)') '}'
        CLOSE(90)
      ELSE
        WRITE(89,'(A)') ','
      ENDIF
C
      WRITE(89,'(A)') '{'
      WRITE(89,'(A)') '  "subroutine": "DQDM_PSILIN",'
      WRITE(89,'(A,I4,A)') '  "wake_idx": ', IWAKE, ','
      WRITE(89,'(A,I4,A)') '  "iw": ', IW, ','
      WRITE(89,'(A,I4,A)') '  "n_airfoil": ', N, ','
      WRITE(89,'(A,E17.10,A)') '  "x": ', XI, ','
      WRITE(89,'(A,E17.10,A)') '  "y": ', YI, ','
      WRITE(89,'(A,E17.10,A)') '  "nx": ', NXI, ','
      WRITE(89,'(A,E17.10,A)') '  "ny": ', NYI, ','
C---- Output first 20 DQDM values (airfoil panels influence on wake Qtan)
      WRITE(89,'(A)') '  "dqdm": ['
      DO 10 J=1, MIN(20, N)
        IF(J.LT.MIN(20, N)) THEN
          WRITE(89,'(A,E17.10,A)') '    ', DQDM(J), ','
        ELSE
          WRITE(89,'(A,E17.10)') '    ', DQDM(J)
        ENDIF
   10 CONTINUE
      WRITE(89,'(A)') '  ]'
      WRITE(89,'(A)') '}'
      RETURN
C
      ENTRY DBGDQDM_CLOSE()
C---- Close the DQDM debug file
      WRITE(89,'(A)') ']}'
      CLOSE(89)
      RETURN
      END


C---- Dump closure function HSL result
      SUBROUTINE DBGHSL(HK, RT, MSQ, HS, HS_HK, HS_RT, HS_MSQ)
      REAL HK, RT, MSQ, HS, HS_HK, HS_RT, HS_MSQ
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      COMMON /XDBGLOCK/ LDBGLOCK
      LOGICAL LDBG, LDBGLOCK
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
      IF(LDBGLOCK) RETURN
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "HSL",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Hk": ', HK, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Rtheta": ', RT, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Msq": ', MSQ, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Hs": ', HS, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Hs_Hk": ', HS_HK, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Hs_Rt": ', HS_RT, ','
      WRITE(LUDBG,'(A,E15.8)') '  "Hs_Msq": ', HS_MSQ
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump closure function CFL result
      SUBROUTINE DBGCFL(HK, RT, MSQ, CF, CF_HK, CF_RT, CF_MSQ)
      REAL HK, RT, MSQ, CF, CF_HK, CF_RT, CF_MSQ
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      COMMON /XDBGLOCK/ LDBGLOCK
      LOGICAL LDBG, LDBGLOCK
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
      IF(LDBGLOCK) RETURN
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "CFL",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Hk": ', HK, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Rtheta": ', RT, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Msq": ', MSQ, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Cf": ', CF, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Cf_Hk": ', CF_HK, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Cf_Rt": ', CF_RT, ','
      WRITE(LUDBG,'(A,E15.8)') '  "Cf_Msq": ', CF_MSQ
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump DAMPL amplification rate
      SUBROUTINE DBGDAMPL(HK, TH, RT, AX, AX_HK, AX_TH, AX_RT)
      REAL HK, TH, RT, AX, AX_HK, AX_TH, AX_RT
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      COMMON /XDBGLOCK/ LDBGLOCK
      LOGICAL LDBG, LDBGLOCK
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
      IF(LDBGLOCK) RETURN
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "DAMPL",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Hk": ', HK, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "theta": ', TH, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Rtheta": ', RT, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Ax": ', AX, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Ax_Hk": ', AX_HK, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Ax_theta": ', AX_TH, ','
      WRITE(LUDBG,'(A,E15.8)') '  "Ax_Rt": ', AX_RT
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump AXSET RMS averaging (transition detection)
      SUBROUTINE DBGAXSET(HK1, T1, RT1, A1, HK2, T2, RT2, A2, ACRIT,
     &                    AX1, AX2, AXA, DAX, AX)
      REAL HK1, T1, RT1, A1, HK2, T2, RT2, A2, ACRIT
      REAL AX1, AX2, AXA, DAX, AX
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "AXSET",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
C---- Station 1 inputs
      WRITE(LUDBG,'(A,E15.8,A)') '  "Hk1": ', HK1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "theta1": ', T1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Rtheta1": ', RT1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "ampl1": ', A1, ','
C---- Station 2 inputs
      WRITE(LUDBG,'(A,E15.8,A)') '  "Hk2": ', HK2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "theta2": ', T2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Rtheta2": ', RT2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "ampl2": ', A2, ','
C---- Critical N
      WRITE(LUDBG,'(A,E15.8,A)') '  "Ncrit": ', ACRIT, ','
C---- Individual amplification rates
      WRITE(LUDBG,'(A,E15.8,A)') '  "Ax1": ', AX1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Ax2": ', AX2, ','
C---- RMS average
      WRITE(LUDBG,'(A,E15.8,A)') '  "Axa_rms": ', AXA, ','
C---- Near-Ncrit correction
      WRITE(LUDBG,'(A,E15.8,A)') '  "Dax": ', DAX, ','
C---- Final output
      WRITE(LUDBG,'(A,E15.8)') '  "Ax_final": ', AX
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump TRCHEK transition check result (call with IS, IBL from caller)
      SUBROUTINE DBGTRCHEK(IS, IBL)
      INTEGER IS, IBL
      INCLUDE 'XBL.INC'
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "TRCHEK",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,I2,A)') '  "side": ', IS, ','
      WRITE(LUDBG,'(A,I4,A)') '  "ibl": ', IBL, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "x1": ', X1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "x2": ', X2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "ampl1": ', AMPL1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "ampl2": ', AMPL2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Hk1": ', HK1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Hk2": ', HK2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Rt1": ', RT1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Rt2": ', RT2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Ncrit": ', AMCRIT, ','
      IF(TRAN) THEN
        WRITE(LUDBG,'(A)') '  "transition": true,'
      ELSE
        WRITE(LUDBG,'(A)') '  "transition": false,'
      ENDIF
      WRITE(LUDBG,'(A,E15.8)') '  "x_transition": ', XT
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump TRCHEK2 iteration details for debugging N-factor evolution
C     Uses ISDBG, IBLDBG from /XDBGCTX/ set by DBGSETCTX before TRCHEK call
      SUBROUTINE DBGTRCHEK_DETAIL(ITAM, AXV, AMPL2V, RESV,
     &                            TRANV, XTV, WF1V, WF2V)
      INTEGER ITAM
      REAL AXV, AMPL2V, RESV, XTV, WF1V, WF2V
      LOGICAL TRANV
      INCLUDE 'XBL.INC'
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      COMMON /XDBGCTX/ ISDBG, IBLDBG
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER ISDBG, IBLDBG
C
      IF(.NOT.LDBG) RETURN
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "TRCHEK2_ITER",'
      WRITE(LUDBG,'(A,I4,A)') '  "global_iter": ', IDBGITER, ','
      WRITE(LUDBG,'(A,I2,A)') '  "side": ', ISDBG, ','
      WRITE(LUDBG,'(A,I4,A)') '  "ibl": ', IBLDBG, ','
      WRITE(LUDBG,'(A,I4,A)') '  "trchek_iter": ', ITAM, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "x1": ', X1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "x2": ', X2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "ampl1": ', AMPL1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "ampl2": ', AMPL2V, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "ax": ', AXV, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "residual": ', RESV, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "wf1": ', WF1V, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "wf2": ', WF2V, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "xt": ', XTV, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Hk1": ', HK1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Hk2": ', HK2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Rt1": ', RT1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Rt2": ', RT2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "T1": ', T1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "T2": ', T2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "U1": ', U1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "U2": ', U2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Ncrit": ', AMCRIT, ','
      IF(TRANV) THEN
        WRITE(LUDBG,'(A)') '  "transition": true'
      ELSE
        WRITE(LUDBG,'(A)') '  "transition": false'
      ENDIF
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump TRCHEK2 final result with convergence info
C     Uses ISDBG, IBLDBG from /XDBGCTX/ set by DBGSETCTX before TRCHEK call
      SUBROUTINE DBGTRCHEK_FINAL(NITER, CONVERGED, AXV, XTV)
      INTEGER NITER
      LOGICAL CONVERGED
      REAL AXV, XTV
      INCLUDE 'XBL.INC'
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      COMMON /XDBGCTX/ ISDBG, IBLDBG
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER ISDBG, IBLDBG
C
      IF(.NOT.LDBG) RETURN
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "TRCHEK2_FINAL",'
      WRITE(LUDBG,'(A,I4,A)') '  "global_iter": ', IDBGITER, ','
      WRITE(LUDBG,'(A,I2,A)') '  "side": ', ISDBG, ','
      WRITE(LUDBG,'(A,I4,A)') '  "ibl": ', IBLDBG, ','
      WRITE(LUDBG,'(A,I4,A)') '  "n_iterations": ', NITER, ','
      IF(CONVERGED) THEN
        WRITE(LUDBG,'(A)') '  "converged": true,'
      ELSE
        WRITE(LUDBG,'(A)') '  "converged": false,'
      ENDIF
      WRITE(LUDBG,'(A,E15.8,A)') '  "x1": ', X1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "x2": ', X2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "ampl1": ', AMPL1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "ampl2_final": ', AMPL2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "ax_final": ', AXV, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "xt_final": ', XTV, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Ncrit": ', AMCRIT, ','
      IF(TRAN) THEN
        WRITE(LUDBG,'(A)') '  "transition": true,'
      ELSE
        WRITE(LUDBG,'(A)') '  "transition": false,'
      ENDIF
      IF(TRFORC) THEN
        WRITE(LUDBG,'(A)') '  "forced": true'
      ELSE
        WRITE(LUDBG,'(A)') '  "forced": false'
      ENDIF
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump BLSOLV Newton solve state
      SUBROUTINE DBGBLSOLV(NSYS)
      INTEGER NSYS
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "BLSOLV",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,I4)') '  "system_size": ', NSYS
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump full SETBL Newton system (VA, VB, VDEL arrays)
C     Called AFTER SETBL builds the global Newton system.
C     VA(3,2,*) - diagonal block matrices [3 equations x 2 vars (ctau/theta)]
C     VB(3,2,*) - sub-diagonal block matrices
C     VM(3,*,*) - mass defect coupling matrix
C     VDEL(3,2,*) - RHS (col 1) and alpha sensitivity (col 2)
      SUBROUTINE DBGSETBLSYSTEM(NSYS, VA, VB, VM, VDEL, IVX, IZX)
      INTEGER NSYS, IVX, IZX
      REAL VA(3,2,IVX), VB(3,2,IVX)
      REAL VM(3,IZX,IVX), VDEL(3,2,IVX)
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER IV, K, J, NOUT
C
      IF(.NOT.LDBG) RETURN
C
C---- Limit output to first 20 stations for large systems
      NOUT = MIN(NSYS, 20)
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "SETBL_SYSTEM",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,I4,A)') '  "nsys": ', NSYS, ','
C
C---- VA blocks (diagonal) - output all stations up to NOUT
      WRITE(LUDBG,'(A)') '  "VA": ['
      DO 10 IV=1, NOUT
        WRITE(LUDBG,'(A)') '    ['
        DO 15 K=1, 3
          IF(K.LT.3) THEN
            WRITE(LUDBG,'(A,E14.7,A,E14.7,A)')
     &        '      [', VA(K,1,IV), ', ', VA(K,2,IV), '],'
          ELSE
            WRITE(LUDBG,'(A,E14.7,A,E14.7,A)')
     &        '      [', VA(K,1,IV), ', ', VA(K,2,IV), ']'
          ENDIF
   15   CONTINUE
        IF(IV.LT.NOUT) THEN
          WRITE(LUDBG,'(A)') '    ],'
        ELSE
          WRITE(LUDBG,'(A)') '    ]'
        ENDIF
   10 CONTINUE
      WRITE(LUDBG,'(A)') '  ],'
C
C---- VB blocks (sub-diagonal) - output all stations up to NOUT
      WRITE(LUDBG,'(A)') '  "VB": ['
      DO 20 IV=1, NOUT
        WRITE(LUDBG,'(A)') '    ['
        DO 25 K=1, 3
          IF(K.LT.3) THEN
            WRITE(LUDBG,'(A,E14.7,A,E14.7,A)')
     &        '      [', VB(K,1,IV), ', ', VB(K,2,IV), '],'
          ELSE
            WRITE(LUDBG,'(A,E14.7,A,E14.7,A)')
     &        '      [', VB(K,1,IV), ', ', VB(K,2,IV), ']'
          ENDIF
   25   CONTINUE
        IF(IV.LT.NOUT) THEN
          WRITE(LUDBG,'(A)') '    ],'
        ELSE
          WRITE(LUDBG,'(A)') '    ]'
        ENDIF
   20 CONTINUE
      WRITE(LUDBG,'(A)') '  ],'
C
C---- VDEL (RHS and alpha sensitivity) - output all stations up to NOUT
      WRITE(LUDBG,'(A)') '  "VDEL": ['
      DO 30 IV=1, NOUT
        IF(IV.LT.NOUT) THEN
          WRITE(LUDBG,'(A,3(E14.7,A),A)')
     &      '    [', (VDEL(K,1,IV), ', ', K=1,2), VDEL(3,1,IV), '],'
        ELSE
          WRITE(LUDBG,'(A,3(E14.7,A),A)')
     &      '    [', (VDEL(K,1,IV), ', ', K=1,2), VDEL(3,1,IV), ']'
        ENDIF
   30 CONTINUE
      WRITE(LUDBG,'(A)') '  ],'
C
C---- VM diagonal (mass coupling on diagonal) - sample first NOUT
      WRITE(LUDBG,'(A)') '  "VM_diagonal": ['
      DO 40 IV=1, NOUT
        IF(IV.LT.NOUT) THEN
          WRITE(LUDBG,'(A,3(E14.7,A),A)')
     &      '    [', (VM(K,IV,IV), ', ', K=1,2), VM(3,IV,IV), '],'
        ELSE
          WRITE(LUDBG,'(A,3(E14.7,A),A)')
     &      '    [', (VM(K,IV,IV), ', ', K=1,2), VM(3,IV,IV), ']'
        ENDIF
   40 CONTINUE
      WRITE(LUDBG,'(A)') '  ],'
C
C---- VM first row sample (coupling from station 1 to all others)
      WRITE(LUDBG,'(A)') '  "VM_row1": ['
      DO 50 J=1, NOUT
        IF(J.LT.NOUT) THEN
          WRITE(LUDBG,'(A,3(E14.7,A),A)')
     &      '    [', (VM(K,J,1), ', ', K=1,2), VM(3,J,1), '],'
        ELSE
          WRITE(LUDBG,'(A,3(E14.7,A),A)')
     &      '    [', (VM(K,J,1), ', ', K=1,2), VM(3,J,1), ']'
        ENDIF
   50 CONTINUE
      WRITE(LUDBG,'(A)') '  ]'
C
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump BLSOLV solution (Newton deltas after solving)
C     VDEL(1:3,1,1:NSYS) contains the solution [dCtau, dTheta, dMass]
      SUBROUTINE DBGBLSOLVSOLUTION(NSYS, VDEL, IVX)
      INTEGER NSYS, IVX
      REAL VDEL(3,2,IVX)
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER IV, K, NOUT
C
      IF(.NOT.LDBG) RETURN
C
C---- Limit output to first 30 stations for large systems
      NOUT = MIN(NSYS, 30)
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "BLSOLV_SOLUTION",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,I4,A)') '  "nsys": ', NSYS, ','
C
C---- Solution deltas [dCtau, dTheta, dMass] for each station
      WRITE(LUDBG,'(A)') '  "deltas": ['
      DO 10 IV=1, NOUT
        IF(IV.LT.NOUT) THEN
          WRITE(LUDBG,'(A,3(E14.7,A),A)')
     &      '    [', (VDEL(K,1,IV), ', ', K=1,2), VDEL(3,1,IV), '],'
        ELSE
          WRITE(LUDBG,'(A,3(E14.7,A),A)')
     &      '    [', (VDEL(K,1,IV), ', ', K=1,2), VDEL(3,1,IV), ']'
        ENDIF
   10 CONTINUE
      WRITE(LUDBG,'(A)') '  ]'
C
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump HKIN result
      SUBROUTINE DBGHKIN(H, MSQ, HK, HK_H, HK_MSQ)
      REAL H, MSQ, HK, HK_H, HK_MSQ
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      COMMON /XDBGLOCK/ LDBGLOCK
      LOGICAL LDBG, LDBGLOCK
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
C---- Skip if locked (nested call from another debug routine)
      IF(LDBGLOCK) RETURN
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "HKIN",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "H": ', H, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Msq": ', MSQ, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Hk": ', HK, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Hk_H": ', HK_H, ','
      WRITE(LUDBG,'(A,E15.8)') '  "Hk_Msq": ', HK_MSQ
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C***********************************************************************
C    INVISCID SOLVER DEBUG ROUTINES
C    These dump intermediate values for validating the Rust implementation
C***********************************************************************


C---- Dump NCALC result (normal vectors)
      SUBROUTINE DBGNCALC(N, NX, NY, S)
      INTEGER N
      REAL NX(N), NY(N), S(N)
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER I
C
      IF(.NOT.LDBG) RETURN
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "NCALC",'
      WRITE(LUDBG,'(A,I4,A)') '  "n": ', N, ','
C---- NX array
      WRITE(LUDBG,'(A)') '  "NX": ['
      DO 10 I=1, N
        IF(I.LT.N) THEN
          WRITE(LUDBG,'(A,E17.10,A)') '    ', NX(I), ','
        ELSE
          WRITE(LUDBG,'(A,E17.10)') '    ', NX(I)
        ENDIF
   10 CONTINUE
      WRITE(LUDBG,'(A)') '  ],'
C---- NY array
      WRITE(LUDBG,'(A)') '  "NY": ['
      DO 20 I=1, N
        IF(I.LT.N) THEN
          WRITE(LUDBG,'(A,E17.10,A)') '    ', NY(I), ','
        ELSE
          WRITE(LUDBG,'(A,E17.10)') '    ', NY(I)
        ENDIF
   20 CONTINUE
      WRITE(LUDBG,'(A)') '  ],'
C---- S array
      WRITE(LUDBG,'(A)') '  "S": ['
      DO 30 I=1, N
        IF(I.LT.N) THEN
          WRITE(LUDBG,'(A,E17.10,A)') '    ', S(I), ','
        ELSE
          WRITE(LUDBG,'(A,E17.10)') '    ', S(I)
        ENDIF
   30 CONTINUE
      WRITE(LUDBG,'(A)') '  ]'
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump APCALC result (panel angles)
      SUBROUTINE DBGAPCALC(N, APANEL)
      INTEGER N
      REAL APANEL(N)
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER I
C
      IF(.NOT.LDBG) RETURN
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "APCALC",'
      WRITE(LUDBG,'(A,I4,A)') '  "n": ', N, ','
C---- APANEL array
      WRITE(LUDBG,'(A)') '  "APANEL": ['
      DO 10 I=1, N
        IF(I.LT.N) THEN
          WRITE(LUDBG,'(A,E17.10,A)') '    ', APANEL(I), ','
        ELSE
          WRITE(LUDBG,'(A,E17.10)') '    ', APANEL(I)
        ENDIF
   10 CONTINUE
      WRITE(LUDBG,'(A)') '  ]'
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump TECALC result (trailing edge geometry)
      SUBROUTINE DBGTECALC(ANTE, ASTE, DSTE, SHARP)
      REAL ANTE, ASTE, DSTE
      LOGICAL SHARP
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "TECALC",'
      WRITE(LUDBG,'(A,E17.10,A)') '  "ANTE": ', ANTE, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "ASTE": ', ASTE, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "DSTE": ', DSTE, ','
      IF(SHARP) THEN
        WRITE(LUDBG,'(A)') '  "SHARP": true'
      ELSE
        WRITE(LUDBG,'(A)') '  "SHARP": false'
      ENDIF
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump PSILIN intermediate values for a panel pair
      SUBROUTINE DBGPSILIN(I, JO, PSIS, PSID, DZDG_JO, DZDG_JP)
      INTEGER I, JO
      REAL PSIS, PSID, DZDG_JO, DZDG_JP
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "PSILIN",'
      WRITE(LUDBG,'(A,I4,A)') '  "i": ', I, ','
      WRITE(LUDBG,'(A,I4,A)') '  "jo": ', JO, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "PSIS": ', PSIS, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "PSID": ', PSID, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "DZDG_JO": ', DZDG_JO, ','
      WRITE(LUDBG,'(A,E17.10)') '  "DZDG_JP": ', DZDG_JP
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump GGCALC result (influence matrix and RHS)
      SUBROUTINE DBGGGCALC(N, AIJ, GAMU1, GAMU2, IZX)
      INTEGER N, IZX
      REAL AIJ(IZX,IZX), GAMU1(N), GAMU2(N)
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER I, J, NOUT
C
      IF(.NOT.LDBG) RETURN
      CALL DBGCOMMA()
C---- Limit output to first 20 entries for large matrices
      NOUT = MIN(N, 20)
C
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "GGCALC",'
      WRITE(LUDBG,'(A,I4,A)') '  "n": ', N, ','
C---- Diagonal entries (sample)
      WRITE(LUDBG,'(A)') '  "AIJ_diagonal": ['
      DO 10 I=1, NOUT
        IF(I.LT.NOUT) THEN
          WRITE(LUDBG,'(A,E17.10,A)') '    ', AIJ(I,I), ','
        ELSE
          WRITE(LUDBG,'(A,E17.10)') '    ', AIJ(I,I)
        ENDIF
   10 CONTINUE
      WRITE(LUDBG,'(A)') '  ],'
C---- First row (sample)
      WRITE(LUDBG,'(A)') '  "AIJ_row1": ['
      DO 20 J=1, NOUT
        IF(J.LT.NOUT) THEN
          WRITE(LUDBG,'(A,E17.10,A)') '    ', AIJ(1,J), ','
        ELSE
          WRITE(LUDBG,'(A,E17.10)') '    ', AIJ(1,J)
        ENDIF
   20 CONTINUE
      WRITE(LUDBG,'(A)') '  ],'
C---- GAMU(:,1) - alpha=0 solution
      WRITE(LUDBG,'(A)') '  "GAMU_0": ['
      DO 30 I=1, NOUT
        IF(I.LT.NOUT) THEN
          WRITE(LUDBG,'(A,E17.10,A)') '    ', GAMU1(I), ','
        ELSE
          WRITE(LUDBG,'(A,E17.10)') '    ', GAMU1(I)
        ENDIF
   30 CONTINUE
      WRITE(LUDBG,'(A)') '  ],'
C---- GAMU(:,2) - alpha=90 solution
      WRITE(LUDBG,'(A)') '  "GAMU_90": ['
      DO 40 I=1, NOUT
        IF(I.LT.NOUT) THEN
          WRITE(LUDBG,'(A,E17.10,A)') '    ', GAMU2(I), ','
        ELSE
          WRITE(LUDBG,'(A,E17.10)') '    ', GAMU2(I)
        ENDIF
   40 CONTINUE
      WRITE(LUDBG,'(A)') '  ]'
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump SPECAL result (combined gamma for given alpha)
      SUBROUTINE DBGSPECAL(N, ALFA, GAM, QINV)
      INTEGER N
      REAL ALFA, GAM(N), QINV(N)
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER I, NOUT
C
      IF(.NOT.LDBG) RETURN
      CALL DBGCOMMA()
      NOUT = MIN(N, 40)
C
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "SPECAL",'
      WRITE(LUDBG,'(A,E17.10,A)') '  "alpha_rad": ', ALFA, ','
      WRITE(LUDBG,'(A,I4,A)') '  "n": ', N, ','
C---- GAM array (sample)
      WRITE(LUDBG,'(A)') '  "GAM": ['
      DO 10 I=1, NOUT
        IF(I.LT.NOUT) THEN
          WRITE(LUDBG,'(A,E17.10,A)') '    ', GAM(I), ','
        ELSE
          WRITE(LUDBG,'(A,E17.10)') '    ', GAM(I)
        ENDIF
   10 CONTINUE
      WRITE(LUDBG,'(A)') '  ],'
C---- QINV array (sample)
      WRITE(LUDBG,'(A)') '  "QINV": ['
      DO 20 I=1, NOUT
        IF(I.LT.NOUT) THEN
          WRITE(LUDBG,'(A,E17.10,A)') '    ', QINV(I), ','
        ELSE
          WRITE(LUDBG,'(A,E17.10)') '    ', QINV(I)
        ENDIF
   20 CONTINUE
      WRITE(LUDBG,'(A)') '  ]'
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump CLCALC result (lift and moment coefficients)
      SUBROUTINE DBGCLCALC(N, CL, CM, CDP, CP)
      INTEGER N
      REAL CL, CM, CDP, CP(N)
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER I, NOUT
C
      IF(.NOT.LDBG) RETURN
      CALL DBGCOMMA()
      NOUT = MIN(N, 40)
C
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "CLCALC",'
      WRITE(LUDBG,'(A,I4,A)') '  "n": ', N, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "CL": ', CL, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "CM": ', CM, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "CDP": ', CDP, ','
C---- CP array (sample)
      WRITE(LUDBG,'(A)') '  "CP": ['
      DO 10 I=1, NOUT
        IF(I.LT.NOUT) THEN
          WRITE(LUDBG,'(A,E17.10,A)') '    ', CP(I), ','
        ELSE
          WRITE(LUDBG,'(A,E17.10)') '    ', CP(I)
        ENDIF
   10 CONTINUE
      WRITE(LUDBG,'(A)') '  ]'
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump STFIND result (stagnation point)
      SUBROUTINE DBGSTFIND(IST, SST, X, Y, S, N)
      INTEGER IST, N
      REAL SST, X(N), Y(N), S(N)
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      REAL XST, YST, DS, T
C
      IF(.NOT.LDBG) RETURN
C
C---- Interpolate position at SST
      DS = S(IST+1) - S(IST)
      IF(DS .GT. 1.0E-10) THEN
        T = (SST - S(IST)) / DS
        XST = X(IST) + T * (X(IST+1) - X(IST))
        YST = Y(IST) + T * (Y(IST+1) - Y(IST))
      ELSE
        XST = X(IST)
        YST = Y(IST)
      ENDIF
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "STFIND",'
      WRITE(LUDBG,'(A,I4,A)') '  "IST": ', IST, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "SST": ', SST, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "XST": ', XST, ','
      WRITE(LUDBG,'(A,E17.10)') '  "YST": ', YST
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump IBLPAN result (BL station to panel mapping)
      SUBROUTINE DBGIBLPAN(IST, NBL1, NBL2, IBLTE1, IBLTE2)
      INTEGER IST, NBL1, NBL2, IBLTE1, IBLTE2
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "IBLPAN",'
      WRITE(LUDBG,'(A,I4,A)') '  "IST": ', IST, ','
      WRITE(LUDBG,'(A,I4,A)') '  "NBL_upper": ', NBL1, ','
      WRITE(LUDBG,'(A,I4,A)') '  "NBL_lower": ', NBL2, ','
      WRITE(LUDBG,'(A,I4,A)') '  "IBLTE_upper": ', IBLTE1, ','
      WRITE(LUDBG,'(A,I4)') '  "IBLTE_lower": ', IBLTE2
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump transition detection result
      SUBROUTINE DBGTRAN(ISIDE, ITRAN, XTRAN, TRFORC)
      INTEGER ISIDE, ITRAN
      REAL XTRAN
      LOGICAL TRFORC
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "TRANSITION",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,I2,A)') '  "side": ', ISIDE, ','
      WRITE(LUDBG,'(A,I4,A)') '  "ITRAN": ', ITRAN, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "x_transition": ', XTRAN, ','
      IF(TRFORC) THEN
        WRITE(LUDBG,'(A)') '  "forced": true'
      ELSE
        WRITE(LUDBG,'(A)') '  "forced": false'
      ENDIF
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump viscous solution final results
      SUBROUTINE DBGVISCFINAL(CL, CD, CM, CDFORM, CDFRIC,
     &                        XTR1, XTR2, XSEP1, XSEP2, NITER, CONV)
      REAL CL, CD, CM, CDFORM, CDFRIC, XTR1, XTR2, XSEP1, XSEP2
      INTEGER NITER
      LOGICAL CONV
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "VISCOUS_FINAL",'
      WRITE(LUDBG,'(A,E17.10,A)') '  "CL": ', CL, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "CD": ', CD, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "CM": ', CM, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "CD_form": ', CDFORM, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "CD_friction": ', CDFRIC, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "x_tr_upper": ', XTR1, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "x_tr_lower": ', XTR2, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "x_sep_upper": ', XSEP1, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "x_sep_lower": ', XSEP2, ','
      WRITE(LUDBG,'(A,I4,A)') '  "iterations": ', NITER, ','
      IF(CONV) THEN
        WRITE(LUDBG,'(A)') '  "converged": true'
      ELSE
        WRITE(LUDBG,'(A)') '  "converged": false'
      ENDIF
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump BL station initial state after IBLPAN/XICALC
      SUBROUTINE DBGBLINIT(ISIDE, IBL, XSSI, UEDG, THETA, DSTAR)
      INTEGER ISIDE, IBL
      REAL XSSI, UEDG, THETA, DSTAR
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "BL_INIT",'
      WRITE(LUDBG,'(A,I2,A)') '  "side": ', ISIDE, ','
      WRITE(LUDBG,'(A,I4,A)') '  "ibl": ', IBL, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "xssi": ', XSSI, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "Ue": ', UEDG, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "theta": ', THETA, ','
      WRITE(LUDBG,'(A,E17.10)') '  "delta_star": ', DSTAR
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump SETBL Newton system coefficients (VM, VA, VB matrices)
      SUBROUTINE DBGSETBL(IS, IBL, IV, NSYS, 
     &                    VA, VB, VM, VDEL, IVX)
      INTEGER IS, IBL, IV, NSYS, IVX
      REAL VA(3,2,IVX), VB(3,2,IVX)
      REAL VM(3,IVX,IVX), VDEL(3,2,IVX)
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER K, J
C
      IF(.NOT.LDBG) RETURN
C---- Only output for first few stations to avoid huge output
      IF(IBL.GT.5) RETURN
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "SETBL",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,I2,A)') '  "side": ', IS, ','
      WRITE(LUDBG,'(A,I4,A)') '  "ibl": ', IBL, ','
      WRITE(LUDBG,'(A,I4,A)') '  "iv": ', IV, ','
      WRITE(LUDBG,'(A,I4,A)') '  "nsys": ', NSYS, ','
C---- VA block (diagonal)
      WRITE(LUDBG,'(A)') '  "VA": ['
      DO 10 K=1, 3
        IF(K.LT.3) THEN
          WRITE(LUDBG,'(A,E14.7,A,E14.7,A)') 
     &      '    [', VA(K,1,IV), ', ', VA(K,2,IV), '],'
        ELSE
          WRITE(LUDBG,'(A,E14.7,A,E14.7,A)') 
     &      '    [', VA(K,1,IV), ', ', VA(K,2,IV), ']'
        ENDIF
   10 CONTINUE
      WRITE(LUDBG,'(A)') '  ],'
C---- VB block (lower diagonal)
      WRITE(LUDBG,'(A)') '  "VB": ['
      DO 20 K=1, 3
        IF(K.LT.3) THEN
          WRITE(LUDBG,'(A,E14.7,A,E14.7,A)') 
     &      '    [', VB(K,1,IV), ', ', VB(K,2,IV), '],'
        ELSE
          WRITE(LUDBG,'(A,E14.7,A,E14.7,A)') 
     &      '    [', VB(K,1,IV), ', ', VB(K,2,IV), ']'
        ENDIF
   20 CONTINUE
      WRITE(LUDBG,'(A)') '  ],'
C---- VDEL (residual and alpha sensitivity)
      WRITE(LUDBG,'(A)') '  "VDEL": ['
      DO 30 K=1, 3
        IF(K.LT.3) THEN
          WRITE(LUDBG,'(A,E14.7,A,E14.7,A)') 
     &      '    [', VDEL(K,1,IV), ', ', VDEL(K,2,IV), '],'
        ELSE
          WRITE(LUDBG,'(A,E14.7,A,E14.7,A)') 
     &      '    [', VDEL(K,1,IV), ', ', VDEL(K,2,IV), ']'
        ENDIF
   30 CONTINUE
      WRITE(LUDBG,'(A)') '  ],'
C---- VM row sample (first 10 columns of VM for this station)
      WRITE(LUDBG,'(A)') '  "VM_sample": ['
      DO 40 K=1, 3
        WRITE(LUDBG,'(A)',ADVANCE='NO') '    ['
        DO 45 J=1, MIN(10, NSYS)
          IF(J.LT.MIN(10,NSYS)) THEN
            WRITE(LUDBG,'(E12.5,A)',ADVANCE='NO') VM(K,J,IV), ', '
          ELSE
            WRITE(LUDBG,'(E12.5)',ADVANCE='NO') VM(K,J,IV)
          ENDIF
   45   CONTINUE
        IF(K.LT.3) THEN
          WRITE(LUDBG,'(A)') '],'
        ELSE
          WRITE(LUDBG,'(A)') ']'
        ENDIF
   40 CONTINUE
      WRITE(LUDBG,'(A)') '  ]'
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump 4x4 Newton system BEFORE GAUSS solve (captures original Jacobian)
C     This is called in MRCHUE right before CALL GAUSS to capture the
C     exact system being solved: VS2 * dx = VSREZ
      SUBROUTINE DBGVS2_BEFORE(IS, IBL, ITBL, VS2, VSREZ)
      INTEGER IS, IBL, ITBL
      REAL VS2(4,5), VSREZ(4)
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER I, J
C
      IF(.NOT.LDBG) RETURN
C---- Limit output to targeted station ranges
      IF(IS.EQ.1) THEN
        IF(IBL.LT.24 .OR. IBL.GT.35) RETURN
      ELSE IF(IS.EQ.2) THEN
        IF(IBL.LT.70 .OR. IBL.GT.75) RETURN
      ELSE
        RETURN
      ENDIF
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "VS2_BEFORE",'
      WRITE(LUDBG,'(A,I2,A)') '  "side": ', IS, ','
      WRITE(LUDBG,'(A,I4,A)') '  "ibl": ', IBL, ','
      WRITE(LUDBG,'(A,I3,A)') '  "newton_iter": ', ITBL, ','
C---- VS2 matrix (4x4 Jacobian, columns 1-4)
      WRITE(LUDBG,'(A)') '  "VS2_4x4": ['
      DO 10 I=1,4
        IF(I.LT.4) THEN
          WRITE(LUDBG,'(A,4(E15.8,A),A)')
     &      '    [', (VS2(I,J), ', ', J=1,3), VS2(I,4), '],'
        ELSE
          WRITE(LUDBG,'(A,4(E15.8,A),A)')
     &      '    [', (VS2(I,J), ', ', J=1,3), VS2(I,4), ']'
        ENDIF
   10 CONTINUE
      WRITE(LUDBG,'(A)') '  ],'
C---- VSREZ (RHS before solve = -residual)
      WRITE(LUDBG,'(A,4(E15.8,A),A)')
     &  '  "VSREZ_rhs": [', (VSREZ(I), ', ', I=1,3), VSREZ(4), ']'
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump MRCHUE mode decision details (direct vs inverse)
      SUBROUTINE DBGMRCHUE_MODE(IS, IBL, ITBL, DIRECT,
     &     HTEST, HKTEST, HMAX, DMAX, RLX, UEI, THI, DSI)
      INTEGER IS, IBL, ITBL
      LOGICAL DIRECT
      REAL HTEST, HKTEST, HMAX, DMAX, RLX, UEI, THI, DSI
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
      IF(IS.EQ.1) THEN
        IF(IBL.LT.24 .OR. IBL.GT.35) RETURN
      ELSE IF(IS.EQ.2) THEN
        IF(IBL.LT.70 .OR. IBL.GT.75) RETURN
      ELSE
        RETURN
      ENDIF
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "MRCHUE_MODE",'
      WRITE(LUDBG,'(A,I2,A)') '  "side": ', IS, ','
      WRITE(LUDBG,'(A,I4,A)') '  "ibl": ', IBL, ','
      WRITE(LUDBG,'(A,I3,A)') '  "newton_iter": ', ITBL, ','
      IF(DIRECT) THEN
        WRITE(LUDBG,'(A)') '  "direct_before": true,'
        WRITE(LUDBG,'(A)') '  "direct_after": true,'
      ELSE
        WRITE(LUDBG,'(A)') '  "direct_before": false,'
        WRITE(LUDBG,'(A)') '  "direct_after": false,'
      ENDIF
      WRITE(LUDBG,'(A,E15.8,A)') '  "Htest": ', HTEST, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Hk_test": ', HKTEST, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Hmax": ', HMAX, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "dmax": ', DMAX, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "rlx": ', RLX, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Ue": ', UEI, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "theta": ', THI, ','
      WRITE(LUDBG,'(A,E15.8)') '  "delta_star": ', DSI
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump VSREZ AFTER GAUSS solve (contains solution vector dx)
C     This is called in MRCHUE right after CALL GAUSS to capture
C     the Newton update: dx = [d_ctau/d_ampl, d_theta, d_dstar, d_ue]
      SUBROUTINE DBGVSREZ_AFTER(IS, IBL, ITBL, VSREZ)
      INTEGER IS, IBL, ITBL
      REAL VSREZ(4)
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER I
C
      IF(.NOT.LDBG) RETURN
C---- Limit output to targeted station ranges
      IF(IS.EQ.1) THEN
        IF(IBL.LT.24 .OR. IBL.GT.35) RETURN
      ELSE IF(IS.EQ.2) THEN
        IF(IBL.LT.70 .OR. IBL.GT.75) RETURN
      ELSE
        RETURN
      ENDIF
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "VSREZ_AFTER",'
      WRITE(LUDBG,'(A,I2,A)') '  "side": ', IS, ','
      WRITE(LUDBG,'(A,I4,A)') '  "ibl": ', IBL, ','
      WRITE(LUDBG,'(A,I3,A)') '  "newton_iter": ', ITBL, ','
C---- VSREZ (solution vector after GAUSS)
      WRITE(LUDBG,'(A,4(E15.8,A),A)')
     &  '  "VSREZ_solution": [', (VSREZ(I), ', ', I=1,3), VSREZ(4), ']'
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump BLDIF primary state (X1/U1/T1/D1/S1 and X2/U2/T2/D2/S2)
      SUBROUTINE DBGBLDIF_STATE(IS, IBL, ITYP,
     &                          X1, U1, T1, D1, S1,
     &                          X2, U2, T2, D2, S2)
      INTEGER IS, IBL, ITYP
      REAL X1, U1, T1, D1, S1
      REAL X2, U2, T2, D2, S2
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      COMMON /XDBGIT/ ITBLDBG
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER ITBLDBG
C
      IF(.NOT.LDBG) RETURN
C---- Limit output to early stations for Newton debugging
      IF(IS.EQ.1) THEN
        IF(IBL.GT.5) RETURN
      ELSE IF(IS.EQ.2) THEN
        IF(IBL.GT.5) RETURN
      ELSE
        RETURN
      ENDIF
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I0,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "BLDIF_STATE",'
      WRITE(LUDBG,'(A,I0,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,I0,A)') '  "side": ', IS, ','
      WRITE(LUDBG,'(A,I0,A)') '  "ibl": ', IBL, ','
      WRITE(LUDBG,'(A,I0,A)') '  "newton_iter": ', ITBLDBG, ','
      WRITE(LUDBG,'(A,I0,A)') '  "flow_type": ', ITYP, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "X1": ', X1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "U1": ', U1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "T1": ', T1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "D1": ', D1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "S1": ', S1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "X2": ', X2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "U2": ', U2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "T2": ', T2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "D2": ', D2, ','
      WRITE(LUDBG,'(A,E15.8)') '  "S2": ', S2
      WRITE(LUDBG,'(A)') '}'
C
      RETURN
      END


C---- Dump shape equation Jacobian intermediates
C     This captures all the Z_ coefficients and derivatives for verifying
C     the shape equation Jacobian computation in bldif
      SUBROUTINE DBGSHAPE(IS, IBL, ITYP,
     &                    Z_HS2, Z_CF2, Z_DI2, Z_T2, Z_U2,
     &                    Z_HCA, Z_HA, Z_UPW,
     &                    HS2_T2, HS2_D2, HS2_U2,
     &                    CF2_T2, CF2_D2, CF2_U2,
     &                    DI2_T2, DI2_D2, DI2_U2, DI2_S2,
     &                    HC2_T2, HC2_D2, HC2_U2,
     &                    H2_T2, H2_D2,
     &                    UPW_T2, UPW_D2, UPW_U2,
     &                    VS2_31, VS2_32, VS2_33, VS2_34,
     &                    CF1, CF2, DI1, DI2, XOT1, XOT2,
     &                    XLOG, CFX_UPW, DIX_UPW, Z_CFX, Z_DIX)
      INTEGER IS, IBL, ITYP
      REAL Z_HS2, Z_CF2, Z_DI2, Z_T2, Z_U2
      REAL Z_HCA, Z_HA, Z_UPW
      REAL HS2_T2, HS2_D2, HS2_U2
      REAL CF2_T2, CF2_D2, CF2_U2
      REAL DI2_T2, DI2_D2, DI2_U2, DI2_S2
      REAL HC2_T2, HC2_D2, HC2_U2
      REAL H2_T2, H2_D2
      REAL UPW_T2, UPW_D2, UPW_U2
      REAL VS2_31, VS2_32, VS2_33, VS2_34
      REAL CF1, CF2, DI1, DI2, XOT1, XOT2
      REAL XLOG, CFX_UPW, DIX_UPW, Z_CFX, Z_DIX
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      COMMON /XDBGIT/ ITBLDBG
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER ITBLDBG
C
      IF(.NOT.LDBG) RETURN
C---- Limit output to targeted station ranges
      IF(IS.EQ.1) THEN
        IF(IBL.LT.24 .OR. IBL.GT.35) RETURN
      ELSE IF(IS.EQ.2) THEN
        IF(IBL.LT.70 .OR. IBL.GT.75) RETURN
      ELSE
        RETURN
      ENDIF
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "SHAPE_JACOBIAN",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,I2,A)') '  "side": ', IS, ','
      WRITE(LUDBG,'(A,I4,A)') '  "ibl": ', IBL, ','
      WRITE(LUDBG,'(A,I3,A)') '  "newton_iter": ', ITBLDBG, ','
      WRITE(LUDBG,'(A,I2,A)') '  "flow_type": ', ITYP, ','
C---- Z coefficients
      WRITE(LUDBG,'(A,E15.8,A)') '  "Z_HS2": ', Z_HS2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Z_CF2": ', Z_CF2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Z_DI2": ', Z_DI2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Z_T2": ', Z_T2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Z_U2": ', Z_U2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Z_HCA": ', Z_HCA, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Z_HA": ', Z_HA, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Z_UPW": ', Z_UPW, ','
C---- Hs derivatives
      WRITE(LUDBG,'(A,E15.8,A)') '  "HS2_T2": ', HS2_T2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "HS2_D2": ', HS2_D2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "HS2_U2": ', HS2_U2, ','
C---- Cf derivatives
      WRITE(LUDBG,'(A,E15.8,A)') '  "CF2_T2": ', CF2_T2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "CF2_D2": ', CF2_D2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "CF2_U2": ', CF2_U2, ','
C---- DI derivatives
      WRITE(LUDBG,'(A,E15.8,A)') '  "DI2_T2": ', DI2_T2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "DI2_D2": ', DI2_D2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "DI2_U2": ', DI2_U2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "DI2_S2": ', DI2_S2, ','
C---- Hc derivatives
      WRITE(LUDBG,'(A,E15.8,A)') '  "HC2_T2": ', HC2_T2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "HC2_D2": ', HC2_D2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "HC2_U2": ', HC2_U2, ','
C---- H derivatives
      WRITE(LUDBG,'(A,E15.8,A)') '  "H2_T2": ', H2_T2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "H2_D2": ', H2_D2, ','
C---- UPW derivatives
      WRITE(LUDBG,'(A,E15.8,A)') '  "UPW_T2": ', UPW_T2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "UPW_D2": ', UPW_D2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "UPW_U2": ', UPW_U2, ','
C---- Final VS2 row 3 values (shape equation Jacobian)
      WRITE(LUDBG,'(A,E15.8,A)') '  "VS2_3_1": ', VS2_31, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "VS2_3_2": ', VS2_32, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "VS2_3_3": ', VS2_33, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "VS2_3_4": ', VS2_34, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "CF1": ', CF1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "CF2": ', CF2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "DI1": ', DI1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "DI2": ', DI2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "XOT1": ', XOT1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "XOT2": ', XOT2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "XLOG": ', XLOG, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "CFX_UPW": ', CFX_UPW, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "DIX_UPW": ', DIX_UPW, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Z_CFX": ', Z_CFX, ','
      WRITE(LUDBG,'(A,E15.8)') '  "Z_DIX": ', Z_DIX
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump momentum equation Jacobian intermediates
C     Captures terms used in VS2(2,2..4) for comparison with Rust
      SUBROUTINE DBGMOM(IS, IBL, ITYP, X1, U1, T1, D1, X2, U2, T2, D2,
     &                  Z_HA, Z_CFM, Z_CF2, Z_T2, Z_U2,
     &                  H2_T2, H2_D2,
     &                  CFM_T2, CFM_D2, CFM_U2,
     &                  CF2_T2, CF2_D2, CF2_U2,
     &                  VS2_22, VS2_23, VS2_24)
      INTEGER IS, IBL, ITYP
      REAL X1, U1, T1, D1, X2, U2, T2, D2
      REAL Z_HA, Z_CFM, Z_CF2, Z_T2, Z_U2
      REAL H2_T2, H2_D2
      REAL CFM_T2, CFM_D2, CFM_U2
      REAL CF2_T2, CF2_D2, CF2_U2
      REAL VS2_22, VS2_23, VS2_24
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      COMMON /XDBGIT/ ITBLDBG
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER ITBLDBG
C
      IF(.NOT.LDBG) RETURN
C---- Limit output to targeted station ranges
      IF(IS.EQ.1) THEN
        IF(IBL.LT.24 .OR. IBL.GT.35) RETURN
      ELSE IF(IS.EQ.2) THEN
        IF(IBL.LT.70 .OR. IBL.GT.75) RETURN
      ELSE
        RETURN
      ENDIF
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "MOM_JACOBIAN",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,I2,A)') '  "side": ', IS, ','
      WRITE(LUDBG,'(A,I4,A)') '  "ibl": ', IBL, ','
      WRITE(LUDBG,'(A,I3,A)') '  "newton_iter": ', ITBLDBG, ','
      WRITE(LUDBG,'(A,I2,A)') '  "flow_type": ', ITYP, ','
C---- primary state used in BLDIF
      WRITE(LUDBG,'(A,E15.8,A)') '  "X1": ', X1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "U1": ', U1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "T1": ', T1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "D1": ', D1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "X2": ', X2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "U2": ', U2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "T2": ', T2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "D2": ', D2, ','
C---- Z coefficients
      WRITE(LUDBG,'(A,E15.8,A)') '  "Z_HA": ', Z_HA, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Z_CFM": ', Z_CFM, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Z_CF2": ', Z_CF2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Z_T2": ', Z_T2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Z_U2": ', Z_U2, ','
C---- H derivatives
      WRITE(LUDBG,'(A,E15.8,A)') '  "H2_T2": ', H2_T2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "H2_D2": ', H2_D2, ','
C---- CFM derivatives
      WRITE(LUDBG,'(A,E15.8,A)') '  "CFM_T2": ', CFM_T2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "CFM_D2": ', CFM_D2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "CFM_U2": ', CFM_U2, ','
C---- CF2 derivatives
      WRITE(LUDBG,'(A,E15.8,A)') '  "CF2_T2": ', CF2_T2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "CF2_D2": ', CF2_D2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "CF2_U2": ', CF2_U2, ','
C---- Final VS2 entries
      WRITE(LUDBG,'(A,E15.8,A)') '  "VS2_2_2": ', VS2_22, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "VS2_2_3": ', VS2_23, ','
      WRITE(LUDBG,'(A,E15.8)') '  "VS2_2_4": ', VS2_24
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump laminar amplification Jacobian intermediates (row 1)
C     Captures terms used in VS2(1,2..4) for comparison with Rust
      SUBROUTINE DBGLAMAX(IS, IBL, ITYP,
     &                    AX, AX_HK2, AX_T2, AX_RT2, AX_A2,
     &                    HK2_T2, HK2_D2, HK2_U2,
     &                    RT2_T2, RT2_U2,
     &                    Z_AX, VS2_12, VS2_13, VS2_14)
      INTEGER IS, IBL, ITYP
      REAL AX, AX_HK2, AX_T2, AX_RT2, AX_A2
      REAL HK2_T2, HK2_D2, HK2_U2
      REAL RT2_T2, RT2_U2
      REAL Z_AX, VS2_12, VS2_13, VS2_14
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      COMMON /XDBGIT/ ITBLDBG
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER ITBLDBG
C
      IF(.NOT.LDBG) RETURN
C---- Only dump lower surface IBL 64-67 to keep output manageable
      IF(IS.NE.2) RETURN
      IF(IBL.LT.64 .OR. IBL.GT.67) RETURN
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "LAMINAR_JACOBIAN",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,I2,A)') '  "side": ', IS, ','
      WRITE(LUDBG,'(A,I4,A)') '  "ibl": ', IBL, ','
      WRITE(LUDBG,'(A,I3,A)') '  "newton_iter": ', ITBLDBG, ','
      WRITE(LUDBG,'(A,I2,A)') '  "flow_type": ', ITYP, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "AX": ', AX, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "AX_HK2": ', AX_HK2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "AX_T2": ', AX_T2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "AX_RT2": ', AX_RT2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "AX_A2": ', AX_A2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "HK2_T2": ', HK2_T2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "HK2_D2": ', HK2_D2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "HK2_U2": ', HK2_U2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "RT2_T2": ', RT2_T2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "RT2_U2": ', RT2_U2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Z_AX": ', Z_AX, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "VS2_1_2": ', VS2_12, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "VS2_1_3": ', VS2_13, ','
      WRITE(LUDBG,'(A,E15.8)') '  "VS2_1_4": ', VS2_14
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C***********************************************************************
C    FULL ARRAY DEBUG ROUTINES FOR INVISCID SOLVER COMPARISON
C    These dump complete arrays for validating the Rust implementation
C***********************************************************************


C---- Dump FULL gamma, velocity, and Cp distributions for given alpha
C     This outputs complete arrays (not samples) for comparison
      SUBROUTINE DBGFULLGAMMA(N, GAM, QINV, CP, CL, ALFA)
      INTEGER N
      REAL GAM(N), QINV(N), CP(N)
      REAL CL, ALFA
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER I
C
      IF(.NOT.LDBG) RETURN
      CALL DBGCOMMA()
C
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "FULLGAMMA",'
      WRITE(LUDBG,'(A,I4,A)') '  "n": ', N, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "alpha_rad": ', ALFA, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "CL": ', CL, ','
C---- Full GAM array
      WRITE(LUDBG,'(A)') '  "GAM": ['
      DO 10 I=1, N
        IF(I.LT.N) THEN
          WRITE(LUDBG,'(A,E17.10,A)') '    ', GAM(I), ','
        ELSE
          WRITE(LUDBG,'(A,E17.10)') '    ', GAM(I)
        ENDIF
   10 CONTINUE
      WRITE(LUDBG,'(A)') '  ],'
C---- Full QINV array
      WRITE(LUDBG,'(A)') '  "QINV": ['
      DO 20 I=1, N
        IF(I.LT.N) THEN
          WRITE(LUDBG,'(A,E17.10,A)') '    ', QINV(I), ','
        ELSE
          WRITE(LUDBG,'(A,E17.10)') '    ', QINV(I)
        ENDIF
   20 CONTINUE
      WRITE(LUDBG,'(A)') '  ],'
C---- Full CP array
      WRITE(LUDBG,'(A)') '  "CP": ['
      DO 30 I=1, N
        IF(I.LT.N) THEN
          WRITE(LUDBG,'(A,E17.10,A)') '    ', CP(I), ','
        ELSE
          WRITE(LUDBG,'(A,E17.10)') '    ', CP(I)
        ENDIF
   30 CONTINUE
      WRITE(LUDBG,'(A)') '  ]'
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump FULL AIC base solutions (GAMU arrays for alpha=0 and alpha=90)
C     This outputs complete arrays (not samples) for comparison
      SUBROUTINE DBGFULLAIC(N, GAMU, IZX)
      INTEGER N, IZX
      REAL GAMU(IZX,2)
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER I
C
      IF(.NOT.LDBG) RETURN
      CALL DBGCOMMA()
C
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "FULLAIC",'
      WRITE(LUDBG,'(A,I4,A)') '  "n": ', N, ','
C---- Full GAMU(:,1) - alpha=0 solution
      WRITE(LUDBG,'(A)') '  "gamu_0": ['
      DO 10 I=1, N
        IF(I.LT.N) THEN
          WRITE(LUDBG,'(A,E17.10,A)') '    ', GAMU(I,1), ','
        ELSE
          WRITE(LUDBG,'(A,E17.10)') '    ', GAMU(I,1)
        ENDIF
   10 CONTINUE
      WRITE(LUDBG,'(A)') '  ],'
C---- Full GAMU(:,2) - alpha=90 solution
      WRITE(LUDBG,'(A)') '  "gamu_90": ['
      DO 20 I=1, N
        IF(I.LT.N) THEN
          WRITE(LUDBG,'(A,E17.10,A)') '    ', GAMU(I,2), ','
        ELSE
          WRITE(LUDBG,'(A,E17.10)') '    ', GAMU(I,2)
        ENDIF
   20 CONTINUE
      WRITE(LUDBG,'(A)') '  ]'
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump FULL DIJ influence matrix for comparison with RustFoil
C     DIJ(I,J) = dUe_i / d(delta_star * Ue)_j
C     Outputs all matrix elements as nested JSON arrays
      SUBROUTINE DBGFULLDIJ(NSYS, DIJ, IZX)
      INTEGER NSYS, IZX
      REAL DIJ(IZX,IZX)
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER I, J
C
      IF(.NOT.LDBG) RETURN
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "FULLDIJ",'
      WRITE(LUDBG,'(A,I4,A)') '  "nsys": ', NSYS, ','
C---- Output flattened DIJ matrix in row-major order
      WRITE(LUDBG,'(A)') '  "dij": ['
      DO 10 I=1, NSYS
        WRITE(LUDBG,'(A)',ADVANCE='NO') '    ['
        DO 20 J=1, NSYS
          IF(J.LT.NSYS) THEN
            WRITE(LUDBG,'(E17.10,A)',ADVANCE='NO') DIJ(I,J), ', '
          ELSE
            WRITE(LUDBG,'(E17.10)',ADVANCE='NO') DIJ(I,J)
          ENDIF
   20   CONTINUE
        IF(I.LT.NSYS) THEN
          WRITE(LUDBG,'(A)') '],'
        ELSE
          WRITE(LUDBG,'(A)') ']'
        ENDIF
   10 CONTINUE
      WRITE(LUDBG,'(A)') '  ]'
      WRITE(LUDBG,'(A)') '}'
C
      RETURN
      END


C---- Dump force calculation detail (CL, CD, CM breakdown + TE quantities)
C     Called after CLCALC and CDCALC to capture force breakdown for comparison
      SUBROUTINE DBGFORCEDETAIL(ITER)
      INTEGER ITER
      INCLUDE 'XFOIL.INC'
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C     Local variables for trailing edge quantities
      REAL THT_UP, THT_LO, DST_UP, DST_LO, UE_UP, UE_LO, H_UP, H_LO
C
      IF(.NOT.LDBG) RETURN
      IF(.NOT.LVISC .OR. .NOT.LBLINI) RETURN
C
C---- Extract trailing edge BL quantities
      THT_UP = THET(IBLTE(1), 1)
      THT_LO = THET(IBLTE(2), 2)
      DST_UP = DSTR(IBLTE(1), 1)
      DST_LO = DSTR(IBLTE(2), 2)
      UE_UP  = UEDG(IBLTE(1), 1)
      UE_LO  = UEDG(IBLTE(2), 2)
C---- Compute shape factors
      IF(THT_UP .GT. 1.0E-12) THEN
        H_UP = DST_UP / THT_UP
      ELSE
        H_UP = 0.0
      ENDIF
      IF(THT_LO .GT. 1.0E-12) THEN
        H_LO = DST_LO / THT_LO
      ELSE
        H_LO = 0.0
      ENDIF
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "FORCE_DETAIL",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', ITER, ','
C---- Force coefficients
      WRITE(LUDBG,'(A,E17.10,A)') '  "CL": ', CL, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "CD_pressure": ', CDP, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "CD_friction": ', CDF, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "CD_total": ', CD, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "CM": ', CM, ','
C---- Upper surface TE quantities
      WRITE(LUDBG,'(A,E17.10,A)') '  "theta_te_upper": ', THT_UP, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "delta_star_te_upper": ',DST_UP,','
      WRITE(LUDBG,'(A,E17.10,A)') '  "ue_te_upper": ', UE_UP, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "h_te_upper": ', H_UP, ','
C---- Lower surface TE quantities
      WRITE(LUDBG,'(A,E17.10,A)') '  "theta_te_lower": ', THT_LO, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "delta_star_te_lower": ',DST_LO,','
      WRITE(LUDBG,'(A,E17.10,A)') '  "ue_te_lower": ', UE_LO, ','
      WRITE(LUDBG,'(A,E17.10)') '  "h_te_lower": ', H_LO
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump CD breakdown for force comparison
C     Outputs friction, pressure, and total CD plus TE BL quantities
C     Called after CDCALC in VISCAL
      SUBROUTINE DBGCDBREAKDOWN(ITER)
      INTEGER ITER
      INCLUDE 'XFOIL.INC'
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C     Local variables for TE quantities
      REAL THT_UP, THT_LO, DST_UP, DST_LO, UE_UP, UE_LO
      REAL H_UP, H_LO, CDF_UP, CDF_LO, CDP_UP, CDP_LO
C
      IF(.NOT.LDBG) RETURN
      IF(.NOT.LVISC .OR. .NOT.LBLINI) RETURN
C
C---- Extract trailing edge BL quantities
      THT_UP = THET(IBLTE(1), 1)
      THT_LO = THET(IBLTE(2), 2)
      DST_UP = DSTR(IBLTE(1), 1)
      DST_LO = DSTR(IBLTE(2), 2)
      UE_UP  = UEDG(IBLTE(1), 1)
      UE_LO  = UEDG(IBLTE(2), 2)
C---- Compute shape factors at TE
      IF(THT_UP .GT. 1.0E-12) THEN
        H_UP = DST_UP / THT_UP
      ELSE
        H_UP = 0.0
      ENDIF
      IF(THT_LO .GT. 1.0E-12) THEN
        H_LO = DST_LO / THT_LO
      ELSE
        H_LO = 0.0
      ENDIF
C
C---- Per-surface friction drag from TAU integration
C     Note: This requires computing the integrals - approximate for now
C     using the fact that CDF is the total integrated value
      CDF_UP = 0.0
      CDF_LO = 0.0
C---- Integrate friction on upper surface (IS=1)
      DO 10 IBL=3, IBLTE(1)
        I = IPAN(IBL,1)
        IM = IPAN(IBL-1,1)
        DX = (X(I) - X(IM))*CA + (Y(I) - Y(IM))*SA
        CDF_UP = CDF_UP + 0.5*(TAU(IBL,1)+TAU(IBL-1,1))*DX
   10 CONTINUE
      CDF_UP = CDF_UP * 2.0/QINF**2
C---- Integrate friction on lower surface (IS=2)
      DO 20 IBL=3, IBLTE(2)
        I = IPAN(IBL,2)
        IM = IPAN(IBL-1,2)
        DX = (X(I) - X(IM))*CA + (Y(I) - Y(IM))*SA
        CDF_LO = CDF_LO + 0.5*(TAU(IBL,2)+TAU(IBL-1,2))*DX
   20 CONTINUE
      CDF_LO = CDF_LO * 2.0/QINF**2
C
C---- Pressure drag per surface (from momentum deficit)
C     CDP = 2 * theta * Ue^((5+H)/2) at TE
      IF(UE_UP .GT. 0.01) THEN
        CDP_UP = 2.0 * THT_UP * UE_UP**((5.0+H_UP)/2.0)
      ELSE
        CDP_UP = 0.0
      ENDIF
      IF(UE_LO .GT. 0.01) THEN
        CDP_LO = 2.0 * THT_LO * UE_LO**((5.0+H_LO)/2.0)
      ELSE
        CDP_LO = 0.0
      ENDIF
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "CD_BREAKDOWN",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', ITER, ','
C---- Total force coefficients
      WRITE(LUDBG,'(A,E17.10,A)') '  "cd_friction": ', CDF, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "cd_pressure": ', CDP, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "cd_total": ', CD, ','
C---- Per-surface friction breakdown
      WRITE(LUDBG,'(A,E17.10,A)') '  "cd_friction_upper": ', CDF_UP, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "cd_friction_lower": ', CDF_LO, ','
C---- Per-surface pressure breakdown (approximate)
      WRITE(LUDBG,'(A,E17.10,A)') '  "cd_pressure_upper": ', CDP_UP, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "cd_pressure_lower": ', CDP_LO, ','
C---- Upper surface TE quantities
      WRITE(LUDBG,'(A,E17.10,A)') '  "theta_te_upper": ', THT_UP, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "delta_star_te_upper": ',DST_UP,','
      WRITE(LUDBG,'(A,E17.10,A)') '  "ue_te_upper": ', UE_UP, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "h_te_upper": ', H_UP, ','
C---- Lower surface TE quantities
      WRITE(LUDBG,'(A,E17.10,A)') '  "theta_te_lower": ', THT_LO, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "delta_star_te_lower": ',DST_LO,','
      WRITE(LUDBG,'(A,E17.10,A)') '  "ue_te_lower": ', UE_LO, ','
      WRITE(LUDBG,'(A,E17.10)') '  "h_te_lower": ', H_LO
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump CL detail for force comparison
C     Outputs CL, CM, CDP, alpha for comparison with RustFoil
C     Called after CLCALC in VISCAL
      SUBROUTINE DBGCLDETAIL(ITER)
      INTEGER ITER
      INCLUDE 'XFOIL.INC'
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C     Local variable for degrees
      REAL ADEG
C
      IF(.NOT.LDBG) RETURN
C
      ADEG = ALFA / DTOR
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "CL_DETAIL",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', ITER, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "cl": ', CL, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "cm": ', CM, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "cdp": ', CDP, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "alpha_rad": ', ALFA, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "alpha_deg": ', ADEG, ','
C---- Also output GAMTE if available (trailing edge circulation)
      WRITE(LUDBG,'(A,E17.10,A)') '  "gamte": ', GAMTE, ','
C---- CL sensitivity to alpha
      WRITE(LUDBG,'(A,E17.10)') '  "cl_alpha": ', CL_ALF
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C***********************************************************************
C    FULL ARRAY DEBUG ROUTINES FOR VISCAL ITERATION COMPARISON
C    These dump complete arrays per iteration for validating RustFoil
C***********************************************************************


C---- Dump FULL BL state arrays per VISCAL iteration
C     Outputs theta, delta_star, Ue, Hk, Cf for ALL stations on both surfaces
      SUBROUTINE DBGFULL_BL_STATE(ITER)
      INTEGER ITER
      INCLUDE 'XFOIL.INC'
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      COMMON /XDBGLOCK/ LDBGLOCK
      LOGICAL LDBG, LDBGLOCK
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER IS, IBL
      REAL MSQ, H, HK, HK_H, HK_M
      REAL HSTINV_L
C
      IF(.NOT.LDBG) RETURN
      IF(.NOT.LBLINI) RETURN
C
C---- Set lock to prevent nested debug calls during array writes
      LDBGLOCK = .TRUE.
C
C---- Compute HSTINV for Hk calculation
      HSTINV_L = GAMM1*(MINF/QINF)**2 / (1.0 + 0.5*GAMM1*MINF**2)
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "FULL_BL_STATE",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', ITER, ','
      WRITE(LUDBG,'(A,I4,A)') '  "nbl_upper": ', NBL(1), ','
      WRITE(LUDBG,'(A,I4,A)') '  "nbl_lower": ', NBL(2), ','
C
C---- Upper surface (IS=1)
      WRITE(LUDBG,'(A)') '  "upper_surface": {'
C---- Theta
      WRITE(LUDBG,'(A)') '    "theta": ['
      DO 10 IBL=2, NBL(1)
        IF(IBL.LT.NBL(1)) THEN
          WRITE(LUDBG,'(A,E17.10,A)') '      ', THET(IBL,1), ','
        ELSE
          WRITE(LUDBG,'(A,E17.10)') '      ', THET(IBL,1)
        ENDIF
   10 CONTINUE
      WRITE(LUDBG,'(A)') '    ],'
C---- Delta_star
      WRITE(LUDBG,'(A)') '    "delta_star": ['
      DO 20 IBL=2, NBL(1)
        IF(IBL.LT.NBL(1)) THEN
          WRITE(LUDBG,'(A,E17.10,A)') '      ', DSTR(IBL,1), ','
        ELSE
          WRITE(LUDBG,'(A,E17.10)') '      ', DSTR(IBL,1)
        ENDIF
   20 CONTINUE
      WRITE(LUDBG,'(A)') '    ],'
C---- Ue
      WRITE(LUDBG,'(A)') '    "Ue": ['
      DO 30 IBL=2, NBL(1)
        IF(IBL.LT.NBL(1)) THEN
          WRITE(LUDBG,'(A,E17.10,A)') '      ', UEDG(IBL,1), ','
        ELSE
          WRITE(LUDBG,'(A,E17.10)') '      ', UEDG(IBL,1)
        ENDIF
   30 CONTINUE
      WRITE(LUDBG,'(A)') '    ],'
C---- Hk (computed)
      WRITE(LUDBG,'(A)') '    "Hk": ['
      DO 40 IBL=2, NBL(1)
        IF(THET(IBL,1).GT.1.0E-12) THEN
          H = DSTR(IBL,1) / THET(IBL,1)
        ELSE
          H = 2.5
        ENDIF
        MSQ = UEDG(IBL,1)**2 * HSTINV_L 
     &      / (GAMM1*(1.0 - 0.5*UEDG(IBL,1)**2*HSTINV_L))
        IF(MSQ.LT.0.0) MSQ = 0.0
        CALL HKIN(H, MSQ, HK, HK_H, HK_M)
        IF(IBL.LT.NBL(1)) THEN
          WRITE(LUDBG,'(A,E17.10,A)') '      ', HK, ','
        ELSE
          WRITE(LUDBG,'(A,E17.10)') '      ', HK
        ENDIF
   40 CONTINUE
      WRITE(LUDBG,'(A)') '    ],'
C---- Cf (using stored TAU)
      WRITE(LUDBG,'(A)') '    "Cf": ['
      DO 50 IBL=2, NBL(1)
        IF(IBL.LT.NBL(1)) THEN
          WRITE(LUDBG,'(A,E17.10,A)') '      ', 
     &      2.0*TAU(IBL,1)/(UEDG(IBL,1)**2+1.0E-20), ','
        ELSE
          WRITE(LUDBG,'(A,E17.10)') '      ', 
     &      2.0*TAU(IBL,1)/(UEDG(IBL,1)**2+1.0E-20)
        ENDIF
   50 CONTINUE
      WRITE(LUDBG,'(A)') '    ],'
C---- Ctau/Ampl
      WRITE(LUDBG,'(A)') '    "ctau": ['
      DO 55 IBL=2, NBL(1)
        IF(IBL.LT.NBL(1)) THEN
          WRITE(LUDBG,'(A,E17.10,A)') '      ', CTAU(IBL,1), ','
        ELSE
          WRITE(LUDBG,'(A,E17.10)') '      ', CTAU(IBL,1)
        ENDIF
   55 CONTINUE
      WRITE(LUDBG,'(A)') '    ]'
      WRITE(LUDBG,'(A)') '  },'
C
C---- Lower surface (IS=2)
      WRITE(LUDBG,'(A)') '  "lower_surface": {'
C---- Theta
      WRITE(LUDBG,'(A)') '    "theta": ['
      DO 110 IBL=2, NBL(2)
        IF(IBL.LT.NBL(2)) THEN
          WRITE(LUDBG,'(A,E17.10,A)') '      ', THET(IBL,2), ','
        ELSE
          WRITE(LUDBG,'(A,E17.10)') '      ', THET(IBL,2)
        ENDIF
  110 CONTINUE
      WRITE(LUDBG,'(A)') '    ],'
C---- Delta_star
      WRITE(LUDBG,'(A)') '    "delta_star": ['
      DO 120 IBL=2, NBL(2)
        IF(IBL.LT.NBL(2)) THEN
          WRITE(LUDBG,'(A,E17.10,A)') '      ', DSTR(IBL,2), ','
        ELSE
          WRITE(LUDBG,'(A,E17.10)') '      ', DSTR(IBL,2)
        ENDIF
  120 CONTINUE
      WRITE(LUDBG,'(A)') '    ],'
C---- Ue
      WRITE(LUDBG,'(A)') '    "Ue": ['
      DO 130 IBL=2, NBL(2)
        IF(IBL.LT.NBL(2)) THEN
          WRITE(LUDBG,'(A,E17.10,A)') '      ', UEDG(IBL,2), ','
        ELSE
          WRITE(LUDBG,'(A,E17.10)') '      ', UEDG(IBL,2)
        ENDIF
  130 CONTINUE
      WRITE(LUDBG,'(A)') '    ],'
C---- Hk (computed)
      WRITE(LUDBG,'(A)') '    "Hk": ['
      DO 140 IBL=2, NBL(2)
        IF(THET(IBL,2).GT.1.0E-12) THEN
          H = DSTR(IBL,2) / THET(IBL,2)
        ELSE
          H = 2.5
        ENDIF
        MSQ = UEDG(IBL,2)**2 * HSTINV_L 
     &      / (GAMM1*(1.0 - 0.5*UEDG(IBL,2)**2*HSTINV_L))
        IF(MSQ.LT.0.0) MSQ = 0.0
        CALL HKIN(H, MSQ, HK, HK_H, HK_M)
        IF(IBL.LT.NBL(2)) THEN
          WRITE(LUDBG,'(A,E17.10,A)') '      ', HK, ','
        ELSE
          WRITE(LUDBG,'(A,E17.10)') '      ', HK
        ENDIF
  140 CONTINUE
      WRITE(LUDBG,'(A)') '    ],'
C---- Cf (using stored TAU)
      WRITE(LUDBG,'(A)') '    "Cf": ['
      DO 150 IBL=2, NBL(2)
        IF(IBL.LT.NBL(2)) THEN
          WRITE(LUDBG,'(A,E17.10,A)') '      ', 
     &      2.0*TAU(IBL,2)/(UEDG(IBL,2)**2+1.0E-20), ','
        ELSE
          WRITE(LUDBG,'(A,E17.10)') '      ', 
     &      2.0*TAU(IBL,2)/(UEDG(IBL,2)**2+1.0E-20)
        ENDIF
  150 CONTINUE
      WRITE(LUDBG,'(A)') '    ],'
C---- Ctau/Ampl
      WRITE(LUDBG,'(A)') '    "ctau": ['
      DO 155 IBL=2, NBL(2)
        IF(IBL.LT.NBL(2)) THEN
          WRITE(LUDBG,'(A,E17.10,A)') '      ', CTAU(IBL,2), ','
        ELSE
          WRITE(LUDBG,'(A,E17.10)') '      ', CTAU(IBL,2)
        ENDIF
  155 CONTINUE
      WRITE(LUDBG,'(A)') '    ]'
      WRITE(LUDBG,'(A)') '  }'
C
      WRITE(LUDBG,'(A)') '}'
C
C---- Release lock
      LDBGLOCK = .FALSE.
      RETURN
      END


C---- Dump FULL gamma array per VISCAL iteration
C     Outputs gamma[i] for all panels
      SUBROUTINE DBGFULL_GAMMA_ITER(ITER)
      INTEGER ITER
      INCLUDE 'XFOIL.INC'
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER I
C
      IF(.NOT.LDBG) RETURN
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "FULL_GAMMA_ITER",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', ITER, ','
      WRITE(LUDBG,'(A,I4,A)') '  "n_panels": ', N, ','
C
C---- Full gamma array
      WRITE(LUDBG,'(A)') '  "gamma": ['
      DO 10 I=1, N
        IF(I.LT.N) THEN
          WRITE(LUDBG,'(A,E17.10,A)') '    ', GAM(I), ','
        ELSE
          WRITE(LUDBG,'(A,E17.10)') '    ', GAM(I)
        ENDIF
   10 CONTINUE
      WRITE(LUDBG,'(A)') '  ]'
C
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump FULL N-factor (amplification ratio) array per VISCAL iteration
C     Outputs CTAU[ibl] for laminar stations (which stores log(ampl ratio))
      SUBROUTINE DBGFULL_NFACTOR(ITER)
      INTEGER ITER
      INCLUDE 'XFOIL.INC'
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER IS, IBL, NLAM
C
      IF(.NOT.LDBG) RETURN
      IF(.NOT.LBLINI) RETURN
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "FULL_NFACTOR",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', ITER, ','
      WRITE(LUDBG,'(A,I4,A)') '  "itran_upper": ', ITRAN(1), ','
      WRITE(LUDBG,'(A,I4,A)') '  "itran_lower": ', ITRAN(2), ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "ncrit_upper": ', ACRIT(1), ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "ncrit_lower": ', ACRIT(2), ','
C
C---- Upper surface N-factor (laminar stations only)
C     CTAU stores log(amplification ratio) for laminar, Ctau for turbulent
      NLAM = MIN(ITRAN(1)-1, NBL(1))
      WRITE(LUDBG,'(A)') '  "upper_ampl": ['
      DO 10 IBL=2, NLAM
        IF(IBL.LT.NLAM) THEN
          WRITE(LUDBG,'(A,E17.10,A)') '    ', CTAU(IBL,1), ','
        ELSE
          WRITE(LUDBG,'(A,E17.10)') '    ', CTAU(IBL,1)
        ENDIF
   10 CONTINUE
      IF(NLAM.LT.2) WRITE(LUDBG,'(A)') '    0.0'
      WRITE(LUDBG,'(A)') '  ],'
C
C---- Lower surface N-factor (laminar stations only)
      NLAM = MIN(ITRAN(2)-1, NBL(2))
      WRITE(LUDBG,'(A)') '  "lower_ampl": ['
      DO 20 IBL=2, NLAM
        IF(IBL.LT.NLAM) THEN
          WRITE(LUDBG,'(A,E17.10,A)') '    ', CTAU(IBL,2), ','
        ELSE
          WRITE(LUDBG,'(A,E17.10)') '    ', CTAU(IBL,2)
        ENDIF
   20 CONTINUE
      IF(NLAM.LT.2) WRITE(LUDBG,'(A)') '    0.0'
      WRITE(LUDBG,'(A)') '  ],'
C
C---- Also output transition x locations
      WRITE(LUDBG,'(A,E17.10,A)') '  "xtr_upper": ', XOCTR(1), ','
      WRITE(LUDBG,'(A,E17.10)') '  "xtr_lower": ', XOCTR(2)
C
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump shear-lag equation Jacobian intermediates (row 1)
C     Captures terms used in VS2(1,2..4) for comparison with Rust
      SUBROUTINE DBGSHEAR(IS, IBL, ITYP,
     &                    Z_UPW, Z_DE2, Z_US2, Z_CQ2, Z_CF2, Z_HK2,
     &                    Z_D2, Z_U2, Z_S2,
     &                    UPW_T2, UPW_D2, UPW_U2,
     &                    DE2_T2, DE2_D2, DE2_U2,
     &                    US2_T2, US2_D2, US2_U2,
     &                    CQ2_T2, CQ2_D2, CQ2_U2,
     &                    CF2_T2, CF2_D2, CF2_U2,
     &                    HK2_T2, HK2_D2, HK2_U2,
     &                    VS2_12, VS2_13, VS2_14)
      INTEGER IS, IBL, ITYP
      REAL Z_UPW, Z_DE2, Z_US2, Z_CQ2, Z_CF2, Z_HK2
      REAL Z_D2, Z_U2, Z_S2
      REAL UPW_T2, UPW_D2, UPW_U2
      REAL DE2_T2, DE2_D2, DE2_U2
      REAL US2_T2, US2_D2, US2_U2
      REAL CQ2_T2, CQ2_D2, CQ2_U2
      REAL CF2_T2, CF2_D2, CF2_U2
      REAL HK2_T2, HK2_D2, HK2_U2
      REAL VS2_12, VS2_13, VS2_14
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      COMMON /XDBGIT/ ITBLDBG
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER ITBLDBG
C
      IF(.NOT.LDBG) RETURN
C---- Only dump lower surface IBL 64-67 to keep output manageable
      IF(IS.NE.2) RETURN
      IF(IBL.LT.64 .OR. IBL.GT.67) RETURN
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "SHEAR_JACOBIAN",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,I2,A)') '  "side": ', IS, ','
      WRITE(LUDBG,'(A,I4,A)') '  "ibl": ', IBL, ','
      WRITE(LUDBG,'(A,I3,A)') '  "newton_iter": ', ITBLDBG, ','
      WRITE(LUDBG,'(A,I2,A)') '  "flow_type": ', ITYP, ','
C---- Z coefficients
      WRITE(LUDBG,'(A,E15.8,A)') '  "Z_UPW": ', Z_UPW, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Z_DE2": ', Z_DE2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Z_US2": ', Z_US2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Z_CQ2": ', Z_CQ2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Z_CF2": ', Z_CF2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Z_HK2": ', Z_HK2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Z_D2": ', Z_D2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Z_U2": ', Z_U2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "Z_S2": ', Z_S2, ','
C---- Derivatives used in VS2(1,2..4)
      WRITE(LUDBG,'(A,E15.8,A)') '  "UPW_T2": ', UPW_T2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "UPW_D2": ', UPW_D2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "UPW_U2": ', UPW_U2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "DE2_T2": ', DE2_T2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "DE2_D2": ', DE2_D2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "DE2_U2": ', DE2_U2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "US2_T2": ', US2_T2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "US2_D2": ', US2_D2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "US2_U2": ', US2_U2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "CQ2_T2": ', CQ2_T2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "CQ2_D2": ', CQ2_D2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "CQ2_U2": ', CQ2_U2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "CF2_T2": ', CF2_T2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "CF2_D2": ', CF2_D2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "CF2_U2": ', CF2_U2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "HK2_T2": ', HK2_T2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "HK2_D2": ', HK2_D2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "HK2_U2": ', HK2_U2, ','
C---- Final VS2 row 1 values
      WRITE(LUDBG,'(A,E15.8,A)') '  "VS2_1_2": ', VS2_12, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "VS2_1_3": ', VS2_13, ','
      WRITE(LUDBG,'(A,E15.8)') '  "VS2_1_4": ', VS2_14
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END

C---- Dump GGCALC (inviscid system) with geometry - alternate version
C     Captures AIJ matrix elements, RHS vectors, and GAMU solutions
      SUBROUTINE DBGGGCALC_GEOM(N, X, Y, NX, NY, AIJ, GAMU)
      INTEGER N
      REAL X(N), Y(N), NX(N), NY(N)
      REAL AIJ(600,600), GAMU(600,2)
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER I, J
C
      IF(.NOT.LDBG) RETURN
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      IDBGCALL = IDBGCALL + 1
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "GGCALC",'
      WRITE(LUDBG,'(A,I4,A)') '  "n_panels": ', N, ','
C
C---- Panel geometry (first 10)
      WRITE(LUDBG,'(A)') '  "panel_x_first10": ['
      DO I=1,MIN(10,N)
        IF(I.LT.MIN(10,N)) THEN
          WRITE(LUDBG,'(4X,E18.10,A)') X(I), ','
        ELSE
          WRITE(LUDBG,'(4X,E18.10)') X(I)
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ],'
C
      WRITE(LUDBG,'(A)') '  "panel_y_first10": ['
      DO I=1,MIN(10,N)
        IF(I.LT.MIN(10,N)) THEN
          WRITE(LUDBG,'(4X,E18.10,A)') Y(I), ','
        ELSE
          WRITE(LUDBG,'(4X,E18.10)') Y(I)
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ],'
C
      WRITE(LUDBG,'(A)') '  "panel_nx_first10": ['
      DO I=1,MIN(10,N)
        IF(I.LT.MIN(10,N)) THEN
          WRITE(LUDBG,'(4X,E18.10,A)') NX(I), ','
        ELSE
          WRITE(LUDBG,'(4X,E18.10)') NX(I)
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ],'
C
      WRITE(LUDBG,'(A)') '  "panel_ny_first10": ['
      DO I=1,MIN(10,N)
        IF(I.LT.MIN(10,N)) THEN
          WRITE(LUDBG,'(4X,E18.10,A)') NY(I), ','
        ELSE
          WRITE(LUDBG,'(4X,E18.10)') NY(I)
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ],'
C
C---- AIJ matrix: first row
      WRITE(LUDBG,'(A)') '  "aij_row0_first10": ['
      DO J=1,MIN(10,N)
        IF(J.LT.MIN(10,N)) THEN
          WRITE(LUDBG,'(4X,E18.10,A)') AIJ(1,J), ','
        ELSE
          WRITE(LUDBG,'(4X,E18.10)') AIJ(1,J)
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ],'
C
C---- AIJ diagonal
      WRITE(LUDBG,'(A)') '  "aij_diagonal_first10": ['
      DO I=1,MIN(10,N)
        IF(I.LT.MIN(10,N)) THEN
          WRITE(LUDBG,'(4X,E18.10,A)') AIJ(I,I), ','
        ELSE
          WRITE(LUDBG,'(4X,E18.10)') AIJ(I,I)
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ],'
C
C---- RHS (stored in GAMU before solve)
      WRITE(LUDBG,'(A)') '  "rhs_alpha0_first10": ['
      DO I=1,MIN(10,N)
        IF(I.LT.MIN(10,N)) THEN
          WRITE(LUDBG,'(4X,E18.10,A)') GAMU(I,1), ','
        ELSE
          WRITE(LUDBG,'(4X,E18.10)') GAMU(I,1)
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ],'
C
      WRITE(LUDBG,'(A)') '  "rhs_alpha90_first10": ['
      DO I=1,MIN(10,N)
        IF(I.LT.MIN(10,N)) THEN
          WRITE(LUDBG,'(4X,E18.10,A)') GAMU(I,2), ','
        ELSE
          WRITE(LUDBG,'(4X,E18.10)') GAMU(I,2)
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ]'
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END

C---- Dump GGCALC solution (after LU solve)
      SUBROUTINE DBGGGCALC_SOLN(N, GAMU)
      INTEGER N
      REAL GAMU(600,2)
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER I
C
      IF(.NOT.LDBG) RETURN
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      IDBGCALL = IDBGCALL + 1
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "GGCALC_SOLN",'
      WRITE(LUDBG,'(A,I4,A)') '  "n_panels": ', N, ','
C
C---- GAMU(:,1) - gamma for alpha=0
      WRITE(LUDBG,'(A)') '  "gamu_alpha0_first10": ['
      DO I=1,MIN(10,N)
        IF(I.LT.MIN(10,N)) THEN
          WRITE(LUDBG,'(4X,E18.10,A)') GAMU(I,1), ','
        ELSE
          WRITE(LUDBG,'(4X,E18.10)') GAMU(I,1)
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ],'
C
      WRITE(LUDBG,'(A)') '  "gamu_alpha0_last10": ['
      DO I=MAX(1,N-9),N
        IF(I.LT.N) THEN
          WRITE(LUDBG,'(4X,E18.10,A)') GAMU(I,1), ','
        ELSE
          WRITE(LUDBG,'(4X,E18.10)') GAMU(I,1)
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ],'
C
C---- GAMU(:,2) - gamma for alpha=90
      WRITE(LUDBG,'(A)') '  "gamu_alpha90_first10": ['
      DO I=1,MIN(10,N)
        IF(I.LT.MIN(10,N)) THEN
          WRITE(LUDBG,'(4X,E18.10,A)') GAMU(I,2), ','
        ELSE
          WRITE(LUDBG,'(4X,E18.10)') GAMU(I,2)
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ],'
C
      WRITE(LUDBG,'(A)') '  "gamu_alpha90_last10": ['
      DO I=MAX(1,N-9),N
        IF(I.LT.N) THEN
          WRITE(LUDBG,'(4X,E18.10,A)') GAMU(I,2), ','
        ELSE
          WRITE(LUDBG,'(4X,E18.10)') GAMU(I,2)
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ],'
C
C---- Kutta check
      WRITE(LUDBG,'(A,E18.10,A)') '  "kutta_alpha0": ',
     &  GAMU(1,1) + GAMU(N,1), ','
      WRITE(LUDBG,'(A,E18.10)') '  "kutta_alpha90": ',
     &  GAMU(1,2) + GAMU(N,2)
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END

C---- Dump final gamma at specific alpha - alternate version with CL
      SUBROUTINE DBGSPECAL_CL(N, ALFA, GAM, CL)
      INTEGER N
      REAL ALFA, CL
      REAL GAM(N)
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER I
C
      IF(.NOT.LDBG) RETURN
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      IDBGCALL = IDBGCALL + 1
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "SPECAL",'
      WRITE(LUDBG,'(A,I4,A)') '  "n_panels": ', N, ','
      WRITE(LUDBG,'(A,E18.10,A)') '  "alpha_rad": ', ALFA, ','
      WRITE(LUDBG,'(A,E18.10,A)') '  "alpha_deg": ', ALFA*57.29578, ','
      WRITE(LUDBG,'(A,E18.10,A)') '  "CL": ', CL, ','
C
C---- Gamma first 10
      WRITE(LUDBG,'(A)') '  "gamma_first10": ['
      DO I=1,MIN(10,N)
        IF(I.LT.MIN(10,N)) THEN
          WRITE(LUDBG,'(4X,E18.10,A)') GAM(I), ','
        ELSE
          WRITE(LUDBG,'(4X,E18.10)') GAM(I)
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ],'
C
C---- Gamma near LE (middle)
      WRITE(LUDBG,'(A)') '  "gamma_mid5": ['
      DO I=N/2-2, N/2+2
        IF(I.LT.N/2+2) THEN
          WRITE(LUDBG,'(4X,E18.10,A)') GAM(I), ','
        ELSE
          WRITE(LUDBG,'(4X,E18.10)') GAM(I)
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ],'
C
C---- Gamma last 10
      WRITE(LUDBG,'(A)') '  "gamma_last10": ['
      DO I=MAX(1,N-9),N
        IF(I.LT.N) THEN
          WRITE(LUDBG,'(4X,E18.10,A)') GAM(I), ','
        ELSE
          WRITE(LUDBG,'(4X,E18.10)') GAM(I)
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ],'
C
C---- Kutta check
      WRITE(LUDBG,'(A,E18.10)') '  "kutta_sum": ', GAM(1) + GAM(N)
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump AIJ influence matrix BEFORE LU decomposition
C     For comparison with RustFoil's influence coefficient matrix
      SUBROUTINE DBGAIJMATRIX(N, AIJ, IZX)
      INTEGER N, IZX
      REAL AIJ(IZX,IZX)
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER I, J, NOUT, ROW
C
      IF(.NOT.LDBG) RETURN
C
C---- Limit output to first 20 entries
      NOUT = MIN(N, 20)
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      IDBGCALL = IDBGCALL + 1
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "AIJ_MATRIX",'
      WRITE(LUDBG,'(A,I4,A)') '  "n_panels": ', N, ','
      WRITE(LUDBG,'(A)') '  "note": "AIJ before LU factorization",'
C
C---- Diagonal entries (first 20)
      WRITE(LUDBG,'(A)') '  "aij_diagonal_20": ['
      DO I=1, NOUT
        IF(I.LT.NOUT) THEN
          WRITE(LUDBG,'(4X,E18.10,A)') AIJ(I,I), ','
        ELSE
          WRITE(LUDBG,'(4X,E18.10)') AIJ(I,I)
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ],'
C
C---- First row (first 20 columns)
      WRITE(LUDBG,'(A)') '  "aij_row1_20cols": ['
      DO J=1, NOUT
        IF(J.LT.NOUT) THEN
          WRITE(LUDBG,'(4X,E18.10,A)') AIJ(1,J), ','
        ELSE
          WRITE(LUDBG,'(4X,E18.10)') AIJ(1,J)
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ],'
C
C---- Last row (first 20 columns)
      WRITE(LUDBG,'(A)') '  "aij_rowN_20cols": ['
      DO J=1, NOUT
        IF(J.LT.NOUT) THEN
          WRITE(LUDBG,'(4X,E18.10,A)') AIJ(N,J), ','
        ELSE
          WRITE(LUDBG,'(4X,E18.10)') AIJ(N,J)
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ],'
C
C---- Row 40 (if exists, first 20 columns)
      ROW = 40
      IF(ROW.LE.N) THEN
        WRITE(LUDBG,'(A)') '  "aij_row40_20cols": ['
        DO J=1, NOUT
          IF(J.LT.NOUT) THEN
            WRITE(LUDBG,'(4X,E18.10,A)') AIJ(ROW,J), ','
          ELSE
            WRITE(LUDBG,'(4X,E18.10)') AIJ(ROW,J)
          ENDIF
        ENDDO
        WRITE(LUDBG,'(A)') '  ],'
      ENDIF
C
C---- Row 80 (if exists, first 20 columns)
      ROW = 80
      IF(ROW.LE.N) THEN
        WRITE(LUDBG,'(A)') '  "aij_row80_20cols": ['
        DO J=1, NOUT
          IF(J.LT.NOUT) THEN
            WRITE(LUDBG,'(4X,E18.10,A)') AIJ(ROW,J), ','
          ELSE
            WRITE(LUDBG,'(4X,E18.10)') AIJ(ROW,J)
          ENDIF
        ENDDO
        WRITE(LUDBG,'(A)') '  ],'
      ENDIF
C
C---- Row 120 (if exists, first 20 columns)
      ROW = 120
      IF(ROW.LE.N) THEN
        WRITE(LUDBG,'(A)') '  "aij_row120_20cols": ['
        DO J=1, NOUT
          IF(J.LT.NOUT) THEN
            WRITE(LUDBG,'(4X,E18.10,A)') AIJ(ROW,J), ','
          ELSE
            WRITE(LUDBG,'(4X,E18.10)') AIJ(ROW,J)
          ENDIF
        ENDDO
        WRITE(LUDBG,'(A)') '  ],'
      ENDIF
C
C---- Kutta condition row (N+1, first 20 columns)
      WRITE(LUDBG,'(A)') '  "aij_kutta_row_20cols": ['
      DO J=1, NOUT
        IF(J.LT.NOUT) THEN
          WRITE(LUDBG,'(4X,E18.10,A)') AIJ(N+1,J), ','
        ELSE
          WRITE(LUDBG,'(4X,E18.10)') AIJ(N+1,J)
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ],'
C
C---- Also dump last few columns of row 1 to check boundary treatment
      WRITE(LUDBG,'(A)') '  "aij_row1_lastcols": ['
      DO J=MAX(1,N-9), N+1
        IF(J.LT.N+1) THEN
          WRITE(LUDBG,'(4X,E18.10,A)') AIJ(1,J), ','
        ELSE
          WRITE(LUDBG,'(4X,E18.10)') AIJ(1,J)
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ],'
C
C---- Stagnation region (around row 80, cols 75-85) if exists
      IF(N.GE.85) THEN
        WRITE(LUDBG,'(A)') '  "aij_stagnation_region": ['
        DO I=78, 82
          DO J=78, 82
            IF(I.NE.82 .OR. J.NE.82) THEN
              WRITE(LUDBG,'(A,I3,A,I3,A,E18.10,A)')
     &          '    {"i": ', I, ', "j": ', J, ', "val": ',
     &          AIJ(I,J), '},'
            ELSE
              WRITE(LUDBG,'(A,I3,A,I3,A,E18.10,A)')
     &          '    {"i": ', I, ', "j": ', J, ', "val": ',
     &          AIJ(I,J), '}'
            ENDIF
          ENDDO
        ENDDO
        WRITE(LUDBG,'(A)') '  ],'
      ENDIF
C
C---- Matrix sum for quick sanity check
      WRITE(LUDBG,'(A,E18.10)') '  "aij_trace_sum": ',
     &  AIJ(1,1) + AIJ(N/2,N/2) + AIJ(N,N)
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- Dump panel geometry (X, Y, NX, NY, S arrays) for RustFoil comparison
      SUBROUTINE DBGPANELGEOM()
      INCLUDE 'XFOIL.INC'
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER I
C
      IF(.NOT.LDBG) RETURN
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "PANEL_GEOM",'
      WRITE(LUDBG,'(A,I4,A)') '  "n": ', N, ','
C
C---- X coordinates - first 10
      WRITE(LUDBG,'(A)') '  "x_first10": ['
      DO I=1,MIN(10,N)
        IF(I.LT.MIN(10,N)) THEN
          WRITE(LUDBG,'(4X,E18.10,A)') X(I), ','
        ELSE
          WRITE(LUDBG,'(4X,E18.10)') X(I)
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ],'
C
C---- X coordinates - last 10
      WRITE(LUDBG,'(A)') '  "x_last10": ['
      DO I=MAX(1,N-9),N
        IF(I.LT.N) THEN
          WRITE(LUDBG,'(4X,E18.10,A)') X(I), ','
        ELSE
          WRITE(LUDBG,'(4X,E18.10)') X(I)
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ],'
C
C---- Y coordinates - first 10
      WRITE(LUDBG,'(A)') '  "y_first10": ['
      DO I=1,MIN(10,N)
        IF(I.LT.MIN(10,N)) THEN
          WRITE(LUDBG,'(4X,E18.10,A)') Y(I), ','
        ELSE
          WRITE(LUDBG,'(4X,E18.10)') Y(I)
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ],'
C
C---- Y coordinates - last 10
      WRITE(LUDBG,'(A)') '  "y_last10": ['
      DO I=MAX(1,N-9),N
        IF(I.LT.N) THEN
          WRITE(LUDBG,'(4X,E18.10,A)') Y(I), ','
        ELSE
          WRITE(LUDBG,'(4X,E18.10)') Y(I)
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ],'
C
C---- NX (normal x-component) - first 10
      WRITE(LUDBG,'(A)') '  "nx_first10": ['
      DO I=1,MIN(10,N)
        IF(I.LT.MIN(10,N)) THEN
          WRITE(LUDBG,'(4X,E18.10,A)') NX(I), ','
        ELSE
          WRITE(LUDBG,'(4X,E18.10)') NX(I)
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ],'
C
C---- NX (normal x-component) - last 10
      WRITE(LUDBG,'(A)') '  "nx_last10": ['
      DO I=MAX(1,N-9),N
        IF(I.LT.N) THEN
          WRITE(LUDBG,'(4X,E18.10,A)') NX(I), ','
        ELSE
          WRITE(LUDBG,'(4X,E18.10)') NX(I)
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ],'
C
C---- NY (normal y-component) - first 10
      WRITE(LUDBG,'(A)') '  "ny_first10": ['
      DO I=1,MIN(10,N)
        IF(I.LT.MIN(10,N)) THEN
          WRITE(LUDBG,'(4X,E18.10,A)') NY(I), ','
        ELSE
          WRITE(LUDBG,'(4X,E18.10)') NY(I)
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ],'
C
C---- NY (normal y-component) - last 10
      WRITE(LUDBG,'(A)') '  "ny_last10": ['
      DO I=MAX(1,N-9),N
        IF(I.LT.N) THEN
          WRITE(LUDBG,'(4X,E18.10,A)') NY(I), ','
        ELSE
          WRITE(LUDBG,'(4X,E18.10)') NY(I)
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ],'
C
C---- S (arc length) - first 10
      WRITE(LUDBG,'(A)') '  "s_first10": ['
      DO I=1,MIN(10,N)
        IF(I.LT.MIN(10,N)) THEN
          WRITE(LUDBG,'(4X,E18.10,A)') S(I), ','
        ELSE
          WRITE(LUDBG,'(4X,E18.10)') S(I)
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ],'
C
C---- S (arc length) - last 10
      WRITE(LUDBG,'(A)') '  "s_last10": ['
      DO I=MAX(1,N-9),N
        IF(I.LT.N) THEN
          WRITE(LUDBG,'(4X,E18.10,A)') S(I), ','
        ELSE
          WRITE(LUDBG,'(4X,E18.10)') S(I)
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ]'
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C=======================================================================
C     NEWTON ITERATION DEBUG SUBROUTINES
C     For detailed comparison with RustFoil Newton solver
C=======================================================================

C---- DBGNEWTONITER: Dump at start of each Newton iteration
C     Note: RMSBL, NSYS are in COMMON via XFOIL.INC
      SUBROUTINE DBGNEWTONITER(ITER)
      INTEGER ITER
      INCLUDE 'XFOIL.INC'
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
C---- Only log iterations 0-5 and 14-20 where divergence happens
      IF(ITER.GT.5 .AND. ITER.LT.14) RETURN
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "NEWTON_ITER_START",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', ITER, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "rmsbl_before": ', RMSBL, ','
      WRITE(LUDBG,'(A,I4,A)') '  "nsys": ', NSYS, ','
      WRITE(LUDBG,'(A,I4,A)') '  "nbl_top": ', NBL(1), ','
      WRITE(LUDBG,'(A,I4,A)') '  "nbl_bot": ', NBL(2), ','
      WRITE(LUDBG,'(A,I4,A)') '  "itran_top": ', ITRAN(1), ','
      WRITE(LUDBG,'(A,I4,A)') '  "itran_bot": ', ITRAN(2), ','
      WRITE(LUDBG,'(A,I4,A)') '  "iblte_top": ', IBLTE(1), ','
      WRITE(LUDBG,'(A,I4)') '  "iblte_bot": ', IBLTE(2)
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- DBGBLSTATE: Dump BL state at key stations
      SUBROUTINE DBGBLSTATE(ITER, IS, IBL)
      INTEGER ITER, IS, IBL
      INCLUDE 'XFOIL.INC'
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
C---- Only log iterations 0-5 and 14-20
      IF(ITER.GT.5 .AND. ITER.LT.14) RETURN
C---- Only log stations 10, 20, 40, 60, 80
      IF(MOD(IBL,20).NE.0 .AND. IBL.NE.10) RETURN
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "BL_STATE",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', ITER, ','
      WRITE(LUDBG,'(A,I2,A)') '  "is": ', IS, ','
      WRITE(LUDBG,'(A,I4,A)') '  "ibl": ', IBL, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "ctau": ', CTAU(IBL,IS), ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "thet": ', THET(IBL,IS), ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "dstr": ', DSTR(IBL,IS), ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "uedg": ', UEDG(IBL,IS), ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "mass": ', MASS(IBL,IS), ','
      WRITE(LUDBG,'(A,E15.8)') '  "xssi": ', XSSI(IBL,IS)
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- DBGBLDIFOUT: Dump after BLSYS returns in SETBL
      SUBROUTINE DBGBLDIFOUT(IBL, IS, VS1, VS2, VSREZ, ITYP)
      INTEGER IBL, IS, ITYP
      REAL VS1(4,5), VS2(4,5), VSREZ(4)
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER I, J
      CHARACTER*12 FLOWTYPE
C
      IF(.NOT.LDBG) RETURN
C---- Only log every 20th station to keep output manageable
      IF(MOD(IBL,20).NE.0 .AND. IBL.NE.10) RETURN
C---- Only log iterations 0-5 and 14-20
      IF(IDBGITER.GT.5 .AND. IDBGITER.LT.14) RETURN
C
C---- Determine flow type string
      IF(ITYP.EQ.1) THEN
        FLOWTYPE = 'LAMINAR'
      ELSE IF(ITYP.EQ.2) THEN
        FLOWTYPE = 'TURBULENT'
      ELSE IF(ITYP.EQ.3) THEN
        FLOWTYPE = 'WAKE'
      ELSE
        FLOWTYPE = 'UNKNOWN'
      ENDIF
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "BLSYS_OUTPUT",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,I2,A)') '  "is": ', IS, ','
      WRITE(LUDBG,'(A,I4,A)') '  "ibl": ', IBL, ','
      WRITE(LUDBG,'(A,A,A)') '  "flow_type": "', TRIM(FLOWTYPE), '",'
C---- VS1 matrix (3x5 for first 3 equations)
      WRITE(LUDBG,'(A)') '  "VS1": ['
      DO 10 I=1,3
        IF(I.LT.3) THEN
          WRITE(LUDBG,'(A,5(E14.7,A),A)')
     &      '    [', (VS1(I,J), ',', J=1,4), VS1(I,5), '],'
        ELSE
          WRITE(LUDBG,'(A,5(E14.7,A),A)')
     &      '    [', (VS1(I,J), ',', J=1,4), VS1(I,5), ']'
        ENDIF
   10 CONTINUE
      WRITE(LUDBG,'(A)') '  ],'
C---- VS2 matrix (3x5)
      WRITE(LUDBG,'(A)') '  "VS2": ['
      DO 20 I=1,3
        IF(I.LT.3) THEN
          WRITE(LUDBG,'(A,5(E14.7,A),A)')
     &      '    [', (VS2(I,J), ',', J=1,4), VS2(I,5), '],'
        ELSE
          WRITE(LUDBG,'(A,5(E14.7,A),A)')
     &      '    [', (VS2(I,J), ',', J=1,4), VS2(I,5), ']'
        ENDIF
   20 CONTINUE
      WRITE(LUDBG,'(A)') '  ],'
C---- Residuals (first 3)
      WRITE(LUDBG,'(A,3(E14.7,A),A)')
     &  '  "VSREZ": [', VSREZ(1), ',', VSREZ(2), ',', VSREZ(3), ']'
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- DBGVMBLOCK: Dump VM near-diagonal entries after assembly
C     Note: NSYS is in COMMON via XFOIL.INC
      SUBROUTINE DBGVMBLOCK(IV)
      INTEGER IV
      INCLUDE 'XFOIL.INC'
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER I, J, JMIN, JMAX
C
      IF(.NOT.LDBG) RETURN
C---- Only log every 20th station
      IF(MOD(IV,20).NE.0 .AND. IV.NE.10) RETURN
C---- Only log iterations 0-5 and 14-20
      IF(IDBGITER.GT.5 .AND. IDBGITER.LT.14) RETURN
C
      JMIN = MAX(1, IV-3)
      JMAX = MIN(NSYS, IV+3)
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "VM_BLOCK",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,I4,A)') '  "iv": ', IV, ','
      WRITE(LUDBG,'(A,I4,A)') '  "nsys": ', NSYS, ','
      WRITE(LUDBG,'(A,I4,A)') '  "jmin": ', JMIN, ','
      WRITE(LUDBG,'(A,I4,A)') '  "jmax": ', JMAX, ','
C---- VM row 1 (near diagonal)
      WRITE(LUDBG,'(A)') '  "VM_row1": ['
      DO J=JMIN,JMAX
        IF(J.LT.JMAX) THEN
          WRITE(LUDBG,'(A,E14.7,A)') '    ', VM(1,J,IV), ','
        ELSE
          WRITE(LUDBG,'(A,E14.7)') '    ', VM(1,J,IV)
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ],'
C---- VM row 2 (near diagonal)
      WRITE(LUDBG,'(A)') '  "VM_row2": ['
      DO J=JMIN,JMAX
        IF(J.LT.JMAX) THEN
          WRITE(LUDBG,'(A,E14.7,A)') '    ', VM(2,J,IV), ','
        ELSE
          WRITE(LUDBG,'(A,E14.7)') '    ', VM(2,J,IV)
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ],'
C---- VM row 3 (near diagonal)
      WRITE(LUDBG,'(A)') '  "VM_row3": ['
      DO J=JMIN,JMAX
        IF(J.LT.JMAX) THEN
          WRITE(LUDBG,'(A,E14.7,A)') '    ', VM(3,J,IV), ','
        ELSE
          WRITE(LUDBG,'(A,E14.7)') '    ', VM(3,J,IV)
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ]'
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- DBGVDEL: Dump VDEL (residual/RHS) after assembly
      SUBROUTINE DBGVDEL(IV, VDEL1, VDEL2, VDEL3)
      INTEGER IV
      REAL VDEL1, VDEL2, VDEL3
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
C---- Only log every 20th station
      IF(MOD(IV,20).NE.0 .AND. IV.NE.10) RETURN
C---- Only log iterations 0-5 and 14-20
      IF(IDBGITER.GT.5 .AND. IDBGITER.LT.14) RETURN
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "VDEL_RESIDUAL",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,I4,A)') '  "iv": ', IV, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "vdel1": ', VDEL1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "vdel2": ', VDEL2, ','
      WRITE(LUDBG,'(A,E15.8)') '  "vdel3": ', VDEL3
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- DBGVSREZ_DUE: Dump VSREZ base residual and DUE/DDS forced changes
C     This helps identify where the 2x factor discrepancy originates
      SUBROUTINE DBGVSREZ_DUE(IV, VSREZ1, VSREZ2, VSREZ3,
     &                        DUE1, DUE2, DDS1, DDS2)
      INTEGER IV
      REAL VSREZ1, VSREZ2, VSREZ3, DUE1, DUE2, DDS1, DDS2
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      COMMON /XDBGLOCK/ LDBGLOCK
      LOGICAL LDBG, LDBGLOCK
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
      IF(LDBGLOCK) RETURN
C---- Only log every 10th station (including station 0)
      IF(MOD(IV,10).NE.0) RETURN
C---- Only log iterations 1-5 (where we care about the 2x factor)
      IF(IDBGITER.LT.1 .OR. IDBGITER.GT.5) RETURN
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "VSREZ_DUE",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,I4,A)') '  "iv": ', IV, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "vsrez1": ', VSREZ1, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "vsrez2": ', VSREZ2, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "vsrez3": ', VSREZ3, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "due1": ', DUE1, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "due2": ', DUE2, ','
      WRITE(LUDBG,'(A,E17.10,A)') '  "dds1": ', DDS1, ','
      WRITE(LUDBG,'(A,E17.10)') '  "dds2": ', DDS2
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- DBGSOLUTION: Dump solution deltas after BLSOLV
C     Note: NSYS is in COMMON via XFOIL.INC
      SUBROUTINE DBGSOLUTION()
      INCLUDE 'XFOIL.INC'
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER IV, NOUT
C
      IF(.NOT.LDBG) RETURN
C---- Only log iterations 0-5 and 14-20
      IF(IDBGITER.GT.5 .AND. IDBGITER.LT.14) RETURN
C
      NOUT = MIN(20, NSYS)
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "BLSOLV_SOLUTION",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,I4,A)') '  "nsys": ', NSYS, ','
      WRITE(LUDBG,'(A,I4,A)') '  "nout": ', NOUT, ','
C---- Solution deltas (first NOUT stations)
      WRITE(LUDBG,'(A)') '  "vdel_solution": ['
      DO IV=1,NOUT
        IF(IV.LT.NOUT) THEN
          WRITE(LUDBG,'(A,I4,A,E14.7,A,E14.7,A,E14.7,A)')
     &      '    {"iv": ', IV, ', "d1": ', VDEL(1,1,IV),
     &      ', "d2": ', VDEL(2,1,IV), ', "d3": ', VDEL(3,1,IV), '},'
        ELSE
          WRITE(LUDBG,'(A,I4,A,E14.7,A,E14.7,A,E14.7,A)')
     &      '    {"iv": ', IV, ', "d1": ', VDEL(1,1,IV),
     &      ', "d2": ', VDEL(2,1,IV), ', "d3": ', VDEL(3,1,IV), '}'
        ENDIF
      ENDDO
      WRITE(LUDBG,'(A)') '  ]'
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- DBGUPDATE_NEWTON: Dump update info at each station in Newton loop
      SUBROUTINE DBGUPDATE_NEWTON(IBL, IS, RLX, DN1, DN2, DN3, DN4,
     &                            DCTAU, DTHET, DDSTR, DUEDG)
      INTEGER IBL, IS
      REAL RLX, DN1, DN2, DN3, DN4, DCTAU, DTHET, DDSTR, DUEDG
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
C---- Only log every 20th station
      IF(MOD(IBL,20).NE.0) RETURN
C---- Only log iterations 0-5 and 14-20
      IF(IDBGITER.GT.5 .AND. IDBGITER.LT.14) RETURN
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "UPDATE_NEWTON",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,I2,A)') '  "is": ', IS, ','
      WRITE(LUDBG,'(A,I4,A)') '  "ibl": ', IBL, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "rlx": ', RLX, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "dn1": ', DN1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "dn2": ', DN2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "dn3": ', DN3, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "dn4": ', DN4, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "dctau": ', DCTAU, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "dthet": ', DTHET, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "ddstr": ', DDSTR, ','
      WRITE(LUDBG,'(A,E15.8)') '  "duedg": ', DUEDG
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- DBGTRANSITION: Dump transition station details
C     Captures BL state at specific stations for comparison with RustFoil
      SUBROUTINE DBGTRANSITION(ITER, IS, IBL, ITYP)
      INTEGER ITER, IS, IBL, ITYP
      INCLUDE 'XFOIL.INC'
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER ILAMINAR
C
      IF(.NOT.LDBG) RETURN
C---- Include transition region (stations 54-70 on lower surface)
      IF(IS.NE.2 .OR. IBL.LT.54 .OR. IBL.GT.70) RETURN
C---- Only log iterations 0-10
      IF(ITER.GT.10) RETURN
C
      IF(IBL.LT.ITRAN(IS)) THEN
        ILAMINAR = 1
      ELSE
        ILAMINAR = 0
      ENDIF
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "TRANSITION_DEBUG",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', ITER, ','
      WRITE(LUDBG,'(A,I2,A)') '  "is": ', IS, ','
      WRITE(LUDBG,'(A,I4,A)') '  "ibl": ', IBL, ','
      WRITE(LUDBG,'(A,I2,A)') '  "ityp": ', ITYP, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "ctau": ', CTAU(IBL,IS), ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "thet": ', THET(IBL,IS), ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "dstr": ', DSTR(IBL,IS), ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "uedg": ', UEDG(IBL,IS), ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "mass": ', MASS(IBL,IS), ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "xssi": ', XSSI(IBL,IS), ','
      WRITE(LUDBG,'(A,I4,A)') '  "itran": ', ITRAN(IS), ','
      WRITE(LUDBG,'(A,I1)') '  "laminar": ', ILAMINAR
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END

C---- DBGVSREZ_TRANS: Log VSREZ at transition stations during global Newton
      SUBROUTINE DBGVSREZ_TRANS(IS, IBL, VSREZ, VS2)
      INTEGER IS, IBL
      REAL VSREZ(4), VS2(4,5)
      INCLUDE 'XFOIL.INC'
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
C---- Only log transition stations (ITRAN-2 to ITRAN+5) on lower surface
      IF(IS.NE.2) RETURN
      IF(IBL.LT.ITRAN(IS)-2 .OR. IBL.GT.ITRAN(IS)+5) RETURN
C---- Only log iterations 0-15
      IF(IDBGITER.GT.15) RETURN
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "VSREZ_TRANSITION",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,I2,A)') '  "is": ', IS, ','
      WRITE(LUDBG,'(A,I4,A)') '  "ibl": ', IBL, ','
      WRITE(LUDBG,'(A,I4,A)') '  "itran": ', ITRAN(IS), ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "vsrez_0": ', VSREZ(1), ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "vsrez_1": ', VSREZ(2), ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "vsrez_2": ', VSREZ(3), ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "vsrez_3": ', VSREZ(4), ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "ctau": ', CTAU(IBL,IS), ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "cq": ', CTQ(IBL,IS), ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "vs2_00": ', VS2(1,1), ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "vs2_01": ', VS2(1,2), ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "vs2_02": ', VS2(1,3), ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "vs2_03": ', VS2(1,4), ','
      WRITE(LUDBG,'(A,E15.8)') '  "vs2_04": ', VS2(1,5)
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C---- DBGREZC_TERMS: Dump individual REZC terms at transition station
      SUBROUTINE DBGREZC_TERMS(IS, IBL, SCC, CQA, SA, ALD, DXI, DEA,
     &                         SLOG, UQ, ULOG, REZC, S1, S2, USA)
      INTEGER IS, IBL
      REAL SCC, CQA, SA, ALD, DXI, DEA, SLOG, UQ, ULOG, REZC, S1, S2
      REAL USA
      INCLUDE 'XFOIL.INC'
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      REAL DUXCON
      PARAMETER (DUXCON = 0.01)
      REAL TERM1, TERM2, TERM3
C
      IF(.NOT.LDBG) RETURN
C---- Only log transition station (ITRAN) on lower surface
      IF(IS.NE.2) RETURN
      IF(IBL.NE.ITRAN(IS)) RETURN
C---- Only log iterations 0-3
      IF(IDBGITER.GT.3) RETURN
C
C---- Compute terms
      TERM1 = SCC*(CQA - SA*ALD)*DXI
      TERM2 = -DEA*2.0*SLOG
      TERM3 = DEA*2.0*(UQ*DXI - ULOG)*DUXCON
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "REZC_TERMS",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,I2,A)') '  "is": ', IS, ','
      WRITE(LUDBG,'(A,I4,A)') '  "ibl": ', IBL, ','
      WRITE(LUDBG,'(A,I4,A)') '  "itran": ', ITRAN(IS), ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "SCC": ', SCC, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "CQA": ', CQA, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "SA": ', SA, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "ALD": ', ALD, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "DXI": ', DXI, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "DEA": ', DEA, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "SLOG": ', SLOG, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "UQ": ', UQ, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "ULOG": ', ULOG, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "S1": ', S1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "S2": ', S2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "USA": ', USA, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "term1": ', TERM1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "term2": ', TERM2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "term3": ', TERM3, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "REZC": ', REZC, ','
      WRITE(LUDBG,'(A,E15.8)') '  "vsrez0": ', -REZC
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END

C---- DBGCONVERGE: Dump convergence info at end of Newton iteration
C     Note: RMSBL, RLX, LVCONV are in COMMON via XFOIL.INC
      SUBROUTINE DBGCONVERGE(ITER, LCONV)
      INTEGER ITER
      LOGICAL LCONV
      INCLUDE 'XFOIL.INC'
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER ICONV
C
      IF(.NOT.LDBG) RETURN
C---- Only log iterations 0-5 and 14-20
      IF(ITER.GT.5 .AND. ITER.LT.14) RETURN
C
      IF(LCONV) THEN
        ICONV = 1
      ELSE
        ICONV = 0
      ENDIF
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "NEWTON_CONVERGE",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', ITER, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "rmsbl": ', RMSBL, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "rlx": ', RLX, ','
      WRITE(LUDBG,'(A,I1,A)') '  "converged": ', ICONV, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "rmxbl": ', RMXBL, ','
      WRITE(LUDBG,'(A,I4,A)') '  "imxbl": ', IMXBL, ','
      WRITE(LUDBG,'(A,I2)') '  "ismxbl": ', ISMXBL
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END

C---- DBGXIULE: Dump stagnation point coupling terms
      SUBROUTINE DBGXIULE(IS, IBL, XI_ULE1, XI_ULE2, SST_GO, SST_GP,
     &                    DULE1, DULE2)
      INTEGER IS, IBL
      REAL XI_ULE1, XI_ULE2, SST_GO, SST_GP, DULE1, DULE2
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
C---- Only log at a few key stations per iteration (every 20th)
      IF(MOD(IBL, 20).NE.0) RETURN
C---- Only log first 3 iterations
      IF(IDBGITER.GT.3) RETURN
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "XI_ULE",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,I2,A)') '  "is": ', IS, ','
      WRITE(LUDBG,'(A,I4,A)') '  "ibl": ', IBL, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "XI_ULE1": ', XI_ULE1, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "XI_ULE2": ', XI_ULE2, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "SST_GO": ', SST_GO, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "SST_GP": ', SST_GP, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "DULE1": ', DULE1, ','
      WRITE(LUDBG,'(A,E15.8)') '  "DULE2": ', DULE2
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END

C---- DBGDUE: Dump DUE values at first few stations
      SUBROUTINE DBGDUE(IS, IBL, UEDG_VAL, USAV_VAL, DUE_VAL, 
     &                   UINV_VAL)
      INTEGER IS, IBL
      REAL UEDG_VAL, USAV_VAL, DUE_VAL, UINV_VAL
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      COMMON /XDBGLOCK/ LDBGLOCK
      LOGICAL LDBG, LDBGLOCK
      INTEGER LUDBG, IDBGCALL, IDBGITER
      REAL MASS_CONTRIB
C
      IF(.NOT.LDBG) RETURN
      IF(LDBGLOCK) RETURN
C---- Log every 10th station, aligned with RustFoil (ibl=11,21,31... -> RF ibl=10,20,30)
      IF(MOD(IBL-1,10).NE.0) RETURN
C---- Only log iterations 1-5
      IF(IDBGITER.LT.1 .OR. IDBGITER.GT.5) RETURN
C
C---- Compute mass contribution = USAV - UINV
      MASS_CONTRIB = USAV_VAL - UINV_VAL
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "DUE",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,I2,A)') '  "is": ', IS, ','
      WRITE(LUDBG,'(A,I4,A)') '  "ibl": ', IBL, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "uedg": ', UEDG_VAL, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "usav": ', USAV_VAL, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "uinv": ', UINV_VAL, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "mass_contrib": ', MASS_CONTRIB, ','
      WRITE(LUDBG,'(A,E15.8)') '  "due": ', DUE_VAL
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C***********************************************************************
C    DBGFULLITER: Full Newton iteration dump for brute-force comparison
C    
C    Dumps comprehensive iteration state for systematic comparison
C    with RustFoil. Called after BLSOLV solution is computed.
C    
C    Includes:
C    - Global iteration state (RMSBL, transition, system size)
C    - Full BL state at every 10th station
C    - Full system matrices (all stations, not sampled)
C    - Solution deltas for all stations
C
C    Note: Uses NSYS, VA, VB, VM, VDEL from XFOIL.INC COMMON blocks
C***********************************************************************
      SUBROUTINE DBGFULLITER(ITER)
      INTEGER ITER
      INCLUDE 'XFOIL.INC'
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER IS, IBL, IV, K, J
      REAL HK_VAL
C
      IF(.NOT.LDBG) RETURN
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "FULL_ITER",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', ITER, ','
C
C---- Global state
      WRITE(LUDBG,'(A,E15.8,A)') '  "rmsbl": ', RMSBL, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "rmxbl": ', RMXBL, ','
      WRITE(LUDBG,'(A,I4,A)') '  "imxbl": ', IMXBL, ','
      WRITE(LUDBG,'(A,I2,A)') '  "ismxbl": ', ISMXBL, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "rlx": ', RLX, ','
      WRITE(LUDBG,'(A,I4,A)') '  "nsys": ', NSYS, ','
      WRITE(LUDBG,'(A,I4,A)') '  "nbl_upper": ', NBL(1), ','
      WRITE(LUDBG,'(A,I4,A)') '  "nbl_lower": ', NBL(2), ','
      WRITE(LUDBG,'(A,I4,A)') '  "iblte_upper": ', IBLTE(1), ','
      WRITE(LUDBG,'(A,I4,A)') '  "iblte_lower": ', IBLTE(2), ','
      WRITE(LUDBG,'(A,I4,A)') '  "itran_upper": ', ITRAN(1), ','
      WRITE(LUDBG,'(A,I4,A)') '  "itran_lower": ', ITRAN(2), ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "xtr_upper": ', XOCTR(1), ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "xtr_lower": ', XOCTR(2), ','
C
C---- Full BL state for upper surface (every 10th station)
      WRITE(LUDBG,'(A)') '  "bl_upper": ['
      DO 100 IBL=2, NBL(1), 10
        IF(IBL.GT.2) WRITE(LUDBG,'(A)') ','
        HK_VAL = DSTR(IBL,1) / THET(IBL,1)
        WRITE(LUDBG,'(A)') '    {'
        WRITE(LUDBG,'(A,I4,A)') '      "ibl": ', IBL, ','
        WRITE(LUDBG,'(A,E15.8,A)') '      "x": ', XSSI(IBL,1), ','
        WRITE(LUDBG,'(A,E15.8,A)') '      "theta": ', THET(IBL,1), ','
        WRITE(LUDBG,'(A,E15.8,A)') '      "dstar": ', DSTR(IBL,1), ','
        WRITE(LUDBG,'(A,E15.8,A)') '      "ue": ', UEDG(IBL,1), ','
        WRITE(LUDBG,'(A,E15.8,A)') '      "mass": ', MASS(IBL,1), ','
        WRITE(LUDBG,'(A,E15.8,A)') '      "ctau": ', CTAU(IBL,1), ','
        WRITE(LUDBG,'(A,E15.8,A)') '      "hk": ', HK_VAL, ','
        WRITE(LUDBG,'(A,E15.8)') '      "cf": ', TAU(IBL,1)
        WRITE(LUDBG,'(A,$)') '    }'
  100 CONTINUE
      WRITE(LUDBG,'(A)') ''
      WRITE(LUDBG,'(A)') '  ],'
C
C---- Full BL state for lower surface (every 10th station)
      WRITE(LUDBG,'(A)') '  "bl_lower": ['
      DO 110 IBL=2, NBL(2), 10
        IF(IBL.GT.2) WRITE(LUDBG,'(A)') ','
        HK_VAL = DSTR(IBL,2) / THET(IBL,2)
        WRITE(LUDBG,'(A)') '    {'
        WRITE(LUDBG,'(A,I4,A)') '      "ibl": ', IBL, ','
        WRITE(LUDBG,'(A,E15.8,A)') '      "x": ', XSSI(IBL,2), ','
        WRITE(LUDBG,'(A,E15.8,A)') '      "theta": ', THET(IBL,2), ','
        WRITE(LUDBG,'(A,E15.8,A)') '      "dstar": ', DSTR(IBL,2), ','
        WRITE(LUDBG,'(A,E15.8,A)') '      "ue": ', UEDG(IBL,2), ','
        WRITE(LUDBG,'(A,E15.8,A)') '      "mass": ', MASS(IBL,2), ','
        WRITE(LUDBG,'(A,E15.8,A)') '      "ctau": ', CTAU(IBL,2), ','
        WRITE(LUDBG,'(A,E15.8,A)') '      "hk": ', HK_VAL, ','
        WRITE(LUDBG,'(A,E15.8)') '      "cf": ', TAU(IBL,2)
        WRITE(LUDBG,'(A,$)') '    }'
  110 CONTINUE
      WRITE(LUDBG,'(A)') ''
      WRITE(LUDBG,'(A)') '  ],'
C
C---- Full VA blocks (all stations)
      WRITE(LUDBG,'(A)') '  "va_full": ['
      DO 200 IV=1, NSYS
        IF(IV.GT.1) WRITE(LUDBG,'(A)') ','
        WRITE(LUDBG,'(A,$)') '    ['
        DO 210 K=1, 3
          IF(K.GT.1) WRITE(LUDBG,'(A,$)') ', '
          WRITE(LUDBG,'(A,E14.7,A,E14.7,A,$)')
     &      '[', VA(K,1,IV), ', ', VA(K,2,IV), ']'
  210   CONTINUE
        WRITE(LUDBG,'(A,$)') ']'
  200 CONTINUE
      WRITE(LUDBG,'(A)') ''
      WRITE(LUDBG,'(A)') '  ],'
C
C---- Full VB blocks (all stations)
      WRITE(LUDBG,'(A)') '  "vb_full": ['
      DO 300 IV=1, NSYS
        IF(IV.GT.1) WRITE(LUDBG,'(A)') ','
        WRITE(LUDBG,'(A,$)') '    ['
        DO 310 K=1, 3
          IF(K.GT.1) WRITE(LUDBG,'(A,$)') ', '
          WRITE(LUDBG,'(A,E14.7,A,E14.7,A,$)')
     &      '[', VB(K,1,IV), ', ', VB(K,2,IV), ']'
  310   CONTINUE
        WRITE(LUDBG,'(A,$)') ']'
  300 CONTINUE
      WRITE(LUDBG,'(A)') ''
      WRITE(LUDBG,'(A)') '  ],'
C
C---- Full VDEL (residuals before solve - first column only)
      WRITE(LUDBG,'(A)') '  "vdel_full": ['
      DO 400 IV=1, NSYS
        IF(IV.GT.1) WRITE(LUDBG,'(A)') ','
        WRITE(LUDBG,'(A,E14.7,A,E14.7,A,E14.7,A,$)')
     &    '    [', VDEL(1,1,IV), ', ', VDEL(2,1,IV), 
     &    ', ', VDEL(3,1,IV), ']'
  400 CONTINUE
      WRITE(LUDBG,'(A)') ''
      WRITE(LUDBG,'(A)') '  ],'
C
C---- Full solution deltas (all stations)
      WRITE(LUDBG,'(A)') '  "deltas_full": ['
      DO 500 IV=1, NSYS
        IF(IV.GT.1) WRITE(LUDBG,'(A)') ','
        WRITE(LUDBG,'(A,E14.7,A,E14.7,A,E14.7,A,$)')
     &    '    [', VDEL(1,1,IV), ', ', VDEL(2,1,IV), 
     &    ', ', VDEL(3,1,IV), ']'
  500 CONTINUE
      WRITE(LUDBG,'(A)') ''
      WRITE(LUDBG,'(A)') '  ],'
C
C---- VM diagonal entries (mass coupling self-influence)
      WRITE(LUDBG,'(A)') '  "vm_diag": ['
      DO 600 IV=1, NSYS
        IF(IV.GT.1) WRITE(LUDBG,'(A)') ','
        WRITE(LUDBG,'(A,E14.7,A,E14.7,A,E14.7,A,$)')
     &    '    [', VM(1,IV,IV), ', ', VM(2,IV,IV), 
     &    ', ', VM(3,IV,IV), ']'
  600 CONTINUE
      WRITE(LUDBG,'(A)') ''
      WRITE(LUDBG,'(A)') '  ],'
C
C---- Force results
      WRITE(LUDBG,'(A,E15.8,A)') '  "cl": ', CL, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "cd": ', CD, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "cm": ', CM, ','
      WRITE(LUDBG,'(A,E15.8,A)') '  "cdf": ', CDF, ','
      WRITE(LUDBG,'(A,E15.8)') '  "cdp": ', CD-CDF
      WRITE(LUDBG,'(A)') '}'
C
      RETURN
      END


C***********************************************************************
C    DBGFULLSYSTEM: Dump complete Newton system (all stations)
C    
C    Enhanced version of DBGSETBLSYSTEM without station limits.
C    Dumps VA, VB, VDEL for ALL NSYS stations.
C***********************************************************************
      SUBROUTINE DBGFULLSYSTEM(NSYS, VA, VB, VM, VDEL, IVX, IZX)
      INTEGER NSYS, IVX, IZX
      REAL VA(3,2,IVX), VB(3,2,IVX)
      REAL VM(3,IZX,IVX), VDEL(3,2,IVX)
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER IV, K
C
      IF(.NOT.LDBG) RETURN
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "FULL_SYSTEM",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,I4,A)') '  "nsys": ', NSYS, ','
C
C---- VA blocks (all stations)
      WRITE(LUDBG,'(A)') '  "VA": ['
      DO 10 IV=1, NSYS
        IF(IV.GT.1) WRITE(LUDBG,'(A)') ','
        WRITE(LUDBG,'(A,$)') '    ['
        DO 15 K=1, 3
          IF(K.GT.1) WRITE(LUDBG,'(A,$)') ', '
          WRITE(LUDBG,'(A,E14.7,A,E14.7,A,$)')
     &      '[', VA(K,1,IV), ', ', VA(K,2,IV), ']'
   15   CONTINUE
        WRITE(LUDBG,'(A,$)') ']'
   10 CONTINUE
      WRITE(LUDBG,'(A)') ''
      WRITE(LUDBG,'(A)') '  ],'
C
C---- VB blocks (all stations)
      WRITE(LUDBG,'(A)') '  "VB": ['
      DO 20 IV=1, NSYS
        IF(IV.GT.1) WRITE(LUDBG,'(A)') ','
        WRITE(LUDBG,'(A,$)') '    ['
        DO 25 K=1, 3
          IF(K.GT.1) WRITE(LUDBG,'(A,$)') ', '
          WRITE(LUDBG,'(A,E14.7,A,E14.7,A,$)')
     &      '[', VB(K,1,IV), ', ', VB(K,2,IV), ']'
   25   CONTINUE
        WRITE(LUDBG,'(A,$)') ']'
   20 CONTINUE
      WRITE(LUDBG,'(A)') ''
      WRITE(LUDBG,'(A)') '  ],'
C
C---- VDEL (all stations)
      WRITE(LUDBG,'(A)') '  "VDEL": ['
      DO 30 IV=1, NSYS
        IF(IV.GT.1) WRITE(LUDBG,'(A)') ','
        WRITE(LUDBG,'(A,E14.7,A,E14.7,A,E14.7,A,$)')
     &    '    [', VDEL(1,1,IV), ', ', VDEL(2,1,IV), 
     &    ', ', VDEL(3,1,IV), ']'
   30 CONTINUE
      WRITE(LUDBG,'(A)') ''
      WRITE(LUDBG,'(A)') '  ]'
C
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END


C***********************************************************************
C    DBGFULLSOLUTION: Dump complete solution vector (all stations)
C    
C    Enhanced version of DBGBLSOLVSOLUTION without station limits.
C***********************************************************************
      SUBROUTINE DBGFULLSOLUTION(NSYS, VDEL, IVX)
      INTEGER NSYS, IVX
      REAL VDEL(3,2,IVX)
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      INTEGER IV
C
      IF(.NOT.LDBG) RETURN
C
      CALL DBGCOMMA()
      WRITE(LUDBG,'(A)') '{'
      WRITE(LUDBG,'(A,I6,A)') '  "call_id": ', IDBGCALL, ','
      WRITE(LUDBG,'(A)') '  "subroutine": "FULL_SOLUTION",'
      WRITE(LUDBG,'(A,I4,A)') '  "iteration": ', IDBGITER, ','
      WRITE(LUDBG,'(A,I4,A)') '  "nsys": ', NSYS, ','
C
C---- Solution deltas [dCtau, dTheta, dMass] for ALL stations
      WRITE(LUDBG,'(A)') '  "deltas": ['
      DO 10 IV=1, NSYS
        IF(IV.GT.1) WRITE(LUDBG,'(A)') ','
        WRITE(LUDBG,'(A,E14.7,A,E14.7,A,E14.7,A,$)')
     &    '    [', VDEL(1,1,IV), ', ', VDEL(2,1,IV), 
     &    ', ', VDEL(3,1,IV), ']'
   10 CONTINUE
      WRITE(LUDBG,'(A)') ''
      WRITE(LUDBG,'(A)') '  ]'
C
      WRITE(LUDBG,'(A)') '}'
      RETURN
      END
