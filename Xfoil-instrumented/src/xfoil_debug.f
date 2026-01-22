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
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
      DATA LDBG /.TRUE./
      DATA LUDBG /77/
      DATA IDBGCALL /0/
      DATA IDBGITER /0/
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
      SUBROUTINE DBGBLDIF(ISIDE, IBL, ITYP, VS1, VS2, VSREZ)
      INTEGER ISIDE, IBL, ITYP
      REAL VS1(4,5), VS2(4,5), VSREZ(4)
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
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
      WRITE(LUDBG,'(A,I2,A)') '  "flow_type": ', ITYP, ','
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


C---- Dump closure function HSL result
      SUBROUTINE DBGHSL(HK, RT, MSQ, HS, HS_HK, HS_RT, HS_MSQ)
      REAL HK, RT, MSQ, HS, HS_HK, HS_RT, HS_MSQ
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
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
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
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
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
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


C---- Dump HKIN result
      SUBROUTINE DBGHKIN(H, MSQ, HK, HK_H, HK_MSQ)
      REAL H, MSQ, HK, HK_H, HK_MSQ
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
      IF(.NOT.LDBG) RETURN
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
