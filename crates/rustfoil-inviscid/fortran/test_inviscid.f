C***********************************************************************
C  test_inviscid.f
C  
C  FORTRAN test harness for inviscid solver validation.
C  Generates NACA 0012 airfoil, runs inviscid analysis, and outputs
C  all intermediate values as JSON for comparison with Rust implementation.
C
C  Usage: ./test_inviscid [geometry|psilin|ggcalc|clcalc|all]
C
C  Links against Xfoil-instrumented object files.
C***********************************************************************

      PROGRAM TEST_INVISCID
      INCLUDE '../../../Xfoil-instrumented/src/XFOIL.INC'
C
      CHARACTER*32 ARG
      LOGICAL DO_GEOM, DO_PSILIN, DO_GGCALC, DO_CLCALC
C
C---- XFOIL debug common block
      COMMON /XDEBUG/ LDBG, LUDBG, IDBGCALL, IDBGITER
      LOGICAL LDBG
      INTEGER LUDBG, IDBGCALL, IDBGITER
C
C---- Initialize debug output
      LDBG = .TRUE.
      LUDBG = 6
      IDBGCALL = 0
      IDBGITER = 0
C
C---- Parse command line argument
      CALL GETARG(1, ARG)
      DO_GEOM = .FALSE.
      DO_PSILIN = .FALSE.
      DO_GGCALC = .FALSE.
      DO_CLCALC = .FALSE.
C
      IF(ARG.EQ.'geometry') THEN
        DO_GEOM = .TRUE.
      ELSE IF(ARG.EQ.'psilin') THEN
        DO_PSILIN = .TRUE.
      ELSE IF(ARG.EQ.'ggcalc') THEN
        DO_GGCALC = .TRUE.
      ELSE IF(ARG.EQ.'clcalc') THEN
        DO_CLCALC = .TRUE.
      ELSE
C---- Default: all
        DO_GEOM = .TRUE.
        DO_PSILIN = .TRUE.
        DO_GGCALC = .TRUE.
        DO_CLCALC = .TRUE.
      ENDIF
C
C---- Initialize XFOIL variables
      CALL INIT_XFOIL
C
C---- Generate NACA 0012 with 160 panels
      CALL GEN_NACA0012(160)
C
C---- Set up panels
      CALL PANGEN(.FALSE.)
C
C---- Start JSON output
      WRITE(*,'(A)') '{'
      WRITE(*,'(A)') '  "test": "inviscid_solver",'
      WRITE(*,'(A,I4,A)') '  "n_panels": ', N, ','
      WRITE(*,'(A)') '  "airfoil": "NACA0012",'
      WRITE(*,'(A)') '  "data": ['
C
C---- Output geometry data
      IF(DO_GEOM) THEN
        CALL OUTPUT_GEOMETRY
      ENDIF
C
C---- Build influence matrix (GGCALC)
      LGAMU = .FALSE.
      CALL GGCALC
      LGAMU = .TRUE.
C
C---- Output PSILIN samples
      IF(DO_PSILIN) THEN
        CALL OUTPUT_PSILIN_SAMPLES
      ENDIF
C
C---- Output GGCALC results  
      IF(DO_GGCALC) THEN
        CALL OUTPUT_GGCALC
      ENDIF
C
C---- Solve at multiple angles of attack
      IF(DO_CLCALC) THEN
        CALL OUTPUT_CLCALC(0.0)
        CALL OUTPUT_CLCALC(4.0)
        CALL OUTPUT_CLCALC(8.0)
      ENDIF
C
C---- Close JSON
      WRITE(*,'(A)') '  ]'
      WRITE(*,'(A)') '}'
C
      END


C***********************************************************************
      SUBROUTINE INIT_XFOIL
C---- Initialize XFOIL state variables
      INCLUDE '../../../Xfoil-instrumented/src/XFOIL.INC'
C
C---- Initialize various flags
      LVISC = .FALSE.
      LALFA = .TRUE.
      LWAKE = .FALSE.
      LPACC = .FALSE.
      LBLINI = .FALSE.
      LIPAN = .FALSE.
      LQAIJ = .FALSE.
      LADIJ = .FALSE.
      LWDIJ = .FALSE.
      LGAMU = .FALSE.
      LQINU = .FALSE.
      SHARP = .FALSE.
C
C---- Freestream conditions
      QINF = 1.0
      MINF = 0.0
      MINF1 = 0.0
      ALFA = 0.0
      COSA = 1.0
      SINA = 0.0
C
C---- Reference point for moment
      XCMREF = 0.25
      YCMREF = 0.0
C
      RETURN
      END


C***********************************************************************
      SUBROUTINE GEN_NACA0012(NPANEL)
C---- Generate NACA 0012 coordinates
      INCLUDE '../../../Xfoil-instrumented/src/XFOIL.INC'
      INTEGER NPANEL
C
      REAL PI, T, BETA, XX, YT
      INTEGER I, NHALF
C
      PI = 4.0*ATAN(1.0)
      T = 0.12
      NHALF = NPANEL/2
C
C---- Upper surface: TE to LE (x decreasing)
      N = 0
      DO 10 I=0, NHALF
        BETA = PI * FLOAT(I) / FLOAT(NHALF)
        XX = 0.5 * (1.0 - COS(BETA))
C---- NACA 4-digit thickness
        YT = 5.0*T*(0.2969*SQRT(XX) - 0.126*XX - 0.3516*XX**2
     &            + 0.2843*XX**3 - 0.1036*XX**4)
        N = N + 1
        X(N) = XX
        Y(N) = YT
   10 CONTINUE
C
C---- Lower surface: LE+1 to TE (x increasing, skip LE)
      DO 20 I=1, NHALF
        BETA = PI * FLOAT(NHALF-I) / FLOAT(NHALF)
        XX = 0.5 * (1.0 - COS(BETA))
        YT = 5.0*T*(0.2969*SQRT(XX) - 0.126*XX - 0.3516*XX**2
     &            + 0.2843*XX**3 - 0.1036*XX**4)
        N = N + 1
        X(N) = XX
        Y(N) = -YT
   20 CONTINUE
C
C---- Set name
      NAME = 'NACA 0012'
C
      RETURN
      END


C***********************************************************************
      SUBROUTINE OUTPUT_GEOMETRY
C---- Output geometry reference data as JSON
      INCLUDE '../../../Xfoil-instrumented/src/XFOIL.INC'
      INTEGER I
C
      WRITE(*,'(A)') '    {'
      WRITE(*,'(A)') '      "type": "geometry",'
      WRITE(*,'(A,I4,A)') '      "n": ', N, ','
C
C---- X coordinates
      WRITE(*,'(A)') '      "X": ['
      DO 10 I=1, N
        IF(I.LT.N) THEN
          WRITE(*,'(A,E17.10,A)') '        ', X(I), ','
        ELSE
          WRITE(*,'(A,E17.10)') '        ', X(I)
        ENDIF
   10 CONTINUE
      WRITE(*,'(A)') '      ],'
C
C---- Y coordinates
      WRITE(*,'(A)') '      "Y": ['
      DO 20 I=1, N
        IF(I.LT.N) THEN
          WRITE(*,'(A,E17.10,A)') '        ', Y(I), ','
        ELSE
          WRITE(*,'(A,E17.10)') '        ', Y(I)
        ENDIF
   20 CONTINUE
      WRITE(*,'(A)') '      ],'
C
C---- Arc length S
      WRITE(*,'(A)') '      "S": ['
      DO 30 I=1, N
        IF(I.LT.N) THEN
          WRITE(*,'(A,E17.10,A)') '        ', S(I), ','
        ELSE
          WRITE(*,'(A,E17.10)') '        ', S(I)
        ENDIF
   30 CONTINUE
      WRITE(*,'(A)') '      ],'
C
C---- NX (normal x-component)
      WRITE(*,'(A)') '      "NX": ['
      DO 40 I=1, N
        IF(I.LT.N) THEN
          WRITE(*,'(A,E17.10,A)') '        ', NX(I), ','
        ELSE
          WRITE(*,'(A,E17.10)') '        ', NX(I)
        ENDIF
   40 CONTINUE
      WRITE(*,'(A)') '      ],'
C
C---- NY (normal y-component)
      WRITE(*,'(A)') '      "NY": ['
      DO 50 I=1, N
        IF(I.LT.N) THEN
          WRITE(*,'(A,E17.10,A)') '        ', NY(I), ','
        ELSE
          WRITE(*,'(A,E17.10)') '        ', NY(I)
        ENDIF
   50 CONTINUE
      WRITE(*,'(A)') '      ],'
C
C---- APANEL (panel angles)
      WRITE(*,'(A)') '      "APANEL": ['
      DO 60 I=1, N
        IF(I.LT.N) THEN
          WRITE(*,'(A,E17.10,A)') '        ', APANEL(I), ','
        ELSE
          WRITE(*,'(A,E17.10)') '        ', APANEL(I)
        ENDIF
   60 CONTINUE
      WRITE(*,'(A)') '      ],'
C
C---- XP (dX/dS)
      WRITE(*,'(A)') '      "XP": ['
      DO 70 I=1, N
        IF(I.LT.N) THEN
          WRITE(*,'(A,E17.10,A)') '        ', XP(I), ','
        ELSE
          WRITE(*,'(A,E17.10)') '        ', XP(I)
        ENDIF
   70 CONTINUE
      WRITE(*,'(A)') '      ],'
C
C---- YP (dY/dS)
      WRITE(*,'(A)') '      "YP": ['
      DO 80 I=1, N
        IF(I.LT.N) THEN
          WRITE(*,'(A,E17.10,A)') '        ', YP(I), ','
        ELSE
          WRITE(*,'(A,E17.10)') '        ', YP(I)
        ENDIF
   80 CONTINUE
      WRITE(*,'(A)') '      ],'
C
C---- TE geometry
      WRITE(*,'(A,E17.10,A)') '      "ANTE": ', ANTE, ','
      WRITE(*,'(A,E17.10,A)') '      "ASTE": ', ASTE, ','
      WRITE(*,'(A,E17.10,A)') '      "DSTE": ', DSTE, ','
      IF(SHARP) THEN
        WRITE(*,'(A)') '      "SHARP": true,'
      ELSE
        WRITE(*,'(A)') '      "SHARP": false,'
      ENDIF
C
C---- LE geometry
      WRITE(*,'(A,E17.10,A)') '      "XLE": ', XLE, ','
      WRITE(*,'(A,E17.10,A)') '      "YLE": ', YLE, ','
      WRITE(*,'(A,E17.10,A)') '      "SLE": ', SLE, ','
C
C---- Chord
      WRITE(*,'(A,E17.10)') '      "CHORD": ', CHORD
C
      WRITE(*,'(A)') '    },'
C
      RETURN
      END


C***********************************************************************
      SUBROUTINE OUTPUT_PSILIN_SAMPLES
C---- Output selected PSILIN influence coefficients
C     We compute DZDG for a few representative (I, JO) pairs
      INCLUDE '../../../Xfoil-instrumented/src/XFOIL.INC'
C
      INTEGER I, JO, JP, K
      REAL SX, SY, DSO, DSIO
      REAL RX1, RY1, RX2, RY2
      REAL X1, X2, YY, RS1, RS2
      REAL G1, G2, T1, T2
      REAL PSIS, PSID, DZDGJO, DZDGJP
      REAL QOPI
C
      QOPI = 0.25 / (4.0*ATAN(1.0))
C
      WRITE(*,'(A)') '    {'
      WRITE(*,'(A)') '      "type": "psilin",'
      WRITE(*,'(A)') '      "samples": ['
C
C---- Sample (I, JO) pairs: regular panels, TE panel, near-LE
      K = 0
C
C---- Loop through selected field points
      DO 100 I=1, N, MAX(N/10,1)
C---- Loop through selected panels
      DO 90 JO=1, N-1, MAX(N/10,1)
        JP = JO + 1
        IF(JP.GT.N) JP = 1
C
C---- Panel geometry
        SX = X(JP) - X(JO)
        SY = Y(JP) - Y(JO)
        DSO = SQRT(SX**2 + SY**2)
        IF(DSO.LT.1.0E-12) GO TO 90
        DSIO = 1.0/DSO
        SX = SX * DSIO
        SY = SY * DSIO
C
C---- Field point vectors
        RX1 = X(I) - X(JO)
        RY1 = Y(I) - Y(JO)
        RX2 = X(I) - X(JP)
        RY2 = Y(I) - Y(JP)
C
C---- Local coordinates
        X1 = SX*RX1 + SY*RY1
        X2 = SX*RX2 + SY*RY2
        YY = SX*RY1 - SY*RX1
C
        RS1 = RX1**2 + RY1**2
        RS2 = RX2**2 + RY2**2
C
C---- Log/atan terms with singularity handling
        IF(I.NE.JO .AND. RS1.GT.1.0E-20) THEN
          G1 = LOG(RS1)
          T1 = ATAN2(X1, YY)
        ELSE
          G1 = 0.0
          T1 = 0.0
        ENDIF
C
        IF(I.NE.JP .AND. RS2.GT.1.0E-20) THEN
          G2 = LOG(RS2)
          T2 = ATAN2(X2, YY)
        ELSE
          G2 = 0.0
          T2 = 0.0
        ENDIF
C
C---- PSIS and PSID
        PSIS = 0.5*X1*G1 - 0.5*X2*G2 + X2 - X1 + YY*(T1-T2)
C
        IF(ABS(X1-X2).GT.1.0E-20) THEN
          PSID = ((X1+X2)*PSIS + 0.5*(RS2*G2 - RS1*G1 + X1**2 - X2**2))
     &           / (X1-X2)
        ELSE
          PSID = 0.0
        ENDIF
C
C---- DZDG contributions
        DZDGJO = QOPI * (PSIS - PSID)
        DZDGJP = QOPI * (PSIS + PSID)
C
C---- Output
        IF(K.GT.0) WRITE(*,'(A)') '        ,'
        K = K + 1
        WRITE(*,'(A)') '        {'
        WRITE(*,'(A,I4,A)') '          "i": ', I, ','
        WRITE(*,'(A,I4,A)') '          "jo": ', JO, ','
        WRITE(*,'(A,E17.10,A)') '          "PSIS": ', PSIS, ','
        WRITE(*,'(A,E17.10,A)') '          "PSID": ', PSID, ','
        WRITE(*,'(A,E17.10,A)') '          "G1": ', G1, ','
        WRITE(*,'(A,E17.10,A)') '          "G2": ', G2, ','
        WRITE(*,'(A,E17.10,A)') '          "T1": ', T1, ','
        WRITE(*,'(A,E17.10,A)') '          "T2": ', T2, ','
        WRITE(*,'(A,E17.10,A)') '          "DZDG_JO": ', DZDGJO, ','
        WRITE(*,'(A,E17.10)') '          "DZDG_JP": ', DZDGJP
        WRITE(*,'(A)') '        }'
C
   90 CONTINUE
  100 CONTINUE
C
      WRITE(*,'(A)') '      ]'
      WRITE(*,'(A)') '    },'
C
      RETURN
      END


C***********************************************************************
      SUBROUTINE OUTPUT_GGCALC
C---- Output GGCALC influence matrix results
      INCLUDE '../../../Xfoil-instrumented/src/XFOIL.INC'
      INTEGER I, J, NOUT
C
      NOUT = MIN(N, 40)
C
      WRITE(*,'(A)') '    {'
      WRITE(*,'(A)') '      "type": "ggcalc",'
      WRITE(*,'(A,I4,A)') '      "n": ', N, ','
C
C---- AIJ diagonal
      WRITE(*,'(A)') '      "AIJ_diagonal": ['
      DO 10 I=1, NOUT
        IF(I.LT.NOUT) THEN
          WRITE(*,'(A,E17.10,A)') '        ', AIJ(I,I), ','
        ELSE
          WRITE(*,'(A,E17.10)') '        ', AIJ(I,I)
        ENDIF
   10 CONTINUE
      WRITE(*,'(A)') '      ],'
C
C---- First column (sample)
      WRITE(*,'(A)') '      "AIJ_col1": ['
      DO 20 I=1, NOUT
        IF(I.LT.NOUT) THEN
          WRITE(*,'(A,E17.10,A)') '        ', AIJ(I,1), ','
        ELSE
          WRITE(*,'(A,E17.10)') '        ', AIJ(I,1)
        ENDIF
   20 CONTINUE
      WRITE(*,'(A)') '      ],'
C
C---- GAMU(:,1) - alpha=0 solution
      WRITE(*,'(A)') '      "GAMU_0": ['
      DO 30 I=1, NOUT
        IF(I.LT.NOUT) THEN
          WRITE(*,'(A,E17.10,A)') '        ', GAMU(I,1), ','
        ELSE
          WRITE(*,'(A,E17.10)') '        ', GAMU(I,1)
        ENDIF
   30 CONTINUE
      WRITE(*,'(A)') '      ],'
C
C---- GAMU(:,2) - alpha=90 solution
      WRITE(*,'(A)') '      "GAMU_90": ['
      DO 40 I=1, NOUT
        IF(I.LT.NOUT) THEN
          WRITE(*,'(A,E17.10,A)') '        ', GAMU(I,2), ','
        ELSE
          WRITE(*,'(A,E17.10)') '        ', GAMU(I,2)
        ENDIF
   40 CONTINUE
      WRITE(*,'(A)') '      ],'
C
C---- Kutta check
      WRITE(*,'(A,E17.10,A)') '      "kutta_0": ', 
     &                        GAMU(1,1)+GAMU(N,1), ','
      WRITE(*,'(A,E17.10)') '      "kutta_90": ',
     &                        GAMU(1,2)+GAMU(N,2)
C
      WRITE(*,'(A)') '    },'
C
      RETURN
      END


C***********************************************************************
      SUBROUTINE OUTPUT_CLCALC(AOADEG)
C---- Solve at given alpha and output results
      INCLUDE '../../../Xfoil-instrumented/src/XFOIL.INC'
      REAL AOADEG
C
      REAL PI
      INTEGER I, NOUT
C
      PI = 4.0*ATAN(1.0)
C
C---- Set angle of attack
      ALFA = AOADEG * PI/180.0
      COSA = COS(ALFA)
      SINA = SIN(ALFA)
C
C---- Combine base solutions (SPECAL/QISET)
      DO 10 I=1, N
        GAM(I) = COSA*GAMU(I,1) + SINA*GAMU(I,2)
        QINV(I) = GAM(I)
   10 CONTINUE
C
C---- Compute Cp
      DO 20 I=1, N
        CPI(I) = 1.0 - (QINV(I)/QINF)**2
   20 CONTINUE
C
C---- Compute CL, CM (CLCALC)
      CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF,
     &            XCMREF,YCMREF,CL,CM,CDP,CL_ALF,CL_MSQ)
C
C---- Find stagnation point
      CALL STFIND
C
      NOUT = MIN(N, 40)
C
      WRITE(*,'(A)') '    {'
      WRITE(*,'(A)') '      "type": "clcalc",'
      WRITE(*,'(A,F8.3,A)') '      "alpha_deg": ', AOADEG, ','
      WRITE(*,'(A,E17.10,A)') '      "alpha_rad": ', ALFA, ','
      WRITE(*,'(A,I4,A)') '      "n": ', N, ','
C
C---- CL, CM
      WRITE(*,'(A,E17.10,A)') '      "CL": ', CL, ','
      WRITE(*,'(A,E17.10,A)') '      "CM": ', CM, ','
      WRITE(*,'(A,E17.10,A)') '      "CDP": ', CDP, ','
      WRITE(*,'(A,E17.10,A)') '      "CL_ALF": ', CL_ALF, ','
C
C---- Stagnation point (IST and SST are in common block)
C     Note: STFIND sets IST in XFOIL.INC via common block
      WRITE(*,'(A,E17.10,A)') '      "SST": ', SST, ','
C
C---- GAM array (sample)
      WRITE(*,'(A)') '      "GAM": ['
      DO 30 I=1, NOUT
        IF(I.LT.NOUT) THEN
          WRITE(*,'(A,E17.10,A)') '        ', GAM(I), ','
        ELSE
          WRITE(*,'(A,E17.10)') '        ', GAM(I)
        ENDIF
   30 CONTINUE
      WRITE(*,'(A)') '      ],'
C
C---- Cp array (sample)
      WRITE(*,'(A)') '      "CP": ['
      DO 40 I=1, NOUT
        IF(I.LT.NOUT) THEN
          WRITE(*,'(A,E17.10,A)') '        ', CPI(I), ','
        ELSE
          WRITE(*,'(A,E17.10)') '        ', CPI(I)
        ENDIF
   40 CONTINUE
      WRITE(*,'(A)') '      ],'
C
C---- QINV array (same as GAM for inviscid)
      WRITE(*,'(A)') '      "QINV": ['
      DO 50 I=1, NOUT
        IF(I.LT.NOUT) THEN
          WRITE(*,'(A,E17.10,A)') '        ', QINV(I), ','
        ELSE
          WRITE(*,'(A,E17.10)') '        ', QINV(I)
        ENDIF
   50 CONTINUE
      WRITE(*,'(A)') '      ]'
C
      WRITE(*,'(A)') '    },'
C
      RETURN
      END
