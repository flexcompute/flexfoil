      PROGRAM TEST_CLOSURES
C     Declare BLPAR variables (COMMON block requires type declarations)
      REAL SCCON, GACON, GBCON, GCCON, DLCON
      REAL CTRCON, CTRCEX, DUXCON, CTCON, CFFAC
      INCLUDE 'BLPAR.INC'
      
C     Initialize constants (from XFOIL's INIT subroutine)
      SCCON  = 5.6
      GACON  = 6.70
      GBCON  = 0.75
      GCCON  = 18.0
      DLCON  = 0.9
      CTCON  = 0.5/(GACON**2 * GBCON)
      CFFAC  = 1.0
      
C     Output JSON object with all tests
      WRITE(*,'(A)') '{'
      CALL TEST_HKIN()
      WRITE(*,'(A)') ','
      CALL TEST_CFL()
      WRITE(*,'(A)') ','
      CALL TEST_CFT()
      WRITE(*,'(A)') '}'
      END
      
      SUBROUTINE TEST_HKIN()
      REAL H, MSQ, HK, HK_H, HK_MSQ
      INTEGER I, J
      
      WRITE(*,'(A)') '"hkin_tests": ['
      DO I = 1, 20
        H = 1.0 + (I-1) * 0.2
        DO J = 1, 20
          MSQ = (J-1) * 0.05
          CALL HKIN(H, MSQ, HK, HK_H, HK_MSQ)
          IF (I.EQ.20 .AND. J.EQ.20) THEN
            WRITE(*,'(A,F8.4,A,F8.4,A,F12.8,A,F12.8,A,F12.8,A)')
     &        '{"h":', H, ',"msq":', MSQ, ',"hk":', HK,
     &        ',"hk_h":', HK_H, ',"hk_msq":', HK_MSQ, '}'
          ELSE
            WRITE(*,'(A,F8.4,A,F8.4,A,F12.8,A,F12.8,A,F12.8,A)')
     &        '{"h":', H, ',"msq":', MSQ, ',"hk":', HK,
     &        ',"hk_h":', HK_H, ',"hk_msq":', HK_MSQ, '},'
          ENDIF
        ENDDO
      ENDDO
      WRITE(*,'(A)') ']'
      END
      
      SUBROUTINE TEST_CFL()
C     Test laminar skin friction coefficient
      REAL HK, RT, MSQ, CF, CF_HK, CF_RT, CF_MSQ
      INTEGER I, J
      
      WRITE(*,'(A)') '"cfl_tests": ['
C     Test grid: HK from 2.0 to 8.0, RT from 100 to 2000
      DO I = 1, 13
        HK = 2.0 + (I-1) * 0.5
        DO J = 1, 10
          RT = 100.0 + (J-1) * 200.0
          MSQ = 0.0
          CALL CFL(HK, RT, MSQ, CF, CF_HK, CF_RT, CF_MSQ)
          IF (I.EQ.13 .AND. J.EQ.10) THEN
            WRITE(*,'(A,F8.4,A,F10.2,A,F14.10,A,F14.10,A,F14.10,A)')
     &        '{"hk":', HK, ',"rt":', RT, ',"cf":', CF,
     &        ',"cf_hk":', CF_HK, ',"cf_rt":', CF_RT, '}'
          ELSE
            WRITE(*,'(A,F8.4,A,F10.2,A,F14.10,A,F14.10,A,F14.10,A)')
     &        '{"hk":', HK, ',"rt":', RT, ',"cf":', CF,
     &        ',"cf_hk":', CF_HK, ',"cf_rt":', CF_RT, '},'
          ENDIF
        ENDDO
      ENDDO
      WRITE(*,'(A)') ']'
      END
      
      SUBROUTINE TEST_CFT()
C     Test turbulent skin friction coefficient
      REAL SCCON, GACON, GBCON, GCCON, DLCON
      REAL CTRCON, CTRCEX, DUXCON, CTCON, CFFAC
      INCLUDE 'BLPAR.INC'
      REAL HK, RT, MSQ, CF, CF_HK, CF_RT, CF_MSQ
      INTEGER I, J, K
      INTEGER TOTAL, COUNT
      
      WRITE(*,'(A)') '"cft_tests": ['
C     Test grid: HK from 1.2 to 3.5, RT from 500 to 20000, MSQ from 0 to 0.5
      TOTAL = 8 * 10 * 6
      COUNT = 0
      DO I = 1, 8
        HK = 1.2 + (I-1) * 0.3
        DO J = 1, 10
          RT = 500.0 + (J-1) * 2000.0
          DO K = 1, 6
            MSQ = (K-1) * 0.1
            CALL CFT(HK, RT, MSQ, CF, CF_HK, CF_RT, CF_MSQ)
            COUNT = COUNT + 1
            IF (COUNT.EQ.TOTAL) THEN
              WRITE(*,'(A,F8.4,A,F10.1,A,F6.3,A,F14.10,A,
     &                  F14.10,A,F14.10,A,F14.10,A)')
     &          '{"hk":', HK, ',"rt":', RT, ',"msq":', MSQ,
     &          ',"cf":', CF, ',"cf_hk":', CF_HK,
     &          ',"cf_rt":', CF_RT, ',"cf_msq":', CF_MSQ, '}'
            ELSE
              WRITE(*,'(A,F8.4,A,F10.1,A,F6.3,A,F14.10,A,
     &                  F14.10,A,F14.10,A,F14.10,A)')
     &          '{"hk":', HK, ',"rt":', RT, ',"msq":', MSQ,
     &          ',"cf":', CF, ',"cf_hk":', CF_HK,
     &          ',"cf_rt":', CF_RT, ',"cf_msq":', CF_MSQ, '},'
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      WRITE(*,'(A)') ']'
      END
