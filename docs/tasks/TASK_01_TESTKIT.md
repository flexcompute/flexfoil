# Task 01: Create rustfoil-testkit Crate

## Objective
Create test infrastructure for comparing Rust implementations against FORTRAN.

## Context
- Workspace: `/Users/harry/flexfoil-boundary-layer/crates/` (Rust crates)
- FORTRAN source: `/Users/harry/flexfoil-boundary-layer/Xfoil/src/`
- gfortran: `/opt/homebrew/bin/gfortran`

## Deliverables

### 1. Create crate structure

```
/Users/harry/flexfoil-boundary-layer/crates/rustfoil-testkit/
├── Cargo.toml
├── src/
│   ├── lib.rs
│   ├── fortran_runner.rs   # Compile & run FORTRAN test programs
│   └── approx.rs           # Float comparison utilities
└── fortran/
    ├── Makefile
    └── test_closures.f     # Test harness for BL closure functions
```

### 2. Cargo.toml
```toml
[package]
name = "rustfoil-testkit"
version = "0.1.0"
edition = "2021"

[dependencies]
serde = { version = "1", features = ["derive"] }
serde_json = "1"
tempfile = "3"

[dev-dependencies]
approx = "0.5"
```

### 3. src/lib.rs
- Export `fortran_runner` and `approx` modules
- Provide `load_reference<T>(path)` function to load JSON test data

### 4. src/fortran_runner.rs
- `compile_fortran(source_files, output)` - invoke gfortran
- `run_fortran_test(executable)` - run and capture stdout as JSON
- Handle XFOIL include paths: `-I/Users/harry/flexfoil-boundary-layer/Xfoil/src`

### 5. src/approx.rs
- `assert_close(a, b, tol, name)` - panics with helpful message
- `relative_error(a, b)` - compute |a-b|/max(|a|,|b|,1e-10)

### 6. fortran/test_closures.f
Test harness that outputs JSON for HKIN function:
```fortran
      PROGRAM TEST_CLOSURES
      IMPLICIT NONE
      INCLUDE 'BLPAR.INC'
      
C     Initialize constants (from XFOIL's INIT)
      SCCON  = 5.6
      GACON  = 6.70
      GBCON  = 0.75
      GCCON  = 18.0
      DLCON  = 0.9
      CTCON  = 0.5/(GACON**2 * GBCON)
      CFFAC  = 1.0
      
      CALL TEST_HKIN()
      END
      
      SUBROUTINE TEST_HKIN()
      REAL H, MSQ, HK, HK_H, HK_MSQ
      INTEGER I, J
      
      WRITE(*,'(A)') '{"hkin_tests": ['
      DO I = 1, 20
        H = 1.0 + (I-1) * 0.2  ! H from 1.0 to 4.8
        DO J = 1, 20
          MSQ = (J-1) * 0.05   ! M² from 0 to 0.95
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
      WRITE(*,'(A)') ']}'
      END
```

### 7. fortran/Makefile
```makefile
FC = gfortran
XFOIL_SRC = /Users/harry/flexfoil-boundary-layer/Xfoil/src
FFLAGS = -I$(XFOIL_SRC) -fallow-argument-mismatch

test_closures: test_closures.f $(XFOIL_SRC)/xblsys.f
	$(FC) $(FFLAGS) -o $@ $^

run: test_closures
	./test_closures > ../testdata/closures_reference.json

clean:
	rm -f test_closures *.o
```

### 8. Update workspace Cargo.toml
Add to `/Users/harry/flexfoil-boundary-layer/Cargo.toml`:
```toml
members = [
    # ... existing members ...
    "crates/rustfoil-testkit",
]
```

## Verification
```bash
cd /Users/harry/flexfoil-boundary-layer/crates/rustfoil-testkit/fortran
make run
# Should produce JSON output
cat ../testdata/closures_reference.json | head -20
```

## Next Task
After completion, proceed to TASK_02_BL_CONSTANTS.md

---

## Documentation Requirements

Also ensure that you update Docusaurus with progress.

Explain what tests were for, what they show, and how they passed/failed/worked and consequences.
