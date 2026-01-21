FC = /opt/homebrew/bin/gfortran
XFOIL_SRC = /Users/harry/flexfoil-boundary-layer/Xfoil/src
FFLAGS = -I$(XFOIL_SRC) -fallow-argument-mismatch

test_closures: test_closures.f $(XFOIL_SRC)/xblsys.f
	$(FC) $(FFLAGS) -o $@ $^

run: test_closures
	./test_closures > ../../testdata/closures_reference.json

clean:
	rm -f test_closures *.o

.PHONY: run clean
