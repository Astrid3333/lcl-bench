FC = gfortran
FFLAGS_REL = -O2 -std=f2008 -Wall -fcheck=bounds -ffp-contract=off
FFLAGS_DBG = -g -O0 -std=f2008 -fcheck=all -fbacktrace
FFLAGS_ASAN = -g -O1 -std=f2008 -fcheck=all -fsanitize=address,undefined -fbacktrace
LIBS = -llapack -lblas
SRC = kinds.f90 fjson.f90 fhandles.f90 fkernels.f90 frunner.f90

all: f_runner

f_runner: $(SRC)
	$(FC) $(FFLAGS_REL) $(SRC) $(LIBS) -o f_runner

f_runner_dbg: $(SRC)
	$(FC) $(FFLAGS_DBG) $(SRC) $(LIBS) -o f_runner_dbg

f_runner_asan: $(SRC)
	$(FC) $(FFLAGS_ASAN) $(SRC) $(LIBS) -o f_runner_asan

clean:
	rm -f *.o *.mod f_runner f_runner_dbg f_runner_asan

.PHONY: all clean
