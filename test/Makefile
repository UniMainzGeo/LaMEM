# Targets:
#  1) Call julia test harness
#  2) Delete basic output (make clean)
#  3) Delete all output (make purge)

# define PETSc installation directories:
# PETSC_DEB = /directory/where/petsc/debug/is/installed
# PETSC_OPT = /directory/where/petsc/optimized/is/installed

# check integer types of both PETSc versions
IS_64BIT_DEB := $(shell grep -c "PETSC_USE_64BIT_INDICES 1" ${PETSC_DEB}/include/petscconf.h)
IS_64BIT_OPT := $(shell grep -c "PETSC_USE_64BIT_INDICES 1" ${PETSC_OPT}/include/petscconf.h)

ifneq ($(IS_64BIT_DEB), $(IS_64BIT_OPT))
   $(error Debug and optimized PETSc version are using different types of integers [32-bit vs. 64-bit])
endif

ifeq ($(IS_64BIT_DEB), 1)
   ARGS = "is64bit"
   BIT  = "64-bit"
else
   ARGS = ""
   BIT  = "32-bit"
endif

.PHONY : test grind

test:
	@echo "===================================================="
	@echo "Running LaMEM Julia testing framework in" ${BIT} "mode"
	@echo "===================================================="
	@julia --project=../. start_tests.jl ${ARGS}

grind:
	@echo "==================================================================="
	@echo "Running LaMEM Julia testing framework under Valgrind in" ${BIT} "mode"
	@echo "(this is a lot slower than 'make test' - be patient)"
	@echo "==================================================================="
	@julia --project=../. start_tests.jl ${ARGS} valgrind
	@echo "========================================================="
	@echo "Checking Valgrind output for leaks & uninitialised values"
	@echo "========================================================="
	@julia --project=../. check_valgrind_report.jl; ret=$$?; \
 	rm -f **/*.xml \
	exit $$ret
