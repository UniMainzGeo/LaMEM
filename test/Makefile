# Targets:
#  1) Call Julia test harness (make test)
#  2) Run tests under Valgrind (make grind)
#  3) Remove remaining files (make clean)

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

# Allow running a subset of tests by passing extra "words" after the target,
# e.g.:
#   make test 01 05 32
#   make test 03-07 11 12-17
# Make would normally try to build "01", "05", etc. as targets of their own,
# so we grab everything besides "test" "grind" "clean" as TEST_SELECTORS and turn
# those words into harmless no-op targets.
TEST_SELECTORS := $(filter-out test grind clean,$(MAKECMDGOALS))
$(eval $(TEST_SELECTORS):;@:)

# Test mode, settable on the command line, e.g.:
#   make test mode=work
#   make test mode=update
# test   (default) - run tests, clean up generated output/work files
# work   - run tests, keep generated files around for inspection
# update - run tests and OVERWRITE the expected (reference) files
mode ?= test
ifeq ($(filter $(mode),test work update),)
   $(error Invalid mode "$(mode)": must be one of "test", "work", "update")
endif

.PHONY : test grind

test:
	@echo "===================================================="
	@echo "Running LaMEM Julia testing framework in" ${BIT} "mode"
	@echo "===================================================="
	@julia --project=../. start_tests.jl ${ARGS} mode=${mode} ${TEST_SELECTORS}

grind:
	@echo "==================================================================="
	@echo "Running LaMEM Julia testing framework under Valgrind in" ${BIT} "mode"
	@echo "(this is a lot slower than 'make test' - be patient)"
	@echo "==================================================================="
	@julia --project=../. start_tests.jl ${ARGS} valgrind ${TEST_SELECTORS}
	@echo "========================================================="
	@echo "Checking Valgrind output for leaks & uninitialised values"
	@echo "========================================================="
	@julia --project=../. check_valgrind_report.jl; ret=$$?; \
 	rm -f **/*.xml \
	exit $$ret
	
clean:
	@rm -f **/*.xml

