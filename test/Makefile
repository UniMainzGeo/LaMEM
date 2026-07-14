# Targets:
#  1) Call Julia test harness, clean up generated files afterwards (make test)
#  2) Call Julia test harness, keep generated files for inspection (make work)
#  3) Call Julia test harness and OVERWRITE the expected/reference files (make update)
#  4) Run tests under Valgrind (make grind)
#  5) Re-check existing Valgrind .xml output from a previous 'make grind' (make report)
#  6) Remove remaining files (make clean)

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
# so we grab everything besides "test" "work" "update" "grind" "report" "clean"
# as TEST_SELECTORS and turn those words into harmless no-op targets.
TEST_SELECTORS := $(filter-out test work update grind report clean,$(MAKECMDGOALS))
$(eval $(TEST_SELECTORS):;@:)

.PHONY : test work update grind report clean

# test   - run tests, clean up generated output/work files
# work   - run tests, keep generated files around for inspection
# update - run tests and OVERWRITE the expected (reference) files
test:
	@echo "===================================================="
	@echo "Running LaMEM Julia testing framework in" ${BIT} "mode"
	@echo "===================================================="
	@julia --project=../. start_tests.jl ${ARGS} mode=test ${TEST_SELECTORS}

work:
	@echo "===================================================="
	@echo "Running LaMEM Julia testing framework in" ${BIT} "mode"
	@echo "===================================================="
	@julia --project=../. start_tests.jl ${ARGS} mode=work ${TEST_SELECTORS}

update:
	@echo "===================================================="
	@echo "Running LaMEM Julia testing framework in" ${BIT} "mode"
	@echo "===================================================="
	@julia --project=../. start_tests.jl ${ARGS} mode=update ${TEST_SELECTORS}

grind:
	@echo "==================================================================="
	@echo "Running LaMEM Julia testing framework under Valgrind in" ${BIT} "mode"
	@echo "(this is a lot slower than 'make test' - be patient)"
	@echo "==================================================================="
	@julia --project=../. start_tests.jl ${ARGS} valgrind ${TEST_SELECTORS}; test_ret=$$?; \
	echo "========================================================="; \
	echo "Checking Valgrind output for leaks & uninitialised values"; \
	echo "========================================================="; \
	julia --project=../. check_valgrind_report.jl; report_ret=$$?; \
	if [ $$test_ret -eq 0 ] && [ $$report_ret -eq 0 ]; then \
		rm -f **/*.xml; \
	else \
		echo "==================================================================="; \
		echo "Test suite exit code: $$test_ret, valgrind report exit code: $$report_ret"; \
		echo "Keeping .xml files in place for inspection"; \
		echo "==================================================================="; \
	fi; \
	exit $$(( $$test_ret != 0 || $$report_ret != 0 ))

# report - re-check existing .xml files from a previous 'make grind' run
#          without re-running the (expensive) valgrind test suite.
#          Leaves the .xml files in place either way.
report:
	@echo "========================================================="
	@echo "Checking Valgrind output for leaks & uninitialised values"
	@echo "========================================================="
	@julia --project=../. check_valgrind_report.jl
	
clean:
	@echo "............................................."
	@echo ".......... Performing full clean ............"
	@echo "............................................."
	@rm -f  **/*.xml
	@rm -f  **/*.pvd
	@rm -f  **/*.out
	@rm -f  **/*.log
	@rm -f  **/*.png
	@rm -f  **/*.vts
	@rm -f  **/*.bin
	@rm -rf **/*Timestep*
	@rm -rf **/*markers*
	@rm -rf **/*restart

