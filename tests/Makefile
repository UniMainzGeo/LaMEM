
# Define compilation mode
# deb - debug 
# opt - optimized (default)

mode = opt

.PHONY : check
check :
ifneq ($(mode), deb)
ifneq ($(mode), opt)
	@echo "............................................."
	@echo "Error! mode="$(mode)" is invalid"
	@echo "Should be either mode=deb or mode=opt"
	@echo "............................................."
	@exit 1
endif
endif
ifeq ($(mode), deb)
	@echo "............................................."
	@echo "...... Running tests in debug mode .........."
	@echo "............................................."
	@ln -sf ../bin/deb/LaMEM ../bin/LaMEM
	@rm -rf Timestep*
	@export PETSC_DIR=${PETSC_DEB}; ./scripts/ExecuteTests.pl
endif
ifeq ($(mode), opt)
	@echo "............................................."
	@echo "..... Running tests in optimized mode ......."
	@echo "............................................."
	@ln -sf ../bin/opt/LaMEM ../bin/LaMEM
	@rm -rf Timestep*
	@export PETSC_DIR=${PETSC_OPT}; ./scripts/ExecuteTests.pl
endif

.PHONY : clean
clean :
	@echo "............................................."
	@echo ".......... Performing full clean ............"
	@echo "............................................."
	@rm -f output/*
	@rm -f log/*
	@rm -f *.info
	@rm -f *.breakpoint
	@rm -f *.out
	@rm -f *.pvd
	@rm -f *.vtk
	@rm -f *.vts
	@rm -f *.dat
	@rm -f *.bin
	@rm -f *.xml
	@rm -rf ./Timestep*
	@rm -rf ./Initial*
	@rm -rf ./Breakpoint
	@rm -f false
	@rm -f ../bin/LaMEM

../output :
	@if [ ! -d $@ ]; then mkdir $@; fi
	
../log :
	@if [ ! -d $@ ]; then mkdir $@; fi