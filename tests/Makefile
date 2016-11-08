# Two functions:
#  1) clean test directory (make clean)
#  2) call python script that runs the tests using the pyTestHarness framework by Dave May & Patrick Sanan

test: 
	@python runLaMEM_Tests.py

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
