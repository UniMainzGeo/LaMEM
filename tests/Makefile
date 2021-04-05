# Targets:
#  1) call python script that runs the tests using the pyTestHarness framework by Dave May & Patrick Sanan
#  2) delete basic output (make clean)
#  3) delete all output (make purge)

# Set environmental variables
export PYTHONUNBUFFERED

# Check whether pythonharness has been downloaded properly

.PHONY : test clean purge

test:
	@echo "Testing whether pyTestHarness exists"
	@test -s ./pythontestharness/lib || { echo "Looks like pyTestHarness does not exist! Installing from git..."; git submodule init; git submodule update; git submodule status;  }

	@echo "Running LaMEM tests using the pythontestharness"
	@echo "You can do this, and have additional options, with:" 
	@echo "$ python runLaMEM_Tests.py"
	@echo "In case problems occur, delete 'pthBatchQueuingSystem.conf', and run the python script directly"
	@python runLaMEM_Tests.py

clean :
	@echo "............................................."
	@echo ".......... Performing full clean ............"
	@echo "............................................."
	@rm -f *.pvd
	@rm -f *.bin
	@rm -f *.xml
	@rm -f *.output
	@rm -f permea.darcy.dat
	@rm -rf ./Timestep*
	@rm -rf ./Strainrate*
	@rm -rf ./restart
	@rm -rf ./markers*
	@rm -rf ./output
	@rm -rf ./RTwav_*
	@find . -name *.pyc -delete

purge : clean
	@rm -rf *.log
	@rm -rf *.output
	@rm -f output/*
	@rm -f log/*
