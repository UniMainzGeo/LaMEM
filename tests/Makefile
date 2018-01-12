# Targets:
#  1) call python script that runs the tests using the pyTestHarness framework by Dave May & Patrick Sanan
#  2) delete basic output (make clean)
#  3) delete all output (make purge)

export PYTHONPATH=$(PWD)/pythontestharness/lib
export PYTHONUNBUFFERED

.PHONY : test clean purge

test:
	@python runLaMEM_Tests.py

clean :
	@echo "............................................."
	@echo ".......... Performing full clean ............"
	@echo "............................................."
	@rm -f *.pvd
	@rm -f *.bin
	@rm -f *.xml
	@rm -rf ./Timestep*
	@rm -rf ./restart
	@rm -rf ./markers*
	@find . -name *.pyc -delete

purge : clean
	@rm -rf *.log
	@rm -f output/*
	@rm -f log/*
