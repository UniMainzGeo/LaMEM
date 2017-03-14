# Targets:
#  1) call python script that runs the tests using the pyTestHarness framework by Dave May & Patrick Sanan
#  2) delete basic output (make clean)
#  3) delete all output (make purge)

.PHONY : test clean purge

test:
	@python runLaMEM_Tests.py

clean :
	@echo "............................................."
	@echo ".......... Performing full clean ............"
	@echo "............................................."
	@rm -f output/*
	@rm -f log/*
	@rm -f *.pvd
	@rm -f *.xml
	@rm -rf ./Timestep*
	@find . -name *.pyc -delete

purge : clean
	@rm -rf ./restart
	@rm -rf ./markers
	@rm -rf *.log
