# Targets:
#  1) Call julia test harness
#  2) Delete basic output (make clean)
#  3) Delete all output (make purge)

.PHONY : test clean purge

test:
	@echo "============================================="
	@echo "Running the LaMEM julia testing framework"
	@echo "============================================="

# compile LaMEM & run test suite:			
	@julia --project=../. start_tests.jl 

clean :
	@echo "............................................."
	@echo ".......... Performing full clean ............"
	@echo "............................................."
	@rm -f *.pvd
	@rm -f *.bin
	@rm -f *.xml
	@rm -f *.output
	@rm -f permea.darcy.dat
	@rm -f ScalingLaw.dat
	@rm -f ScalingLaw_Test.dat
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
	@rm -f **/*.png
	@rm -f **/*.log
	