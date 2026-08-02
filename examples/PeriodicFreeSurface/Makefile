
.PHONY : test clean purge

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
	
