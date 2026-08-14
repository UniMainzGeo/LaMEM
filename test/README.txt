==========================================================================================
LaMEM TESTING
==========================================================================================

This directory contains tests to verify the correct functionality of LaMEM. 
We use the build-in testing framework of Julia (provided through the Test package) for this

The general idea is that we run existing examples and compare certain keywords within the output
file to that in an existing log file.

1) Running tests

You can run all tests with:

julia --project=../. start_tests.jl 

This compiles LaMEM (deb & opt) and runs the test-suite

Doing "make test" does the same.

You can also do it from within julia with:

julia --project=../.

Within julia change to the package manager & start the tests:
julia> ]
(LaMEM_C) pkg> test

2) Adding new tests to LaMEM

Adding new tests is reasonably straightforward, and consists of several steps:

a) Create a new test directory.

	For every type of test (i.e.: every test that has its own LaMEM input file), you should create a
	new input directory. In general, we name tests as:
	
	t01_FB1_Direct
	
	where txx is the number of the test (please number them consecutively), followed by something
	that explains the meaning of the test (FallingBlock with direct solvers in this case).

	Use t01, t02, ... prefixes instead of t1, t2, ...
	Use t01, t02, ... only for directories, not for files (literally all files must be free of this prefix)

b) Put the relevant LaMEM input file (*.dat) in the new test directory. 
	If you need to create a more complicated input geometry, you might also have to create a new Julia input file. 
	Have a look at 
	
	./t03_SubductionGMGinput/CreateMarkers_SubductionVEP_parallel.jl 
	
	for an example.
	
	Make sure the test does not run long. Limit the resolution by 64x64x64 in 3D (or equivalent).
	Do not use more than 30 time steps, unless it is a 1D test, or heat diffusion problem, or a very small domain.
	
	Never use continue_on_fail = 1 option!
	
	All tests should always converge to a given tolerance.
	
c) Add the test to "runtests.jl". 

	Example:

	if should_run_test("t09_PhaseDiagrams")
	@testset "t09_PhaseDiagrams" begin
	    cd(test_dir)
	    dir = "t09_PhaseDiagrams";
	    
	    ParamFile = "FallingBlock_PhaseDiagrams.dat";

	    keywords = ("|Div|_inf", "|mRes|_2")
	    acc      = ((rtol=1e-5, atol=1e-7), (rtol=1e-3, atol=1e-4))

	    # Perform tests
	    @test perform_lamem_test(dir,ParamFile, "FallingBlock_PhaseDiagrams",
	                             args="-mfmax 0.15",
	                             keywords=keywords, accuracy=acc, cores=2, mpiexec=mpiexec,
	                             create_expected_file=update_expected, clean_dir=clean_files)
	end
	end

	Provide expected file name without ".expected" extension. This is is added automatically.
	Do not use postfix -p1, -p2, ... indicating number of MPI processes (this is always broken).
	
	In some cases you may have to first generate a setup (see "t03_Subduction") 
	or you want to compare the results with those of an analytical solution and/or create plots ("t13_Rheology0D" or "t14_1DStrengthEnvelope")
	
	Sometimes you want to permanently set clean_dir=false, e.g. when you want to compare the results with closed-form solution ("t13_Rheology0D").
	But later you anyway should clean the entire test directory like this:
	
	if clean_files
		clean_test_directory(dir)
	end
	
	or like this:
	
	if clean_files
		clean_directory(dir)
	end
	
d) Create the "expected" file for your new test.

	This can be done from the command-line (provided your test directory is t09 as in the example)
	
	make update 09
	
	If you want to review the failure reason for the existing test you can use different Makefile target:

	make work 09
   
	This will keep output files for inspection. Once you figure our what is wrong and fix the issue, you can again regenerate the expected files as described above.

e) Make sure that all works fine by running the full test-suite again

f) Please delete all unused files before committing

   Take care of the file extensions!

   Make sure that tests only write files with the following extensions (outside Timestep, marker, or restart directories):
   .xml .png .pvd .vts .out .log .png .bin
   
   Never use these extensions for input files (e.g. use .dat for main text input and .poly for geomIO binary files)
   
   Use .out for all textual diagnostic output like main console output, permeability info, or adjoint diagnostics.

   The following extensions will be automatically deleted by test framework .pvd .vts .out .log .png .bin
   
   Successful Valgrind runs will delete .xml files. Use make clean to delete .xml files when errors are detected.
   
   .png files will only be removed by make clean.

g) Commit to LaMEM 
	Push your new tests to the LaMEM repository (including the changes to runtests and the required input/expected files).
	and check that it works on other machines as well.

h) If you do not have writing rights to LaMEM: 
	- fork the code
	- create a new branch for your changes
	- push the changes there
	- create a pull request to the main branch

As a general remark, it is very important that all aspects of your work are being tested as this is
the only way that we can guarantee that things are still running, once we upgrade the code or make changes. 

Questions/remarks?
kaus@uni-mainz.de
popov@uni-mainz.de
