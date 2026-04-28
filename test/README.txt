==========================================================================================
LaMEM TESTING
==========================================================================================

This directory contains tests to verify the correct functionality of LaMEM. 
We use the build-in testing framework of julia (provided through the Test package) for this

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

	WARNING! USE t01, t02, ... INSTEAD OF t1, t2, ...
	WARNING! USE t01, t02, ... ONLY FOR DIRECTORIES, NOT FOR FILES (LITERALLY ALL FILES MUST BE FREE OF THIS PREFIX)

b) Put the relevant LaMEM input file (*.dat) in the new test directory. 
	If you need to create a more complicated input geometry, you might also have to create a new julia input file. 
	Have a look at 
	
	./t03_SubductionGMGinput/CreateMarkers_SubductionVEP_parallel.jl 
	
	for an example.
	
c) Add the test to "runtests.jl". 

	The main routine to do the tests is 

		perform_lamem_test
    
	It has quite a few options, which you can read with the help information:
	
	julia>?perform_lamem_test	

	The way most LaMEM tests work is that we store the LaMEM output file in a file and 
	compare certain keywords with an "expected" file generated at an earlier stage.
	In the simplest case, this would be:

 	dir = "t01_FB1_Direct";
    ParamFile = "FallingBlock_mono_PenaltyDirect.dat";
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-8), (rtol=1e-5,atol=1e-8), (rtol=1e-4,atol=1e-8));
    
    # Perform test:
    @test perform_lamem_test(dir,ParamFile,"FB1_a_Direct_opt", 
                            keywords=keywords, accuracy=acc, cores=1, opt=true)

	WARNING! PROVIDE EXPECTED FILE NAME WITHOUT EXTENSION .expected IS ADDED INTERNALLY
	WARNING! DO NOT USE POSTFIX -p1, -p2, ... INDICATING NUMBER OF MPI PROCESSES (THIS IS ALWAYS BROKEN)
	
	In some cases you may have to first generate a setup (see "t03_Subduction") 
	or you want to compare the results with those of an analytical solution and/or create plots ("t13_Rheology0D" or "t14_1DStrengthEnvelope")


d) Create the "expected" file for your new test. 
	This can be done from the command-line. 
	Assume that you are in the test directory, and that new test is in the directory "./t45_new" and is called "t45_new_test":

	julia> using LaMEM_C, Test, GeophysicalModelGenerator, LaMEM.IO, CairoMakie
	julia> include("test_utils.jl")
	julia> keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    julia> acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-5, atol=1e-11), (rtol=2e-4,atol=1e-10));

	julia> perform_lamem_test(dir,"1D_VP.dat","t14_1D_VP_Direct_opt",
                            keywords=keywords, accuracy=acc, cores=1, opt=true,
							create_expected_file=refresh_expected, clean_dir=clean_files)
	
	WARNING! PROVIDE EXPECTED FILE NAME WITHOUT EXTENSION .expected IS ADDED INTERNALLY
		
	For designing the test you can set the flags as follows:
	
	refresh_expected = true
	clean_files      = false

	Make sure to revert the flags after completing the test design:
	
	refresh_expected = false
	clean_files      = true
	
	Sometimes you want to permanently set clean_dir=false, e.g. when you want to compare the results with closed-form solution ("t13_Rheology0D").
	But later you anyway should clean the entire test directory like this:
	
	if clean_files
		clean_test_directory(dir)
	end
	
	or like this:
	
	if clean_files
    	clean_directory(dir)
    end
    
    Here clean_files flag has the same meaning as above. You can use it to keep the files while designing/debugging the test.
    
    In case that test is failing you can set the flag clean_files=false and compare .out file against .expected
    
    e.g.: t14_1D_VP_Direct_opt.out vs. t14_1D_VP_Direct_opt.expected
    
	Once you figure our what is wrong and fix the issue, you can again regenerate the expected file as described above.

e) Tests to make sure that it works by running the full test-suite again

f) PLEASE DELETE ALL UNUSED FILES BEFORE COMMITTING!

g) Commit to LaMEM 
	Push your new tests to the LaMEM repository (including the changes to runtests and the required input/expected files) 
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
