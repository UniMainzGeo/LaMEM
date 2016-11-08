LaMEM TESTING
==========================================================================================

This directory contains tests to verify the correct functionality of LaMEM. We use the
PythonTessHarness (pyTH) framework of Dave May & Patrick Sanan for this (see
https://bitbucket.org/dmay/pythontestharness), which contains a number of python routines that
largely simplify running such tests.

The general idea is that we run existing examples and compare certain keywords within the output
file to that in an existing log file.

1) Running tests

You can run all tests with:

python runLaMEM_Tests.py

doing "make test" does the same.

If you only want to run a single test (e.g., unit_FB2_a), do that with:

python runLaMEM_Tests.py -t unit_FB2_a

Note that the default option assumes that you do not use a queueing system and that you have mpiexec
available. If that is not the case, delete the file pthBatchQueuingSystem.conf and follow the two
configuration steps.


2) Adding new tests to LaMEM

Adding new tests is reasonably straightforward, and consists of several steps:

a) Create a new test directory

	For every type of test (i.e.: every test that has its own LaMEM input file), you should create a
	new input directory. In general, we name tests as:

	t1_FB1_Direct

	where t? is the number of the test (please number them consecutively), followed by something
	that explains the meaning of the test (FallingBlock with direct solvers in this case).

b) Put the relevant files in the new test directory In most case you will need a new LaMEM input
	file, so copy that here.

c) Update the python script that runs the test An example of test file is
	/t1_FB1_Direct/test_1_FB1.py Copy this file in your new directory and change the name to reflect the
	test directory that you are adding. Within this test file, you can create several subtests (e.g.,
	running the code on 1 core, or on multiple cores, using either debug or optimized LaMEM versions).
	Each of these subtests should be called test_a, test_b etc.

d) Create the "expected" file for each of the sub-tests. Doing this requires you to run LaMEM once,
	to create the correct output. For example, test_c in test_1_FB1.py has the following options:

	ranks = 2 launch = '../bin/opt/LaMEM -ParamFile ./t1_FB1_Direct/FallingBlock_mono_PenaltyDirect.dat -jp_pc_factor_mat_solver_package mumps'
	expected_file = 't1_FB1_Direct/FB1_Direct_2cores.expected'

  	The file "FB1_Direct_2cores.expected" in the directory ${LaMEM}/tests/t1_FB1_Direct contains
  	the correct output. We generated this file by running, from the /tests/ directory:

  	mpiexec -n 2 ../bin/opt/LaMEM -ParamFile ./t1_FB1_Direct/FallingBlock_mono_PenaltyDirect.dat -jp_pc_factor_mat_solver_package mumps > t1_FB1_Direct/FB1_Direct_2cores.expected

e) Add the tests to the runLaMEM.py Again this is done in two steps:

	 * Add the directory with your new test by adding a new line here:
	 	sys.path.append(os.path.join(os.environ['PWD'], 't1_FB1_Direct'))
	 	sys.path.append(os.path.join(os.environ['PWD'], 't2_FB2_MG'))

	 * Import the test: import test_1_FB1 as FB1

	 	here, the file "test_1_FB1.py" is the name of the new test you just added. Make sure that
	 	this test file only occurs once in all directories! FB1 is how this test file will be named
	 	later in the script (see next item), so choose a unique name.

	 * 	Register all sub-tests that are contained in the file test_1_FB1.py 
	 	registeredTests = [FB1.test_a(), FB1.test_b(), FB1.test_c(), FB1.test_d()]

	 	Simply list all sub-tests here that are in your file, so they can be executed later.


f) Tests to make sure that it works Run the test suite again, to ensure that your new tests are
	incorporated correctly. Note that you don't have to run the full test-suite for this (see above).

g) Commit to LaMEM Push your new tests to the LaMEM repository and check that it works on other
	machines as well.

As a general remark, it is very important that all aspects of your work are being tested as this is
the only way that we can guarantee that things are still running, once we upgrade the code or make
other changes. As each tests is a full python script, you can also create more complicated test
cases, which create partitioning files, run MATLAB if that is available and create






