#!/usr/bin/perl

use File::Path;

# Declare the subroutines
sub trim($);
sub ltrim($);
sub rtrim($);
sub get_mpiexec;
sub strip_equals;


sub execute;
sub PerformComparison;
sub check_line_for_pointer;
sub Comparison_literal_diff;

sub Comparison_literal_diff_pl;

sub Comparison_tolerance_diff_pl;
sub compare_lines_using_tolerance_on_numbers;
sub check_for_scientific_notation;
sub check_for_scientific_notation_E;
sub is_num;

sub Comparison_catch_exception_pl;
sub check_for_error_message;

# ~~~~~~~~~~~~~~~~~~~~~~~ MAIN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

$debug = 'false';

if ($ARGV[0] eq 'clean' ) {
    if( -d 'log' ) {
        if( $debug eq 'true' ) { print "Removing all log files \n"; }    
        `rm -f log/*`;
    }
    if( -d 'output' ) {
        if( $debug eq 'true' ) { print "Removing all output files \n"; }
        `rm -f output/*`;
 }
} else {
    execute;
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


##############################################################################
##  Executes a testing procedure within a directory.
##  Each x.test file contains the execute command of the test to conduct
##  and the type of test which should be performed. See PerformComparison()
##  for list of test allowed to be performed.
##############################################################################
sub execute {

printf( "=======================================================================================================\n");

# User must define which mpi to use
# Examine ${PETSC_DIR}/bmake/${PETSC_ARCH}/petsconf and look for MPIEXEC
$path_to_mpiexec = '';
$mpi_exec = '';
($path_to_mpiexec,$mpi_exec) = get_mpiexec();
if( $mpi_exec eq 'mpirun' ) {
	$MPI_LAUNCH = $path_to_mpiexec . " -np";
}
elsif( $mpi_exec eq 'mpiexec' ) {
	$MPI_LAUNCH = $path_to_mpiexec . " -n";
}
else {
	$error_mesg = "Don't know how to specfiy processor flag for MPIEXEC = " . $path_to_mpiexec . " as set in $ENV{'PETSC_DIR'}/bmake/$ENV{'PETSC_ARCH'}/petscconf \n";
	$error_mesg = $error_mesg . "You should add and appropriate test into testing/ExecuteTests.pl \n";
	die( $error_mesg );
}
if( $debug eq 'true' ) { printf("ExecuteTests-.pl is using $MPI_LAUNCH to run tests \n" ); }

# get directory name
use Cwd;
$path = getcwd;

$tests_performed = 0;
$tests_passed = 0;
$tests_failed = 0;

my $found_EXEC;
my $found_TEST_TYPE;


my $single_test = 0;

# scan arguments and look for --test=XXX
foreach $_arg_i (@ARGV) {
  if( $_arg_i =~ m/\-\-test\=/ ) {
    @data = split(/=/, $_arg_i);
    @files = $data[1];
    $single_test = 1;
  }
}
if( $single_test == 0 ) {
  # Scan directory and get all files with the extension .test
  # Test name is the prefix of each .test
  @files = <*.test>;
}
print "Testing the following files:\n";
foreach $file (@files) {
    print "    $file \n";
}

foreach $file (@files) {

    $found_EXEC = 'false';
    $found_TEST_TYPE = 'false';
  
    # remove directories that might spoil output
    #    rmtree('Timestep*');
    `rm -rf Timestep*`;
    
  
    # look for last instance of the fullstop "."
    $file_L = length $file;
    $k = 0;
    # split into array of chars and search for the dot
    my @chars = split //,$file;
    foreach $c (@chars) {
      if( $c eq "." ) {
        $dot_loc = $k;
      }
      $k++;
    }

    # strip off the .test extension or whatever happens to be AFTER the last fullstop
    $test_name = substr( $file, 0,  -($file_L - $dot_loc) ); # Remove 5 characters from file name
    
    # Open file and get the execute command and diff type
    open( LAUNCH_FILE, "$file") || die "Can't open '$file': $!\n";
    @testContents = <LAUNCH_FILE>;

    foreach $string (@testContents) {
        @values = split('=', $string);

        $t2 = trim( @values[0] ); # remove white space on left and right

        # check if exec
        if( $t2 eq "EXEC" ) {
            $found_EXEC = 'true';
            $exec = trim(@values[1]);     # remove white space on left/right
            chomp( $exec );     # strip line break from end of executable so that stdout is properly captured
        #    printf "EXEC = $exec \n";
        }

        # check if nproc
        if( $t2 eq "NPROCS" ) {
            $num_procs = trim(@values[1]);     # remove white space on left/right
            chomp( $num_procs );     # strip line break from end of executable so that stdout is properly captured
        #    printf "NPROCS = $num_procs \n";
        }

        # check if diff type
        if( $t2 eq "TEST_TYPE" ) {
            $found_TEST_TYPE = 'true';
            $test_type = trim(@values[1]);     # remove white space on left/right
            chomp( $test_type );     # strip line break from end of executable so that stdout is properly captured
        #    printf "TEST_TYPE = $test_type \n";
        }
    }
    close( LAUNCH_FILE );

    my $L;
    $L = length( $test_type );
    if( $L == 0 ) {
        die "You need to specify a particular test to compare the ouput!!. Try TEST_TYPE = literal_diff \n";
    }

    if( $found_TEST_TYPE eq 'false' ) {
        die "ERROR: Did not find essential parameter TEST_TYPE in $file";
    }
    if( $found_EXEC eq 'false' ) {
        die "ERROR: Did not find essential parameter EXEC in $file";
    }
    


    # Launch executable and send output to output/$test_name.output
#    printf "  Execute: $MPI_LAUNCH $num_procs $exec > output/$test_name.output \n";
      if( $num_procs == 1 ) {
        system( "$path/$exec > output/$test_name.output" );    
      }
      else {
        system( "$MPI_LAUNCH $num_procs $path/$exec > output/$test_name.output" );    
      }
    
    # Execute the appropriate diff operation and send result to log/$test_name.diff
    # diff will return 0 on success, 1 if failed
    # If no difference occurred, report PASSED else report FAILED
    $diff_status = PerformComparison( $test_type, "output/$test_name.output", "expected/$test_name.expected", "log/$test_name.diff"  );


    # Report status of diff
    if( $diff_status == 0 ) {
        printf "$test_name  :\t PASSED\n";
        $tests_passed = $tests_passed + 1;
    } else {
        printf "$test_name  :\t FAILED\n";
        $tests_failed = $tests_failed + 1;
    }
    
    $tests_performed = $tests_performed + 1;

}


# generate a truncated path name which is printed to the screen
my @split_path = split( /\//, $path );
my $N  =@split_path;
my $pname = "";
for( $i=$N-3; $i<$N; $i++ ) {
	$pname = $pname . "/" . $split_path[$i];
}

#printf "Test summary: [$path]    >>>>>>>>>>>>>>>>>>>    ";
printf "Test summary: [$pname]    >>>>>>>>>>>>>>>>>>>    ";
if( $tests_failed == 0 ) {
	printf "  Passed ( $tests_passed / $tests_performed ) \n";
} else {
	printf "  Failed ( $tests_failed / $tests_performed ) \n";
}

}





## ==========================================================================================================
## Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

## ==========================================================================================================
## Left trim function to remove leading whitespace
sub ltrim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	return $string;
}

## ==========================================================================================================
## Right trim function to remove trailing whitespace
sub rtrim($)
{
	my $string = shift;
	$string =~ s/\s+$//;
	return $string;
}

## =============================================================================
##  Examines the file ${PETSC_DIR}/bmake/${PETSC_ARCH}/petscconf and looks 
##  for MPIEXEC. Returns the full path to binary and the binary name
##
sub get_mpiexec {
	my $file_name = $ENV{'PETSC_DIR'} . "/" . $ENV{'PETSC_ARCH'} . "/conf/petscvariables";	
	my $n = 0;
	my $FULLPATH_MPIEXEC = '';
	my $MPIEXEC = '';
	my $first = '';
	
	open( DAT, $file_name ) || die( "Could not open file datafile ($file_name) ! \n" );
	@raw_data=<DAT>;

	foreach $line (@raw_data) {
		chomp( $line );
		my $tline = trim( $line );
		my @blist = split( /\=/, $tline );

	        # Use the first part to determine if the line corresponds to MPIEXEC = 
	        $first = trim( $blist[0] );

	        if( $first eq 'MPIEXEC' ) {
#			print "1)   $first \n";
	
	            # Everything including the /'s, expect the first part, MPIEXEC = , corresponds to the full path to mpiexec
	            my $t1 = strip_equals( $tline );  # take entire line and remove equals sign
	            $t1 =~ s/MPIEXEC//;   # remove word MPIEXEC from line
		    $FULLPATH_MPIEXEC = trim( $t1 );

        	    # The last part of the list will be the executable. ie. /usr/local/mpich/bin/mpirun
		    my @clist = split( /\//, $t1 );
        	    my $N = @clist;
        	    $t1 = $clist[$N-1];
        	    $MPIEXEC = trim( $t1 );
        	    break;
       		}

		$n = $n + 1;
	 }

#	print "2) mpiexec = $MPIEXEC \n ";
#	print "3) path = $FULLPATH_MPIEXEC \n ";

    if( $FULLPATH_MPIEXEC eq '' ) {
    	my $error_mesg = "Failed to locate the variable MPIEXEC from within " .  "$ENV{'PETSC_DIR'}/bmake/$ENV{'PETSC_ARCH'}/petscconf";
        die( $error_mesg );
    }

#    printf "PATHTO_MPIEXEC = $FULLPATH_MPIEXEC \n";
#    printf "MPIEXEC = $MPIEXEC \n";
	
	close( DAT );

    return ( $FULLPATH_MPIEXEC, $MPIEXEC );

}

## ==========================================================================================================
##
sub strip_equals
{
	my $string = shift;

    my $t1 = trim($string);

	$string =~ s/=//g;

    my $t1 = trim($string);

	return $t1;
}


# Returns 0 if passes, anything else indicates comparison test failed
# Comparisions supported
#   literal_diff        <unix diff>
#   literal_diff_pl     <perl version of diff which will ignore pointers>
#   tolerance_diff_pl [tolerance]   <perl diff which till use a compare numbers and only error if the abs(diff) < tol
#   catch_exception_pl [assert_statement_to_search_for] <will look for a string in the output>   
#
#   TEST_TYPE = test_name,arg1,arg2,arg3
sub PerformComparison 
{
    my $test_type = $_[0];
    my $job_output = $_[1];
    my $expected_output = $_[2];
    my $log_file = $_[3];
 
    my @supported = ( 'literal_diff', 'tolerance_diff_pl','catch_exception_pl' );
 
    chomp( $test_type );
    my @list = split( /,/, $test_type );


#    print "PerformComparison \n" ;
#    print "  Between    $job_output     :      $expected_output \n";
#    print "  Using: @list \n" ;    

    my $status = 1;
    if( $list[0] eq "literal_diff" ) {
        $status = Comparison_literal_diff( $job_output, $expected_output, $log_file );
    }
    elsif( $list[0] eq "catch_exception_pl" ) {
        $status = Comparison_catch_exception_pl( $job_output, $list[1], $log_file );
    }
    elsif ( $list[0] eq "tolerance_diff_pl" ) {
        $status = Comparison_tolerance_diff_pl( $job_output, $expected_output, $log_file, $list[1] );
    }
    else {
        print "You requested TEST_TYPE = $test_type \n";
        print "Test types supported include: \n";
        print "    ";
        foreach $i (@supported) {
            print "$i  ";
        }
        print "\n";
        die "";
    }
    
    return $status;
}

# ==============================================================================
sub Comparison_literal_diff
{
    my $job_output = $_[0];
    my $expected_output = $_[1];
    my $log = $_[2];

    my $diff_status = system( "diff $job_output $expected_output > $log" );
    $diff_status = $diff_status/256;
    ## Must divide by 256 to get actual exit value.
    ## See:
    ##  http://www.cs.cmu.edu/People/rgs/pl-exp-sys.html

    return $diff_status;
}

sub check_line_for_pointer
{
    my $line = $_[0];
    return( $line =~ m/0x/ );
}

# ==============================================================================
sub Comparison_literal_diff_pl
{

}


# ==============================================================================
sub Comparison_tolerance_diff_pl
{
    my $job_output = $_[0];
    my $exp_output = $_[1];
    my $log = $_[2];
    my $tol = $_[3];


    my $line_stat = 'pass';
    my $com_log;
    my $stat = 0;

    # Check tolerance is provided
    if( $tol eq '' ) {
        print "ERROR(tolerance_diff_pl): You must provide a tolerance \n";
        die ;
    }

    open( DAT, $job_output ) || die( "Could not open file job output ($job_output) ! \n" );
    my @data=<DAT>;

    open( EDAT, $exp_output ) || die( "Could not open expected file ($exp_output) ! \n" );
    my @edata=<EDAT>;

    open( LDAT, ">$log" ) || die( "Could not open log file ($log) ! \n" );

	my $N_lines_expected_file = @edata;
    my $lc_expected = 0;
	my $lc_output = 0;
	foreach $my_line ( @data ) {

		# OUTPUT_FILE: check line does not start with a hash
        my $first_char    = substr($my_line, 0, 1);
        if( $first_char eq '#' ) {
			$lc_output++;
	        next;
        }
        #  OUTPUT_FILE: check lines are either both contain pointers or both do NOT contain pointers
        my $a_has_ptr = check_line_for_pointer( $my_line );
        if( $a_has_ptr ) {
			$lc_output++;
			next;
		}
		# OUTPUT FILE: check line is not blank
		if( $my_line =~ /^$/ ) {
			$lc_output++;
			next;
		}
        
        # OUTPUT FILE: check line is not blank
        if( $my_line =~ /^$/ ) {
            $lc_output++;
            next;
        }
        
        # OUTPUT FILE: check line does not contain (sec)
        my $substr = '(sec)';
        if (index($my_line, $substr) != -1) {
            $lc_output++;
            next;
        }

        # OUTPUT FILE: check line does not contain Compiled
        my $substr = 'Compiled';
        if (index($my_line, $substr) != -1) {
            $lc_output++;
            next;
        }



		# find the next line in the expected file to check
		my $start = $lc_expected;
		my $II  = 0;
		for( $II=$start; $II < $N_lines_expected_file; $II++ ) {
	        $exp_line = $edata[ $II ];
			# check for comment line
    	    my $e_first_char = substr($exp_line, 0, 1);
        	if( $e_first_char eq '#' ) {
            	next;
        	}
			# EXP_FILE: check for pointer
	        my $b_has_ptr = check_line_for_pointer( $exp_line );
			if( $b_has_ptr ) {
				next;
			}
			# EXP_FILE FILE: check line is not blank
			if( $exp_line =~ /^$/ ) {
				next;
			}
            # EXP_FILE FILE: check line does not contain (sec)
            my $substr = '(sec)';
            if (index($exp_line, $substr) != -1) {
                $lc_output++;
                next;
            }
            # EXP_FILE FILE: check line does not contain Compiled
            my $substr = 'Compiled';
            if (index($exp_line, $substr) != -1) {
                $lc_output++;
                next;
            }
            
			$lc_expected = $II;
			last;
		}
        
		
        chomp( $my_line );
        chomp( $exp_line );
            ($line_stat,$com_log) = compare_lines_using_tolerance_on_numbers( $my_line, $exp_line, $tol );
			my $output_line_num = $lc_output + 1;
			my $exp_line_num      = $lc_expected + 1;
            if( $line_stat eq 'fail' ) {
                $stat = 1;
            
                print LDAT "comparing line\tOUTPUT   [$output_line_num] $my_line \n";
                print LDAT "\t\tEXPECTED [$exp_line_num] $exp_line \n";
                print LDAT "    status: $com_log \n";
            }
			else {
                #   print LDAT "comparing line\tOUTPUT   [$output_line_num] $my_line \n";
                #   print LDAT "\t\tEXPECTED [$exp_line_num] $exp_line \n";
                #   print LDAT "    status: passed \n\n";
			}
			
			$lc_expected++;
			$lc_output++;
			
    }

    close( DAT );
    close( EDAT );
    close( LDAT );

    return $stat;    
}


# ==============================================================================
#    Break line into seperate strings
#       1) Check each string for (float,scientific_notation,word)
sub compare_lines_using_tolerance_on_numbers
{
    my $line_test     = $_[0];
    my $line_expected = $_[1];
    my $tol = $_[2];

    my $status;
    my $r = -1;
    my $line_status = 'pass';
    
    
    
    # break into seperate strings
    chomp( $line_test );
    chomp( $line_expected );
    my $tt = trim( $line_test );
    my $te = trim( $line_expected );
#    $tline_test = $tt;
#    $line_expected = $te;

    
#    my @s_line_test = split( / /, $line_test);
#    my @s_line_expected = split( / /, $line_expected );
    my @s_line_test = split( / /, $tt);
    my @s_line_expected = split( / /, $te );
    
    # check if same length
    my $N1 = @s_line_test;
    my $N2 = @s_line_expected;
    if( $N1 != $N2 ) {
    #     print "    fail: Diff length  (@s_line_test) : (@s_line_expected) \n";
        $status = $status . "  fail: lines different length\n";
        $line_status = 'fail';
        return ($line_status,$status);
    }
    

    # For each word
    for( $i=0; $i<$N1; $i++ ) {
        $test = $s_line_test[$i];
        $exp = $s_line_expected[$i];
    
        # check if both are numbers
        $a = is_num( $test );
        $b = is_num( $exp );
        if( $a == 0 && $b == 0 ) { # both not numbers, try a literal string comparison
            if( $test eq $exp ) {
            #    $status = $status . "\n" . "  pass[$i]: strings ($test,$exp) match";
            }
            else {
                $status = $status . "  fail(word $i): strings ($test,$exp) differ\n";
                $line_status = 'fail';
            }
        }
        elsif( $a == 0 || $b == 0 ) {   # one number one not number
            print "    fail[$i]: trying to compare a number with a string ($test,$exp)\n";
            print " EXPECTED: $line_expected \n";
            print " MINE:         $line_test \n";
            $status =  $status ."  fail(word $i): trying to compare a number with a string ($test,$exp)\n";
            $line_status = 'fail';
        }
        else {
            $r = $test - $exp;
            $r = abs( $r );
            if( $r > $tol ) {
            #    print "    Failed ($test:$exp) : Res = $r > $tol \n";
                 $status = $status . "  fail(word $i): residual=$r > tolerance=$tol, for ($test:$exp)\n";
                 $line_status = 'fail';
            }
            else {
            #    print "    Pass ($test:$exp) : Res = $r < $tol \n";
                $status = $status . "  pass(word $i): residual ($r) < tolerance ($tol)\n";
            }
        }
    }
    
#    print "    $line_status \n";
    return ($line_status,$status);
}


sub check_for_scientific_notation
{
    $line = $_[0];
    $rm_dot_test = $line;
    $rm_dot_test =~ s/\.//g;


    if( $rm_dot_test =~ m/\d+(\+|\-){1}\d+/ ) {
#        print "Found sci notation $line \n";
         return 1;
    }
}

sub check_for_scientific_notation_E
{
    $line = $_[0];
    $rm_dot_test = $line;
    $rm_dot_test =~ s/\.//g;


    if( $rm_dot_test =~ m/\d+(E|e)(\+|\-)\d+/ ) {
#        print "Found sci notation_E $line \n";
        return 1;
    }
}


# Assumes (by default $v is decimal number, i.e. only contains digits
# Checks if $v is scientific of the form i)  xxxxx(E/e)(+/-)xxxxx
#                                            ii) xxxxxxx(+/-)xxxxx
# If $v is not scientific and contains any charcters, it just be a string

sub is_num 
{
    $v = $_[0];
    $r = 0;
    
    # strip out decimal point
    $v =~ s/\.//g;

    if( check_for_scientific_notation($v) ) {
        $r = 1;
    }
    elsif( check_for_scientific_notation_E($v) ) {
        $r = 1;
    }
    elsif( $v =~ m/\D/ ) {  # If it is not scientific AND contains any non-digit characters, its a string.
        $r = 0;
    }
    else {      # Must be a regular decimal number. We could check if it only contains digits here as well.
        $r = 1;
    }
#   if( $r ) { print "Found number = $_[0] \n"; }
    
    return $r;
}













# ==============================================================================
sub Comparison_catch_exception_pl
{
    my $job_output = $_[0];
    my $key_phrase = $_[1];
    my $log = $_[2];

    my $stat = 1;
    my $check = 'fail';
    
    # Check key_phrase is given
    my $L = length( $key_phrase );
    if( $L == 0 ) {
        print "ERROR: (catch_exception_pl) - You need to specify a key_phrase. \n";
        print "  e.g. TEST_TYPE = catch_exception_pl, key_phrase \n";
        die;
    }
    
    open( DAT, $job_output ) || die( "Could not open job output file ($job_output) ! \n" );
    my @data=<DAT>;


    # Look for keyphrase
    foreach $line (@data) {
        
        $check = check_for_error_message( $line, $key_phrase );
        if( $check eq 'pass' ) {
            $stat = 0;
            last;
        }
    }

    close( DAT );
    
    return $stat;
}

sub check_for_error_message
{
    my $line_test     = $_[0];
    my $line_expected = $_[1];

    my $status;

    # strip out any whitespace from start / end of test line
    chomp( $line_test );
    my $tline = trim( $line_test );

   chomp( $line_expected );
    my $tline_e = trim( $line_expected );

#    print "comparing    $tline     with      $line_expected \n";

    if( $tline =~ m/$tline_e/ ) {
        $status = 'pass';
    }
    else {
        $status = 'fail';
    }

    return ($status);
}

