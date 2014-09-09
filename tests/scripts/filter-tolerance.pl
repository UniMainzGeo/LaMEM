#!/usr/bin/perl

sub check_for_tolerance; # (pass,res), (fail,res), (skip,-1)
sub check_for_scientific_notation;
sub check_for_scientific_notation_E;
sub is_num;
sub trim;
sub check_for_tolerance_2;
sub compare_lines_using_tolerance_on_numbers;   # This is the one to go with. Does rel tolerance on numbers and diff on strings


print "\n\nTESTING: scientific notation \n";

check_for_scientific_notation( '1.0e-5' );
check_for_scientific_notation( '1.0E+9' );
check_for_scientific_notation( '1.0Ex' );
check_for_scientific_notation( '10.9E+312412' );
check_for_scientific_notation( '10.9e312412' );
check_for_scientific_notation( '10.9e--312412' );
check_for_scientific_notation( '10.9ee312412' );

check_for_scientific_notation( '1.0-5' );
check_for_scientific_notation( '1.0+9' );
check_for_scientific_notation( '1.0x' );
check_for_scientific_notation( '10.9+312412' );
check_for_scientific_notation( '10.9312412' );
check_for_scientific_notation( '10.9--312412' );
check_for_scientific_notation( '10.9++312412' );


check_for_scientific_notation_E( '1.0e-5' );
check_for_scientific_notation_E( '1.0E+9' );
check_for_scientific_notation_E( '1.0Ex' );
check_for_scientific_notation_E( '10.9E+312412' );
check_for_scientific_notation_E( '10.9e312412' );
check_for_scientific_notation_E( '10.9e--312412' );
check_for_scientific_notation_E( '10.9ee312412' );

check_for_scientific_notation_E( '1.0-5' );
check_for_scientific_notation_E( '1.0+9' );
check_for_scientific_notation_E( '1.0x' );
check_for_scientific_notation_E( '10.9+312412' );
check_for_scientific_notation_E( '10.9312412' );
check_for_scientific_notation_E( '10.9--312412' );
check_for_scientific_notation_E( '10.9++312412' );



is_num( '1.0-5' );
is_num( '1.0+9' );
is_num( '1.0x' );
is_num( '10.9+312412' );
is_num( '10.9312412' );
is_num( '10.9--312412' );
is_num( '10.9++312412' );

is_num( '1.0E-5' );
is_num( '1.0e+9' );
is_num( '1.0x' );
is_num( '10.9e+312412' );
is_num( '10.9312412' );
is_num( '10.9E--312412' );
is_num( '10.9e++312412' );



$tol = 1.0e-3;


print "\n\nTESTING: check_for_tolerance \n";

($stat,$r) = check_for_tolerance( '0.00012', '2359.454', $tol );
print "Status: $stat : $r > $tol \n";

($stat,$r) = check_for_tolerance( '0.00012', '0.00', $tol );
print "Status: $stat : $r > $tol \n";

($stat,$r) = check_for_tolerance( '0.0r0012', '0.v00', $tol );
print "Status: $stat : $r > $tol \n";

($stat,$r) = check_for_tolerance( '12.999', '13', $tol );
print "Status: $stat : $r > $tol \n";

($stat,$r) = check_for_tolerance( '12.999 11.999', '13 12', $tol );
print "Status: $stat : $r > $tol \n";


($stat,$r) = check_for_tolerance( '1.0e-5', '1.0e-3', $tol );
print "Status: $stat : $r > $tol \n";




print "\n\nTESTING: check_for_tolerance_2 \n";
($line_stat,$stat) = check_for_tolerance_2( '0.00012', '2359.454', $tol );
print "Line status : $line_stat -- Status: $stat \n";

($line_stat,$stat) = check_for_tolerance_2( '0.00012', '0.00', $tol );
print "Line status : $line_stat -- Status: $stat \n";

($line_stat,$stat) = check_for_tolerance_2( '0.0r0012', '0.v00', $tol );
print "Line status : $line_stat -- Status: $stat \n";

($line_stat,$stat) = check_for_tolerance_2( '12.999', '13', $tol );
print "Line status : $line_stat -- Status: $stat \n";

($line_stat,$stat) = check_for_tolerance_2( '12.999 11.999', '13 12', $tol );
print "Line status : $line_stat -- Status: $stat \n";

($line_stat,$stat) = check_for_tolerance_2( '1.0-5', '1.0-3', $tol );
print "Line status : $line_stat -- Status: $stat \n";


($line_stat,$stat) = check_for_tolerance_2( '12.999 11.999 a b c 3.300001', '13 12 a d e 3.3', $tol );
print "Line status : $line_stat -- Status: $stat \n";

($line_stat,$stat) = check_for_tolerance_2( '12.999 11.999 a b c 3.300001', '13 12 a d e j', $tol );
print "Line status : $line_stat -- Status: $stat \n";


print "\n\nTESTING: compare_lines_using_tolerance_on_numbers \n";
($line_stat,$stat) = compare_lines_using_tolerance_on_numbers( '12.999 11.999 a b c 3.300001', '13 12 a d e 3.3', $tol );
print "Line status : $line_stat -- Status: $stat \n";

($line_stat,$stat) = compare_lines_using_tolerance_on_numbers( '12.999 11.999 a b c 3.300001', '13 12 a d e j', $tol );
print "Line status : $line_stat -- Status: $stat \n";

($line_stat,$stat) = compare_lines_using_tolerance_on_numbers( 
        '12.999 11.999 a d e j', 
        '13 12 a d e j', 
        $tol );
print "Line status : $line_stat -- Status: $stat \n";

($line_stat,$stat) = compare_lines_using_tolerance_on_numbers( 
        '12.999 11.999 a d e j', 
        '13 12 a d e j', 
        $tol );
print "Line status : $line_stat -- Status: $stat \n";

print "\n\nTest gmres with cg \n";
    # gmres         #cg
($line_stat,$stat) = compare_lines_using_tolerance_on_numbers(  '0 KSP Residual norm 3.811398540713e+00',          '0 KSP Residual norm 3.811398540713e+00', $tol );
print "Line status : $line_stat -- Status: $stat \n";
($line_stat,$stat) = compare_lines_using_tolerance_on_numbers(  '1 KSP Residual norm 1.444731352039e+00',          '1 KSP Residual norm 1.445199464584e+00', $tol ); 
print "Line status : $line_stat -- Status: $stat \n";
($line_stat,$stat) = compare_lines_using_tolerance_on_numbers(  '2 KSP Residual norm 4.712840878823e-01',          '2 KSP Residual norm 4.752140927010e-01' , $tol );
print "Line status : $line_stat -- Status: $stat \n";
($line_stat,$stat) = compare_lines_using_tolerance_on_numbers(  '3 KSP Residual norm 6.943414362981e-02',          '3 KSP Residual norm 6.954360130519e-02' , $tol );
print "Line status : $line_stat -- Status: $stat \n";
($line_stat,$stat) = compare_lines_using_tolerance_on_numbers(  '4 KSP Residual norm 9.502667693381e-03',          '4 KSP Residual norm 9.655213731885e-03' , $tol );
print "Line status : $line_stat -- Status: $stat \n";
($line_stat,$stat) = compare_lines_using_tolerance_on_numbers(  '5 KSP Residual norm 7.905767061099e-04',          '5 KSP Residual norm 7.984952761226e-04' , $tol );
print "Line status : $line_stat -- Status: $stat \n";
($line_stat,$stat) = compare_lines_using_tolerance_on_numbers(  '6 KSP Residual norm 1.108383257088e-04',          '6 KSP Residual norm 1.115808341851e-04' , $tol );
print "Line status : $line_stat -- Status: $stat \n";




# ==============================================================================

sub check_for_tolerance
{
    $line_test     = $_[0];
    $line_expected = $_[1];
    $tol = $_[2];

    my $status;
    my $r;
    
    # strip out decimal point.    
    $rm_dot_test = $line_test;
    $rm_dot_exp = $line_expected;
    $rm_dot_test =~ s/\.//g;
    $rm_dot_exp =~ s/\.//g;
#    print "$rm_dot_test and $rm_dot_exp \n";
   # if the test line is contains an characters don't do a numeric test
    if( $rm_dot_test =~ m/\D/ || $rm_dot_exp =~ m/\D/ ) {
        $status = 'skip';
        $r = -1;
        
        return ($status,$r);
    }
    # We should allow for scientific notation XXXe+XXX and XXXe-XXX



    $r = $line_test - $line_expected;
    $r = abs( $r );
    if( $r > $tol ) {
    #    print "Failed: Res = $r > $tol \n";
        $status = 'fail';
    }
    else {
    #    print "Pass: Res = $r < $tol \n";
        $status = 'pass';
    }

    return ($status,$r);
    
}



# ==============================================================================
#    Break line into seperate strings
#       1) Check each string for (float,scientific_notation,word)
sub check_for_tolerance_2
{
    $line_test     = $_[0];
    $line_expected = $_[1];
    $tol = $_[2];

    my $status;
    my $r = -1;
    my $line_status = 'pass';
    
    # break into seperate strings
    @s_line_test = split( / /, $line_test);
    @s_line_expected = split( / /, $line_expected );
    
    # check if same length
    $N1 = @s_line_test;
    $N2 = @s_line_expected;
    if( $N1 != $N2 ) {
        $status = $status . "\n" . "fail: lines different length";
        $line_status = 'fail';
        return ($line_status,$status);
    }
    
    
    for( $i=0; $i<$N1; $i++ ) {
        # strip out decimal point from each string
        $test = $s_line_test[$i];
        $exp = $s_line_expected[$i];
    
        # check if either contains any characters
        $a = is_num( $test );
        $b = is_num( $exp );
        if( $a == 0 && $b == 0 ) {
            $status = $status . "\n" . "skip[$i]: both strings ($test,$exp)";
        }
        elsif( $a == 0 || $b == 0 ) {
            $status =  $status . "\n" ."fail[$i]: trying to compare a number with a string ($test,$exp)";
            $line_status = 'fail';
        }
        else {
            $r = $test - $exp;
            $r = abs( $r );
            if( $r > $tol ) {
            #    print "Failed: Res = $r > $tol \n";
                 $status = $status . "\n" . "fail[$i]: residual ($r) > tolerance ($tol)";
            }
            else {
            #    print "Pass: Res = $r < $tol \n";
                $status = $status . "\n" . "pass[$i]";
            }
        }
    }
    
    return ($line_status,$status);
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
    $tline_test = $tt;
    $line_expected = $te;

    
    my @s_line_test = split( / /, $line_test);
    my @s_line_expected = split( / /, $line_expected );
    
    # check if same length
    my $N1 = @s_line_test;
    my $N2 = @s_line_expected;
    if( $N1 != $N2 ) {
        $status = $status . "\n" . "fail: lines different length";
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
                $status = $status . "\n" . "pass[$i]: strings ($test,$exp) match";
            }
            else {
                $status = $status . "\n" . "fail[$i]: strings ($test,$exp) differ";
                $line_status = 'fail';
            }
        }
        elsif( $a == 0 || $b == 0 ) {   # one number one not number
            $status =  $status . "\n" ."fail[$i]: trying to compare a number with a string ($test,$exp)";
            $line_status = 'fail';
        }
        else {
            $r = $test - $exp;
            $r = abs( $r );
            if( $r > $tol ) {
            #    print "Failed: Res = $r > $tol \n";
                 $status = $status . "\n" . "fail[$i]: residual ($r) > tolerance ($tol)";
                 $line_status = 'fail';
            }
            else {
            #    print "Pass: Res = $r < $tol \n";
                $status = $status . "\n" . "pass[$i]: residual ($r) < tolerance ($tol)";
            }
        }
    }
    
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


sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
