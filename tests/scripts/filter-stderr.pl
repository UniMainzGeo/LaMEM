#!/usr/bin/perl

sub check_for_error_message; # returns (pass or fail)
sub trim;

#
# Given a file from a job which crashed, we wish to examine it for a particular
# assert statement, expected or error message
#


# As a test, read in file and scan contents line by line

my $error_message = 'Matrix ARG symbolic LU!';
chomp( $error_message );
my $terror_message = trim( $error_message );

$data_file = 'error.out';
open( DAT, $data_file ) || die( "Could not open file ($data_file) ! \n" );
@raw_data=<DAT>;

my $N = 0;
foreach $line (@raw_data) {

    $status = check_for_error_message( $line, $terror_message );
    if( $status eq 'pass' ) {
        print "status=pass - Found \"$terror_message\" at line $N \n";
    }
    $N++;
}
close( DAT );




sub check_for_error_message
{
    my $line_test     = $_[0];
    my $line_expected = $_[1];

    my $status;

    # strip out any whitespace from start / end of test line
    chomp( $line_test );
    my $tline = trim( $line_test );

#    print "comparing    $tline     with      $line_expected \n";

    if( $tline =~ m/$line_expected/ ) {
        $status = 'pass';
    }
    else {
        $status = 'fail';
    }

    return ($status);
}


sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
