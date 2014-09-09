#!/usr/bin/perl

sub check_line_for_pointer; # 32/64 bit
sub check_line_for_pointer32;
sub check_line_for_pointer32_64;

$test32 = '&i  = 0x0';
if( check_line_for_pointer( $test32 ) ) {
    print "Found ptr $test32 \n";
}
else {
    print "This is not a ptr $test32 \n";
}

$test32 = '&i  = 0xbffff944';
if( check_line_for_pointer( $test32 ) ) {
    print "Found ptr $test32 \n";
}
else {
    print "This is not a ptr $test32 \n";
}

$test32 = '&i  = cvbffff944';
if( check_line_for_pointer( $test32 ) ) {
    print "Found ptr $test32 \n";
}
else {
    print "This is not a ptr $test32 \n";
}

$test32 = '&i  = 0xbffff944bffff94';
if( check_line_for_pointer( $test32 ) ) {
    print "Found ptr $test32 \n";
}
else {
    print "This is not a ptr $test32 \n";
}


$test32 = '&i  = 0xbffff944bffff94';
if( check_line_for_pointer( $test32 ) ) {
    print "Found ptr $test32 \n";
}
else {
    print "This is not a ptr $test32 \n";
}



$test = '&i  = 0xbffff944bffff9453';
if( check_line_for_pointer32( $test ) ) {
    print "Found 32ptr $test \n";
}
else {
    print "This is not a 32ptr $test \n";
}
$test = '&i  = 0xbffff944 sdg';
if( check_line_for_pointer32( $test ) ) {
    print "Found 32ptr $test \n";
}
else {
    print "This is not a 32ptr $test \n";
}



print "\n\n";

####################################################

$test64 = '&i  = 0xbffff944bffff944';
if( check_line_for_pointer( $test64 ) ) {
    print "Found ptr $test64 \n";
}
else {
    print "This is not a ptr $test64 \n";
}

$test64 = '&i  = cvbffff944bffff944';
if( check_line_for_pointer( $test64 ) ) {
    print "Found 64ptr $test64 \n";
}
else {
    print "This is not a 64ptr $test64 \n";
}

$test64 = '&i  = 0xbffff944bffff94';
if( check_line_for_pointer( $test64 ) ) {
    print "Found 64ptr $test64 \n";
}
else {
    print "This is not a 64ptr $test64 \n";
}

print "\n\n";

#####################################################

$test = '0xbffff944';
if( check_line_for_pointer32_64( $test ) ) {
    print "Found ptr $test \n";
}
else {
    print "This is not a ptr $test \n";
}


$test = '0xbffff944bffff944';
if( check_line_for_pointer32_64( $test ) ) {
    print "Found ptr $test \n";
}
else {
    print "This is not a ptr $test \n";
}


$test = '&i  = 0xbffff944 sdg';
if( check_line_for_pointer32_64( $test ) ) {
    print "Found ptr $test \n";
}
else {
    print "This is not a ptr $test \n";
}


$test = '&i  = 0xbffff944bffff944';
if( check_line_for_pointer32_64( $test ) ) {
    print "Found ptr $test \n";
}
else {
    print "This is not a ptr $test \n";
}

$test = '&i  = 0xbffff94 sdg';
if( check_line_for_pointer32_64( $test ) ) {
    print "Found ptr $test \n";
}
else {
    print "This is not a ptr $test \n";
}


$test = '&i  = 0xbffff944bffff944sfdsg';
if( check_line_for_pointer32_64( $test ) ) {
    print "Found ptr $test \n";
}
else {
    print "This is not a ptr $test \n";
}





sub check_line_for_pointer
{
    $line = $_[0];
    return( $line =~ m/0x/ );
}

# ==
#  Should do check for 8 chars and white space AND 16 chars and white space
#  Thus catching both 32bit and 64 bit pointers
sub check_line_for_pointer32
{
    $line = $_[0];
    return( $line =~ m/0x.{8}\s/ );
}


# ==
#  Should do check for 8 chars and white space AND 16 chars and white space
#  Thus catching both 32bit and 64 bit pointers
#
#  Arggh, this one is fucked too cause we need to know if pointer is at end of line
#  Just stick with 0x anywhere to indicate a pointer - KISS

sub check_line_for_pointer32_64
{
    $line = $_[0];
    
    # for cases where pointer is last item on line
    $N = length( $line );
    if( $N == 10 || $N == 18 ) { # 0x + 8 ot 0x + 16
        return check_line_for_pointer( $line );
    }
    
    if( $line =~ m/0x.{8}\s/ || $line =~ m/0x.{16}\s/ ) {
        return 1;
    }
}



