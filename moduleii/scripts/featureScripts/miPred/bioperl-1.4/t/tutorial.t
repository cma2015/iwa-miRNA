#-*-Perl-*- mode
# $Id: tutorial.t,v 1.1 2003/03/06 18:20:59 jason Exp $

BEGIN {
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    use vars qw($NUMTESTS);
    $NUMTESTS = 21;
    plan tests => $NUMTESTS;
    @ARGV = (-1);
    require 'bptutorial.pl';
}

END {
    unlink 'bptutorial.out';
}

# run the first 21 tests
for my $test ( 1..21 ) {
    ok(&run_examples($test));
}


