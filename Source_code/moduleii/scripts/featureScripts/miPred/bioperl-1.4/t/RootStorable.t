# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id: RootStorable.t,v 1.2 2003/09/05 18:00:37 heikki Exp $

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;

BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't', '.';
    }
    use Test;
    plan tests => 34;

}

$| = 1;

use Bio::Root::Storable;

foreach my $mode( "BINARY", "ASCII" ){
    if( $mode eq "ASCII" ){
        $Bio::Root::Storable::BINARY = 0;
    }

    #------------------------------
    # Test the easy bits that don't need file IO
    my $obj = Bio::Root::Storable->new();
    ok defined($obj) && $obj->isa('Bio::Root::Storable');

    eval { $obj->throw('Testing throw') };
    ok $@ =~ /Testing throw/;   # 'throw failed';

    $obj->{_test}  = "_TEST";   # Provide test attributes
    $obj->{__test} = "__TEST";  # 

    my $state = $obj->serialise;
    ok length($state) > 0;

    my $clone = $obj->clone;
    ok defined($clone) and $clone->isa('Bio::Root::Storable');
    ok $clone->{_test} eq "_TEST" && $clone->{__test}  eq "__TEST";

    #------------------------------
    # Test standard file IO 
    my $file = $obj->store;
    ok $file && -f $obj->statefile;

    my $retrieved = Bio::Root::Storable->retrieve( $file );
    ok defined($retrieved) && $retrieved->isa('Bio::Root::Storable');
    ok $retrieved->{_test} eq "_TEST" && ! exists $retrieved->{__test};

    my $skel = $obj->new_retrievable;
    ok defined($skel) && $skel->isa('Bio::Root::Storable');
    ok ! exists $skel->{_test} && ! exists $skel->{__test};
    ok $skel->retrievable;

    $skel->retrieve;
    ok ! $skel->retrievable;
    ok $skel->{_test} eq "_TEST" && ! exists $skel->{__test};

    my $obj2 = Bio::Root::Storable->new();
    $obj2->template('TEST_XXXXXX');
    $obj2->suffix('.state');
    my $file2 = $obj2->store;
    ok $file2 =~ /TEST_\w{6}?\.state$/ and -f $file2;

    #------------------------------
    # Test recursive file IO
    $obj->{_test_lazy} = $obj2;
    $obj->store;
    my $retrieved2 = Bio::Root::Storable->retrieve( $obj->token );
    ok $retrieved2->{_test_lazy} && $retrieved2->{_test_lazy}->retrievable;

    #------------------------------
    # Clean up
    # Should only be 2 object files; all others were clones in one way or another
    $obj->remove;
    ok ! -f $obj->statefile;
    $obj2->remove;
    ok ! -f $obj2->statefile;
}

1;
