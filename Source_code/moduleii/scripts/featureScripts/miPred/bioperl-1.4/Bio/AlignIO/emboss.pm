# $Id: emboss.pm,v 1.11 2002/10/22 07:45:10 lapp Exp $
#
# BioPerl module for Bio::AlignIO::emboss
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::AlignIO::emboss - Parse EMBOSS alignment output (from applications water and needle)

=head1 SYNOPSIS

    # do not use the object directly
    use Bio::AlignIO;
    # read in an alignment from the EMBOSS program water
    my $in = new Bio::AlignIO(-format => 'emboss',
			      -file   => 'seq.water');
    while( my $aln = $in->next_aln ) {
	# do something with the alignment
    }

=head1 DESCRIPTION

This object handles parsing and writing pairwise sequence alignments
from the EMBOSS suite.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::AlignIO::emboss;
use vars qw(@ISA $EMBOSSTitleLen $EMBOSSLineLen);
use strict;

use Bio::AlignIO;
use Bio::LocatableSeq;

@ISA = qw(Bio::AlignIO );

BEGIN { 
    $EMBOSSTitleLen    = 13;
    $EMBOSSLineLen     = 50;
}

sub _initialize {
    my($self,@args) = @_;
    $self->SUPER::_initialize(@args);
    $self->{'_type'} = undef;
}

=head2 next_aln

 Title   : next_aln
 Usage   : $aln = $stream->next_aln()
 Function: returns the next alignment in the stream.
 Returns : L<Bio::Align::AlignI> object - returns 0 on end of file
	    or on error
 Args    : NONE

=cut

sub next_aln {
    my ($self) = @_;
    my $seenbegin = 0;
    my %data = ( 'seq1' => { 
		     'start'=> undef,
		     'end'=> undef,		
		     'name' => '',
		     'data' => '' },
		 'seq2' => { 
		     'start'=> undef,
		     'end'=> undef,
		     'name' => '',
		     'data' => '' },
		 'align' => '',
		 'type'  => $self->{'_type'},  # to restore type from 
		                                     # previous aln if possible
		 );
    my %names;
    while( defined($_ = $self->_readline) ) {
	next if( /^\#?\s+$/ || /^\#+\s*$/ );
	if( /^\#(\=|\-)+\s*$/) {
	    last if( $seenbegin);
	} elsif( /(Local|Global):\s*(\S+)\s+vs\s+(\S+)/ ||
		 /^\#\s+Program:\s+(\S+)/ )
	{
	    my ($name1,$name2) = ($2,$3);
	    if( ! defined $name1 ) { # Handle EMBOSS 2.2.X
		$data{'type'} = $1;
		$name1 = $name2 = '';
	    } else { 
		$data{'type'} = $1 eq 'Local' ? 'water' : 'needle';
	    }	    
	    $data{'seq1'}->{'name'} = $name1;
	    $data{'seq2'}->{'name'} = $name2;

	    $self->{'_type'} = $data{'type'};

	} elsif( /Score:\s+(\S+)/ ) {
	    $data{'score'} = $1;		
	} elsif( /^\#\s+(1|2):\s+(\S+)/ && !  $data{"seq$1"}->{'name'} ) {
	    my $nm = $2;
	    $nm = substr($nm,0,$EMBOSSTitleLen); # emboss has a max seq length
	    if( $names{$nm} ) {
		$nm .= "-". $names{$nm};
	    }
	    $names{$nm}++;
	    $data{"seq$1"}->{'name'} = $nm;	
	} elsif( $data{'seq1'}->{'name'} &&
		 /^$data{'seq1'}->{'name'}/ ) {	    
	    my $count = 0;
	    $seenbegin = 1;
	    my @current;
	    while( defined ($_) ) {
		my $align_other = '';
		my $delayed;		
		if($count == 0 || $count == 2 ) {
		    my @l = split;
		    my ($seq,$align,$start,$end);
		    if( $count == 2 && $data{'seq2'}->{'name'} eq '' ) {
			# weird boundary condition 
			($start,$align,$end) = @l;
		    } elsif( @l == 3 ) {
			$align = '';
			($seq,$start,$end) = @l
		    } else { 
			($seq,$start,$align,$end) = @l;
 		    }

		    my $seqname = sprintf("seq%d", ($count == 0) ? '1' : '2'); 
		    $data{$seqname}->{'data'} .= $align;
		    $data{$seqname}->{'start'} ||= $start;
		    $data{$seqname}->{'end'} = $end;
		    $current[$count] = [ $start,$align || ''];
		} else { 
		    s/^\s+//;
		    s/\s+$//;
		    $data{'align'} .= $_;
		}

	      BOTTOM:
		last if( $count++ == 2);
		$_ = $self->_readline();
	    }

	    if( $data{'type'} eq 'needle' ) {
		# which ever one is shorter we want to bring it up to 
		# length.  Man this stinks.
		my ($s1,$s2) =  ($data{'seq1'}, $data{'seq2'});
		
		my $d = length($current[0]->[1]) - length($current[2]->[1]);
		if( $d < 0 ) { # s1 is smaller, need to add some
		    # compare the starting points for this alignment line
		    if( $current[0]->[0] <= 1 && $current[2]->[0] > 1) { 
			$s1->{'data'} = ('-' x abs($d)) . $s1->{'data'};
			$data{'align'} = (' 'x abs($d)).$data{'align'};
		    } else { 
			$s1->{'data'} .= '-' x abs($d);
			$data{'align'} .= ' 'x abs($d);
		    }
		} elsif( $d > 0) { # s2 is smaller, need to add some  
		    if( $current[2]->[0] <= 1 && $current[0]->[0] > 1) { 
			$s2->{'data'} = ('-' x abs($d)) . $s2->{'data'};
			$data{'align'} = (' 'x abs($d)).$data{'align'};
		    } else { 
			$s2->{'data'} .= '-' x abs($d);
			$data{'align'} .= ' 'x abs($d);
		    }
		}
	    }
	    
	}
    }
    return undef unless $seenbegin;
    my $aln =  Bio::SimpleAlign->new(-verbose => $self->verbose(),
				     -source => "EMBOSS-".$data{'type'});
    
    foreach my $seqname ( qw(seq1 seq2) ) { 
	return undef unless ( defined $data{$seqname} );	
	$data{$seqname}->{'name'} ||= $seqname;
	my $seq = new Bio::LocatableSeq('-seq' => $data{$seqname}->{'data'},
					'-id'  => $data{$seqname}->{'name'},
					'-start'=> $data{$seqname}->{'start'},
					'-end' => $data{$seqname}->{'end'},
					);
	$aln->add_seq($seq);
    }
    return $aln;
}

=head2 write_aln

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: writes the $aln object into the stream in emboss format
 Returns : 1 for success and 0 for error
 Args    : L<Bio::Align::AlignI> object


=cut

sub write_aln {
    my ($self,@aln) = @_;

    $self->throw("Sorry: writing emboss output is not currently available! \n");
}

1;
