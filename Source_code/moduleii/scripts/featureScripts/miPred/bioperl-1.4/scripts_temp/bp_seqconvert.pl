#!/usr/bin/perl perl
# $Id: seqconvert.PLS,v 1.2 2003/04/25 13:35:41 bosborne Exp $

use strict;
use Getopt::Long;
use Bio::SeqIO;

my $help;
my $from=undef;
my $to=undef;

### please add to this list (see the modules under Bio/SeqIO):
my @known_formats=
  qw(gcg fasta ace raw fastq phd pir scf swiss genbank locuslink
     embl game qual bsml tab raw abi chado alf ctf exp ztr pln);

my $script=substr($0, 1+rindex($0,'/'));
my $usage="Usage:

  $script --from in-format --to out-format < file.in-format > file.out-format

Known formats:\n  " . join(' ', @known_formats) . "\n\n";

die $usage unless
  &GetOptions( 'from:s'   => \$from,
               'to:s'     => \$to,
               'h|help'   => \$help
	     )
  && !$help &&  $from && $to
  && grep($from eq $_, @known_formats)
  && grep($to eq $_, @known_formats);

my $in  = Bio::SeqIO->newFh(-fh => \*STDIN , '-format' => $from);
my $out = Bio::SeqIO->newFh(-fh=> \*STDOUT, '-format' => $to);

print $out $_ while <$in>;


__END__

=head1 NAME

seqconvert - generic BioPerl sequence format converter

=head1 SYNOPSIS

  seqconvert --from in-format --to out-format < file.in-format > file.out-format
  # or
  seqconvert -f in-format -t out-format < file.in-format > file.out-format

=head1 DESCRIPTION

This script gives command line interface to BioPerl Bio::SeqIO.

=head1 SEE ALSO

L<Bio::SeqIO>

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

=head1 AUTHOR - Philip Lijnzaad

Email p.lijnzaad@med.uu.nl

=cut
