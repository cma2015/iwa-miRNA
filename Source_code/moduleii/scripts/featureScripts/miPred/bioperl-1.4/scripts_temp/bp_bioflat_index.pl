#!/usr/bin/perl -w
#$Id: bioflat_index.PLS,v 1.1 2003/11/20 12:39:02 bosborne Exp $

use strict;
use Bio::Root::Root;
use Bio::Root::IO;
use Bio::DB::Flat;
use Getopt::Long;
use File::Path qw(mkpath rmtree);

use constant USAGE => <<END;
Usage: $0 <options> file1 file2 file3...

Create or update a biological sequence database indexed with the
Bio::DB::Flat indexing scheme.  The arguments are a list of flat files
containing the sequence information to be indexed.

 Options:
     --create              Create or reinitialize the index.  If not specified,
                           the index must already exist.

     --format   <format>   The format of the sequence files.  Must be one
                           of "genbank", "swissprot", "embl" or "fasta".

     --location <path>     Path to the directory in which the index files
                           are stored.

     --dbname <name>       The symbolic name of the database to be created.

     --indextype <type>    Type of index to create.  Either "bdb" or "flat".
                           "binarysearch" is the same as "flat".

Options can be abbreviated.  For example, use -i for --indextype.

The following environment variables will be used as defaults if the 
corresponding options are not provided:

     OBDA_FORMAT      format of sequence file
     OBDA_LOCATION    path to directory in which index files are stored
     OBDA_DBNAME      name of database
     OBDA_INDEX       type of index to create

END
;

my ($CREATE,$FORMAT,$LOCATION,$DBNAME,$INDEXTYPE);

GetOptions( 'create'      => \$CREATE,
	    'format:s'    => \$FORMAT,
	    'location:s'  => \$LOCATION,
	    'dbname:s'    => \$DBNAME,
	    'indextype:s' => \$INDEXTYPE ) or die USAGE;

$FORMAT    = $ENV{OBDA_FORMAT}    unless defined $FORMAT;
$LOCATION  = $ENV{OBDA_LOCATION}  unless defined $LOCATION;
$DBNAME    = $ENV{OBDA_DBNAME}    unless defined $DBNAME;
$INDEXTYPE = $ENV{OBDA_INDEXTYPE} unless defined $INDEXTYPE;

my $root = 'Bio::Root::Root';
my $io   = 'Bio::Root::IO';

# confirm that database directory is there
defined $LOCATION or $root->throw("please provide a base directory with the --location option");
-d $LOCATION or $root->throw("$LOCATION is not a valid directory; use --create to create a new index");

defined $DBNAME or $root->throw("please provide a database name with the --dbname option");

defined $FORMAT or $root->throw("please specify the format for the input files with the --format option");

defined $INDEXTYPE or ($INDEXTYPE='flat'
		       && $root->warn('setting index type to "flat", use the --indextype option to override'));

# Confirm that database is there and that --create flag is sensible.
my $path = $io->catfile($LOCATION,$DBNAME,'config.dat');
if (-e $path) {
  if ($CREATE) {
    $root->warn("existing index detected; deleting.");
    rmtree($io->catfile($LOCATION,$DBNAME),1,1);
  } else {
    $root->warn("existing index detected; ignoring --indextype and --format options.");
    undef $INDEXTYPE;
  }
}
elsif (!$CREATE) {
  $root->throw("Cannot find database config file at location $path; use --create to create a new index");
}

# open for writing/updating
my $db = Bio::DB::Flat->new(-directory  => $LOCATION,
			    -dbname     => $DBNAME,
			    $INDEXTYPE ? (
					  -index      => $INDEXTYPE
					  )
			               : (),
			    -write_flag => 1,
			    -format     => $FORMAT) or $root->throw("can't create Bio::DB::Flat object");
my $entries = $db->build_index(@ARGV);

print STDERR "(Re)indexed $entries entries.\n ";
1;
