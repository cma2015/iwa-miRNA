## --Makefile.PL for Bioperl core module distribution--##
## ----------------------------------------------------##
## $Id: Makefile.PL,v 1.69 2003/12/15 11:21:47 heikki Exp $
##----------------------------------------------------

#
# Basically this is a standard Makefile with a hash of
# needed packages which are tested at the start
#
# If you want help with this Makefile get in touch with
# the developers by sending a message to bioperl-l@bioperl.org.



# All Bioperl dependencies on external modules are defined here
%packages = (
	     'HTTP::Request::Common'    => '0.0/Part of the LWP package/GenBank+GenPept sequence retrieval, remote http Blast jobs/Bio::DB::*,Bio::Tools::Run::RemoteBlast',
	     'LWP::UserAgent'           => '0.0/Part of the LWP package/GenBank+GenPept sequence retrieval, remote http Blast jobs/Bio::DB::*,Bio::Tools::Run::RemoteBlast',
	     'HTML::Entities'           => '0.0/To parse POST request URL formas/Remote Analysis POST submissions/Bio::Tools::Run::Analysis::*',
	     'Ace'                      => '0.0/Aceperl/access of ACeDB database/Bio::DB::Ace',
	     'IO::Scalar'               => '0.0/IO handle to read or write to a scalar/remote http Blast jobs/Bio::Tools::Blast::Run::Webblast',
	     'IO::String'               => '0.0/IO handle to read or write to a string/GenBank+GenPept sequence retrieval, Variation code/Bio::DB::*,Bio::Variation::*,Bio::Tools::Blast::Run::Webblast, Bio::Index::Blast',
	     'XML::Parser'              => '0.0/Parsing of XML documents/Bio::Variation code, GAME parser/Bio::SeqIO::game,Bio::Variation::*',
	     'XML::Writer'              => '0.0/Parsing + writing of XML documents/Bio::Variation code, GAME parser/Bio::SeqIO::game,Bio::Variation::*',
	     'XML::Parser::PerlSAX'              => '0.0/Parsing of XML documents/Bio::Variation code, GAME parser/Bio::SeqIO::game,Bio::Variation::*',
	     'XML::Twig'         	=> '0.0/Available on CPAN/parsing of XML documents/Module Bio::Variation::IO::xml.pm',
	     'File::Temp'               => '0.0/Temporary File creation/Bio::Root::IO uses if present/Bio::DB::WebDBSeqI, Bio::Seq::LargePrimarySeq',
	     'SOAP::Lite'               => '0.0/SOAP protocol/XEMBL Services (also Bibliographic queries in Biblio::) /Bio::DB::XEMBLService',
	     'HTML::Parser'             => '0.0/HTML parsing/HTML parsing of GDB page/Bio::DB::GDB',
	     'DBD::mysql'               => '0.0/Mysql driver/loading and querying of Mysql-based GFF feature databases/Bio::DB::GFF',
	     'GD'                       => '1.3/Graphical Drawing Toolkit/Rendering Sequences and Features/Bio::Graphics',
	     'Text::Shellwords'         => '0.0/Execute shell commands/test scripts/Bio::Graphics',
	     'SVG'                      => '2.26/Generate optional SVG output/Bio::Graphics',
	     'GD::SVG'                  => '0.0/Generate optional SVG output/Bio::Graphics',
	     'Storable'                 => '0.0/Persistent object storage and retrieval/Storing sequence objects in local file cache/Bio::DB::FileCache',
	     'Graph::Directed'          => '0.0/Generic Graph data stucture and algorithms/Ontology Engine implementation for the GO parser/Bio::Ontology::SimpleOntologyEngine',
# Bio::Microarray is not in core bioperl yet
#	     'enum'                     => '0.0/C style enumerated types and bitmask flags in Perl/Bio::Microarray::Affymetrix::*',

# rare function let's not scare people
# 	     'srsperl.pm'               => '0.0/Sequence Retrieval System (SRS) perl API/alternative way of retrieving sequences/Bio::LiveSeq::IO::SRS.pm'
	     );

sub check_package {
    my ($name,$str) = @_;
    my ($ver,$desc,$expl,$module) = split( /\// , $str);
    if( !eval "require($name)") {
	print "External Module $name, $desc,\n is not installed on this computer.\n  The $module in Bioperl needs it for $expl\n\n";
	return 1;
    }

    return 0;
}

#
# Generate sub testing package system
#
#

# make mini test targets for playing around with
# in testing...

opendir(DIR,"t") || warn "No test directory. Weird!";
@files = readdir(DIR);
shift @files;
shift @files;
foreach my $file (@files) {
    $file =~ /(\w+)\.t/ || next;
    $file = $1;
    $line = "test_$file :: pure_all\n\tPERL_DL_NONLAZY=1 \$(FULLPERL) -I\$(INST_ARCHLIB) -I\$(INST_LIB) -I\$(PERL_ARCHLIB) -I\$(PERL_LIB) -e \'use Test::Harness qw(&runtests \$\$verbose); \$\$verbose=\$(TEST_VERBOSE); runtests \@ARGV;\' t/$file.t\n";
    push(@tline,$line);
    push(@targets,$file);
    $tset  .= "$file \\\n";
}

$tline = join('',@tline);
$tline = "show_tests :\n\t\@echo \'type make test_<subtest> to run\'\n\t\@echo '$tset'\n$tline\n";

print STDERR "Generated sub tests. go make show_tests to see available subtests\n";

#
# Prompt user for BioDBGFF stuff
#

sub prompt_for_biodbgff {
  if(  eval { require DBI } && grep {/mysql/} DBI->available_drivers) {
    print "\n** Bio::DB::GFF Tests section** \n";
    my $proceed = prompt("Do you want to run the BioDBGFF tests (requires a mysql, postgres or oracle database) y/n",'n');
    if( $proceed =~ /^[yY]/) {
      my $cfg = {dbd_driver => 'mysql'};
      GFFQueryDb($cfg,{'prompt'=>1,'verbose'=>1});
      if (open T,">t/do_biodbgff.tests") {
	foreach (keys %$cfg) {
	  print T $_,"\t",$cfg->{$_},"\n";
	}
	close T;
      }
    } else {
      if( -e "t/do_biodbgff.tests" ) {
	unlink "t/do_biodbgff.tests";
      }
    }
  }
}

sub MY::postamble {
  $tline . "\n\n";
}

sub get_scripts_to_install {
  print "\n*** Script Install Section ****\n";
  print "\nBioperl comes with a number of useful scripts which you may wish to install.\n";
  my $prompt_mode = prompt("Install [a]ll Bioperl scripts, [n]one, or choose groups [i]nteractively?",'a');
  return if $prompt_mode =~ /^[nN]/;
  my @scripts;
  my $script_dir = './scripts';
  my $dest_dir    = './scripts_temp';
  rmtree ($dest_dir) if -e $dest_dir;
  mkdir($dest_dir) or die "Can't mkdir($dest_dir): $!";
  opendir(F,$script_dir) or die "Can't opendir $script_dir: $!";
  while (my $file_or_dir = readdir(F)) {
    next if $file_or_dir =~ /^\./;
    next if $file_or_dir eq 'CVS';
    next unless -d "$script_dir/$file_or_dir";
    next unless prompt_to_install($script_dir,$file_or_dir,$prompt_mode);
    my @files = install_contents("$script_dir/$file_or_dir",$dest_dir);
    push @scripts,@files;
  }
  closedir F;
  return @scripts;
}

sub prompt_to_install {
  my $prefix = shift;
  my $f      = shift;
  my $prompt_mode = shift;
  print "\n* Script Directory $f *\n";
  if (-e "$prefix/$f/TAG" && (my $g = IO::File->new("$prefix/$f/TAG"))) {
    print while <$g>;
  }
  return 1 if $prompt_mode =~ /^[aA]/;  # all
  my $result = prompt("Install scripts in $f? [y]es [n]o [q]uit",'n');
  die if $result =~ /^[qQ]/;
  return $result =~ /^[yY]/;
}

sub install_contents {
  my $dir  = shift;
  my $dest = shift;
  my $bangline = $Config{startperl};

  my @files;
  opendir (D,$dir) or die "Can't open $dir: $!\n";
  while (my $script = readdir(D)) {
    next unless $script =~ /\.PLS$/;
    my $in  = IO::File->new("$dir/$script")    or die "Can't open $dir/$script: $!";
    $script =~ s/\.PLS$/\.pl/;                   # change from .PLS to .pl
    $script =~ s/^/bp_/ unless $script =~ /^bp/; # add the "bp" prefix
    print "\tActivating $script....\n";
    unlink "$dest/$script" if -e "$dest/$script";
    my $out = IO::File->new(">$dest/$script")  or die "Can't open $dest/$script: $!";
    my $doneit;
    while (<$in>) {
      next if $doneit;
      if (s/^\#\!\S+/$bangline/) {
	$doneit++;
      }
    } continue {
      print $out $_;
    }
    close $in;
    close $out;
    chmod 0555,"$dest/$script" or die "Can't change mode of $script to ",MODE,": $!";
    push @files,"$dest/$script";
  }
  closedir D;
  return @files;
}


############################################################################
#
#   Name:    GFFQueryDb
#
#   Purpose: Query user for settings for running the Bio::DB::GFF test suite
#
#   Inputs:  $cfg - Config hash ref
#
#   Returns: Nothing; creates
#            $cfg->{$driver}->{test_(db|host|dsn|user|pass)}
#
############################################################################

sub GFFQueryDb {
  my $cfg     = shift;
  my $options = shift || {};
  my $db     = $cfg->{'description'} || $cfg->{'dbd_driver'};
  my $driver = $cfg->{'dbd_driver'};

  my $prompt = $options->{'prompt'};

  $driver  = prompt("Which database driver should I use for testing BioDBGFF (mysql,Pg,Oracle)?",
		    $driver) if $prompt;

  my $d = lc $driver;
  $cfg->{'dbd_driver'} = $driver;

  my $test_db =  exists($options->{"$d-test-db"}) ?
    $options->{"$d-test-db"} : ($cfg->{'test_db'} || 'test');

  $test_db = prompt
    ("Which database should I use for testing the $driver driver?",
     $test_db) if $prompt;

  my $test_host = exists($options->{"$d-test-host"}) ?
    $options->{"$d-test-host"} : ($cfg->{'test_host'} || 'localhost');
  $test_host = prompt
    ("On which host is database $test_db running (hostname, ip address or host:port)", $test_host) if $prompt;

  my($test_user, $test_pass);

  $test_user = exists($options->{"$d-test-user"}) ?
    $options->{"$d-test-user"} : ($cfg->{'test_user'} || "undef");
  $test_user = prompt
    ("User name for connecting to database $test_db?", $test_user)
      if $prompt;
  $test_user = undef if $test_user eq 'undef';

  $test_pass = exists($options->{"$d-test-pass"}) ?
    $options->{"$d-test-pass"} : ($cfg->{'test_pass'} || "undef");
  $test_pass = prompt
    ("Password for connecting to database $test_db?", $test_pass)
      if $prompt;
  $test_pass = undef if $test_pass eq 'undef';

  $cfg->{'test_db'}   = $test_db;
  $cfg->{'test_host'} = $test_host;
  $cfg->{'test_user'} = $test_user;
  $cfg->{'test_pass'} = $test_pass;
  if ($test_host eq 'undef'  ||  $test_host eq 'localhost') {
    $test_host = '';
  }
  if ($d eq 'pg') {
    $test_dsn = "dbi:$driver:dbname=$test_db"; 
  } else {
    $test_dsn = "dbi:$driver:database=$test_db";
  }
  $cfg->{test_dsn} = $test_dsn;
  if ($test_host) {
    $cfg->{'test_dsn'} .= ";host=$test_host";
  }


  if ($options->{'verbose'}) {
    local $^W=0;  # some uninitialized variable warning coming through
    print("Driver $driver is using the following settings for tests:\n",
	  "    Database $test_db\n",
	  "    Host     $test_host\n",
	  "    DSN      $test_dsn\n",
	  "    User     $test_user\n",
	  "    Password $test_pass\n");
  }
}

				           
#
# Talk to Ewan (<birney@ebi.ac.uk>) or bioperl-l@bioperl.org for more info.
#

# Let the code begin...

require 5.005;
use ExtUtils::MakeMaker ;
use File::Path 'rmtree';
use IO::File;
use Config;

my @scripts_to_install = eval {get_scripts_to_install()};
prompt_for_biodbgff();


$do_autoload_finesse = 0;
$NAME     = 'Bio';
$DISTNAME = "bioperl";

if( $do_autoload_finesse == 1 ) {
    fatal("Autoload finessing does not work at the moment!");
}



@BACKUP_FILES = qw(./Seq.pm.bak);

#
# This does an inplace edit on Seq.pm
#

if( $do_autoload_finesse ) {

    print "Bioperl module installation\nPreparing modules for AUTOLOADing...";
    `perl -p -i.bak -e 's/#__END__/__END__/' Seq.pm`;
    print "done.\n\n";
}

$error = 0;
print "\n";
foreach $name ( keys %packages ) {
    if( &check_package($name,$packages{$name}) == 1 ) {
	$error = 1;
    }
}

if( $error == 1 ) {
    print <<QQ_ERROR_QQ;

Information:

   There are some external packages and perl modules, listed above, which 
   bioperl uses. This only effects the functionality which is listed above:
   the rest of bioperl will work fine, which includes nearly all of the
   core packages.

   The installation of these external packages is very simple. You
   can read more about bioperl external dependencies in the INSTALL
   file or at:

   http://bioperl.org/Core/Latest/INSTALL

   Enjoy the rest of bioperl, which you can use after going 'make install'

QQ_ERROR_QQ
}

WriteMakefile(
    NAME	=> $NAME,
    DISTNAME    => $DISTNAME,
    VERSION_FROM => "Bio/Root/Version.pm",
    'dist'      => {  COMPRESS     => 'gzip -9f',
                      SUFFIX       => '.gz', 
                      DIST_DEFAULT => 'all tardist',
                   },
    ($] ge '5.005')  ? (
     'AUTHOR'   => 'Bioperl Team (bioperl-l@bioperl.org)',
     'ABSTRACT' => 'Bioinformatics Toolkit',
    ) : (),
    'realclean' => {  FILES        => join(' ',@BACKUP_FILES) },
    DIR         => [  ], # directories with other Makefile.PL's in them
    EXE_FILES   => \@scripts_to_install,
    PL_FILES	=> { 'doc/makedoc.PL' => ['bioscripts.pod','biodatabases.pod','biodesign.pod','bioperl.pod'] },
    PREREQ_PM   => {
	'IO::String'       => 0,
    	'File::Temp'       => 0,
	'File::Spec'       => 0,
	'HTML::Entities'   => 0,
	'DB_File'          => 0,
	'IO::Scalar'       => 0,
    },
);


sub MY::manifypods {
  my $self = shift;
  #print STDERR "In manifypods moment\n";
  if( 1 ) {
    return "\nmanifypods : pure_all\n\t$self->{NOECHO}\$(NOOP)\n"
  }
  else {
    return $self->SUPER::manifypods(@_);
  }
}

sub MY::clean {
  package MY;
  return shift->SUPER::clean(@_) . <<END;
	-rm -rf scripts_temp
END
}
