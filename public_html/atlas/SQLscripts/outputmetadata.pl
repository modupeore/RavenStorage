#!/usr/bin/perl

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MANUAL FOR outputmetadata.pl

=pod

=head1 NAME

$0 -- Save all the metadata information from the database : transcriptatlas

=head1 SYNOPSIS

outputmetadata.pl [--help] [--manual]

=head1 OPTIONS

=over 3

=item B<-1, -a, in>

List of libraries numbers separated by comma (e.g 123,975).  (Required)



=back

=head1 DEPENDENCIES

Requires the following Perl libraries (all standard in most Perl installs).
   DBI
   Getopt::Long
   Pod::Usage

=head1 AUTHOR

Written by Modupe Adetunji, 
Animal and Food Sciences Department, University of Delaware.

=head1 REPORTING BUGS

Report bugs to amodupe@udel.edu

=head1 COPYRIGHT

Copyright 2015 MOA.  

=cut

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - U S E R  V A R I A B L E S- - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# CODE FOR genelistquery.pl

use strict;
use DBI;
use Getopt::Long;
use Pod::Usage;

#ARGUMENTS
my($specifics,$output,$metadata, $sequence,$help,$manual);
GetOptions(
				"1|a|in|in1|list=s"	=>	\$specifics,
				"2|b|out|output=s"        =>      \$output,
                        "m"        =>      \$metadata,
                        "z"        =>      \$sequence,
				"h|help"        =>      \$help,
				"man|manual"	=>      \$manual );

# VALIDATE ARGS
pod2usage( -verbose => 2 )  if ($manual);
pod2usage( -verbose => 1 )  if ($help);


# DATABASE ATTRIBUTES
my $dsn = 'dbi:mysql:transcriptatlas';
my $user = 'frnakenstein';
my $passwd = 'maryshelley';

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - G L O B A L  V A R I A B L E S- - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# DATABASE VARIABLES
my ($dbh, $sth, $syntax, @row);

# CONNECT TO THE DATABASE
#open (OUT,">/home/modupe/genelist.txt");
$dbh = DBI->connect($dsn, $user, $passwd) or die "Connection Error: $DBI::errstr\n";

# OPENING OUTPUT FILE
open (OUT, ">$output")  or die "it aint opening";


#SPECIFYING LIBRARIES OF INTEREST
my @headers = split("\,", $specifics);

if ($metadata) {
  # HEADER print out
  print OUT "library_id\tbird_id\tspecies\tline\ttissue\tmethod\t";
  print OUT "index\tchip_result\tscientist\tdate\tnotes\n";

  $syntax = "select * from bird_libraries where library_id in 
        ($specifics);";
  $sth = $dbh->prepare($syntax);
  $sth->execute or die "SQL Error: $DBI::errstr\n";
  
  while (@row = $sth->fetchrow_array() ) { 
        foreach my $real (0..$#row-1){
          print OUT "$row[$real]\t";
        }
        print OUT "$row[$#row]\n";
  }
}
if ($sequence) {
  #HEADER print
  print OUT "Library id\tLine\tSpecies\tTissue\tTotal reads\tMapped reads\tGenes\tIsoforms\tVariants\tSNPs\tINDELs\tNotes\n";
  
  $syntax = "select v.library_id,v.line,v.species, v.tissue, t.total_reads, v.mapped_reads, v.genes, v.isoforms,
              v.total_VARIANTS VARIANTS,v.total_SNPs SNPs, v.total_INDELs INDELs, v.notes from transcripts_summary as t
              join vw_libraryinfo as v on t.library_id = v.library_id where v.library_id in 
              ($specifics);";
  $sth = $dbh->prepare($syntax);
  $sth->execute or die "SQL Error: $DBI::errstr\n";
  
  while (@row = $sth->fetchrow_array() ) { 
        foreach my $real (0..$#row-1){
          print OUT "$row[$real]\t";
        }
        print OUT "$row[$#row]\n";
  }
}
# DISCONNECT FROM THE DATABASE
$dbh->disconnect();
close(OUT);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - -T H E  E N D - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
exit;
