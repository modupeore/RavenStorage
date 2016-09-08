#!/usr/bin/perl

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MANUAL FOR outputgenequery.pl

=pod

=head1 NAME

$0 -- Extract genes information from specified library_id in the mysql database : transcriptatlas

=head1 SYNOPSIS

outputgenequery.pl -1 <libraryid> -2 <output> [--help] [--manual]

 
=head1 OPTIONS

=over 3

=item B<-1, -a, in>

List of libraries numbers separated by comma (e.g 123,975).  (Required)

=item B<-2, --b, -out>

Output file name  (Required)

=item B<-h, --help>

Displays the usage message.  (Optional) 

=item B<-man, --manual>

Displays full manual.  (Optional) 

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
my($specifics,$output1,$help,$manual);
GetOptions(
				"1|a|in|in1|list=s"	=>	\$specifics,
				"2|b|out1|output1=s"        =>      \$output1,
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
$dbh = DBI->connect($dsn, $user, $passwd) or die "Connection Error: $DBI::errstr\n";


#EMPTY SPACES COUNTER
my $no = 0;
#HASH TABLES
my (%CHROM, %FPKM, %POSITION, %REALPOST);
my ($realstart, $realstop);
# OPENING OUTPUT FILE
open (OUT, ">$output1");

#SPECIFYING LIBRARIES OF INTEREST
my @headers = split("\,", $specifics);

# HEADER print out
print OUT "GENE\tCHROM\t";
foreach my $name (0..$#headers-1){
	print OUT "library_$headers[$name]\t";
}
print OUT "library_$headers[$#headers]\n";
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - M A I N  W O R K F L O W - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#TABLE COLUMNS
$syntax = "select gene_short_name, fpkm, 
	library_id, chrom_no, chrom_start, chrom_stop from genes_fpkm where library_id in 
	($specifics) ORDER BY gene_id desc;";
$sth = $dbh->prepare($syntax);
$sth->execute or die "SQL Error: $DBI::errstr\n";
while (my ($gene_id, $fpkm, $library_id, $chrom, $start, $stop) = $sth->fetchrow_array() ) {
	$FPKM{$gene_id}{$library_id} = $fpkm;
	$CHROM{$gene_id} = $chrom;
	$POSITION{$gene_id}{$library_id} = "$start|$stop";
}
# DISCONNECT FROM THE DATABASE
$dbh->disconnect();

foreach my $genest (sort keys %POSITION) {
	if ($genest =~ /^[0-9a-zA-Z]/){
		my $status = "nothing";
		my @newstartarray; my @newstoparray;
		foreach my $libest (sort keys % {$POSITION{$genest}} ){
			my @newposition = split('\|',$POSITION{$genest}{$libest},2);  
			my $status = "nothing";
			
			if ($newposition[0] > $newposition[1]) {
				$status = "reverse";
			}
			elsif ($newposition[0] < $newposition[1]) {
				$status = "forward";
			}
			push @newstartarray, $newposition[0];
			push @newstoparray, $newposition[1];
			
			if ($status == "forward"){
				$realstart = (sort {$a <=> $b} @newstartarray)[0];
				$realstop = (sort {$b <=> $a} @newstoparray)[0];
				
			}
			elsif ($status == "reverse"){
				$realstart = (sort {$b <=> $a} @newstartarray)[0];
				$realstop = (sort {$a <=> $b} @newstoparray)[0];
			}
			else { die "Something is wrong\n"; }
			
		}
		$REALPOST{$genest} = "$realstart|$realstop";
		
	}
}
foreach my $genename (sort keys %FPKM){
	if ($genename =~ /^[0-9a-zA-Z]/){
		my ($newrealstart,$newrealstop) = split('\|',$REALPOST{$genename},2);
		print OUT $genename."\t".$CHROM{$genename}."\:".$newrealstart."\-".$newrealstop."\t";
		foreach my $lib (0..$#headers-1){ 
			if (exists $FPKM{$genename}{$headers[$lib]}){
				print OUT "$FPKM{$genename}{$headers[$lib]}\t";
			}
			else {
				print OUT "0\t";
			}
		}
		if (exists $FPKM{$genename}{$headers[$#headers]}){
			print OUT "$FPKM{$genename}{$headers[$#headers]}\n";
		}
		else {
			print OUT "0\n";
		} 
	}
}
# DISCONNECT FROM THE DATABASE
close(OUT);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - -T H E  E N D - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
exit;
