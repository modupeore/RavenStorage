#!/usr/bin/perl

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MANUAL FOR g...

=pod

=head1 NAME

$0 -- ... : transcriptatlas

=head1 SYNOPSIS

 [--help] [--manual]

 
=head1 OPTIONS

=over 3

=item B<-h, --help>

Displays the usage message.  (Optional) 

=item B<-man, --manual>

Displays full manual.  (Optional) 

=back

=head1 DEPENDENCIES

Requires the following Perl libraries (all standard in most Perl installs).
   DBI
   DBD::mysql
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
my($help,$manual);
GetOptions(
                        "h|help"	=>      \$help,
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
my ($dbh, $sth, $syntax);

# CONNECT TO THE DATABASE
#open (OUT,">/home/modupe/genelist.txt");
$dbh = DBI->connect($dsn, $user, $passwd) or die "Connection Error: $DBI::errstr\n";

#HASH TABLES
my (%SNP, %GEN, %VAR, %IND);
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - M A I N  W O R K F L O W - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#TABLE COLUMNS
$syntax = "select species Species, format(sum(genes),0) Genes, format(sum(total_VARIANTS ),0) Variants from vw_libraryinfo group by species";

$sth = $dbh->prepare($syntax);
$sth->execute or die "SQL Error: $DBI::errstr\n";
 
while (my ($species, $genes, $variants) = $sth->fetchrow_array() ) {
    $GEN{$species} = $genes;
    $VAR{$species} = $variants;
}
print "<table class=\"summary\">
        <tr>
          <th class=\"summary\">Species</th>
          <th class=\"summary\">Genes</th>
          <th class=\"summary\">Variants</th>
        </tr>\n";
foreach my $first (sort {$a cmp $b} keys %GEN){
  print "<tr><td class=\"summary\"><b>$first</b></td>
            <td class=\"summary\">$GEN{$first}</td>
            <td class=\"summary\">$VAR{$first}</td></tr>\n";
}
#Final Row
$syntax = "select format(sum(genes),0) Genes, format(sum(total_VARIANTS ),0) Variants from vw_libraryinfo";
$sth = $dbh->prepare($syntax);
$sth->execute or die "SQL Error: $DBI::errstr\n";
while (my ($genes, $variants) = $sth->fetchrow_array() ) {
 print "<tr><th class=\"summary\"><b>Total</b></td>
            <td class=\"summary\"><b>$genes</b></td>
            <td class=\"summary\"><b>$variants</b></td></tr>\n"; 
}
print "</table>\n";

# DISCONNECT FROM THE DATABASE
$sth->finish();
$dbh->disconnect();

# DISCONNECT FROM THE DATABASE
$sth->finish();
$dbh->disconnect();

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - -T H E  E N D - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
exit;
