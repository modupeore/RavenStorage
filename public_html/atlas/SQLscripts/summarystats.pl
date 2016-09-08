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
my (%PRO, %KEY);
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - M A I N  W O R K F L O W - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#TABLE COLUMNS
$syntax = "select a.species Species, format(count(a.species),0) Recorded, format(count(b.species),0) Processed
            from bird_libraries a left outer join vw_libraryinfo b on a.library_id = b.library_id
            group by a.species";

$sth = $dbh->prepare($syntax);
$sth->execute or die "SQL Error: $DBI::errstr\n";
 
while (my ($species, $recorded, $processed) = $sth->fetchrow_array() ) {
    $KEY{$species} = $recorded;
    $PRO{$species} = $processed;
}
print "<table class=\"summary\">
        <tr>
          <th class=\"summary\">Species</th>
          <th class=\"summary\">Recorded</th>
          <th class=\"summary\">Processed</th>
        </tr>\n";
foreach my $first (sort {$a cmp $b} keys %KEY){
  print "<tr><td class=\"summary\"><b>$first</b></td>
            <td class=\"summary\">$KEY{$first}</td>
            <td class=\"summary\">$PRO{$first}</td></tr>\n";
}
#Final Row
$syntax = "select format(count(a.species),0) Recorded, format(count(b.species),0) Processed from bird_libraries a left outer join vw_libraryinfo b
            on a.library_id = b.library_id";
$sth = $dbh->prepare($syntax);
$sth->execute or die "SQL Error: $DBI::errstr\n";
while (my ($recorded, $processed) = $sth->fetchrow_array() ) {
 print "<tr><th class=\"summary\"><b>Total</b></td>
            <td class=\"summary\"><b>$recorded</b></td>
            <td class=\"summary\"><b>$processed</b></td></tr>\n"; 
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
