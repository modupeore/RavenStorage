#!/usr/bin/perl

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MANUAL FOR outputgeneslist.pl

=pod

=head1 NAME

$0 -- Extract fpkm information of specified genes based on tissue and line in the mysql database : transcriptatlas

=head1 SYNOPSIS

outputgeneslist.pl -1 <genes> -2 <tissue> -3 <species> [--help] [--manual]

 
=head1 OPTIONS

=over 3

=item B<-1, -a, in1, gene>

Gene name (e.g ASC).  (Required)

=item B<-2, -b, in2, tissue>

Tissues of interest separated by commas (e.g liver,kidney).  (Required)

=item B<-3, -c, -in3, species>

Species (e.g gallus).  (Required)

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
my($gene,$tissue,$species, $output1,$specs,$help,$manual);
GetOptions(
			"1|a|in1|gene=s"	=>	\$gene,
                        "2|b|in2|tissue=s"	=>	\$tissue,
                        "3|c|in3|species=s"	=>	\$species,
                        "h|help"		=>      \$help,
			"man|manual"		=>      \$manual );

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

#HASH TABLES
my (%CHROM, %FPKM, %GENES);

# OPENING OUTPUT FILE


# HEADER print out

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - M A I N  W O R K F L O W - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#TABLE COLUMNS
#gene_short_name, chrom_no, chrom_start, chrom_stop, fpkm, library_id
if ($tissue =~ /\,/) {
  my @tissue = split("\,",$tissue);
  foreach my $ftissue (@tissue){
    if ($species =~ /gallus/) {
      $syntax = "call usp_genedgenedtissue(\"".$gene."\",\"".$ftissue."\",\"".$species."\")";
      $sth = $dbh->prepare($syntax);
      $sth->execute or die "SQL Error: $DBI::errstr\n";
      
      while (my ($genename, $line, $max, $avg, $min) = $sth->fetchrow_array() ) { 
        $GENES{$genename}{$line} = "$max|$avg|$min";
      }
      GETINFO($ftissue);
    }
    else {
      $syntax = "call usp_genedgenedtissueless(\"".$gene."\",\"".$ftissue."\",\"".$species."\")";
      $sth = $dbh->prepare($syntax);
      $sth->execute or die "SQL Error: $DBI::errstr\n";
      
      while (my ($genename, $line, $max, $avg, $min) = $sth->fetchrow_array() ) { 
        $GENES{$genename}{$line} = "$max|$avg|$min";
      }
      GETINFO($ftissue);
    }
  }
}
else {
  if ($species =~ /gallus/) {
    $syntax = "call usp_genedgenedtissue(\"".$gene."\",\"".$tissue."\",\"".$species."\")";
    $sth = $dbh->prepare($syntax);
    $sth->execute or die "SQL Error: $DBI::errstr\n";
    
    while (my ($genename, $line, $max, $avg, $min) = $sth->fetchrow_array() ) { 
      $GENES{$genename}{$line} = "$max|$avg|$min";
    }
      GETINFO($tissue);
  }
  else {
    $syntax = "call usp_genedgenedtissueless(\"".$gene."\",\"".$tissue."\",\"".$species."\")";
    $sth = $dbh->prepare($syntax);
    $sth->execute or die "SQL Error: $DBI::errstr\n";
    
    while (my ($genename, $line, $max, $avg, $min) = $sth->fetchrow_array() ) { 
      $GENES{$genename}{$line} = "$max|$avg|$min";
    }
      GETINFO($tissue);
  }
}


# DISCONNECT FROM THE DATABASE
$sth->finish();
$dbh->disconnect();

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - -S U B R O U T I N E S- - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub GETINFO {
  my $tisue = uc $_[0];
  print "<xx>$tisue<xx>";
    print "<table class=\"gened\">
            <tr>
              <th class=\"gened\">Line</th>
              <th class=\"gened\">Maximum Fpkm</th>
              <th class=\"gened\">Average Fpkm</th>
              <th class=\"gened\">Minimum Fpkm</th>
            </tr>\n";
      foreach my $a (keys %GENES){
        print "<tr>
              <th class=\"geneds\" colspan=100%>$a</th></tr>\n";
        foreach my $b (sort keys % {$GENES{$a} }){
          my @all = split('\|', $GENES{$a}{$b}, 3);
          print "<tr><td class=\"gened\"><b>$b</b></td>
                <td class=\"gened\">$all[0]</td>
                <td class=\"gened\">$all[1]</td>
                <td class=\"gened\">$all[2]</td></tr>\n";
        }    
      }
      print "</table>\n";
      undef %GENES;
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - -T H E  E N D - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
exit;

