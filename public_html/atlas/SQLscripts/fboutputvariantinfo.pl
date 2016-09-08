#!/usr/bin/perl

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MANUAL FOR genelistquery.pl

=pod

=head1 NAME

$0 -- Extract genes information from specified library_id in the mysql database : transcriptatlas

=head1 SYNOPSIS

genelistquery.pl [--help] [--manual]

 
=head1 OPTIONS

=over 3

=item B<-1, -a, in>

List of libraries numbers separated by comma (e.g 123,975).  (Required)

=item B<-2, --b, -out>

Output1 .  (Optional)

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
use Getopt::Long;
use Pod::Usage;

#ARGUMENTS
my($gene,$begin,$species,$end,$chrom, $output,$help,$manual);
GetOptions(
			"1|g|in1|gene=s"	=>	\$gene,
                        "2|b|in2|begin=s"	=>	\$begin,
                        "3|e|in3|end=s"		=>	\$end,
                        "4|s|in4|species=s"	=>	\$species,
                        "5|c|in5|chrom=s"	=>	\$chrom,
                        "6|o|out|output=s"	=>	\$output,
                        "h|help"		=>      \$help,
			"man|manual"		=>      \$manual );

# VALIDATE ARGS
pod2usage( -verbose => 2 )  if ($manual);
pod2usage( -verbose => 1 )  if ($help);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - G L O B A L  V A R I A B L E S- - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#DIRECTORY WHERE THE FILES ARE STORED
my $chickenpath = "/home/modupe/public_html/FBTranscriptAtlas/chickenvariants";
my $mousepath = "/home/modupe/public_html/FBTranscriptAtlas/mousevariants";
my $alligatorpath = "/home/modupe/public_html/FBTranscriptAtlas/alligatorvariants";

#OUTPUT FILE (temporary and permanent)
my $tempoutput = substr($output,0,-4);
my $downloadoutput = substr($output,0,-4).".txt";

#QUERYING OPTIONS
my ($ibis, $syntax);

#HASH TABLES
my (%CHROM, %FPKM, %GENES);

# OPENING OUTPUT FILE
open(OUT, '>', $output);
open(OUTDOWN, '>', $downloadoutput);

#print headers
print OUT "<table class=\"gened\">
      <tr>
        <th class=\"gened\">Chrom</th>
        <th class=\"gened\">Position</th>
        <th class=\"gened\">Ref</th>
        <th class=\"gened\">Alt</th>
        <th class=\"gened\">Class</th>
        <th class=\"gened\">Annotation</th>
        <th class=\"gened\">Gene Name</th>
        <th class=\"gened\">dbSNP</th>
        <th class=\"gened\">library_ids</th>
      </tr>\n";
print OUTDOWN "Line\tChrom\tPosition\tRef\tAlt\tClass\tAnnotation\tGene Name\tdbSNP\tlibrary_ids\n";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - M A I N  W O R K F L O W - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#TABLE COLUMNS
if ($species =~ /gallus/) {
  $ibis = "/home/modupe/.bin/bin/ibis -d $chickenpath -q \"";
  if ($chrom) {
    $syntax = "select line,chrom,position,ref,alt,class,consequence,
                genename,dbsnp,group_concat(library) where chrom = \'$chrom\' and position between $begin and $end\" -v -o ";
  }
  elsif($gene) {
    $syntax = "select line, chrom,position, ref, alt, class, consequence,
                genename, dbsnp,group_concat(library) where genename = \'$gene\' order by library\" -v -o ";
  }
  print "$ibis$syntax$tempoutput\n\n";
  `$ibis$syntax$tempoutput`;
  
  open(IN,'<',$tempoutput);
  while (<IN>){
    chomp;
    my ($line, $chrom, $position, $ref, $alt, $class, $ann, $genename, $dbsnp, $library) = split(/\, /, $_, 10);
    #removing the quotation marks from the words
    $line = uc(substr($line,1,-1));$chrom = substr($chrom,1,-1);
    $ref = substr($ref,1,-1); $alt = substr($alt,1,-1);
    $class = substr($class,1,-1); $ann = substr($ann,1,-1);
    $library = substr($library,1,-1); my $newlib = undef;
    foreach (sort {$a <=> $b} split(",",$library)) { $newlib .= $_.", "};
    $newlib = substr($newlib,0,-2);
    $genename = substr($genename,1,-1); $dbsnp = substr($dbsnp,1,-1);
    if ($genename =~ /NULL/) { $genename = "-"; }
    if ($dbsnp =~ /NULL/) { $dbsnp = " "; }
    
    #storing into a hash table.
    $GENES{$line}{$chrom}{$position} = "$ref|$alt|$class|$ann|$genename|$dbsnp|$newlib";
  }
  close (IN); `rm -rf $tempoutput`;
  
  foreach my $a (sort keys %GENES){
    print OUT "<tr><th class=\"geneds\" colspan=100%>$a</th></tr>\n";
    foreach my $b (sort keys % {$GENES{$a} }){
      foreach my $c (sort keys % {$GENES{$a}{$b} }){
        my @all = split('\|', $GENES{$a}{$b}{$c}, 7);
        print OUTDOWN "$a\t$b\t$c\t";
        print OUT "<tr><td class=\"gened\"><b>$b</b></td><td class=\"gened\"><b>$c</b></td>";
        
        foreach my $ii (0..$#all-1){
          print OUTDOWN "$all[$ii]\t";
          print OUT "<td class=\"gened\">$all[$ii]</td>";
        }
        print OUTDOWN "$all[$#all]\n";
        my $finalcol = $all[$#all]; $finalcol =~ s/\s+//g;
        print OUT "<td class=\"gened\"><a href=\"metadata.php?libs=$finalcol\" target=\"_blank\">$all[$#all]</a></td></tr>\n";
      }
    }
  }
  print OUT "</table>\n";
}
else {
  if ($species =~ /mus_musculus/) {$ibis = "ibis -d $mousepath -q \"";}
  elsif ($species =~ /alligator/) {$ibis = "ibis -d $alligatorpath -q \"";}
  if ($chrom) {
    $syntax = "select chrom,position,ref,alt,class,consequence,
                genename,dbsnp,group_concat(library) where chrom = \'$chrom\' and position between $begin and $end\" -v -o ";
  }
  elsif($gene) {
    $syntax = "select chrom,position,ref,alt,class,consequence,
                genename, dbsnp,group_concat(library) where genename = \'$gene\'\" -v -o ";
  }
  print "$ibis$syntax$tempoutput\n\n";
  `$ibis$syntax$tempoutput`;
  
  open(IN,'<',$tempoutput);
  while (<IN>){
    chomp;
    my ($chrom, $position, $ref, $alt, $class, $ann, $genename, $dbsnp) = split(/\, /, $_, 8);
    #removing the quotation marks from the words
    $chrom = substr($chrom,1,-1);$ref = substr($ref,1,-1);
    $alt = substr($alt,1,-1);$class = substr($class,1,-1);
    $ann = substr($ann,1,-1);$genename = substr($genename,1,-1);
    $dbsnp = substr($dbsnp,1,-1);
    if ($genename =~ /NULL/) { $genename = "-"; }
    if ($dbsnp =~ /NULL/) { $dbsnp = " "; }
    
    #storing into a hash table.
    $GENES{$chrom}{$position} = "$ref|$alt|$class|$ann|$genename|$dbsnp";
  }
  close (IN); `rm -rf $tempoutput`;
 
  foreach my $a (keys %GENES){
    foreach my $b (sort keys % {$GENES{$a} }){
      my @all = split('\|', $GENES{$a}{$b}, 6);
      print OUT "<tr><td class=\"gened\"><b>$a</b></td><td class=\"gened\"><b>$b</b></td>";
      print "<tr><td class=\"gened\"><b>$a</b></td><td class=\"gened\"><b>$b</b></td>";
      print OUTDOWN "$a\t$b\t";
      foreach my $ii (0..$#all-1){
        print OUTDOWN "$all[$ii]\t";
        print OUT "<td class=\"gened\">$all[$ii]</td>";
      }
      print OUTDOWN "$all[$#all]\n";
      print OUT "<td class=\"gened\">$all[$#all]</td></tr>\n";
    }
  }
  print OUT "</table>\n";
}

close (OUT); close (OUTDOWN);
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - -T H E  E N D - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
exit;

