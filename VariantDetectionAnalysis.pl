#!/usr/bin/perl

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MANUAL FOR VariantdetectionAnalysis.pl

=pod

=head1 NAME

$0 -- Comprehensive pipeline : Variant detection using PICARD, GATK and produces and output folder.

=head1 SYNOPSIS

VariantdetectionAnalysis.pl -a configfile [--help] [--manual]

=head1 DESCRIPTION

Accepts all folders from frnakenstein output.
 
=head1 OPTIONS

=over 3

=item B<-a, --config>=FILE

Configuration file (a template can be found @ .. ).  (Required)

=item B<-h, --help>

Displays the usage message.  (Optional) 

=item B<-man, --manual>

Displays full manual.  (Optional) 

=back

=head1 DEPENDENCIES

Requires the following Perl libraries (all standard in most Perl installs).
   Getopt::Long
   Pod::Usage

=head1 AUTHOR

Written by Modupe Adetunji, 
Center for Bioinformatics and Computational Biology Core Facility, University of Delaware.

=head1 REPORTING BUGS

Report bugs to amodupe@udel.edu

=head1 COPYRIGHT

Copyright 2015 MOA.  
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.  
This is free software: you are free to change and redistribute it.  
There is NO WARRANTY, to the extent permitted by law.  

Please acknowledge author and affiliation in published work arising from this script's usage
=cut

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - U S E R  V A R I A B L E S- - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# CODE FOR 
use strict;
use File::Basename;
use Getopt::Long;
use Time::localtime;
use Pod::Usage;
use Time::Piece;
use File::stat;
use DateTime;

 
#ARGUMENTS
my($help,$manual,$config,%CONFIGURE);
GetOptions (
                                "a|config=s"    =>      \$config,
                                "h|help"        =>      \$help,
                                "man|manual"	=>      \$manual );

# VALIDATE ARGS
pod2usage( -verbose => 2 )  if ($manual);
pod2usage( -verbose => 1 )  if ($help);

#file path for configuration file
open(CONFIG, "<", $config) or die "Configuration File \"$config\" can't be found\nTERMINATED!\n";
while (<CONFIG>){
  chomp;
  if ($_ =~ /\=/){
    my @ALL = split"\="; $ALL[1] =~ s/\[.+\]//g; #removing comments 
    foreach my $i (0..$#ALL) {$ALL[$i] =~ s/^\s+|\s+$//g;} #removing whitespace
    $CONFIGURE{$ALL[0]} = $ALL[1];
  }
}close CONFIG;
#INDEPENDENT PROGRAMS TO RUN
my $PICARDDIR =$CONFIGURE{"PICARDDIR"}."/picard.jar";
my $GATKDIR=$CONFIGURE{"GATKDIR"}."/GenomeAnalysisTK.jar";
my $VEP=$CONFIGURE{"VEP"} . "/scripts/variant_effect_predictor/variant_effect_predictor.pl";
my $SNPdat=$CONFIGURE{"SNPDAT"};
my $SNPEFF=$CONFIGURE{"SNPEFF"};
my $REF = $CONFIGURE{"REFERENCE"};
my $ANN = $CONFIGURE{"GFF"};
my $outputfolder = $CONFIGURE{"FOLDER"};
my $specie = $CONFIGURE{"ORGANISM"};
my $variantsname = "$outputfolder/all_variants.vcf";
my ($email, $notify); #notification options
#EMAILS
my $subject = 'VAP-'.`date +%m-%d-%y_%T`; chomp $subject; my $notification = $outputfolder."/".$subject.'.log';
if (length($CONFIGURE{"SUBJECT"}) >= 1) { $subject = $CONFIGURE{"SUBJECT"}; }
if (length ($CONFIGURE{"EMAIL"}) >= 1) {
  $email = $CONFIGURE{"EMAIL"};
  $notify = "yes";
}

#OUTPUT FOLDER
`mkdir $outputfolder`;
unless ($CONFIGURE{"TYPEOFDATA"} =~ /^(DNA|RNA)$/i) {
  die "TYPEOFDATA can only be DNA/RNA\t Sorry!\n";
}
#ANNOTATION TOOL OPTIONS
my ($vepformat);
if ($CONFIGURE{"TOOL"} =~ /vep/i) {
  if ($CONFIGURE{"VEPFORMAT"}) {
    $vepformat= "--".lc($CONFIGURE{"VEPFORMAT"})."--everything on --terms ensembl";
  }
  else {
    $vepformat = "--database --terms ensembl"; 
  }
}
elsif ($CONFIGURE{"TOOL"} =~ /snpeff/i) {
  #
}

#PROCESSING
if ($CONFIGURE{"TYPEOFINPUT"} =~ /^BAM$/i) {
  if ($notify) { NOTIFICATION("Processing Variants from Bam File");  }
  VARIANTS();
}
elsif ($CONFIGURE{"TYPEOFINPUT"} =~ /^FASTQ$/i) {
  ASSEMBLY();
}
else {
  die "TYPEOFINPUT can only be = BAM/FASTQ\t Sorry!\n";
}

sub ASSEMBLY {
  #`$executemouse`;
}
sub FILTERING {
  my $input = $_[1];
  my $wkdir = $_[0];
  unless(open(FILE,$input)){
    print "File \'$input\' doesn't exist\n";
    exit;
  }
  my $out = fileparse($input, qr/(\.vcf)?$/);
  my $output = "all_filtered.vcf";
  open(OUT,">$wkdir/$output");

  my @file = <FILE>; chomp @file; close (FILE);
  foreach my $chr (@file){
    unless ($chr =~ /^\#/){
      my @chrdetails = split('\t', $chr);
      my $chrIwant = $chrdetails[7];
      my @morechrsplit = split(';', $chrIwant);
      foreach my $Imptchr (@morechrsplit){
        if ($Imptchr =~ m/^DP/) {
          my @addchrsplit = split('=', $Imptchr);
          if ($addchrsplit[1] > $CONFIGURE{"DP"}){print OUT "$chr\n";}
        }
      }
    }
    else {
      print OUT "$chr\n";
    }
  }
  close (OUT);
}
sub VARIANTS{
  my $bamfile = $CONFIGURE{"FILENAME"}; 
  my $DICT = $outputfolder."/".fileparse($REF, qr/(\..*)?$/).".dict";

  #CREATE DICTIONARY for gatk
  `java -jar $PICARDDIR CreateSequenceDictionary R=$REF O=$DICT`;
  if ($notify) { NOTIFICATION("Sequence Dictionary complete");  }
  
  #SORT BAM
  `java -jar $PICARDDIR SortSam INPUT=$bamfile OUTPUT=$outputfolder/aln_sorted.bam SO=coordinate`;
  if ($notify) { NOTIFICATION("Sort Bam complete");  }      
  #ADDREADGROUPS
  my $addreadgroup = "java -jar $PICARDDIR AddOrReplaceReadGroups INPUT=$outputfolder/aln_sorted.bam OUTPUT=$outputfolder/aln_sorted_add.bam SO=coordinate RGID=LAbel RGLB=Label RGPL=illumina RGPU=Label RGSM=Label";
  `$addreadgroup`;
  if ($notify) { NOTIFICATION("Add read groups complete");  }
  
  #MARKDUPLICATES
  my $markduplicates = "java -jar $PICARDDIR MarkDuplicates INPUT=$outputfolder/aln_sorted_add.bam OUTPUT=$outputfolder/aln_sorted_mdup.bam M=$outputfolder/aln_sorted_mdup.metrics CREATE_INDEX=true";
  `$markduplicates`;
  if ($notify) { NOTIFICATION("Mark duplicates complete");  }
  
  if ($CONFIGURE{"TYPEOFDATA"} =~ /^RNA$/i) {
    #SPLIT&TRIM
    my $splittrim = "java -jar $GATKDIR -T SplitNCigarReads -R $REF -I $outputfolder/aln_sorted_mdup.bam -o $outputfolder/aln_sorted_split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 --filter_reads_with_N_cigar";
    `$splittrim`;
    if ($notify) { NOTIFICATION("SplitNCigars complete");  }
    
    #GATK
    my $gatk = "java -jar $GATKDIR -T HaplotypeCaller -R $REF -I $outputfolder/aln_sorted_split.bam -o $variantsname";
    `$gatk`;
    if ($notify) { NOTIFICATION("Haplotype caller complete");  }
  }
  else {
    my $gatk = "java -jar $GATKDIR -T HaplotypeCaller -R $REF -I $outputfolder/aln_sorted_mdup.bam -o $variantsname";
    `$gatk`;
    if ($notify) { NOTIFICATION("Haplotype caller complete");  }
  }
  #perl to select DP > 5 & get header information
  if ($CONFIGURE{"FILTER"} =~ /YES/){
    FILTERING($outputfolder, "$outputfolder/$variantsname");
    if ($notify) { NOTIFICATION("Fitering complete");  }
    $variantsname = "$outputfolder/all_filtered.vcf";
  }
  ANNOTATION($CONFIGURE{"TOOL"});
}
sub ANNOTATION {
  if ($_[0] =~ /vep/i){
    #ANNOTATIONS : running VEP
    if (length $specie > 1 ){
      my $veptxt = "perl $VEP -i $variantsname --fork 24 --species $specie $vepformat -o $outputfolder/all_VEP.txt";
      `$veptxt`;
      if ($notify) { NOTIFICATION("Vep \"TXT\" version complete");  }
      my $vepvcf = "perl $VEP -i $variantsname --fork 24 --species $specie  --dir /home/modupeore17/.vep/ --cache --vcf --merged --everything on --terms ensembl --output_file $outputfolder/all_VEP.vcf";
      `$vepvcf`;
      if ($notify) { NOTIFICATION("Vep \"VCF\" version complete");  }
    }
  }
  elsif ($_[0] =~ /snpeff/i){
    #nothing yet
  }
}
sub NOTIFICATION {
  open (NOTE, ">>$notification");
  print NOTE "Subject: ". $subject .": Status is : $_[0]\n";
  system "sendmail $email < $notification";
  close NOTE;
  system "rm -rf $notification";
}

#cufflinks -p 24 -g ~/.big_ten/chicken/chicken.gff -o cufflinks_out $tophat/accepted_hits.bam";