#!/usr/bin/perl

open(INPUT,"<", $ARGV[0]) or die "Didnt open";
#open(OUT, ">", $ARGV[1]);
open(GEN,"<GeneticCode.txt") or die "Didn't work";

%AminoAcid = '';
while(<GEN>){
	chomp;
	@aaa = split("\t",$_,3);
	$AminoAcid{$aaa[0]} = $aaa[2];
}
close (GEN);
$status = 1; %SParent = ''; %SSample = '';
while(<INPUT>){
	chomp;
	if($status == 1){
		$header = $_;
		$status =2;
	}
	elsif($status == 2){
		$parent = $_;
		$status = 3;
	}
	elsif($status == 3) {
		$sample = $_;
		$status = 4;
	}
	if($status == 4){
		$SParent{$header} = $parent;
		print "$header\t$parent\t$sample\n";
		$SSample{$header} = $sample;
		$status =1;
	}
}
close(INPUT);
#convert to aminoacid
%AParent = ''; %ASample = ''; $aaseq=undef;
foreach $codons (keys %SParent){
	$seqlength = length($SParent{$codons});
	for ($leng = 0; $leng <= $seqlength;$leng=$leng+3){
		$aaseq .= $AminoAcid{substr($SParent{$codons},$leng,3)};
	}
	$AParent{$codons} = $aaseq;$aaseq=undef;
} $aaseq=undef;
foreach $codons (keys %SSample){
        $seqlength = length($SSample{$codons});
        for ($leng = 0; $leng <= $seqlength;$leng=$leng+3){
                $aaseq .= $AminoAcid{substr($SSample{$codons},$leng,3)};
        }
        $ASample{$codons} = $aaseq;$aaseq=undef;
}


foreach $keys (keys %AParent){
	if (length($keys) > 1) {
        print "$keys\t$AParent{$keys}\t$ASample{$keys}\n";
	}
}
#nucleotide sites

#foreach $names (keys %SParent){
#	@nuc = split('',$SParent{$names});
	
