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
#ref nucleotide expected sites
%NParent = ''; 
@nucleotides = qw(A G C T);
foreach $sequence (keys %SParent){
	if (length($sequence) > 1) {
	@allcod = undef;
	$seqcodon = length($AParent{$sequence});
	@allAA = split('',$AParent{$sequence});
	for ($leng = 0; $leng <= ($seqcodon*3);$leng=$leng+3){
		push @allcod,(substr($SParent{$sequence},$leng,3));
	}
	@allcod = @allcod[ 1 .. $#allcod-1];
	$finalN = 0;
	foreach $eachcodon (0..$#allcod){
		$NON=0;
		@allnuc = split('',$allcod[$eachcodon]); print "allnuc $allcod[$eachcodon]\n";
		foreach $step (0..2){
			foreach $nucchange (@nucleotides) {
				if ($step == 0){ 
					$newcodon = $nucchange.$allnuc[1].$allnuc[2];
					unless ($AminoAcid{$newcodon} =~ $allAA[$eachcodon]){$NON++;} 
                        }
                        elsif($step==1){
					$newcodon = $allnuc[0].$nucchange.$allnuc[2];
					unless ($AminoAcid{$newcodon} =~ $allAA[$eachcodon]){$NON++;} 
				}
				elsif($step==2){
					$newcodon = $allnuc[0].$allnuc[1].$nucchange;
					unless ($AminoAcid{$newcodon} =~ $allAA[$eachcodon]){$NON++;} 
				}
			}
		}
		$capitalN = $NON/3; print "$NON\n";
		$finalN = $finalN+$capitalN;
	}
	$finalS = (3*$seqcodon - $finalN);
	$NParent{$sequence} = "$finalN|$finalS";
	print "$sequence\t$finalN\t$finalS\n";
	}
}

#ref sample observed sites
%NSample = '';
foreach $sequence (keys %SParent){
	if (length($sequence) > 1) {
	@allcod = undef;
	$seqcodon = length($AParent{$sequence});
	@allAA = split('',$AParent{$sequence});
	for ($leng = 0; $leng <= ($seqcodon*3);$leng=$leng+3){
		push @allcod,(substr($SParent{$sequence},$leng,3));
	}
	@allcod = @allcod[ 1 .. $#allcod-1];
	$finalN = 0;
	foreach $eachcodon (0..$#allcod){
		$NON=0;
		@allnuc = split('',$allcod[$eachcodon]); print "allnuc $allcod[$eachcodon]\n";
		foreach $step (0..2){
			foreach $nucchange (@nucleotides) {
				if ($step == 0){ 
					$newcodon = $nucchange.$allnuc[1].$allnuc[2];
					unless ($AminoAcid{$newcodon} =~ $allAA[$eachcodon]){$NON++;} 
                        }
                        elsif($step==1){
					$newcodon = $allnuc[0].$nucchange.$allnuc[2];
					unless ($AminoAcid{$newcodon} =~ $allAA[$eachcodon]){$NON++;} 
				}
				elsif($step==2){
					$newcodon = $allnuc[0].$allnuc[1].$nucchange;
					unless ($AminoAcid{$newcodon} =~ $allAA[$eachcodon]){$NON++;} 
				}
			}
		}
		$capitalN = $NON/3; print "$NON\n";
		$finalN = $finalN+$capitalN;
	}
	$finalS = (3*$seqcodon - $finalN);
	$NParent{$sequence} = "$finalN|$finalS";
	print "$sequence\t$finalN\t$finalS\n";
	}
}
	
