#!/bin/env perl
$|=1;
use strict;
my $version = '240630';

my($file, $L, $outfile, $outfilenorm) = @ARGV;

# This script computes genome fingerprints for one or more samples, all in the same VCF file.
# First parameter: the vcf/bcf with multiple samples.
# Second parameter (optional): the comma-delimited list of fingerprint lengths to compute. [120]
# Third parameter (optional): the file where raw output should be saved. [DFM_multiVCF_out]
# Fourth parameter (optional): the file where normalized output should be saved. [DFM_multiVCF_outn]

unless (length($file) && -e $file) {
	print "Usage: computeDMF-multiVCF-asTable.pl inputVCF [fingerprintLength] [outputFile]\n";
	print "       fingerprintLength defaults to 120; outputFiles default to DMF_multiVCF_out and DMF_multiVCF_outn\n";
	print "Examples: computeDMF-multiVCF-asTable.pl VCF\n";
	print "          computeDMF-multiVCF-asTable.pl BCF 0 desiredOutputFileRaw desiredOutputFileNorm\n";
	exit;
}

$L ||= 120;
my $tooCloseCutoff = 20;
$outfile ||= "DMF_multiVCF_out";
$outfilenorm ||= "DMF_multiVCF_outn";
my @keys = qw/ACAC ACAG ACAT ACCA ACCG ACCT ACGA ACGC ACGT ACTA ACTC ACTG AGAC AGAG AGAT AGCA AGCG AGCT AGGA AGGC AGGT AGTA AGTC AGTG ATAC ATAG ATAT ATCA ATCG ATCT ATGA ATGC ATGT ATTA ATTC ATTG CAAC CAAG CAAT CACA CACG CACT CAGA CAGC CAGT CATA CATC CATG CGAC CGAG CGAT CGCA CGCG CGCT CGGA CGGC CGGT CGTA CGTC CGTG CTAC CTAG CTAT CTCA CTCG CTCT CTGA CTGC CTGT CTTA CTTC CTTG GAAC GAAG GAAT GACA GACG GACT GAGA GAGC GAGT GATA GATC GATG GCAC GCAG GCAT GCCA GCCG GCCT GCGA GCGC GCGT GCTA GCTC GCTG GTAC GTAG GTAT GTCA GTCG GTCT GTGA GTGC GTGT GTTA GTTC GTTG TAAC TAAG TAAT TACA TACG TACT TAGA TAGC TAGT TATA TATC TATG TCAC TCAG TCAT TCCA TCCG TCCT TCGA TCGC TCGT TCTA TCTC TCTG TGAC TGAG TGAT TGCA TGCG TGCT TGGA TGGC TGGT TGTA TGTC TGTG/;

my $cat;
if ($file =~ /\.gz$/) {
	$cat = 'gunzip -c';
} elsif ($file =~ /\.bcf$/) {
	$cat = 'bcftools view';
} elsif ($file =~ /\.vcf$/) {
	$cat = 'cat';
} elsif ($file =~ /\.bz2$/) {
	$cat = 'bzcat';
}
	
open F, "$cat $file |";
my @samples;
while (<F>) {
	next if /^##/;
	if (/^#CHROM/) {
		chomp;
		(undef, undef, undef, undef, undef, undef, undef, undef, undef, @samples) = split /\t/;
		last;
	}
}

my($prevChrom, @prevPos, @prevKey, @close, @count, @binary, @snvPairs);
while (<F>) {
	chomp;
	my($chrom, $pos, $rsid, $ref, $alts, $qual, $filter, $infostring, $format, @obs) = split /\t/;
	$chrom =~ s/^chr//i;
	next unless $chrom =~ /^\d+$/; # retain only autosomes
	if ($chrom ne $prevChrom) {
		$prevChrom = $chrom;
		@prevPos = @prevKey = ();
	}
	
	next if $rsid =~ /CNV/;
	next if $alts eq '<NON_REF>';
	next unless $ref =~ /^[ACGT]$/io;
	$alts =~ s/<NON_REF>,//;
	$alts =~ s/,<NON_REF>//;
	next unless $alts =~ /^[ACGT]$/io;
	
	my $key = uc "$ref$alts";
	foreach my $i (0..$#obs) {
		my($gt) = split /:/, $obs[$i], 2;
		next if $gt =~ /0.0/ || $gt =~ /\./;
		if ($prevKey[$i]) {
			my $d = $pos-$prevPos[$i]-1;
			next if $d<0;
			my $pairKey = $prevKey[$i].$key;
			$count[$i]{$pairKey}[$d % $L]++ unless $d<$tooCloseCutoff;
		}
		$prevPos[$i] = $pos;
		$prevKey[$i] = $key;
		$snvPairs[$i]++;
	}
}
close F;

open  OUTF, ">$outfile";
print OUTF "#source\t$file\n";
print OUTF "#tooCloseCutoff\t$tooCloseCutoff\n";
print OUTF "#vectorLengths\t$L\n";

open  OUTFN, ">$outfilenorm";
print OUTFN "#source\t$file\n";
print OUTFN "#tooCloseCutoff\t$tooCloseCutoff\n";
print OUTFN "#vectorLengths\t$L\n";

foreach my $i (0..$#samples) {
	my $sample = $samples[$i];
	
	# Collect and save raw fingerprint.
	my @vector;
	foreach my $key (@keys) {
		push @vector, map {$count[$i]{$key}[$_] || 0} (0..$L-1);
	}
	print OUTF join("\t", $sample, @vector), "\n";
	
	# Normalize fingerprint per reduced distance.
	foreach my $col (0..$L-1) {
		my @v = ();
		push @v, $count[$i]{$_}[$col] foreach @keys;
		my($avg, $std) = avgstd(\@v);
		$std ||= 1;
		$count[$i]{$_}[$col] = ($count[$i]{$_}[$col]-$avg)/$std foreach @keys;
	}
	# Normalize fingerprint per SNV pair key.
	foreach my $sig (@keys) {
		my($avg, $std) = avgstd($count[$i]{$sig});
		$std ||= 1;
		$_ = ($_-$avg)/$std foreach @{$count[$i]{$sig}};
	}
	
	# Collect and save normalized fingerprint.
	@vector = ();
	foreach my $key (@keys) {
		push @vector, map {sprintf("%.3f", $_)} @{$count[$i]{$key}};
	}
	print OUTFN join("\t", $sample, @vector), "\n";
}
close OUTFN;
close OUTF;

###
sub avgstd {
	my($values) = @_;
	my($sum, $devsqsum);

	my $n = scalar @$values;
	return unless $n>1;
	$sum += $_ foreach @$values;
	my $avg = $sum / $n;
	$devsqsum += ($_-$avg)**2 foreach @$values;
	my $std = sqrt($devsqsum/($n-1));
	return $avg, $std;
}
