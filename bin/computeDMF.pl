#!/bin/env perl
use strict;
my $version = '181217';
####
#
# This software computes a genome fingerprint for a single personal genome.
# The method is described in:
#    Glusman G, Mauldin DE, Hood LE, Robinson M. Ultrafast Comparison of Personal
#    Genomes via Precomputed Genome Fingerprints. Front Genet. 2017 Sep 26;8:136. doi:
#    10.3389/fgene.2017.00136. eCollection 2017. PubMed PMID: 29018478; PubMed Central
#    PMCID: PMC5623000.
# 
# Copyright 2017 by Gustavo Glusman, Institute for Systems Biology, Seattle, WA, USA.
# It is provided by the Institute for Systems Biology as open source software,
# free for non-commercial use.
# See the accompanying LICENSE file for information about the governing license.
#
####
#
# Accepted input formats include VCF, BCF, RCF (ISB's range call format), and Complete Genomics' var and masterVar formats.
# The first parameter is an 'id' for the job; output files will use this as base. You can include in it a path to where you want the output files to be located.
# The second parameter is the input file (the genome as VCF, RCF, var or masterVar).
# The third (optional) parameter is the format of the input file: 'vcf', 'bcf', 'rcf', 'var' or 'masterVar'. Defaults to 'vcf'.
# The fourth (optional) parameter is the fingerprint size. Multiple sizes can be specified, comma-delimited. Defaults to including several sizes.
# The fifth (optional) parameter is the distance between consecutive SNVs that are considered 'too close'. Default is 20.
# The sixth (optional) parameter is a bed file specifying regions of interest to be included in the analysis. For example, one could specify the definition of exome segments to compute an exome-compatible fingerprint from whole-genome data. This is available only for VCF and RCF input.
#
####
#
# Examples of usage:
#   computeDMF.pl myGenome vcfs/myGenome.vcf.gz
#   computeDMF.pl fingerprints/anotherGenome vcfs/aGenome.vcf.gz vcf 5,20,120 20 exomeRegions.bed
#
####

my($id, $file, $format, $L, $C, $bedmask) = @ARGV;
my @vls = (5, 7, 11, 13, 17, 19, 20, 40, 50, 80, 100, 120, 200);
@vls = split /,/, $L if $L;	
$C ||= 20;
$format ||= 'vcf';
my @keys = qw/ACAC ACAG ACAT ACCA ACCG ACCT ACGA ACGC ACGT ACTA ACTC ACTG AGAC AGAG AGAT AGCA AGCG AGCT AGGA AGGC AGGT AGTA AGTC AGTG ATAC ATAG ATAT ATCA ATCG ATCT ATGA ATGC ATGT ATTA ATTC ATTG CAAC CAAG CAAT CACA CACG CACT CAGA CAGC CAGT CATA CATC CATG CGAC CGAG CGAT CGCA CGCG CGCT CGGA CGGC CGGT CGTA CGTC CGTG CTAC CTAG CTAT CTCA CTCG CTCT CTGA CTGC CTGT CTTA CTTC CTTG GAAC GAAG GAAT GACA GACG GACT GAGA GAGC GAGT GATA GATC GATG GCAC GCAG GCAT GCCA GCCG GCCT GCGA GCGC GCGT GCTA GCTC GCTG GTAC GTAG GTAT GTCA GTCG GTCT GTGA GTGC GTGT GTTA GTTC GTTG TAAC TAAG TAAT TACA TACG TACT TAGA TAGC TAGT TATA TATC TATG TCAC TCAG TCAT TCCA TCCG TCCT TCGA TCGC TCGT TCTA TCTC TCTG TGAC TGAG TGAT TGCA TGCG TGCT TGGA TGGC TGGT TGTA TGTC TGTG/;

# Sanitize inputs.
die if $file =~ /[\s\;]/;
die if $bedmask =~ /[\s\;]/;

# Preparation.
my $rawFileExt = 'out';
my $closeFileExt = 'out.close';
my $normFileExt = 'outn';

my($prevChrom, $prevStart, $prevKey, %count, %close, %binary, $snvPairs);
my $cat = 'cat';
if ($file =~ /\.gz$/) {
	$cat = 'gunzip -c';
} elsif ($file =~ /\.bz2$/) {
	$cat = 'bzcat';
} elsif ($file =~ /\.bcf$/) {
	$cat = 'bcftools view';
}
my %filter = (
	'vcf' => "grep -v ^\# | grep -v \"0[\/\|]0\"",
	'bcf' => "grep -v ^\# | grep -v \"0[\/\|]0\"",
	'rcf' => 'grep -v ^\#',
	'masterVar' => 'grep -v no-call',
	'var' => 'grep -v ref | grep -v no-call');
my %fields = (
	'vcf' => '1,2,4,5', 'bcf' => '1,2,4,5', 'rcf' => '1,2,5,6',
	'masterVar' => '3,4,8,9,10', 'var' => '4,5,8,9');

# Process input file.
if ($bedmask && ($format eq 'vcf' || $format eq 'rcf')) {
	open INF, "bedtools intersect -a $file -b $bedmask | $filter{$format} | cut -f$fields{$format} |";
} else {
	open INF, "$cat $file | $filter{$format} | cut -f$fields{$format} |";
}

while (<INF>) {
	chomp;
	my($chrom, $start, $ref, $var, $othervar) = split /\t/;

	# Focus the analysis on autosomes only, excluding sex chromosomes, mitochondrial chromosome, alternative haplotypes, etc.
	next unless $chrom =~ /^(chr)?\d+$/;
	
	# Interpret othervar statement from masterVar.
	$var = $othervar if $othervar && $var eq $ref;
	
	# Pay attention only to SNVs.
	next unless $var =~ /^[ACGT]$/i && $ref =~ /^[ACGT]$/i && uc $var ne uc $ref;
	
	# Compute the key for the current SNV.
	my $key = uc "$ref$var";
	if ($chrom eq $prevChrom) {
		# Compute the distance. The -1 is to shift to base zero for the modulo function.
		my $d = $start-$prevStart-1;
		next if $d<0;
		
		# Compute the key for the pair.
		my $pairKey = $prevKey.$key;
		
		# Add to table, segregating by $C, by pair key and by reduced distance.
		if ($d<$C) {
			$close{$pairKey}[$d]++;
		} else {
			$binary{$pairKey}[$d % 2]++;
			$count{$_}{$pairKey}[$d % $_]++ foreach @vls;
		}
		$snvPairs++;
	}
	
	# Store info on current SNV for next round.
	$prevChrom = $chrom;
	$prevStart = $start;
	$prevKey = $key;
}
close INF;

my $bin = join("", map {$binary{$_}[1]>$binary{$_}[0] ? 1 : 0} @keys);

my @headers = (
	['#software-version', $version],
	['#source', $file],
	['#format', $format],
	['#SNVpairs', $snvPairs],
	['#vectorLengths', join("\t", @vls)],
	['#tooCloseCutoff', $C],
	['#created', `date`],
	);
my $header = join("\n", map {join("\t", @{$_})} @headers);

# Output main fingerprint table.
open OUTF, ">$id.$rawFileExt";
print OUTF $header;
print OUTF "#binary\t$bin\n";
foreach my $vl (@vls) {
	foreach my $key (@keys) {
		print OUTF join("\t", $vl, $key, map {$count{$vl}{$key}[$_] || 0} (0..$vl-1)), "\n";
	}
}
close OUTF;

# Output secondary (short distance) fingerprint table.
open OUTF, ">$id.$closeFileExt";
print OUTF $header;
foreach my $key (@keys) {
	print OUTF join("\t", $key, map {$close{$key}[$_] || 0} (0..$C-1)), "\n";
}
close OUTF;

# Normalize the fingerprints.
foreach my $vl (@vls) {
	# Normalize fingerprint per reduced distance.
	foreach my $col (0..$vl-1) {
		my @v = ();
		push @v, $count{$vl}{$_}[$col] foreach @keys;
		my($avg, $std) = avgstd(\@v);
		$std ||= 1;
		$count{$vl}{$_}[$col] = ($count{$vl}{$_}[$col]-$avg)/$std foreach @keys;
	}
	# Normalize fingerprint per SNV pair key.
	foreach my $sig (@keys) {
		my($avg, $std) = avgstd($count{$vl}{$sig});
		$std ||= 1;
		$_ = ($_-$avg)/$std foreach @{$count{$vl}{$sig}};
	}
}

# Output normalized fingerprint table.
open OUTF, ">$id.$normFileExt";
print OUTF $header;
print OUTF "#binary\t$bin\n";
foreach my $vl (@vls) {
	foreach my $key (@keys) {
		print OUTF join("\t", $vl, $key, map {sprintf("%.3f", $_)} @{$count{$vl}{$key}}), "\n";
	}
}
close OUTF;

# Compress output files.
`gzip -f $id.$rawFileExt; gzip -f $id.$closeFileExt; gzip -f $id.$normFileExt`;


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
