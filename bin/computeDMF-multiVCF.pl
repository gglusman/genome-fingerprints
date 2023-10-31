#!/bin/env perl
$|=1;
use strict;

my($file, $vls, $outdir) = @ARGV;

# This script computes genome fingerprints for one or more samples, all in the same VCF file.
# First parameter: the vcf/bcf with multiple samples.
# Second parameter (optional): the comma-delimited list of fingerprint lengths to compute. [120]
# Third parameter (optional): the directory where output should be saved. [DFM_multiVCF_outdir]

unless (length($file) && -e $file) {
	print "Usage: computeDMF-multiVCF.pl inputVCF [fingerprintLengths] [outputDirectory]\n";
	print "       fingerprintLengths defaults to 120; outputDirectory defaults to DMF_multiVCF_outdir\n";
	print "Examples: computeDMF-multiVCF.pl VCF\n";
	print "          computeDMF-multiVCF.pl VCF 20,100\n";
	print "          computeDMF-multiVCF.pl BCF 0 desiredOutputDir\n";
	exit;
}

my @vls = (120);
@vls = split /,/, $vls if $vls;
my $tooCloseCutoff = 20;
my $computeBinary = 0;
my $computeClose = 0;
$outdir ||= "DMF_multiVCF_outdir";
mkdir $outdir, 0700;
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
		#print "Found ", scalar @samples, " samples\n";
		last;
	}
}

my($prevChrom, @prevPos, @prevKey, @close, @count, @binary, @snvPairs);
while (<F>) {
	exit if -e 'stop';
	chomp;
	my($chrom, $pos, $rsid, $ref, $alts, $qual, $filter, $infostring, $format, @obs) = split /\t/;
	next if $chrom =~ /\D/;
	#$chrom = "chr$chrom" unless $chrom =~ /^chr/;
	if ($chrom ne $prevChrom) {
		$prevChrom = $chrom;
		@prevPos = @prevKey = ();
		my $now = `date`;
		chomp($now);
		print "$now\tstarting\t$chrom\n";
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
			if ($d<$tooCloseCutoff) {
				$close[$i]{$pairKey}[$d]++ if $computeClose;
			} else {
				$binary[$i]{$pairKey}[$d % 2]++ if $computeBinary;
				foreach my $vl (@vls) {
					$count[$i]{$vl}{$pairKey}[$d % $vl]++;
				}
			}
		}
		$prevPos[$i] = $pos;
		$prevKey[$i] = $key;
		$snvPairs[$i]++;
	}
}
close F;

foreach my $i (0..$#samples) {
	my $sample = $samples[$i];
	$sample =~ s/[^a-z0-9\-\.]+//i; #sanitize sample names before they become part of filepaths!!
	$sample ||= "_sample_$i";
	
	if ($computeBinary) {
		open OUTF, ">$outdir/$sample.binary";
		print OUTF "#source\t$file\n";
		print OUTF "#tooCloseCutoff\t$tooCloseCutoff\n";
		print OUTF "#SNVpairs\t$snvPairs[$i]\n";
		foreach my $key (@keys) {
			print OUTF join("\t", $key, map {$binary[$i]{$key}[$_] || 0} (0,1)), "\n";
		}
		close OUTF;
	}
	
	open OUTF, ">$outdir/$sample.out";
	print OUTF "#source\t$file\n";
	print OUTF "#tooCloseCutoff\t$tooCloseCutoff\n";
	print OUTF "#SNVpairs\t$snvPairs[$i]\n";
	foreach my $vl (@vls) {
		foreach my $key (@keys) {
			print OUTF join("\t", $vl, $key, map {$count[$i]{$vl}{$key}[$_] || 0} (0..$vl-1)), "\n";
		}
	}
	close OUTF;
	
	if ($computeClose) {
		open OUTF, ">$outdir/$sample.close";
		print OUTF "#source\t$file\n";
		print OUTF "#tooCloseCutoff\t$tooCloseCutoff\n";
		print OUTF "#SNVpairs\t$snvPairs[$i]\n";
		foreach my $key (@keys) {
			print OUTF join("\t", $key, map {$close[$i]{$key}[$_] || 0} (0..$tooCloseCutoff-1)), "\n";
		}
		close OUTF;
	}
	
	# Normalize the fingerprints.
	foreach my $vl (@vls) {
		# Normalize fingerprint per reduced distance.
		foreach my $col (0..$vl-1) {
			my @v = ();
			push @v, $count[$i]{$vl}{$_}[$col] foreach @keys;
			my($avg, $std) = avgstd(\@v);
			$std ||= 1;
			$count[$i]{$vl}{$_}[$col] = ($count[$i]{$vl}{$_}[$col]-$avg)/$std foreach @keys;
		}
		# Normalize fingerprint per SNV pair key.
		foreach my $sig (@keys) {
			my($avg, $std) = avgstd($count[$i]{$vl}{$sig});
			$std ||= 1;
			$_ = ($_-$avg)/$std foreach @{$count[$i]{$vl}{$sig}};
		}
	}
	
	# Output normalized fingerprint table.
	open OUTF, ">$outdir/$sample.outn";
	print OUTF "#source\t$file\n";
	print OUTF "#tooCloseCutoff\t$tooCloseCutoff\n";
	print OUTF "#SNVpairs\t$snvPairs[$i]\n";
	foreach my $vl (@vls) {
		foreach my $key (@keys) {
			print OUTF join("\t", $vl, $key, map {sprintf("%.3f", $_)} @{$count[$i]{$vl}{$key}}), "\n";
		}
	}
	close OUTF;
	`gzip -f $outdir/$sample.*`;
}

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
