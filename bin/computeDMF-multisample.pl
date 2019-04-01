#!/usr/bin/env perl
$|=1;
use strict;

my($dir, $vls, $out) = @ARGV;

# First parameter: the directory where the vcfs/bcfs are located. [cwd]
# Second parameter (optional): the comma-delimited list of fingerprint lengths to compute. [120]
# Third parameter (optional): the directory where output should be saved. [DFM_multisample_outdir]

unless (length($dir) && -e $dir) {
	print "Usage: computeDMF-multisample.pl inputDirectory [fingerprintLengths] [outputDirectory]\n";
	print "       fingerprintLengths defaults to 120; outputDirectory defaults to DMF_multisample_outdir\n";
	print "Examples: computeDMF-multisample.pl dirWithVCFs\n";
	print "          computeDMF-multisample.pl dirWithVCFs 20,100\n";
	print "          computeDMF-multisample.pl dirWithBCFs 0 desiredOutputDir\n";
	exit;
}

my @vls = (120);
@vls = split /,/, $vls if $vls;
my $tooCloseCutoff = 20;
my $computeBinary = 0;
my $computeClose = 0;
my $numberOfPartitions = 100;
$out ||= "DMF_multisample_outdir";
mkdir $out, 0700;
my @keys = qw/ACAC ACAG ACAT ACCA ACCG ACCT ACGA ACGC ACGT ACTA ACTC ACTG AGAC AGAG AGAT AGCA AGCG AGCT AGGA AGGC AGGT AGTA AGTC AGTG ATAC ATAG ATAT ATCA ATCG ATCT ATGA ATGC ATGT ATTA ATTC ATTG CAAC CAAG CAAT CACA CACG CACT CAGA CAGC CAGT CATA CATC CATG CGAC CGAG CGAT CGCA CGCG CGCT CGGA CGGC CGGT CGTA CGTC CGTG CTAC CTAG CTAT CTCA CTCG CTCT CTGA CTGC CTGT CTTA CTTC CTTG GAAC GAAG GAAT GACA GACG GACT GAGA GAGC GAGT GATA GATC GATG GCAC GCAG GCAT GCCA GCCG GCCT GCGA GCGC GCGT GCTA GCTC GCTG GTAC GTAG GTAT GTCA GTCG GTCT GTGA GTGC GTGT GTTA GTTC GTTG TAAC TAAG TAAT TACA TACG TACT TAGA TAGC TAGT TATA TATC TATG TCAC TCAG TCAT TCCA TCCG TCCT TCGA TCGC TCGT TCTA TCTC TCTG TGAC TGAG TGAT TGCA TGCG TGCT TGGA TGGC TGGT TGTA TGTC TGTG/;

FILE: foreach my $file (slicedirlist($dir, "[bv]cf")) {
	my $cat;
	if ($file =~ /\.gz$/) {
		$cat = 'gunzip -c';
	} elsif ($file =~ /\.bcf$/) {
		$cat = 'bcftools view';
	} elsif ($file =~ /\.vcf$/) {
		$cat = 'cat';
	} elsif ($file =~ /\.bz2$/) {
		$cat = 'bzcat';
	} else {
		print "skipping file: $file\n";
		next;
	}
	
	open F, "$cat $dir/$file |";
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
	
	my($currentChromosome, $outdir, $lines, @prevPos, @prevKey, @close, @count, @binary, @snvPairs);
	while (<F>) {
		chomp;
		my($chrom, $pos, $rsid, $ref, $alt, $qual, $filter, $infostring, $format, @obs) = split /\t/;
		$chrom = "chr$chrom" unless $chrom =~ /^chr/;
		if (-e "$out/$chrom") {
			if ($chrom ne $currentChromosome) {
				close F;
				next FILE;
			}
		} elsif ($currentChromosome) {
			die "Unexpected chromosome $chrom when working on $currentChromosome\n";
		} else {
			$outdir = "$out/$chrom";
			mkdir $outdir, 0700;
			$currentChromosome = $chrom;
			print "(", scalar @samples, " samples) $chrom";
		}
		
		next if $rsid =~ /CNV/;
		next unless $alt =~ /^[ACGT]$/io && $ref =~ /^[ACGT]$/io;
		my $key = uc "$ref$alt";
		
		foreach my $i (0..$#obs) {
			my $gt = $obs[$i];
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
		
		$lines++;
		#last if $lines>10000;
		print "." unless $lines % 10000;
		print " " unless $lines % 1000000;
	}
	close F;
	
	foreach my $i (0..$#samples) {
		my $sample = $samples[$i];
		my $partition = partition($sample);
		mkdir "$outdir/$partition", 0700;
		
		if ($computeBinary) {
			open OUTF, ">$outdir/$partition/$sample.binary";
			print OUTF "#source\t$file\n";
			print OUTF "#tooCloseCutoff\t$tooCloseCutoff\n";
			print OUTF "#SNVpairs\t$snvPairs[$i]\n";
			foreach my $key (@keys) {
				print OUTF join("\t", $key, map {$binary[$i]{$key}[$_] || 0} (0,1)), "\n";
			}
			close OUTF;
		}
		
		open OUTF, ">$outdir/$partition/$sample.out";
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
			open OUTF, ">$outdir/$partition/$sample.close";
			print OUTF "#source\t$file\n";
			print OUTF "#tooCloseCutoff\t$tooCloseCutoff\n";
			print OUTF "#SNVpairs\t$snvPairs[$i]\n";
			foreach my $key (@keys) {
				print OUTF join("\t", $key, map {$close[$i]{$key}[$_] || 0} (0..$tooCloseCutoff-1)), "\n";
			}
			close OUTF;
		}
		
		`gzip -f $outdir/$partition/$sample.*`;
		print "o" unless ($i+1) % 1000;
		print " " unless ($i+1) % 10000;
	}
	print "\n";
}

sub slicedirlist {
	my($dir, $pat) = @_;
	opendir (DIR, $dir);
	my @files = grep /$pat/, readdir DIR;
	closedir DIR;
	return @files;
}

sub partition {
	my($name) = @_;
	
	my $v = $name;
	$v =~ s/\D//g;
	unless (length($v)) {
		$v += ord() foreach split //, $name;
	}
	$v = $v % $numberOfPartitions;
	$v = "0$v" if length($v)<2;
	return $v;
}

