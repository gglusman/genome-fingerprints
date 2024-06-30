#!/bin/env perl
$|=1;
use strict;
my $version = '240630';

my($file, $L, $outfile) = @ARGV;

# This script counts het sites for one or more samples, all in the same VCF file.
# First parameter: the vcf/bcf with multiple samples.
# Second parameter (optional): the file where output should be saved. [count_sites.out]

unless (length($file) && -e $file) {
	print "Usage: countSites.pl inputVCF [outputFile]\n";
	print "Examples: countSites.pl VCF\n";
	print "          countSites.pl BCF 0 desiredOutputFile\n";
	exit;
}

my $tooCloseCutoff = 20;
$outfile ||= "count_sites.out";

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

my($prevChrom, @count);
while (<F>) {
	chomp;
	my($chrom, $pos, $rsid, $ref, $alts, $qual, $filter, $infostring, $format, @obs) = split /\t/;
	$chrom =~ s/^chr//i;
	next unless $chrom =~ /^\d+$/; # retain only autosomes
	
	next if $rsid =~ /CNV/;
	next if $alts eq '<NON_REF>';
	next unless $ref =~ /^[ACGT]$/io;
	$alts =~ s/<NON_REF>,//;
	$alts =~ s/,<NON_REF>//;
	next unless $alts =~ /^[ACGT]$/io;
	
	foreach my $i (0..$#obs) {
		my($gt) = split /:/, $obs[$i], 2;
		next if $gt =~ /\./;
		if ($gt =~ /0.0/) {
			$count[$i][0]++;
		} elsif ($gt =~ /1.1/) {
			$count[$i][2]++;
		} elsif (substr($gt,0,1) ne substr($gt,2,1)) {
			$count[$i][1]++;
		} else {
			$count[$i][3]++;
		}
	}
}
close F;

open  OUTF, ">$outfile";
print OUTF "#source\t$file\n";
print OUTF "#tooCloseCutoff\t$tooCloseCutoff\n";

foreach my $i (0..$#samples) {
	print OUTF join("\t", $samples[$i], @{$count[$i]}), "\n";
}
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
