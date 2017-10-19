#!/bin/env perl
use strict;
my $version = '171019';
####
#
# This software compares two genome fingerprints.
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
# The first and second parameters are the files holding the genome fingerprints to be compared.
# The third (optional) parameter is the fingerprint size to be used. Multiple sizes can be specified, comma-delimited. Defaults to comparing all possible fingerprint sizes, depending on what sizes are available for each of the genomes.
#
# The output displays the number of SNV pairs used in each of the two genomes ('q_pairs' and 't_pairs' for 'query' and 'target'), the similarity between their binary fingerprints, and their Spearman correlation for each fingerprint size considered.
#
####
#
# Examples of usage:
#   compareDMFs.pl fingerprints/genome1.outn.gz fingerprints/genome2.outn.gz
#   compareDMFs.pl fingerprints/genome1.outn.gz fingerprints/genome2.outn.gz 5,20,120
#
####

my($f0, $f1, $L) = @ARGV;
my %todo;
$todo{$_}++ foreach (split /,/, $L);

my($dmf0, $bin0, $pairs0) = readDMF($f0);
my($dmf1, $bin1, $pairs1) = readDMF($f1);
my $folded0 = foldDMF($dmf0);
my $folded1 = foldDMF($dmf1);
my $bincorr = compareBinaries($bin0, $bin1);

my @h = qw/q_pairs t_pairs binary/;
my @f = ($pairs0 || 'NA', $pairs1 || 'NA', $bin0 && $bin1 ? sprintf("%.4f", $bincorr) : 'NA');

foreach my $vl (sort {$a<=>$b} keys %$dmf0) {
	next unless defined $dmf1->{$vl};
	next if %todo && !$todo{$vl};
	push @h, "L=$vl";
	my $spearman = correlation(ranks($folded0->{$vl}), ranks($folded1->{$vl}));
	push @f, sprintf("%.4f", $spearman);
}

print join("\t", @h), "\n";
print join("\t", @f), "\n";


###
sub readDMF {
	my($file) = @_;
	my(%v, $binary, $pairs);
	
	if ($file =~ /\.gz$/) {
		open F, "gunzip -c $file |";
	} else {
		open F, $file;
	}
	while (<F>) {
		chomp;
		my($vl, $key, @v) = split /\t/;
		if ($vl =~ /^#/) {
			if ($vl eq '#binary') {
				$binary = $key;
			} elsif ($vl eq '#SNVpairs') {
				$pairs = $key;
			}
			next;
		}
		if ($vl =~ /^[ACGT]+$/) {
			#old-style fingerprint
			unshift @v, $key;
			$key = $vl;
			$vl = scalar @v;
		}
		next if %todo && !$todo{$vl};
		$v{$vl}{$key} = \@v;
	}
	close F;
	return \%v, $binary, $pairs;
}

sub foldDMF {
	my($v) = @_;
	my %f;
	
	foreach my $vl (keys %$v) {
		push @{$f{$vl}}, @{$v->{$vl}{$_}} foreach sort keys %{$v->{$vl}};
	}
	return \%f;
}

sub correlation {
	my($set1ref, $set2ref, $avg1, $std1, $avg2, $std2) = @_;
	my($n, $i, $corr);

	$n = scalar @{$set1ref};
	return unless $n && ($n == scalar @{$set2ref});
	($avg1, $std1) = avgstd($set1ref) unless $std1;
	($avg2, $std2) = avgstd($set2ref) unless $std2;
	foreach $i (0..$n-1) {
		$corr += ($$set1ref[$i]-$avg1)*($$set2ref[$i]-$avg2);
	}
	return $corr/($n-1)/$std1/$std2;
}

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

sub ranks {
	my($v) = @_;
	
	my @ranks;
	my $i;
	foreach (sort {$v->[$a] <=> $v->[$b]} (0..$#$v)) { $ranks[$_] = $i++ }
	return \@ranks;
}

sub compareBinaries {
	my($b0, $b1) = @_;
	return unless $b0 && $b1;
	my $c;
	my @v0 = split //, $b0;
	my @v1 = split //, $b1;
	foreach my $i (0..$#v0) {
		$c++ if $v0[$i] == $v1[$i];
	}
	return ($c/scalar @v0)**2;
}

