#!/bin/env perl
use strict;
my $version = '170719';
####
#
# This software collects genome fingerprints into a table for input into PCA and the like.
# The method is described in:
#    Glusman G, Mauldin DE, Hood LE, Robinson M. Ultrafast Comparison of Personal
#    Genomes via Precomputed Genome Fingerprints. Front Genet. 2017 Sep 26;8:136. doi:
#    10.3389/fgene.2017.00136. eCollection 2017. PubMed PMID: 29018478; PubMed Central
#    PMCID: PMC5623000.
# 
# Copyright 2017 by Gwenlyn Glusman, Institute for Systems Biology, Seattle, WA, USA.
# It is provided by the Institute for Systems Biology as open source software,
# free for non-commercial use.
# See the accompanying LICENSE file for information about the governing license.
#
####
#
# The first parameter is the filename base for the database to be created. The file extension .dmf is added to this filename base.
# The second parameter is the fingerprint size to be used.
# The remaining parameters are the files holding the genome fingerprints to be collected, or lists of such files.
# If the database already existed, it is overwritten!
#
####
#
# Example of usage:
#   collectDMFs.pl dmfDB 120 fingerprints/*.outn.gz
#   collectDMFs.pl dmfDB 120 @listOfFiles fingerprints/*.outn.gz
#     --> dmfDB.dmf
#
####

my($outbase, $L, @files) = @ARGV;
die "The second parameter must be a number.\n" unless $L+1-1;

my %done;
open OUTI, ">$outbase.dmf";

my $done;
my @todo;
foreach my $file (@files) {
	next if $done{$file};
	if ($file =~ /^\@(.+)/ && -s $1) {
		open LST, $1;
		while (<LST>) {
			chomp;
			($_) = split /\t/;
			push @todo, $_ if !$done{$_} && -s $_;
		}
		close LST;
	} else {
		push @todo, $file;
	}
}

my $simple = simplifyNames(@todo);

foreach my $file (@todo) {
	next if $done{$simple->{$file}};
	my($dmf, undef, $pairs) = readDMF($file, $L);
	unless (defined $dmf->{$L}) {
		print "Warning: $file lacks a fingerprint of length $L\n";
		next;
	}
	
	my $folded = foldDMF($dmf);
	$folded = $folded->{$L};
	print join("\t", $file, $pairs), "\n";
	unless (144*$L == scalar @$folded) {
		die "Warning: $file yields an incorrect number of values for L=$L\n";
	}
	
	print OUTI join("\t", $simple->{$file}, $pairs, @$folded), "\n";;
}
close OUTI;


###
sub readDMF {
	my($file, $L) = @_;
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
		$v{$vl}{$key} = \@v if $vl==$L;
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

sub simplifyNames {
	my(@names) = @_;
	my %simple;
	my @dict;
	foreach my $name (@names) {
		my @parts = split /[\/\.]/, $name;
		foreach my $i (0..$#parts) {
			$dict[$i]{$parts[$i]}++;
		}
	}
	my $start;
	my $end = $#dict;
	while (1==scalar keys %{$dict[$start]}) { $start++ }
	while (1==scalar keys %{$dict[$end]})   { $end-- }
	return {} if $start>$end;
	#print join("\t", scalar @names, $start, $end), "\n";
	#print join("\t", map {scalar keys %{$dict[$_]}} (0..$#dict)), "\n";
	foreach my $name (@names) {
		my @parts = split /[\/\.]/, $name;
		$simple{$name} = join(".", @parts[$start..$end]), "\n";
		#print join("\t", $simple{$name}, $name), "\n";
	}
	return \%simple;
}


