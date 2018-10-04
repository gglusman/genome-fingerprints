#!/bin/env perl
use strict;
my $version = '170719';
####
#
# This software collects genome fingerprints into a database for efficient searching.
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
# The first parameter is the filename base for the database to be created or modified.
# The second parameter is the fingerprint size to be used.
# The remaining parameters are the files holding the genome fingerprints to be collected, or lists of such files.
# If the database already existed, only new fingerprints are added to it.
#
# The output consists of two files: *.fp, a condensed encoding of the fingerprints; and *.id, a companion file with metadata (seek positions in the *.fp file, and input filenames).
#
####
#
# Example of usage:
#   serializeDMFs.pl dmfDB 120 fingerprints/*.outn.gz
#   serializeDMFs.pl dmfDB 120 @listOfFiles fingerprints/*.outn.gz
#     --> dmfDB.fp
#     --> dmfDB.id
#
####

my($outbase, $L, @files) = @ARGV;
die "The second parameter must be a number.\n" unless $L+1-1;

my %done;
if (-e "$outbase.id") {
	open OUTI, "$outbase.id";
	while (<OUTI>) {
		chomp;
		if (/^#/) {
			die "Incompatible fingerprint sizes\n" if /^#L\t(\d+)/ && $1 != $L;
			next;
		}
		my(undef, $name) = split /\t/;
		$done{$name}++;
	}
	close OUTI;
	open OUTI, ">>$outbase.id";
	open OUTF, ">>$outbase.fp";
} else {
	open OUTI, ">$outbase.id";
	print OUTI join("\t", "#created", `date`);
	print OUTI join("\t", "#version", $version), "\n";
	print OUTI join("\t", "#L", $L), "\n";
	
	open OUTF, ">$outbase.fp";
	print OUTF pack('i', 144*$L);
}

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
	my($dmf) = readDMF($file, $L);
	unless (defined $dmf->{$L}) {
		print "Warning: $file lacks a fingerprint of length $L\n";
		next;
	}
	
	my $folded = foldDMF($dmf);
	my $ranked = ranks($folded->{$L});
	unless (144*$L == scalar @$ranked) {
		print "Warning: $file yields an incorrect number of values for L=$L\n";
		next;
	}
	
	print OUTI join("\t", tell OUTF, $simple->{$file} || 'NA', $file), "\n";
	print OUTF join('', map {pack('i', $_)} @$ranked);
	$done++;
}
close OUTF;
close OUTI;

$done ||= "zero";
print "Added $done fingerprints to $outbase.fp\n";

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

sub ranks {
	my($v) = @_;
	
	my @ranks;
	my $i;
	foreach (sort {$v->[$a] <=> $v->[$b]} (0..$#$v)) { $ranks[$_] = $i++ }
	return \@ranks;
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


