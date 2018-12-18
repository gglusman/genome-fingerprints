#!/usr/bin/env perl
$|=1;
use strict;

my($dir, $out) = @ARGV;

my @vls = (120);
my $tooCloseCutoff = 20;
my $computeBinary = 0;
my $computeClose = 0;
my $numberOfPartitions = 100;
$out ||= "DMF_multisample_outdir";
mkdir $out, 0700;

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
			print "Found ", scalar @samples, " samples\n";
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
			print $chrom;
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
			foreach my $key (sort keys %{$binary[$i]}) {
				print OUTF join("\t", $key, @{$binary[$i]{$key}}), "\n";
			}
			close OUTF;
		}
		
		open OUTF, ">$outdir/$partition/$sample.out";
		print OUTF "#source\t$file\n";
		print OUTF "#tooCloseCutoff\t$tooCloseCutoff\n";
		print OUTF "#SNVpairs\t$snvPairs[$i]\n";
		foreach my $vl (@vls) {
			foreach my $key (sort keys %{$count[$i]{$vl}}) {
				print OUTF join("\t", $vl, $key, @{$count[$i]{$vl}{$key}}), "\n";
			}
		}
		close OUTF;
		
		if ($computeClose) {
			open OUTF, ">$outdir/$partition/$sample.close";
			print OUTF "#source\t$file\n";
			print OUTF "#tooCloseCutoff\t$tooCloseCutoff\n";
			print OUTF "#SNVpairs\t$snvPairs[$i]\n";
			foreach my $key (sort keys %{$close[$i]}) {
				print OUTF join("\t", $key, @{$close[$i]{$key}}), "\n";
			}
			close OUTF;
		}
		
		`gzip -f $outdir/$partition/$sample.*`;
		print "o" unless $i % 1000;
		print " " unless $i % 10000;
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

