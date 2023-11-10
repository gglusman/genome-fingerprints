#!/bin/env perl
$|=1;
use strict;
use FindBin qw($Bin);

my($set, $vls, $outdir, @todo) = @ARGV;
die "No output directory specified.\n" unless $outdir;
die "It is not recommended to use the same input directory also for the fingerprinting output.\n" if $set eq $outdir;
@todo = (1..22) unless @todo;
my @vls = (200);
@vls = split /,/, $vls if $vls;
my $tooCloseCutoff = 20;

my $tgd = "1000Genomes";

my %warned;
my($popRegions, $allRegions, $allPops) = readRegions("$Bin/../data/populationCodes.txt");
my $sampleinfo = readIndivs("$Bin/../data/20130606_1000genomes.ped.txt");
mkdir $outdir, 0755;
my $file;

foreach my $c (@todo){ #, 'X', 'Y') {
	my $chrom = "chr$c";
	next if -e "$outdir/$chrom";
	mkdir "$outdir/$chrom", 0755;
	
	if ($set eq 'TGP37') {
		#TGP37:
		$file = "ALL.$chrom.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz";
		$file = "ALL.chrX.phase3_shapeit2_mvncall_integrated.20130502.genotypes.vcf.gz" if $c eq 'X';
		$file = "ALL.chrY.phase3_integrated.20130502.genotypes.vcf.gz" if $c eq 'Y';
	} elsif ($set eq 'TGP37r') {
		#TGP37r
		$file = "ALL.$chrom.phase3_shapeit2_mvncall_integrated_v5_related_samples.20130502.genotypes.vcf.gz";
	} elsif ($set eq 'TGP37p1') {
		#TGP37p1
		$file = "ALL.$chrom.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz";
	} elsif ($set eq 'TGP38') {
		#TGP38
		$file = "ALL.$chrom.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz";
	} elsif ($set eq 'TGP38r') {
		#TGP38r
		$file = "related_samples/ALL.$chrom.shapeit2_integrated_snvindels_v2a_related_samples_27022019.GRCh38.phased.vcf.gz";
	} elsif ($set eq 'TGP38S') {
		#TGP38S
		$file = "ALL.$chrom.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz";
	} elsif ($set eq 'TGP38Sr') {
		#TGP38Sr
		$file = "related_samples/ALL.$chrom.shapeit2_integrated_v1a_related_samples.GRCh38.20181129.phased.vcf.gz";
	} elsif ($set eq 'TGP38L') {
		#TGP38L
		$file = "ALL.${chrom}_GRCh38.genotypes.20170504.vcf.gz";
	} elsif ($set eq 'TGP38H') {
		#TGP38H
		$file = "CCDG_13607_B01_GRM_WGS_2019-02-19_${chrom}.recalibrated_variants.vcf.gz";
	} elsif ($set eq 'TGP38N') {
		#TGP38N
		#$file = "TGP38N/CCDG_14151_B01_GRM_WGS_2020-08-05_${chrom}.filtered.shapeit2-duohmm-phased.vcf.gz";
		$file = "${chrom}.vcf.gz";
	} elsif ($set eq 'TGP38Nn') {
		#normalized version of TGP38N
		$file = "${chrom}.vcf.gz";
	} elsif ($set eq 'TGP38Nr') {
		#TGP38Nr
		$file = "${chrom}.vcf.gz";
	} elsif ($set eq 'TGP38C') {
		#TGP38C
		$file = "phase3.${chrom}.GRCh38.GT.crossmap.vcf.gz";
	} elsif ($file = "TGP38Cn") {
		#normalized version of TGP38C
		$file = "phase3.${chrom}.GRCh38.GT.crossmap.vcf.gz";
	}
	
	if (-e "$set/$file") {
		process("$set/$file", "$outdir/$chrom");
	} elsif (-e "$tgd/$set/$file") {
		process("$tgd/$set/$file", "$outdir/$chrom");
	} else {
		print "Cannot find $file\n";
	}
}




sub process {
	my($file, $outdir) = @_;
	
	#print $file, "\t";
	my $lines;
	my(@samples, @p, @r);
	open F, "gunzip -c $file |";
	while (<F>) {
		next if /^##/;
		if (/^#CHROM/) {
			chomp;
			(undef, undef, undef, undef, undef, undef, undef, undef, undef, @samples) = split /\t/;
			foreach my $sample (@samples) {
				unless (defined $sampleinfo->{$sample}) {
					print "Warning: no info for $sample\n" unless $warned{$sample};
					$warned{$sample}++;
					push @p, "?";
					push @r, "?";
				} else {
					push @p, $sampleinfo->{$sample}{'pop'};
					push @r, $sampleinfo->{$sample}{'region'};
				}
				
			}
			last;
		}
	}
	
	my(@prevPos, @prevKey, @close, @count, @binary, @snvPairs, @hetCount);
	
	while (<F>) {
		chomp;
		my($c, $pos, $rsid, $ref, $alt, $qual, $filter, $infostring, $format, @obs) = split /\t/;
		next if $rsid =~ /CNV/;
		next unless $alt =~ /^[ACGT]$/io && $ref =~ /^[ACGT]$/io;
		my $key = uc "$ref$alt";
		
		foreach my $i (0..$#obs) {
			my($gt) = split /:/, $obs[$i]; # needed for phase1
			next if $gt =~ /0.0/ || $gt =~ /\./;
			if ($prevKey[$i]) {
				my $d = $pos-$prevPos[$i]-1;
				next if $d<0;
				my $pairKey = $prevKey[$i].$key;
				if ($d<$tooCloseCutoff) {
					$close[$i]{$pairKey}[$d]++;
				} else {
					$binary[$i]{$pairKey}[$d % 2]++;
					foreach my $vl (@vls) {
						$count[$i]{$vl}{$pairKey}[$d % $vl]++;
					}
				}
			}
			$prevPos[$i] = $pos;
			$prevKey[$i] = $key;
			$snvPairs[$i]++;
			$hetCount[$i]++ if substr($gt,0,1) ne substr($gt,2,1);
		}
		
		$lines++;
		#last if $lines>10000;
		#print "." unless $lines % 100000;
		#print " " unless $lines % 1000000;
	}
	close F;
	
	foreach my $i (0..$#samples) {
		open OUTF, ">$outdir/$samples[$i].binary";
		print OUTF "#source\t$file\n";
		print OUTF "#tooCloseCutoff\t$tooCloseCutoff\n";
		print OUTF "#SNVpairs\t$snvPairs[$i]\n";
		print OUTF "#hetCount\t$hetCount[$i]\n";
		foreach my $key (sort keys %{$binary[$i]}) {
			print OUTF join("\t", $key, @{$binary[$i]{$key}}), "\n";
		}
		close OUTF;
		
		open OUTF, ">$outdir/$samples[$i].out";
		print OUTF "#source\t$file\n";
		print OUTF "#tooCloseCutoff\t$tooCloseCutoff\n";
		print OUTF "#SNVpairs\t$snvPairs[$i]\n";
		print OUTF "#hetCount\t$hetCount[$i]\n";
		foreach my $vl (@vls) {
			foreach my $key (sort keys %{$count[$i]{$vl}}) {
				print OUTF join("\t", $vl, $key, @{$count[$i]{$vl}{$key}}), "\n";
			}
		}
		close OUTF;
		
		open OUTF, ">$outdir/$samples[$i].close";
		print OUTF "#source\t$file\n";
		print OUTF "#tooCloseCutoff\t$tooCloseCutoff\n";
		print OUTF "#SNVpairs\t$snvPairs[$i]\n";
		print OUTF "#hetCount\t$hetCount[$i]\n";
		foreach my $key (sort keys %{$close[$i]}) {
			print OUTF join("\t", $key, @{$close[$i]{$key}}), "\n";
		}
		close OUTF;
		delete $count[$i];
	}
	#print "\n";
	`gzip -f $outdir/*`;
}


######
sub readRegions {
	my($file) = @_;
	my %regions;
	my %allRegions;
	my %allPops;
	
	open F, $file;
	while (<F>) {
		my($region, $pop) = split /\t/;
		$regions{$pop} = $region;
		$allRegions{$region}++;
		$allPops{$pop}++;
	}
	close F;
	
	return \%regions, \%allRegions, \%allPops;
}

sub readIndivs {
	my($file) = @_;
	my %info;
	my %warned;
	
	open F, $file;
	$_ = <F>;
	while (<F>) {
		my($fam, $id, $father, $mother, $sex, undef, $pop) = split /\t/;
		unless ($popRegions->{$pop}) {
			print "Warning: no region for $pop\n" unless $warned{$pop};
			$warned{$pop}++;
		}
		$info{$id} = {
			'fam' => $fam,
			'father' => $father,
			'mother' => $mother,
			'sex' => $sex,
			'pop' => $pop,
			'region' => $popRegions->{$pop} || 'undefined',
		};
	}
	close F;
	return \%info;
}


