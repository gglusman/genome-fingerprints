#!/bin/env perl
$|=1;
use strict;
use FindBin qw($Bin);

my($data) = @ARGV;
die "Can't find $data\n" unless -e $data;

my $tgd = "/15TB_3/KaviarData/sources/1000Genomes";
my %warned;
my($popRegions, $allRegions, $allPops) = readRegions("$Bin/../data/populationCodes.txt");
my $sampleinfo = readIndivs("$Bin/../data/20130606_1000genomes.ped.txt");

my $outdir = "autosomal";
my @files = fulldirlist("$data/chr1");
mkdir "$data/$outdir", 0755;
my @sets = qw/close out binary/;


foreach my $file (@files) {
	my($ind, $set) = split /\./, $file;
	next if -e "$data/$outdir/$ind.out" || -e "$data/$outdir/$ind.out.gz";
	open TOUCH, ">$data/$outdir/$ind.out";
	close TOUCH;
	
	my(@count, @pairs, @hets);
	foreach my $chrom (fulldirlist($data)) {
		next unless $chrom =~ /^chr\d/;
		foreach my $s (0..$#sets) {
			my $set = $sets[$s];
			open F, "gunzip -c $data/$chrom/$ind.$set.gz |";
			while (<F>) {
				chomp;
				if (/^#/) {
					$pairs[$s] += $1 if /^#SNVpairs\t(\d+)/;
					$hets[$s] += $1 if /^#hetCount\t(\d+)/;
					next;
				}
				my($sig, @v);
				my $vl = 1;
				if ($s==1) {
					($vl, $sig, @v)  = split /\t/;
				} else {
					($sig, @v) = split /\t/;
				}
				foreach my $i (0..$#v) { $count[$s]{$vl}{$sig}[$i] += $v[$i] }
			}
			close F;
		}
	}

	foreach my $s (0..2) {
		my $set = $sets[$s];
		open OUTF, ">$data/$outdir/$ind.$set";
		print OUTF join("\t", "#SNVpairs", $pairs[$s]), "\n";
		print OUTF join("\t", "#hetCount", $hets[$s]), "\n" if $hets[$s];
		foreach my $vl (sort {$a<=>$b} keys %{$count[$s]}) {
			foreach my $sig (sort keys %{$count[$s]{$vl}}) {
				if ($s==1) {
					print OUTF join("\t", $vl, $sig, @{$count[$s]{$vl}{$sig}}), "\n";
				} else {
					print OUTF join("\t", $sig, @{$count[$s]{$vl}{$sig}}), "\n";
				}
			}
		}
		close OUTF;
		`gzip -f $data/$outdir/$ind.$set`;
	}
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

sub fulldirlist {
	my($dir) = @_;
	opendir (DIR, $dir);
	my @files = grep /^[^.]/, readdir DIR;
	closedir DIR;
	return @files;
}
