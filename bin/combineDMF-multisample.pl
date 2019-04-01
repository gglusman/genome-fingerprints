#!/usr/bin/env perl
$|=1;
use strict;

my($data, $outdir) = @ARGV;

$outdir ||= "autosomal.norm";
$outdir = "$data.norm" if $data eq $outdir;
mkdir $outdir, 0700;
my @sets = qw/close out binary/;
@sets = qw/out/;

my(%found, %chroms, %part);
foreach my $chrom (fulldirlist($data)) {
	next if $chrom !~ /chr\d/;
	foreach my $part (fulldirlist("$data/$chrom")) {
		foreach my $file (fulldirlist("$data/$chrom/$part")) {
			if ($file =~ /(.+)\.(out|close|binary)\.gz$/) {
				$found{$2}{$1}{$chrom} = "$data/$chrom/$part/$file";
				$chroms{$chrom}++ if $2 eq 'out';
				$part{$1} = $part;
			}
		}
	}
}
my %foundPerChrom;
foreach my $chrom (keys %chroms) {
	$foundPerChrom{$chroms{$chrom}}++;
}
my @sorted = sort {$foundPerChrom{$b}<=>$foundPerChrom{$a}} keys %foundPerChrom;
foreach my $chrom (keys %chroms) {
	if ($chroms{$chrom} != $sorted[0]) {
		print "Unexpected sample count for $chrom: $chroms{$chrom}, will be skipped\n";
		delete $chroms{$chrom};
	} elsif ($chrom !~ /chr\d/) {
		print "$chrom doesn't seem to be an autosome, will be skipped\n";
		delete $chroms{$chrom};
	}
}
my @chroms = sort keys %chroms;
my $nchroms = scalar @chroms;
my $chromString = join(",", @chroms);

foreach my $set (@sets) {
	mkdir "$outdir/$set", 0700;
	foreach my $id (keys %{$found{$set}}) {
		my $part = $part{$id};
		mkdir "$outdir/$set/$part", 0700;
		next if -e "$outdir/$set/$part/$id.$set.gz" || -e "$outdir/$set/$part/$id.$set";
                open OUTF, ">$outdir/$set/$part/$id.$set";
		close OUTF;
		if (join(",", sort keys %{$found{$set}{$id}}) != $chromString) {
			print "Unexpected set of chromosomes for $id: ", join(",", sort keys %{$found{$set}{$id}}), ", will be skipped\n";
			next;
		}
		my(%count, $pairs);
		foreach my $chrom (@chroms) {
			open F, "gunzip -c $found{$set}{$id}{$chrom} |";
			while (<F>) {
				chomp;
				if (/^#/) {
					$pairs += $1 if /^#SNVpairs\t(\d+)/;
					next;
				}
				my($sig, @v);
				my $vl = 1;
				if ($set eq 'out') {
					($vl, $sig, @v)  = split /\t/;
				} else {
					($sig, @v) = split /\t/;
				}
				foreach my $i (0..$#v) { $count{$vl}{$sig}[$i] += $v[$i] }
			}
			close F;
		}
		my @vls = sort {$a<=>$b} keys %count;
		if ($set eq 'out') {
			my @sigs = sort keys %{$count{$vls[0]}};
			foreach my $vl (@vls) {
				## normalize per position
				foreach my $col (0..$vl-1) {
					my @v = ();
					foreach my $sig (@sigs) { push @v, $count{$vl}{$sig}[$col] }
					my($avg, $std) = avgstd(\@v);
					foreach my $sig (@sigs) { $count{$vl}{$sig}[$col] = ($count{$vl}{$sig}[$col]-$avg)/$std }
				}
				## normalize per signature
				foreach my $sig (@sigs) {
					my($avg, $std) = avgstd($count{$vl}{$sig});
					foreach (@{$count{$vl}{$sig}}) { $_ = ($_-$avg)/$std }
				}
			}
		}
		
		mkdir "$outdir/$set/$part", 0700;
		open OUTF, ">$outdir/$set/$part/$id.$set";
		print OUTF join("\t", "#SNVpairs", $pairs), "\n";
		print OUTF join("\t", "#vectorLengths", @vls), "\n";
		foreach my $vl (sort {$a<=>$b} keys %count) {
			foreach my $sig (sort keys %{$count{$vl}}) {
				if ($set eq 'out') {
					print OUTF join("\t", $vl, $sig, map {sprintf("%.3f", $_)} @{$count{$vl}{$sig}}), "\n";
				} else {
					print OUTF join("\t", $sig, @{$count{$vl}{$sig}}), "\n";
				}
			}
		}
		close OUTF;
		`gzip -f $outdir/$set/$part/$id.$set`;
		print join("\t", $id, $pairs, "$outdir/$set/$part/$id.$set.gz"), "\n";
	}
}

######

sub fulldirlist {
	my($dir) = @_;
	opendir (DIR, $dir);
	my @files = grep /^[^.]/, readdir DIR;
	closedir DIR;
	return @files;
}

sub avgstd {
	my($values) = @_;
	my($sum, $devsqsum);

	my $n = scalar @$values;
	return unless $n>1;
	foreach (@$values) { $sum += $_ }
	my $avg = $sum / $n;
	foreach (@$values) { $devsqsum += ($_-$avg)**2 }
	my $std = sqrt($devsqsum/($n-1));
	return $avg, $std;
}
