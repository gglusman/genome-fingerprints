#!/bin/env perl
$|=1;
use strict;

my $indir = shift @ARGV;
die "Usage: normalizeDMFsInDir.pl <input directory>\n" unless -e $indir;
my $outdir = "$indir.norm";
mkdir $outdir, 0755;

foreach my $file (fulldirlist($indir)) {
	## identify name from filename
	my($name, $part) = split /\./, $file;
	next unless $part eq 'out' && $name;
	#print $name;
	
	next if -e "$outdir/$name.out.gz";
	my %count;
	
	open OUTF, ">$outdir/$name.out";
	open F, "gunzip -c $indir/$file |";
	while (<F>) {
		if (/^#/) {
			print OUTF;
			next;
		}
		chomp;
		my($vl, $sig, @v) = split /\t/;
		$count{$vl}{$sig} = \@v;
	}
	close F;
	
	my $binary = computeBinary("$indir/$name.binary.gz");
	print OUTF "#binary\t$binary\n";
	
	my @vls = sort {$a<=>$b} keys %count;
	my @sigs = sort keys %{$count{$vls[0]}};
	my $nsig = scalar @sigs;
	
	print OUTF join("\t", "#vectorLengths", @vls), "\n";
	
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
		
		foreach my $key (@sigs) {
			print OUTF join("\t", $vl, $key, map {sprintf("%.3f", $_)} @{$count{$vl}{$key}}), "\n";
		}
	}
	close OUTF;
	`gzip -f $outdir/$name.out`;
}

###
sub computeBinary {
	my($file) = @_;
	
	my $binary;
	
	open F, "gunzip -c $file |";
	while (<F>) {
		next if /^#/;
		chomp;
		my($sig, @v) = split /\t/;
		$binary .= ($v[1]>$v[0] || 0);
	}
	close F;
	return $binary;
}



###
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


sub fulldirlist {
	my($dir) = @_;
	opendir (DIR, $dir);
	my @files = grep /^[^.]/, readdir DIR;
	closedir DIR;
	return @files;
}

