#!/usr/bin/perl -w
use strict;

foreach my $file (glob "*.contigs")	{
	next unless (-s "$file.m8.out");

	open IN, "< $file.m8.out" or die "Can't open $file.m8.out!\n";
	my $hit = <IN>;
	my @hit = split /\s+/, $hit;
	close (IN);

	my $seq;
	open IN, "< $file" or die "Can't open $file!\n";
	while (<IN>)	{
		next if (/>/);
		s/\s//g;
		$seq .= $_;
	}
	close (IN);

	my @qual;
	open IN, "< $file.qual" or die "Can't open $file.qual!\n";
	while (<IN>)	{
		next if (/>/);
		my @q = split /\s+/;
		push @qual, @q;
	}
	close (IN);
	
	$seq = substr($seq, $hit[6] - 1, $hit[7] - $hit[6] + 1);
	@qual = @qual[($hit[6]-1)..($hit[7]-1)];

	my $D40 = 0;
	foreach my $q (@qual)	{
		if ($q < 20)	{
			$D40 += 10*(40-$q);
		}	elsif ($q < 30)	{
			$D40 += 2*(40-$q);
		}	elsif ($q < 40)	{
			$D40 += 40-$q;
		}
	}

	print "$file\t$D40\n";

	if ($file =~ /\-(\d+)/)	{
		my $id = $1;
		if ($id <= 40)	{
			open OUT, "> $file.qual.plot" or die "Can't write to $file.qual.plot!\n";
			print OUT "rank\tpos\tqual\n";
			my %qual;
			map { $qual{$_} = $qual[$_] } (0..$#qual);
			my @i = sort {$qual{$a}<=>$qual{$b}} keys %qual;
			foreach my $i (0..$#i)	{
				print OUT ($i + 1)."\t$i[$i]\t$qual[$i[$i]]\n";
			}
			close (OUT);

			open OUT, "> tmp.r" or die "Can't write to tmp.r!\n";
			print OUT "library(ggplot2)\n";
			print OUT "d <- read.table(\"$file.qual.plot\", header=TRUE)\n";
			print OUT "p <- ggplot(d, aes(rank, qual)) + geom_point()\n";
			print OUT "png(\"$file.qual.png\")\n";
			print OUT "print(p)\n";
			print OUT "dev.off()\n";
			close (OUT);

			system "Rscript tmp.r";
			system "rm tmp.r";
		}
	}
}

exit;
