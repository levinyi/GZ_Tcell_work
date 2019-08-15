#!/usr/bin/perl -w
use strict;
use Math::Round;

my $ori_count_file = $ARGV[0];
my $cutoff = $ARGV[1];

my %a; my %b;
open IN, "< $ori_count_file" or die "Can't open $ori_count_file!\n";
while (<IN>)	{
	my @count = split /\s+/;
	$count[0] =~ s/^10XD1\.//;
	$count[1] =~ s/^10XD1\.//;
	next if ($count[2] < $cutoff);
	$a{$count[0]}->{$count[1]} = $count[2];
	$b{$count[1]}->{$count[0]} = $count[2];
}
close (IN);

my %a1; my %a2;
foreach my $i (keys %a)	{
	if (defined($a{$i}->{$i}))	{
		$a1{$i} = $a{$i}->{$i};
	}	else	{
		foreach my $j (keys %{$a{$i}})	{
			$a2{$i} += $a{$i}->{$j};
		}
	}
}

my @a = sort {$a1{$b}<=>$a1{$a}} keys %a1;
my @b = reverse @a;
push @a, sort {$a2{$b}<=>$a2{$a}} keys %a2;

my %b1; my %b2;
foreach my $j (keys %b)	{
	if (defined($b{$j}->{$j}))	{
		$b1{$j} = $b{$j}->{$j};
	}	else	{
		foreach my $i (keys %{$b{$j}})	{
			$b2{$j} += $b{$j}->{$i};
		}
	}
}

unshift @b, sort {$b2{$a}<=>$b2{$b}} keys %b2;

open OUT, "> count.txt.cut$cutoff" or die "Can't write to count.txt.cut$cutoff!\n";
print OUT "TRA\tTRB\tRd\tlog2Rd\n";
foreach my $i (0..$#a)	{
	foreach my $j (reverse 0..$#b)	{
		if (defined($a{$a[$i]}->{$b[$j]}))	{
			print OUT "a$a[$i]\tb$b[$j]\t$a{$a[$i]}->{$b[$j]}\t".nearest(0.0001, log($a{$a[$i]}->{$b[$j]})/log(2))."\n";
		}
	}
}
close (OUT);

open OUT, "> order_a.txt.cut$cutoff" or die "Can't write to order_a.txt.cut$cutoff!\n";
print OUT "TRA\n";
foreach my $i (0..$#a)	{
	print OUT "a$a[$i]\n";
}
close (OUT);

open OUT, "> order_b.txt.cut$cutoff" or die "Can't write to order_b.txt.cut$cutoff!\n";
print OUT "TRB\n";
foreach my $j (0..$#b)	{
	print OUT "b$b[$j]\n";
}
close (OUT);


exit;