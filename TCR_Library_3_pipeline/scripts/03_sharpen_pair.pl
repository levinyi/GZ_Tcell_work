#!/usr/bin/perl -w
use strict;
use Math::Round;

my $prefix = $ARGV[0];
my $min_pair_abun = $ARGV[1];
my $min_gene_spec = $ARGV[2];

my $pwd_path = `pwd`;
chomp($pwd_path);
my $result_path = "$pwd_path/../result/";
my %pair;
open IN, "< $result_path/$prefix.pairs" or die "Can't open $result_path/$prefix.pairs!\n";
while (<IN>)	{
	my @pair = split /\s+/;
	$pair{$pair[1]}->{$pair[2]}++;
}
close (IN);

my %freq_a; my %freq_b; my %total_a; my %total_b;
foreach my $a (keys %pair)	{
	foreach my $b (keys %{$pair{$a}})	{
		if ($pair{$a}->{$b} >= $min_pair_abun)	{
			$freq_a{$a}->{$b} = $pair{$a}->{$b};
			$freq_b{$b}->{$a} = $pair{$a}->{$b};
			$total_a{$a} += $pair{$a}->{$b};
			$total_b{$b} += $pair{$a}->{$b};
		}
	}
}

my %uniq_a;
&filter_gene_spec(\%freq_a, \%total_a, $min_gene_spec, \%uniq_a);
my %uniq_b;
&filter_gene_spec(\%freq_b, \%total_b, $min_gene_spec, \%uniq_b);

my %freq_uniq;
my $total = &freq_uniq(\%freq_a, \%uniq_a, \%uniq_b, \%freq_uniq);

my %m1;
&sharpen(\%freq_uniq, \%m1, "M1", $prefix);
&report(\%m1, "M1", $prefix, \%freq_uniq, $total);

my %m2;
&sharpen(\%m1, \%m2, "M2", $prefix);
&report(\%m2, "M2", $prefix, \%freq_uniq, $total);

exit;

sub filter_gene_spec	{
	(my $freq, my $total, my $min_gene_spec, my $uniq) = @_;

	foreach my $tcr (keys %{$freq})	{
		my @tcr = sort {$freq->{$tcr}->{$b}<=>$freq->{$tcr}->{$a}} keys %{$freq->{$tcr}};
		my $count = 0;
		foreach my $i (0..2)	{
			if ($i <= $#tcr)	{
				$count += $freq->{$tcr}->{$tcr[$i]};
			}
		}
		if ($count/$total->{$tcr} >= $min_gene_spec)	{
			$uniq->{$tcr} = 1;
		}
	}
}

sub freq_uniq	{
	(my $freq_a, my $uniq_a, my $uniq_b, my $freq_uniq) = @_;

	my $total = 0;
	foreach my $a (keys %{$freq_a})	{
		next unless (defined($uniq_a->{$a}));
		foreach my $b (keys %{$freq_a->{$a}})	{
			next unless (defined($uniq_b->{$b}));
			$freq_uniq->{$a}->{$b} = $freq_a->{$a}->{$b};
			$total += $freq_a->{$a}->{$b};
		}
	}

	return $total;
}

sub sharpen	{
	(my $i, my $o, my $info, my $prefix) = @_;
	
	my %total_a; my %total_b;
	foreach my $a (keys %{$i})	{
		foreach my $b (keys %{$i->{$a}})	{
			$total_a{$a} += $i->{$a}->{$b};
			$total_b{$b} += $i->{$a}->{$b};
		}
	}

	open OUT, "> $result_path/$prefix.pairs.$info.score" or die "Can't write to $result_path/$prefix.pairs.$info.score!\n";
	foreach my $a (keys %{$i})	{
		foreach my $b (keys %{$i->{$a}})	{
			$o->{$a}->{$b} = $i->{$a}->{$b}**2/($total_a{$a}*$total_b{$b});
			print OUT "$a\t$b\t$o->{$a}->{$b}\n";
		}
	}
	close (OUT);

	return;
}

sub report	{
	(my $m, my $info, my $prefix, my $freq, my $total) = @_;

	my %v;
	foreach my $a (keys %{$m})	{
		foreach my $b (keys %{$m->{$a}})	{
			my $v = int(100*$m->{$a}->{$b});
			my $pair = join "\t", ($a, $b);
			$v{$v}->{$pair} = 1;
		}
	}

	open OUT, ">> $result_path/$prefix.pairs.sharpen.data" or die "Can't write to $result_path/$prefix.pairs.sharpen.data!\n";
	my $count = 0;
	for (my $t = 100; $t >= 0; $t--)	{
		if (defined($v{$t}))	{
			my @pair = keys %{$v{$t}};
			$count += scalar @pair;
		}
		print OUT "$prefix\t$info\ttotal\t".($t/100)."\t$count\n";
	}

	$count = 0;
	for (my $t = 100; $t >= 0; $t--)	{
		if (defined($v{$t}))	{
			foreach my $pair (keys %{$v{$t}})	{
				my @pair = split /\t/, $pair;
				$count += $freq->{$pair[0]}->{$pair[1]};
			}
		}
		print OUT "$prefix\t$info\tfrac\t".($t/100)."\t".nearest(0.0001, $count/$total)."\n";
	}
	close (OUT);

	return;
}
