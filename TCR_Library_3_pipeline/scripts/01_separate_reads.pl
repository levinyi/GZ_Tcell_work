#!/usr/bin/perl -w
use strict;
use File::Basename;

my $list = $ARGV[0];	# a list of samples
my $suffix = $ARGV[1];	# R1.p1 or R2
my $type = $ARGV[2];	# TRA or TRB
my $prefix = $ARGV[3];	# e.g. 16

my %ori;
my $list_path = dirname($list);

open LIST, "< $list" or die "Can't open $list!\n";
while (my $id = <LIST>)	{
	$id =~ s/\s//g;
	open FQ, "gunzip -c $list_path/../data/$id\_$suffix.fq.gz |" or die "Can't open $list_path/../data/$id\_$suffix.fq.gz!\n";
	while (<FQ>)	{
		if (/^\@(\S+)/)	{
			$ori{$1} = $id;
			<FQ>;
			<FQ>;
			<FQ>;
		}
	}
	close (FQ);
}
close (LIST);

print "finished reading $suffix reads!\n";

my %clone;
#open CLONE, "16_$suffix.mixcr.out.clonotypes.$type.txt" or die "Can't open 16_$suffix.mixcr.out.clonotypes.$type.txt!\n";
open CLONE, "$list_path/$prefix\_$suffix.mixcr.out.clonotypes.$type.txt" or die "Can't open $list_path/$prefix\_$suffix.mixcr.out.clonotypes.$type.txt!\n";
<CLONE>;
while (<CLONE>)	{
	next unless (/\t$type/);
	my @clone = split /\t/;
	open FQ, "gunzip -c $list_path/$prefix\_$suffix/$prefix\_$suffix.$clone[0].fastq.gz |" or die "Can't open $list_path/$prefix\_$suffix/$prefix\_$suffix.$clone[0].fastq.gz!\n";
	while (<FQ>)	{
		if (/^\@(\S+)/)	{
			my $read = $1;
			$clone{$ori{$read}}->{$read} = $clone[0];
			<FQ>;
			<FQ>;
			<FQ>;
		}
	}
	close (FQ);
}
close (CLONE);

print "finished reading $type clones!\n";

foreach my $id (sort keys %clone)	{
	open OUT, "> $list_path/$id\_$suffix.$type.reads" or die "Can't write to $list_path/$id\_$suffix.$type.reads!\n";
	foreach my $read (keys %{$clone{$id}})	{
		print OUT "$read\t$clone{$id}->{$read}\n";
	}
	close (OUT);
}

exit;
