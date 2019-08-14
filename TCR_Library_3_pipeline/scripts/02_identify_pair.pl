#!/usr/bin/perl -w
use strict;
use File::Basename;

my $list = $ARGV[0];
my $type = $ARGV[1];   # B_A or A_B

my $list_path = dirname($list);

my %ABtype;
if ($type eq 'B_A'){
	$ABtype{'TRB'} = 'R1';
	$ABtype{'TRA'} = 'R2';
}
elsif ($type eq 'A_B'){
	$ABtype{'TRA'} = 'R1';
	$ABtype{'TRB'} = 'R2';
}else{
	print "wrong $type!\n";
	exit 1;
}
	
open LIST, "< $list" or die "Can't open $list!\n";
while (my $id = <LIST>)	{
	$id =~ s/\s//g;

	my %TRA;
	open TRA, "< $list_path/$id\_$ABtype{'TRA'}.TRA.reads" or die "Can't open $list_path/$id\_$ABtype{'TRA'}.TRA.reads!\n";
	while (<TRA>)	{
		my @read = split /\s+/;
		$TRA{$read[0]} = $read[1];
	}
	close (TRA);

	my %TRB;
	open TRB, "< $list_path/$id\_$ABtype{'TRB'}.TRB.reads" or die "Can't open $list_path/$id\_$ABtype{'TRB'}.TRB.reads!\n";
	while (<TRB>)	{
		my @read = split /\s+/;
		next unless (defined($TRA{$read[0]}));
		$TRB{$read[0]} = $read[1];
	}
	close (TRB);

	my %pair;
	open OUT, "> $list_path/$id.pairs" or die "Can't write to $list_path/$id.pairs!\n";
	foreach my $read (keys %TRB)	{
		print OUT "$read\t$TRA{$read}\t$TRB{$read}\n";
		$pair{$TRA{$read}}->{$TRB{$read}}++;
	}
	close (OUT);

	open OUT, "> $list_path/$id.pairs.freq" or die "Can't write to $list_path/$id.pairs.freq!\n";
	foreach my $TRA (sort keys %pair)	{
		foreach my $TRB (sort keys %{$pair{$TRA}})	{
			print OUT "$TRA\t$TRB\t$pair{$TRA}->{$TRB}\n";
		}
	}
	close (OUT);
	
}
close (LIST);

exit;
