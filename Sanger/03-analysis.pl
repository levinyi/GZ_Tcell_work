#!/usr/bin/perl -w
use strict;

my $index = $ARGV[0];

my @f = glob "*.ab1";

my %p;
foreach my $f (@f)	{
	my $p;
	if ($f =~ /\((.*?)\)/)	{
		$p = $1;
	}
	
	$f =~ s/\.ab1//;
=item
	$f =~ s/\(/\\\(/;
	$f =~ s/\)/\\\)/;
	$f =~ s/\[/\\\[/;
	$f =~ s/\]/\\\]/;
=cut
	push @{$p{$p}}, $f;
}

my $phred = "!\"#\$%&'()*+,-./0123456789:;<=>?\@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";
my @phred = split //, $phred;
my %phred;
map { $phred{$phred[$_]} = $_ } (0..$#phred);

foreach my $p (sort keys %p)	{
	next if ($#{$p{$p}} == 0);
	next if (-s "$p.seq");

	&convert($p, $p{$p}->[0], "F", \%phred);	
	&convert($p, $p{$p}->[1], "R", \%phred);

	system "cap3 $p.seq > $p.seq.cap.log 2>&1";

	system "blastn -db $index -query $p.seq.cap.contigs -out $p.seq.cap.contigs.m8.out -outfmt 6 -evalue 1e-5";
	system "blastn -db $index -query $p.seq.cap.contigs -out $p.seq.cap.contigs.out -evalue 1e-5";
}

exit;

sub convert {
	(my $id, my $file, my $end, my $phred) = @_;

	open IN, "< $file.trim.fq" or die "Can't open $file.trim.fq!\n";
	open OUT1, ">> $id.seq" or die "Can't open $id.seq!\n";
	open OUT2, ">> $id.seq.qual" or die "Can't open $id.seq.qual!\n";
	while (my $l = <IN>)	{
		#$l =~ s/^@//;
		#print OUT1 ">$l";
		#print OUT2 ">$l";
		print OUT1 ">$id.$end\n";
		print OUT2 ">$id.$end\n";
		$l = <IN>;
		print OUT1 $l;
		$l = <IN>;
		$l = <IN>;
		$l =~ s/\s//g;
		my @l = split //, $l;
		my @q = map { $phred->{$_} } @l;
		print OUT2 join " ", @q;
		print OUT2 "\n";
	}
	close (IN);
	close (OUT1);
	close (OUT2);
}

exit;
