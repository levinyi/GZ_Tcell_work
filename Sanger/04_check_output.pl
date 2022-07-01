#!/usr/bin/perl -w
use strict;

my $index = $ARGV[0];

my %index; my $id;
open IN, "< $index" or die "Can't open $index!\n";
while (<IN>)	{
	if (/>(\S+)/)	{
		$id = $1;
	}	else	{
		s/\s//g;
		$index{$id} .= $_;
	}
}
close (IN);

map { $index{$_} = length($index{$_}) } keys %index;

print "#qry\t#ref\t#ref_len\t#aln_len\t#ident\t#mismatch\t#indel\t#ref_vs_aln\t#D40 (smaller is better)\n";
foreach my $file (glob "*.seq.cap.contigs")	{
	$file =~ s/\.seq.cap.contigs//;
	
	unless (-s "$file.seq.cap.contigs.m8.out")	{
		print "$file\tNA\n";
		next;
	}
	
	my @hit;
	open IN, "< $file.seq.cap.contigs.m8.out" or die "Can't open $file.seq.cap.contigs.m8.out!\n";
	while (<IN>)	{
		@hit = split /\s+/;
		print "$file\t$hit[1]\t$index{$hit[1]}\t$hit[3]\t$hit[2]\t$hit[4]\t$hit[5]";
		my @range = sort {$a<=>$b}  @hit[8..9];
		if ($hit[2] eq "100.000" && $range[0] eq 1 && $range[1] eq $index{$hit[1]} && $hit[3] eq $index{$hit[1]})	{
			print "\tTRUE";
		}	else	{
			print "\tFALSE";
		}
		last;
	}
	close (IN);

	my $seq;
	open IN, "< $file.seq.cap.contigs" or die "Can't open $file.seq.cap.contigs!\n";
	<IN>;
	while (<IN>)	{
		if (/>/)	{
			print STDERR "$file has more than 1 assembled seq!\n";
			last;
		}
		s/\s//g;
		$seq .= $_;
	}
	close (IN);

	my @qual;
	open IN, "< $file.seq.cap.contigs.qual" or die "Can't open $file.seq.cap.contigs.qual!\n";
	<IN>;
	while (<IN>)	{
		last if (/>/);
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

	print "\t$D40\n";
}

exit;
#!/usr/bin/perl -w
use strict;

foreach my $file (glob "*.contigs")	{
}

exit;
