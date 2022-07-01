#!/usr/bin/perl -w
use strict;

open OUT, "> tmp.py" or die "Can't write to tmp.py!\n";
print OUT "from Bio import SeqIO\n";

foreach my $file (glob "*.ab1")	{
	$file =~ s/\.ab1//;
	
	print OUT "records = SeqIO.parse(\"$file.ab1\", \"abi\")\n";
	print OUT "count = SeqIO.write(records, \"$file.fq\", \"fastq\")\n";
	print OUT "records = SeqIO.parse(\"$file.ab1\", \"abi\")\n";
	print OUT "count = SeqIO.write(records, \"$file.fa\", \"fasta\")\n";
}

close (OUT);

print "Run \"python tmp.py\" under python3 environment!\n";

exit;
