#!/usr/bin/perl -w
use strict;

foreach my $file (glob "*.ab1")	{
	$file =~ s/\.ab1//;

	$file =~ s/\(/\\\(/g;
	$file =~ s/\)/\\\)/g;
	$file =~ s/\[/\\\[/g;
	$file =~ s/\]/\\\]/g;

	system "java -jar ~/usr/bin/trimmomatic/trimmomatic-0.39.jar SE -threads 1 -phred33 $file.fq $file.trim.fq TRAILING:40";
}

exit;
