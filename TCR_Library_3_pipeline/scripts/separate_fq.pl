#!/usr/bin/perl -w
use strict;

my $format = $ARGV[0]; # log file format: cutadapt.log
my $out_path = $ARGV[1];

foreach my $file (glob "$out_path/*.$format")	{
    $file =~ s/\.$format//;
    print $file,"\n";
    print "read to $file.$format\n";
    print "write to $file.p1.fq\n";
    print "write to $file.p2.fq\n";

    open OUT1, "> $file.p1.fq" or die "Can't open $file.p1.fq!\n";
    open OUT2, "> $file.p2.fq" or die "Can't open $file.p2.fq!\n";
    open IN, "< $file.$format" or die "Can't open $file.$format!\n";
    while (<IN>)    {
        my @rd = split /\s+/;
        if ($rd[2] eq "-1") {
            next;
            print OUT2 "\@$rd[0]\n$rd[3]\n+\n$rd[4]\n";
        }   elsif ($rd[3] eq "0")   {
            print OUT2 "\@$rd[0]\n$rd[6]\n+\n$rd[9]\n";
        }   else    {
            if ($rd[4] eq "451")    {
                print OUT1 "\@$rd[0]\n$rd[5]\n+\n$rd[8]\n";
            }   else    {
                print OUT1 "\@$rd[0]\n$rd[5]\n+\n$rd[9]\n";
                print OUT2 "\@$rd[0]\n$rd[7]\n+\n$rd[11]\n";
            }
        }
    }
    close (IN);
    close (OUT1);
    close (OUT2);
    print "finish\n";
}

exit;
