#!/usr/bin/perl -w
use strict;
use Math::Round;

my $TCR_stats = $ARGV[0];
my $InfiniTIL_data = $ARGV[1];
my $selected_neoTCR = $ARGV[2];	# selected neoTCR csv
my $type = $ARGV[3];		# TCR4 or TCR8

unless ($type eq "TCR4" || $type eq "TCR8") { die "type must be one of TCR4 or TCR8\n"; }
print "clonotype\tCB\tUMAP_1\tUMAP_2\tNeo".$type."\tselect\n";
my @type;
if ($type eq "TCR4")	{
	@type = (3, -2);
}	else	{
	@type = (5, -1);
}

my %TCR;
&csv($TCR_stats, \%TCR, 0);
my %data;
&csv($InfiniTIL_data, \%data, 0);
my %neoTCR;
&csv($selected_neoTCR, \%neoTCR, 1);

foreach my $bc (keys %TCR)	{
	if ($TCR{$bc}->[-3] eq "NA" || $TCR{$bc}->[-2] eq "NA" || $TCR{$bc}->[-1] eq "NA")	{
		delete $TCR{$bc};
	}
}

foreach my $cl (keys %neoTCR)	{
	foreach my $cb (sort keys %TCR)	{
		if ($neoTCR{$cl}->[$type[0]] eq $cb)	{
			print "$TCR{$cb}->[-3]\t";
			print join "\t", @{$data{$cb}}[0..2];
			print "\t$TCR{$cb}->[$type[1]]\tCB\n";
			delete $TCR{$cb};
		}	elsif ($cl eq $TCR{$cb}->[-3])	{
			print "$TCR{$cb}->[-3]\t";
			print join "\t", @{$data{$cb}}[0..2];
			print "\t$TCR{$cb}->[$type[1]]\tCL\n";
			delete $TCR{$cb};
		}
	}
	
}

foreach my $cb (sort keys %TCR)	{
	if (defined($data{$cb}))	{
		print "$TCR{$cb}->[-3]\t";
		print join "\t", @{$data{$cb}}[0..2];
		print "\t$TCR{$cb}->[$type[1]]\tNA\n";
	}	else	{
		print STDERR "$cb\n";
	}
}

exit;

sub csv	{
	(my $file, my $out, my $i) = @_;

	open IN, "< $file" or die "Can't open $file!\n";
	<IN>;
	while (my $l = <IN>)	{
		chomp($l);
		my @l;
		while ($l =~ /,/)	{
			if ($l =~ /^(\S)/)	{
				if ($1 eq "\"")	{
					$l =~ s/^\"(\S+?)\"//;
					push @l, $1;
				}	else	{
					$l =~ s/^([^\,]+)//;
					push @l, $1;
				}
				$l =~ s/^,//;
			}
		}
		push @l, $l;
		$out->{$l[$i]} = \@l;
	}
	close (IN);

	return;
}
