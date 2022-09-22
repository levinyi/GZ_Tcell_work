#!/usr/bin/perl -w
use strict;
use Math::Round;

my $TCR_stats_barcodes = $ARGV[0];	# tcr_stats_barcodes.csv
my $InfiniTIL_data = $ARGV[1];		# InfiniTIL_data.csv
my $selected_neoTCR = $ARGV[2];		# selected neoTCR csv

print "clonotype\tCB\tUMAP_1\tUMAP_2\tTR_1.0\tselect\n";

my %TCR;
&tcr_stats($TCR_stats_barcodes, \%TCR);
my %data;
&infinitil($InfiniTIL_data, \%data);
my %neoTCR;
&selected($selected_neoTCR, \%neoTCR);

foreach my $bc (keys %TCR)	{
	if ($TCR{$bc}->{"TR_1.0_AUCscore"} eq "NA")	{
		delete $TCR{$bc};
	}
}

foreach my $cl (keys %neoTCR)	{
	foreach my $cb (sort keys %TCR)	{
		if ($neoTCR{$cl}->{"TR_1.0_maxCB"} eq $cb)	{
			print $TCR{$cb}->{"original_clone"}."\t".$data{$cb}->{"barcode"}."\t".$data{$cb}->{"UMAP_1"}."\t".$data{$cb}->{"UMAP_2"}."\t".$TCR{$cb}->{"TR_1.0_AUCscore"}."\tCB\n";
			delete $TCR{$cb};
		}	elsif ($cl eq $TCR{$cb}->{"original_clone"})	{
			print $TCR{$cb}->{"original_clone"}."\t".$data{$cb}->{"barcode"}."\t".$data{$cb}->{"UMAP_1"}."\t".$data{$cb}->{"UMAP_2"}."\t".$TCR{$cb}->{"TR_1.0_AUCscore"}."\tCL\n";
			delete $TCR{$cb};
		}
	}
}

foreach my $cb (sort keys %TCR)	{
	if (defined($data{$cb}))	{
		print $TCR{$cb}->{"original_clone"}."\t".$data{$cb}->{"barcode"}."\t".$data{$cb}->{"UMAP_1"}."\t".$data{$cb}->{"UMAP_2"}."\t".$TCR{$cb}->{"TR_1.0_AUCscore"}."\tNA\n";
	}	else	{
		print STDERR "$cb\n";
	}
}

exit;

sub parse_line	{
	my $l = $_[0];

	my @l;
	
	chomp($l);
	while ($l =~ /,/)	{
		if ($l =~ /^(\S)/)	{
			if ($1 eq "\"")	{
				if ($l =~ s/^\"\",/,/)	{
					push @l, " ";
				}	elsif ($l =~ s/^\"(\S+?)\"//)	{
					push @l, $1;
				}
			}	else	{
				$l =~ s/^([^\,]+)//;
				push @l, $1;
			}
			$l =~ s/^,//;
		}
	}
	$l =~ s/\"//g;
	push @l, $l;

	return @l;
}

sub tcr_stats	{
	(my $file, my $tcr) = @_;

	open IN, "< $file" or die "Can't open $file!\n";
	my $header = <IN>;
	my @header = &parse_line($header);
	my %header;
	map { $header{$header[$_]} = $_ } (1..$#header);
	while (my $l = <IN>)	{
		my @l = &parse_line($l);
		map { $tcr->{$l[$header{"barcode"}]}->{$header[$_]} = $l[$_] } (1..$#header);
	}
	close (IN);

	return;
}

sub infinitil	{
	(my $file, my $data) = @_;

	open IN, "< $file" or die "Can't open $file!\n";
	my $header = <IN>;
	my @header = &parse_line($header);
	$header[0] = "barcode";
	my %header;
	map { $header{$header[$_]} = $_ } (0..$#header);
	while (my $l = <IN>)	{
		my @l = &parse_line($l);
		map { $data->{$l[$header{"barcode"}]}->{$header[$_]} = $l[$_] } (0..$#header);
	}
	close (IN);

	return;
}

sub selected	{
	(my $file, my $neoTCR) = @_;

	open IN, "< $file" or die "Can't open $file!\n";
	my $header = <IN>;
	my @header = &parse_line($header);
	my %header;
	map { $header{$header[$_]} = $_ } (0..$#header);
	while (my $l = <IN>)	{
		my @l = &parse_line($l);
		map { $neoTCR->{$l[$header{"clonotype"}]}->{$header[$_]} = $l[$_] } (0..$#header);
	}
	close (IN);

	return;
}
