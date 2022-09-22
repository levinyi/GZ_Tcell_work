#!/usr/bin/perl -w
use strict;
use Math::Round;

my $TCR_stats_barcodes = $ARGV[0];	# TCR_stats_barcodes.csv
my $merged_TCR_table = $ARGV[1];	# merged_tcr_table.csv
my $InfiniTIL_data = $ARGV[2];	# InfiniTIL_data.csv 

my %TCR;
&tcr_stats($TCR_stats_barcodes, \%TCR);

my %data;
&infinitil($InfiniTIL_data, \%data);

foreach my $bc (keys %data)	{
	unless (defined($TCR{$bc}))	{
		delete $data{$bc};
	}
}

my %table;
&tcr_table($merged_TCR_table, \%table);

my %ID_bc; my %ID_clonotype; 
foreach my $bc (keys %data)	{
	my $col = "original_clone";
	if ($data{$bc}->{"ID:CD8"} eq "TRUE")	{
		push @{$ID_bc{"CD8"}}, $bc;
		$ID_clonotype{$TCR{$bc}->{$col}}->{"CD8"}++;
		$ID_clonotype{$TCR{$bc}->{$col}}->{"total"}++;
	}	elsif ($data{$bc}->{"ID:CD4"} eq "TRUE")	{
		push @{$ID_bc{"CD4"}}, $bc;
		$ID_clonotype{$TCR{$bc}->{$col}}->{"CD4"}++;
		$ID_clonotype{$TCR{$bc}->{$col}}->{"total"}++;
	}	elsif ($data{$bc}->{"ID:CD4_Treg"} eq "TRUE")	{
		push @{$ID_bc{"CD4_Treg"}}, $bc;
		$ID_clonotype{$TCR{$bc}->{$col}}->{"CD4_Treg"}++;
		$ID_clonotype{$TCR{$bc}->{$col}}->{"total"}++;
	}
}

foreach my $type ("CD8", "CD4", "CD4_Treg")	{
	my @bc = @{$ID_bc{$type}};
	my %neoTCR;
	&neoTCR(\%TCR, \@bc, \%neoTCR);
	open OUT, "> $type\_Cluster_neoTCR.csv" or die "Can't wrtie to $type\_Cluster_neoTCR.csv!\n";
	print OUT "clonotype\tTR_1.0_maxAUC\tTR_1.0_maxCB\tCD4\tCD8A\tFOXP3\tGAPDH\tCD8 (%)\tCD4 (%)\tCD4_Treg (%)\tCB_total\t#of_synTCR\n";
	foreach my $clonotype (sort keys %neoTCR)	{
		print OUT "$clonotype\t".$neoTCR{$clonotype}->{"TR_1.0"}->{"TR_1.0_AUCscore"}."\t".$neoTCR{$clonotype}->{"TR_1.0"}->{"barcode"};
		my $bc = $neoTCR{$clonotype}->{"TR_1.0"}->{"barcode"};
		foreach my $gene ("CD4", "CD8A", "FOXP3", "GAPDH")	{
			print OUT "\t".$data{$bc}->{$gene};
		}
		unless (defined($ID_clonotype{$clonotype}->{"CD8"}))	{$ID_clonotype{$clonotype}->{"CD8"} = 0;}
		unless (defined($ID_clonotype{$clonotype}->{"CD4"}))	{$ID_clonotype{$clonotype}->{"CD4"} = 0;}
		unless (defined($ID_clonotype{$clonotype}->{"CD4_Treg"}))	{$ID_clonotype{$clonotype}->{"CD4_Treg"} = 0;}
		#print OUT "\t".$ID_clonotype{$clonotype}->{"CD8"}." / ".$ID_clonotype{$clonotype}->{"total"};
		#print OUT "\t".$ID_clonotype{$clonotype}->{"CD4"}." / ".$ID_clonotype{$clonotype}->{"total"};
		#print OUT "\t".$ID_clonotype{$clonotype}->{"CD4_Treg"}." / ".$ID_clonotype{$clonotype}->{"total"}."\n";
		print OUT "\t".nearest(0.01, 100*$ID_clonotype{$clonotype}->{"CD8"}/$ID_clonotype{$clonotype}->{"total"});
		print OUT "\t".nearest(0.01, 100*$ID_clonotype{$clonotype}->{"CD4"}/$ID_clonotype{$clonotype}->{"total"});
		print OUT "\t".nearest(0.01, 100*$ID_clonotype{$clonotype}->{"CD4_Treg"}/$ID_clonotype{$clonotype}->{"total"});
		print OUT "\t".$ID_clonotype{$clonotype}->{'total'}."\t$table{$clonotype}\n";
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

sub tcr_table	{
	(my $file, my $table) = @_;

	open IN, "< $file" or die "Can't open $file!\n";
	my $header = <IN>;
	my @header = &parse_line($header);
	my %header;
	map { $header{$header[$_]} = $_ } (0..$#header);
	while (my $l = <IN>)	{
		my @l = &parse_line($l);
		my @tcr = split /\;/, $l[$header{"cdr3s_aa"}];
		my @count = (0, 0);
		foreach my $tcr (@tcr)	{
			if ($tcr =~ /TRA:/)	{
				$count[0]++;
			}	elsif ($tcr =~ /TRB:/)	{
				$count[1]++;
			}
		}
		if ($count[0] == 0 || $count[1] == 0)	{
			$table->{$l[$header{"raw_clonotype_id"}]} = 0;
		}	elsif ($count[0] == 1 && $count[1] == 1)	{
			$table->{$l[$header{"raw_clonotype_id"}]} = 1;
		}	elsif ($count[0] == 2 && $count[1] == 2)	{
			$table->{$l[$header{"raw_clonotype_id"}]} = 4;
		}	else	{
			$table->{$l[$header{"raw_clonotype_id"}]} = 2;
		}
	}
	close (IN);

	return;
}

sub neoTCR	{
	(my $TCR, my $bc, my $neoTCR) = @_;

	foreach my $id (sort @{$bc})	{
		my $tcr = $TCR->{$id};
	        $tcr->{"TR_1.0_AUCscore"} = ($tcr->{"TR_1.0_AUCscore"} eq "NA") ? -1 : $tcr->{"TR_1.0_AUCscore"};
		my $col = "original_clone";
	        if (defined($neoTCR->{$tcr->{$col}})) {
	                if ($neoTCR->{$tcr->{$col}}->{"TR_1.0"}->{"TR_1.0_AUCscore"} < $tcr->{"TR_1.0_AUCscore"}) {
	                        $neoTCR->{$tcr->{$col}}->{"TR_1.0"} = $tcr;
	                }
	        }       else    {
	                $neoTCR->{$tcr->{$col}}->{"TR_1.0"} = $tcr;
	        }
	}
	close (IN);
	
	return;
}
