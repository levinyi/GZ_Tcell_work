#!/usr/bin/perl -w
use strict;
use Math::Round;

my $TCR_stats = $ARGV[0];	# TCR_stats_barcodes.csv
my $TCR_table = $ARGV[1];	# merged_tcr_table.csv
my $data_for_plot = $ARGV[2];	# InfiniTIL_data.csv 

my %TCR;
&csv($TCR_stats, \%TCR);
my %data;
&csv($data_for_plot, \%data);

foreach my $bc (keys %data)	{
	unless (defined($TCR{$bc}))	{
		delete $data{$bc};
	}
}

my %table;
&tcr_table($TCR_table, \%table);

my %ID_bc; my %ID_clonotype; 
foreach my $bc (keys %data)	{
	if ($data{$bc}->[-3] eq "TRUE")	{
		push @{$ID_bc{"CD8"}}, $bc;
		$ID_clonotype{$TCR{$bc}->[-3]}->{'CD8'}++;
		$ID_clonotype{$TCR{$bc}->[-3]}->{'total'}++;
	}	elsif ($data{$bc}->[-2] eq "TRUE")	{
		push @{$ID_bc{"CD4"}}, $bc;
		$ID_clonotype{$TCR{$bc}->[-3]}->{'CD4'}++;
		$ID_clonotype{$TCR{$bc}->[-3]}->{'total'}++;
	}	elsif ($data{$bc}->[-1] eq "TRUE")	{
		push @{$ID_bc{"CD4_Treg"}}, $bc;
		$ID_clonotype{$TCR{$bc}->[-3]}->{'CD4_Treg'}++;
		$ID_clonotype{$TCR{$bc}->[-3]}->{'total'}++;
	}
}

foreach my $type ("CD8", "CD4", "CD4_Treg")	{
	my @bc = @{$ID_bc{$type}};
	my %neoTCR;
	&neoTCR(\%TCR, \@bc, \%neoTCR);
	open OUT, "> $type\_Cluster_neoTCR.csv" or die "Can't wrtie to $type\_Cluster_neoTCR.csv!\n";
	print OUT "clonotype\tneoTCR4_maxAUC\tneoTCR4_maxCB\tneoTCR8_maxAUC\tneoTCR8_maxCB\tCD4\tCD8A\tFOXP3\tGAPDH\tCD8 (%)\tCD4 (%)\tCD4_Treg (%)\tCB_total\t#of_synTCR\n";
	foreach my $clonotype (sort keys %neoTCR)	{
		print OUT "$clonotype\t".$neoTCR{$clonotype}->{"neoTCR4"}->[-2]."\t".$neoTCR{$clonotype}->{"neoTCR4"}->[0]."\t".$neoTCR{$clonotype}->{"neoTCR8"}->[-1]."\t".$neoTCR{$clonotype}->{"neoTCR8"}->[0]."\t";
		my $bc = ($type eq "CD8") ? $neoTCR{$clonotype}->{"neoTCR8"}->[0] : $neoTCR{$clonotype}->{"neoTCR4"}->[0];
		print OUT join "\t", @{$data{$bc}}[3..6];
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

sub csv	{
	(my $file, my $out) = @_;

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
		$out->{$l[0]} = \@l;
	}
	close (IN);

	return;
}

sub tcr_table	{
	(my $file, my $out) = @_;

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
		my @tcr = split /\;/, $l[-2];
		my @count = (0, 0);
		foreach my $tcr (@tcr)	{
			if ($tcr =~ /TRA:/)	{
				$count[0]++;
			}	elsif ($tcr =~ /TRB:/)	{
				$count[1]++;
			}
		}
		if ($count[0] == 0 || $count[1] == 0)	{
			$out->{$l[-6]} = 0;
		}	elsif ($count[0] == 1 && $count[1] == 1)	{
			$out->{$l[-6]} = 1;
		}	elsif ($count[0] == 2 && $count[1] == 2)	{
			$out->{$l[-6]} = 4;
		}	else	{
			$out->{$l[-6]} = 2;
		}
	}
	close (IN);

	return;
}

sub neoTCR	{
	(my $TCR, my $bc, my $neoTCR) = @_;

	foreach my $id (sort @{$bc})	{
		my $tcr = $TCR->{$id};
	        $tcr->[-2] = ($tcr->[-2] eq "NA") ? -1 : $tcr->[-2];
	        $tcr->[-1] = ($tcr->[-1] eq "NA") ? -1 : $tcr->[-1];
	        if (defined($neoTCR->{$tcr->[-3]})) {
	                if ($neoTCR->{$tcr->[-3]}->{'neoTCR4'}->[-2] < $tcr->[-2]) {
	                        $neoTCR->{$tcr->[-3]}->{'neoTCR4'} = $tcr;
	                }
	                if ($neoTCR->{$tcr->[-3]}->{'neoTCR8'}->[-1] < $tcr->[-1]) {
	                        $neoTCR->{$tcr->[-3]}->{'neoTCR8'} = $tcr;
	                }
	        }       else    {
	                $neoTCR->{$tcr->[-3]}->{'neoTCR4'} = $tcr;
	                $neoTCR->{$tcr->[-3]}->{'neoTCR8'} = $tcr;
	        }
	}
	close (IN);
	
	return;
}
