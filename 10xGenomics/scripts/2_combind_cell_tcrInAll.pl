#!/usr/bin/perl
#combind cell lable and the VDJ clone types

open IN1,"<./clonotypes.csv";
open IN2,"<./filtered_contig_annotations.csv";
open IN3,"<../GB001003E2L1_RNAseq/filtered_feature_bc_matrix/barcodes.tsv";
my @cells;
while (<IN3>){
	chomp $_;
	push @cells,$_;
	}

my %clonotype;
while(<IN1>){
	my @it=split ',',$_;
	$clonotype{$it[0]}->[1]=$it[1];
	$clonotype{$it[0]}->[2]=$it[2];
	$clonotype{$it[0]}->[3]=$it[3];
	}

my $last_cell;
my %cell_clt;
while(<IN2>){
	next if /^barcode\t/;
	@it=split ',',$_;
	my $cell=$it[0];
	next if $it[3] eq "False";#high_confidence
	next if $it[12] eq "None";#cdr3
	next if $it[16] eq "None";#raw_clonotype_id
	
	$cell_clt{$cell}->{'clonotype'}=$it[16] unless exists $cell_clt{$cell}->{'clonotype'};
	$cell_clt{$cell}->{'clonotype_freq'}=$clonotype{$it[16]}->[1] unless exists $cell_clt{$cell}->{'clonotype_freq'};
	$cell_clt{$cell}->{'clonotype_prop'}=$clonotype{$it[16]}->[2] unless exists $cell_clt{$cell}->{'clonotype_prop'};
	if($it[5] eq 'TRA'){
		$cell_clt{$cell}->{'TRA_V'}.='_'.$it[6];
		$cell_clt{$cell}->{'TRA_V_most'}=$it[6] unless exists $cell_clt{$cell}->{'TRA_V_most'};
		$cell_clt{$cell}->{'TRA_J'}.='_'.$it[8];
		$cell_clt{$cell}->{'TRA_J_most'}=$it[8] unless exists $cell_clt{$cell}->{'TRA_J_most'};
		$cell_clt{$cell}->{'TRA_CDR3'}.='_'.$it[12];
		$cell_clt{$cell}->{'TRA_CDR3_most'}=$it[12] unless exists $cell_clt{$cell}->{'TRA_CDR3_most'};
		$cell_clt{$cell}->{'TRA_umi'}.='_'.$it[15];
		$cell_clt{$cell}->{'TRA_umi_most'}=$it[15] unless exists $cell_clt{$cell}->{'TRA_umi_most'};
		$cell_clt{$cell}->{'TRA_n'}+=1;
		}
	if($it[5] eq 'TRB'){
		$cell_clt{$cell}->{'TRB_V'}.='_'.$it[6];
		$cell_clt{$cell}->{'TRB_V_most'}=$it[6] unless exists $cell_clt{$cell}->{'TRB_V_most'};
		$cell_clt{$cell}->{'TRB_J'}.='_'.$it[8];
		$cell_clt{$cell}->{'TRB_J_most'}=$it[8] unless exists $cell_clt{$cell}->{'TRB_J_most'};
		$cell_clt{$cell}->{'TRB_CDR3'}.='_'.$it[12];
		$cell_clt{$cell}->{'TRB_CDR3_most'}=$it[12] unless exists $cell_clt{$cell}->{'TRB_CDR3_most'};
		$cell_clt{$cell}->{'TRB_umi'}.='_'.$it[15];
		$cell_clt{$cell}->{'TRB_umi_most'}=$it[15] unless exists $cell_clt{$cell}->{'TRB_umi_most'};
		$cell_clt{$cell}->{'TRB_n'}+=1;
		}
	}

open OUT,">./cell_clonotype2.csv";
print OUT "barcode,clonotype,clonotype_freq,clonotype_prop,";
print OUT "TRAV,TRAV_most,TRAJ,TRAJ_most,TRA_CDR3,TRA_CDR3_most,TRA_n,TRA_umi,TRA_umi_most,";
print OUT "TRBV,TRBV_most,TRBJ,TRBJ_most,TRB_CDR3,TRB_CDR3_most,TRB_n,TRB_umi,TRB_umi_most\n";
foreach my $cel (@cells){
	print OUT "$cel,$cell_clt{$cel}->{'clonotype'},$cell_clt{$cel}->{'clonotype_freq'},$cell_clt{$cel}->{'clonotype_prop'},";
	print OUT "$cell_clt{$cel}->{'TRA_V'},";
	print OUT "$cell_clt{$cel}->{'TRA_V_most'},";
	print OUT "$cell_clt{$cel}->{'TRA_J'},";
	print OUT "$cell_clt{$cel}->{'TRA_J_most'},";
	print OUT "$cell_clt{$cel}->{'TRA_CDR3'},";
	print OUT "$cell_clt{$cel}->{'TRA_CDR3_most'},";
	print OUT "$cell_clt{$cel}->{'TRA_n'},";
	print OUT "$cell_clt{$cel}->{'TRA_umi'},";
	print OUT "$cell_clt{$cel}->{'TRA_umi_most'},";

	print OUT "$cell_clt{$cel}->{'TRB_V'},";
	print OUT "$cell_clt{$cel}->{'TRB_V_most'},";
	print OUT "$cell_clt{$cel}->{'TRB_J'},";
	print OUT "$cell_clt{$cel}->{'TRB_J_most'},";
	print OUT "$cell_clt{$cel}->{'TRB_CDR3'},";
	print OUT "$cell_clt{$cel}->{'TRB_CDR3_most'},";
	print OUT "$cell_clt{$cel}->{'TRB_n'},";
	print OUT "$cell_clt{$cel}->{'TRB_umi'},";
	print OUT "$cell_clt{$cel}->{'TRB_umi_most'}";

	print OUT "\n";

	}

