#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
	print "# Script to evaluate changes in the ratio of spacer hits for different types of spacers
# Arguments :
# run
";
	die "\n";
}

my $spacer_type_file="../../../Analyses/Spacer_database/Common_spacer_list.tsv";
my $spacer_eco_file="../../../Analyses/Spacer_database/spacer_to_eco-clean.tsv";

my $hit_type_file_vr="../../../Analyses/Target_IMGVR_IMGPR/n_votu_by_cluster.tsv";
my $hit_type_file_pr="../../../Analyses/Target_IMGVR_IMGPR/n_ptu_by_cluster_imgpr.tsv";

my $main_spacer_table="../../../Analyses/Spacer_database/spacer_clusters_and_metadata.tsv";

my $out_file="hit_counts_stats_clean.tsv";

my $i=0;
my %store_spacer;
print "### Reading $spacer_eco_file\n";
open my $tsv,"<",$spacer_eco_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "cluster_id"){next;}
	$store_spacer{$tab[0]}{"eco"}=$tab[1];
    $i++;
    if ($i % 1000000 == 0){print ".. $i ..\n";}
}
close $tsv;

print "### Reading $hit_type_file_vr\n";
$i=0;
open my $tsv,"<",$hit_type_file_vr;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "cluster_id"){next;}
    $store_spacer{$tab[0]}{"n_votu"}=$tab[1];
    $i++;
    if ($i % 1000000 == 0){print ".. $i ..\n";}
}
close $tsv;

print "### Reading $hit_type_file_pr\n";
$i=0;
open my $tsv,"<",$hit_type_file_pr;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "cluster_id"){next;}
    $store_spacer{$tab[0]}{"n_ptu"}=$tab[1];
    $i++;
    if ($i % 100000 == 0){print ".. $i ..\n";}
}
close $tsv;

$i=0;
print "### Reading $spacer_type_file\n";
open my $tsv,"<",$spacer_type_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "spacer"){next;}
    my @t_sets=split(/\|/,$tab[3]);
    foreach my $set (@t_sets){
	    my @t=split(",",$set);
	    $store_spacer{$tab[0]}{"common"}{$t[0]}=1;
    }
    $i++;
    if ($i % 100000 == 0){print ".. $i ..\n";}
}
close $tsv;

$i=0;
print "### Computing the stats\n";
my %tmp=();
my %stats;
open my $tsv,"<",$main_spacer_table;
while(<$tsv>){
	chomp($_);
	%tmp=();
	my @tab=split("\t",$_);
	if ($tab[0] eq "cluster_id"){next;}
	$tmp{"hit_type"}="none";
	$tmp{"sp_type_alpha"}="rare";
	$tmp{"temperate"}="no";
	if (defined($store_spacer{$tab[0]})){
		if (defined($store_spacer{$tab[0]}{"n_votu"})){
			$tmp{"hit_type"}="vr";
			if (defined($store_spacer{$tab[0]}{"n_ptu"})){
				$tmp{"hit_type"}="vr_and_pr";
			}
		}
		elsif(defined($store_spacer{$tab[0]}{"n_ptu"})){
			$tmp{"hit_type"}="pr";
		}
		if (defined($store_spacer{$tab[0]}{"common"})){
			if (defined($store_spacer{$tab[0]}{"common"}{$tab[1]})){
				$tmp{"sp_type_alpha"}="common";
			}
		}
	}
	$tmp{"sp_type_beta"}="single_sample";
	if ($tab[2]>1){
		$tmp{"sp_type_beta"}="multi_sample";
	}
    $stats{"sp_alpha"}{$tmp{"sp_type_alpha"}}{$tmp{"hit_type"}}++;
    $stats{"sp_alpha_and_beta"}{$tmp{"sp_type_alpha"}."_".$tmp{"sp_type_beta"}}{$tmp{"hit_type"}}++;
    if ($tab[6] ne "Unknown"){
        $stats{"eco"}{$store_spacer{$tab[0]}{"eco"}}{$tmp{"hit_type"}}++;
    }
	$i++;
    if ($i % 1000000 == 0){print ".. $i ..\n";}
}
close $tsv;

my @types=("none","vr","pr","vr_and_pr");

open my $s1,">",$out_file;
print $s1 "category\tvalue\ttype\tcount\n";
foreach my $cat (sort keys %stats){
    foreach my $value (sort keys %{$stats{$cat}}){
        foreach my $type (@types){
            if (!defined($stats{$cat}{$value}{$type})){$stats{$cat}{$value}{$type}=0;}
            print $s1 $cat."\t".$value."\t".$type."\t".$stats{$cat}{$value}{$type}."\n";
        }
    }
}
close $s1;