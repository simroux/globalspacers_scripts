#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
use Custom::Utils;
my $h=0;
my $tag_plasmid=0;
my $tag_multitaxa=0;
GetOptions ('help' => \$h, 'h' => \$h, 'p'=>\$tag_plasmid, 'm'=>\$tag_multitaxa); # , 'i=s'=>\$in_file, 'r=s'=>\$out_dir_root);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to get the table of additional spacer hits for cases with substantial targeting, to be loaded later in duckdb
# Arguments :
# -p: for plasmid
# -m: for multitaxa
# run
";
    die "\n";
}

## Note - only HQ viruses, complete plasmids, or viruses flagged as targeted by multitaxa (and not already in HQ)
my $infile_cluster="../../../Data/Spacer_db/All_spacers_info_filtered_clusters-Jul19-24.tsv";

## For viruses
my $in_coverage="../Target_coverage/uvig_hq_coverage_by_spacer.tsv";
my $in_hit_file="../results_nr_imgvrhq/nr_spacers_hq_vs_IMGVR_db_all.matches_clean.tsv.gz";
my $out_file="../../../Data/Spacer_db/Additional_hits_imgvr-Dec13.tsv"; 

## For plasmids
if ($tag_plasmid==1){
    $in_coverage="../Target_coverage/plasmid_complete_coverage_by_spacer.tsv";
    $in_hit_file="../results_nr_imgprhq/nr_spacers_hq_vs_IMGPR_db_all.matches_clean.tsv.gz";
    $out_file="../../../Data/Spacer_db/Additional_hits_imgpr-Dec13.tsv";
}

## For mlutitaxa uvigs
if ($tag_multitaxa==1){
    $in_coverage="../../../Data/Additional_data/List_multiclass_uvigs_nohq.txt";
    $in_hit_file="../results_nr_imgprhq/nr_spacers_hq_vs_IMGPR_db_all.matches_clean.tsv.gz";
    $out_file="../../../Data/Spacer_db/Additional_hits_imgvr-multitaxa-Dec17.tsv"; 
}


my %check;
if ($in_coverage eq "../../../Data/Additional_data/List_multiclass_uvigs_nohq.txt"){ ## Special for multitaxa, because we get a list instead of a coverage file we filter from
    print "Reading $in_coverage\n";
    open my $tsv,"<",$in_coverage;
    while(<$tsv>){
        chomp($_);
        my @tab=split("\t",$_);
        if ($tab[0] eq "uvig"){next;}
        $check{$tab[0]}=1;
    }
    close $tsv;
}
else{
    print "Reading $in_coverage\n";
    open my $tsv,"<",$in_coverage;
    while(<$tsv>){
        chomp($_);
        my @tab=split("\t",$_);
        ### Taking everything with 10 spacers or more, 200bp or more, i.e. nothing "low"
        if ($tab[5]>=10 && $tab[3]>=200){
            $check{$tab[0]}=1;
        }
        elsif($in_lie eq ""){
            $check{$tab[0]}=1;
        }
    }
    close $tsv;
}

my %check_cluster;
my $i=0;
print "Reading $infile_cluster\n";
open my $tsv,"<",$infile_cluster;
while(<$tsv>){
	chomp($_);
	my @tab=split("\t",$_);
	if ($tab[0] eq "cluster_id"){next;}
	$check_cluster{$tab[0]}=1;
	$i++;
	if ($i % 1000000 == 0){print "\t $i ... \n";}
}
close $tsv;

open my $s1,">",$out_file;
$i=0;
print $s1 "cluster_id\ttarget_id\thit_start\thit_end\thit_strand\tn_mismatches\tCIGAR\tMD\tspacer\tprotospacer\tupstream\tdownstream\n";
print "Processing $in_hit_file ...\n";
open my $tsv,"gunzip -c $in_hit_file |";
$i=0;
while(<$tsv>){
	chomp($_);
	my @tab=split("\t",$_);
	if ($tab[0] eq "Spacer id"){next;}
       if ($check{$tab[1]}==1){ ## UViG of interest
           if ($tab[12] eq "mismatches"){ ## The only flag that prevented this hit to be taken in was the fact that there were too many ismatches
               if ($check_cluster{$tab[0]}==1){
                   print $s1 join("\t",@tab[0..11])."\n";
                    $i++;
                    if ($i % 500000 == 0){print "\t $i ... \n";}
                }
                else{
                    # print "$tab[0] is not interesting anymore because not in the final db\n";
                }
            }
        }
	}
	close $tsv;
}
close $s1;
