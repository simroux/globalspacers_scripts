#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
use Custom::Utils;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to prepare the input file for the figure showing how broadly each unique spacer hits viruses
# Arguments :
# run
";
    die "\n";
}

my $out_file_cross="multihit_votu_vs_distance_sample_hqonly.tsv";
my $out_file_full="multihit_votu_totaldistrib_hqonly.tsv";
my $in_full="n_votu_by_cluster_hqonly.tsv";
my $in_distance="../../../Data/Additional_data/selected_subset_spacers_fordist_distance.tsv"; ## Note: this file was generated using the script "get_distrib_dist_by_cluster.pl" in ../../../Analyses/Target_IMGVR_IMGPR/

my %stats;
print "Reading $in_full\n";
open my $tsv,"<",$in_full;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "cluster_id"){next;}
    my $cat=&get_cat($tab[1]);
    $stats{$cat}++;
    if ($cat eq "1000p"){
        print "We see $cat vOTU hits !?!?! --- raw data is $_\n";
    }
}
close $tsv;

my @ordered_cat=("1","[2-5]","[6-10]","[11-50]","[51-100]","[101-1000]");
open my $s1,">",$out_file_full;
print $s1 "category\tcount\n";
foreach my $cat (@ordered_cat){
    if (!defined($stats{$cat})){$stats{$cat}=0;}
    print $s1 $cat."\t".$stats{$cat}."\n";
}
close $s1;

my %cross;
open my $s2,">",$out_file_cross;
print $s2 "cluster\tn_hq_uvig\tmax_dist\n";
print "Reading $in_distance\n";
open my $tsv,"<",$in_distance;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "sp_cluster_id"){next;}
    if ($tab[3]>0){ ## At least one hq uvig target
        if ($tab[3]==1){ ## Only one, distance is 0
            print $s2 $tab[0]."\t".$tab[3]."\t0\n";
        }
        else{
            print $s2 $tab[0]."\t".$tab[3]."\t".$tab[5]."\n";
        }
    }
}
close $tsv;
close $s2;



sub get_cat(){
    my $n=$_[0];
    my $cat="NA";
    if ($n==1){$cat="1";}
    elsif($n<=5){$cat="[2-5]";}
    elsif($n<=10){$cat="[6-10]";}
    elsif($n<=50){$cat="[11-50]";}
    elsif($n<=100){$cat="[51-100]";}
    elsif($n<=1000){$cat="[101-1000]";}
    else{$cat="1000p";}
    return $cat;
}
