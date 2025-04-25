#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use File::Basename;
use Cwd;
my $h=0;
my $plasmid=0;
my $tag_multitaxa=0;
GetOptions ('help' => \$h, 'h' => \$h, 'p'=>\$plasmid, 'm'=>\$tag_multitaxa);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
	print "# Script to add info about common vs rare and unique vs multi-votu hits to the coverage of hq info
# Arguments :
# -p: for plasmid
";
	die "\n";
}

my $common_spacers="../../Spacer_database/Common_spacer_list.tsv";
my $sp_to_votu="../n_votu_by_cluster.tsv";
my $vr_info="../../../Data/Spacer_db//IMGVR_sequence_information_Oct17.tsv";
my $cl_to_array="All_spacers_vr_hq.tsv";
my $hit_file="All_hits_vr_hq.tsv";

my $in_file="uvig_hq_coverage_by_spacer.tsv";
my $out_file="uvig_hq_coverage_by_spacer-with_sp_info.tsv";

if ($plasmid==1){
    $vr_info="../../../Data/Spacer_db/IMGPR_sequence_information_Aug26.tsv";
    $cl_to_array="All_spacers_pr_complete.tsv";
    $hit_file="All_hits_pr_complete.tsv";
    $in_file="plasmid_complete_coverage_by_spacer.tsv";
    $out_file="plasmid_complete_coverage_by_spacer-with_sp_info.tsv";
}
if ($tag_multitaxa==1){
    $cl_to_array="All_spacers_vr_additional_multitaxa.tsv";
    $hit_file="All_hits_vr_additional_multitaxa.tsv";
    $in_file="uvigadditional_multitaxa_coverage_by_spacer.tsv";
    $out_file="uvigadditional_multitaxa_coverage_by_spacer-with_sp_info.tsv";
}

my %info_sp;
print "Reading $sp_to_votu .. \n";
open my $tsv,"<",$sp_to_votu;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "cluster_id"){next;}
    if ($tab[1]>1){$info_sp{$tab[0]}{"multivotu"}=1;}   
}
close $tsv;

print "Reading $common_spacers .. \n";
open my $tsv,"<",$common_spacers;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "spacer"){next;}
    my @t=split(/\|/,$tab[3]);
    foreach my $cell (@t){
        my @t2=split(",",$cell);
        $info_sp{$tab[0]}{"common"}{$t2[0]}=1;
    }
    
}
close $tsv;

print "Reading uvig to votu (or plasmid to ptu) link from $vr_info\n";
my %uvig_to_votu;
open my $tsv,"<",$vr_info;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "uvig"){next;}
    $uvig_to_votu{$tab[0]}=$tab[2];
}
close $tsv;

print "Loading cluster to array link\n";
my %cl_to_sp;
my $i=0;
open my $tsv,"<",$cl_to_array;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    $cl_to_sp{$tab[0]}{$tab[6]}=1;
    $i++;
    if ($i % 1000000 == 0){
        print "... $i ...\n";
    }
}
close $tsv;

my %store_counts;
$i=0;
print "Counting hit per uvig - array combination\n";
open my $tsv,"<",$hit_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "cluster_id"){next;}
    foreach my $array (keys %{$cl_to_sp{$tab[0]}}){
        if (defined($info_sp{$tab[0]}{"multivotu"})){
            $store_counts{$tab[1]}{$array}{"multivotu"}++;
        }
        else{
            $store_counts{$tab[1]}{$array}{"unique"}++;
        }
        if (defined($info_sp{$tab[0]}{"common"}{$array})){
            $store_counts{$tab[1]}{$array}{"common"}++;
        }
        else{
            $store_counts{$tab[1]}{$array}{"rare"}++;
        }
    }
    $i++;
    if ($i % 100000 == 0){
        print "... $i ...\n";
    }
}
close $tsv;


print "Preparing output file $out_file\n";
open my $s1,">",$out_file;
open my $tsv,"<",$in_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "uvig"){
        print $s1 $_."\tn_hit_common\tn_hit_rare\tn_hit_multivotu\tn_hit_unique\tvotu\ttargeting_type\n";
    }
    elsif($plasmid==1 && $tab[0] eq ""){
        print $s1 $_."\tn_hit_common\tn_hit_rare\tn_hit_multivotu\tn_hit_unique\tptu\ttargeting_type\n";
    }
    else{
        my $hash=$store_counts{$tab[0]}{$tab[1]};
        if (!defined($$hash{"multivotu"})){$$hash{"multivotu"}=0;}
        if (!defined($$hash{"unique"})){$$hash{"unique"}=0;}
        if (!defined($$hash{"common"})){$$hash{"common"}=0;}
        if (!defined($$hash{"rare"})){$$hash{"rare"}=0;}
        my $type="NA";
        if ($tab[3]<200 || $tab[5]<10){ ## Covering less than 200bp of the virus or less than 10 spacers
            $type="low"
        }
        elsif($tab[4]>=10 || $tab[3]>=1000){ ## Covering more than 10% of the genome or more than 1kb
            $type="high";
        }
        else{
            $type="medium";
        }
        print $s1 $_."\t".$$hash{"common"}."\t".$$hash{"rare"}."\t".$$hash{"multivotu"}."\t".$$hash{"unique"}."\t".$uvig_to_votu{$tab[0]}."\t".$type."\n";
    }
}
close $tsv;
close $s1;
