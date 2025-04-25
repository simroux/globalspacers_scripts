#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
my $panel="";
GetOptions ('help' => \$h, 'h' => \$h, 'p=s'=>\$panel);
if ($h==1 || $panel eq ""){ 
    print "# Script to get the input file for the coverage visualization of viruses targeted by distant taxa
# Arguments :
# -p: c or d depending on the panel we want
";
    die "\n";
}

my $in_dir="";
if ($panel eq "c"){
    $in_dir="panel_c/";
}
elsif ($panel eq "d"){
    $in_dir="panel_d/";
}
else{die("I don't know this panel\n");}
if (!(-d $in_dir)){
    &run_cmd("mkdir $in_dir");
}

my $in_net="../../Main/Fig_5/".$in_dir."/selected_uvigs_coverage_net-edges.tsv";
my $info_file="../../../Analyses/Target_Custom_datasets/All_spacers_vs_selected_viruses/nr_spacers_hq_vs_selected_uvigs_db_all_hits_spacer_info.tsv";
my $hit_file="../../../Analyses/Target_Custom_datasets/All_spacers_vs_selected_viruses/nr_spacers_hq_vs_selected_uvigs_db_all_hits.tsv";
my $out_file=$in_dir."/coverage_for_fig.tsv";

my %check;
my %check_host;
print "Reading $in_net ..\n";
open my $tsv,"<",$in_net;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "node_1"){next;}
    $check{$tab[0]}{$tab[1]}=1;
    $check_host{$tab[1]}=1;
    print "Checking $tab[0] / $tab[1]\n";
}
close $tsv;

open my $s1,">",$out_file;
print $s1 "virus\tspacer\tstart\tend\tn_mis\ttaxon\n";
my %check_spacer;
print "Reading $info_file..\n";
my $i=0;
open my $tsv,"<",$info_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "spid"){next;}
    ### Note: we ignore arrays that are without taxonomy or with only contig taxonomy
    if ($tab[25]=~/^Genome/){
        my $clean=&clean($tab[26]);
        if ($check_host{$clean}==1){
            $check_spacer{$tab[0]}{$clean}=1;
        }
    }
    $i++;
    if ($i % 10000 == 0){
        print "... $i ...\n";
    }
}
close $tsv;
print "Now reading $hit_file ...\n";
open my $tsv,"<",$hit_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "Spacer id"){next;}
    my $virus=$tab[1];
    if ($tab[1]=~/\|/){
        my @t=split(/\|/,$tab[1]);
        $tab[1]=$t[0];
    }
    if (defined($check_spacer{$tab[0]})){
        if (defined($check{$tab[1]})){
            foreach my $taxon (keys %{$check_spacer{$tab[0]}}){
                if (defined($check{$tab[1]}{$taxon})){
                    print $s1 $tab[1]."\t".$tab[0]."\t".$tab[2]."\t".$tab[3]."\t".$tab[5]."\t".$taxon."\n";
                }
            }
        }
    }
    
}
close $tsv;
close $s1;



sub clean(){
    my @t=split(";",$_[0]);
    $t[1]=~s/p__//;
    $t[5]=~s/g__//;
    my $clean_tax=$t[1].";".$t[5];
    return $clean_tax;
}
