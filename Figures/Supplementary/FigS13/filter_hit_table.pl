#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
use Custom::Utils;
my $h=0;
my $in_file;
my $selected_array;
GetOptions ('help' => \$h, 'h' => \$h, 'i=s'=>\$in_file, 'a=s'=>\$selected_array);
if ($h==1 || $in_file eq "" || $selected_array eq ""){ # If asked for help or did not set up any argument
    print "# Script to get simplified tables for figure
# Arguments :
# -i: input file (e.g. IMGVR_UViG_3300007356_000002_spacer_hits_full.tsv)
# -a: repeat of interest (e.g. Ac_12820)
";
    die "\n";
}

## Get count by sample and get headers while we are at it
open my $tsv,"<",$in_file;
my %header;
my %count;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "cluster_id"){
        for (my $i=0;$i<=$#tab;$i++){$header{$tab[$i]}=$i;}
    }
    else{
        my $array=$tab[$header{"crispr_array"}];
        my $sample=$tab[$header{"library"}];
        if ($array eq $selected_array){
            $count{$sample}++;
        }
    }
}
close $tsv;



my @selected_columns=("target_id","crispr_array","library","cluster_id","hit_start","hit_end","n_mismatches","n_spacers","is_common","n_votu","individual_sample_cover");
open my $s1,">",$out_file;
open my $tsv,"<",$in_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "cluster_id"){
        print $s1 join("\t",@selected_columns)."\t"."hits_in_sample\n";
    }
    else{
        my $array=$tab[$header{"crispr_array"}];
        my $sample=$tab[$header{"library"}];
        if ($array eq $selected_array){
            print $array."\t".$sample."\n";
            my $line="";
            foreach my $col (@selected_columns){
                $line.=$tab[$header{$col}]."\t";
            }
            $line.=$count{$sample};
            print $s1 $line."\n";
        }
    }
}
close $tsv;
close $s1;