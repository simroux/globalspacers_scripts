#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h); 
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to compile numbers for figure about known host vs spacer hits
# Arguments :
# run
";
    die "\n";
}

my $in_file="../../../Analyses/Target_IMGVR_IMGPR/Known_hosts/Consistency_known_host_with_hits.tsv";
my $out_file="Summary_known_host_vs_hits.tsv";

my @categories=("1-9","10-and-more");

my %stats;
open my $tsv,"<",$in_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "uvig"){next;}
    my $cat_hits=&get_cat_hit($tab[2]);
    $stats{$cat_hits}{$tab[9]}{$tab[11]}{$tab[8]}++;
}
close $tsv;

open my $s1,">",$out_file;
print $s1 "n_hits\tref_host_source\thit_host_source\tmatch\tcounts\n";
foreach my $cat (@categories){
    foreach my $ref_host (sort keys %{$stats{$cat}}){
        foreach my $hit_host (sort keys %{$stats{$cat}{$ref_host}}){
            foreach my $match (sort keys %{$stats{$cat}{$ref_host}{$hit_host}}){
                my $rematch="";
                if ($match eq "consistent"){
                    $rematch="known_host_genus";
                }
                elsif($match eq "consistent_before_genus"){
                    $rematch="unclassified_genus_but_consistent_taxonomy";
                }
                elsif($match eq "inconsistent"){
                    $rematch="other_taxon";
                }
                print $s1 $cat."\t".$ref_host."\t".$hit_host."\t".$rematch."\t".$stats{$cat}{$ref_host}{$hit_host}{$match}."\n";
            }
        }
    }
}
close $s1;



sub get_cat_hit(){
    my $n=$_[0];
    my $res="NA";
    if ($n>=1 && $n<10){
        $res="1-9";
    }
    elsif($n>=10){
        $res="10-and-more";
    }
    return $res;
}

