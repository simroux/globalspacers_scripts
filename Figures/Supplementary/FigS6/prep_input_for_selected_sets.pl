#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h); 
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to get a nice table describing the repeats we used for the beta diversity analysis
# Arguments :
# none
";
    die "\n";
}

my $min_spacer=10; ## min spacer per set to consider
my $min_maxcover=20; ## min max cover spacer to consider per set
my $min_samples=10; ## min sample per array to consider - this time to 10, for the subsample, to get a nice curve

my $taxa_list="list_taxa_for_betadiv.tsv";
my $in_file_set_size="../../../Analyses/Spacer_database/spacer_sets_alphadiv.tsv";
my $file_array_info="../../../Data/Spacer_db/Array_info_filtered_for_db-Apr23-26.tsv";
my $file_sample_info="../../../Data/Spacer_db/Runs_to_ecosystem_and_sequencing_and_study_for_db-Jul28-24.tsv";

my $out_file="input_description_betadiv.tsv";

my $n=0;
my %check_taxon;
open my $tsv,"<",$taxa_list;
while(<$tsv>){
    chomp($_);
    $n++;
    my $ori=$_;
    $_=~s/;Other$//;
    $check_taxon{$_}=$ori;
}
close $tsv;

my %info_array;
open my $tsv,"<",$file_array_info;
print ".. reading $file_array_info ..\n";
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "repeat_cluster"){next;}
    $info_array{$tab[0]}{"type"}=$tab[1];
    my $tax="";
    my @t=split(";",$tab[3]);
    if (defined($check_taxon{$t[0]})){$tax=$t[0];}
    for (my $i=0;$i<=$#t;$i++){
        my $test=join(";",@t[0..$i]);
        if (defined($check_taxon{$test})){
            $tax=$check_taxon{$test};
        }
    }
    if ($tax eq ""){
        $tax="NA";
        if ($tab[3] ne "NA" && !($tab[3]=~/Viruses/)){die("\n");}
    }
    $info_array{$tab[0]}{"taxo"}=$tax;
}
close $tsv;

my %info_sample;
open my $tsv,"<",$file_sample_info;
print ".. reading $file_sample_info ..\n";
while(<$tsv>){
      chomp($_);
      my @tab=split("\t",$_);
      if ($tab[0] eq "library"){next;}
      $info_sample{$tab[0]}{"eco"}=$tab[3];
      $info_sample{$tab[1]}{"eco"}=$tab[3];
}
close $tsv;

my %selected;
my %info_set;
open my $tsv,"<",$in_file_set_size;
print ".. reading $in_file_set_size ..\n";
my $i=0;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "array"){next;}
    if ($tab[2]>=$min_spacer && $tab[4]>=$min_maxcover){
        $selected{$tab[0]}{$tab[1]}=$tab[2];
        $info_set{$tab[0]}{$tab[1]}{"max_cover"}=$tab[4];
    }
    $i++;
    if ($i % 1000000 == 0){print " .. $i\n";}
}
close $tsv;


my %check_array;

my %tmp;
open my $s1,">",$out_file;
print $s1 "repeat\ttaxonomy\ttype\tecosystem\tn_samples\tn_spacers\n";
foreach my $array (keys %selected){
    %tmp=();
    my $test=scalar(keys %{$selected{$array}});
    if ($test>=$min_samples){
        $check_array{$array}=1;
        foreach my $sample (sort keys %{$selected{$array}}){
            $tmp{$info_sample{$sample}{"eco"}}{"n_sample"}++;
            $tmp{$info_sample{$sample}{"eco"}}{"n_spacers"}+=$selected{$array}{$sample};
        }
    }
    foreach my $eco (sort keys %tmp){
        print $s1 $array."\t".$info_array{$array}{"taxo"}."\t".$info_array{$array}{"type"}."\t".$eco."\t".$tmp{$eco}{"n_sample"}."\t".$tmp{$eco}{"n_spacers"}."\n";
    }
}
close $s1;

my $test=scalar(keys %check_array);
print "Total relevant array -> $test\n";