#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
my %translate_level=("d"=>"domain","p"=>"phylum","c"=>"class","o"=>"order","f"=>"family","g"=>"genus","s"=>"species");
GetOptions ('help' => \$h, 'h' => \$h); # , 'i=s'=>\$in_file, 'r=s'=>\$out_dir_root);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
	print "# Script to prepare a file to show the distribution of spacer coverage for a subsample of spacer sets
# Arguments :
# none
";
	die "\n";
}

## For the plot we will take 100 spacer sets (all from different arrays) with at least 1000 total spacers (i.e. total coverage) and a maximum coverage of 20x
## This is done separately for each ecosystem, except for Other, where we take the 75 that qualify
## Total should be: 1,175 (100 * 11, plus 75)
my $min_cumulated_cover=1000;
my $min_max_cover=20;

my $to_subsample=100;
my $to_subsample_spacers=1000;

my $file_array_info="Data/Spacer_db/Array_info_filtered_for_db-Apr12-24.tsv";
my $file_sample_info="Data/Spacer_db/Runs_to_ecosystem_and_sequencing_and_study_for_db-Jul28-24.tsv";

my $in_alpha_base="Analyses/Spacer_database/spacer_sets_alphadiv.tsv";
my $in_file_sets="Analyses/Spacer_database/spacer_clusters_and_metadata_for_alphadiv.tsv";
my $in_file_sets_bis="Analyses/Spacer_database/spacer_clusters_and_metadata_for_alphadiv_md10.tsv";
my $out_file="spacer_sets_alphadiv_forplot_subsamples.tsv";

my %info_array;
open my $tsv,"<",$file_array_info;
print "Reading $file_array_info\n";
while(<$tsv>){
      chomp($_);
      my @tab=split("\t",$_);
      if ($tab[0] eq "repeat_cluster"){next;}
      $info_array{$tab[0]}{"type"}=$tab[1];
      $info_array{$tab[0]}{"taxo_source"}=$tab[2];
      $info_array{$tab[0]}{"taxo_genus"}=$tab[6];
      if ($tab[3] eq "NA"){$info_array{$tab[0]}{"taxo_lvl"}="none";}
      else{$info_array{$tab[0]}{"taxo_lvl"}=&guess_level($tab[3]);}
}
close $tsv;

my %info_sample;
open my $tsv,"<",$file_sample_info;
print "Reading $file_sample_info\n";
while(<$tsv>){
      chomp($_);
      my @tab=split("\t",$_);
      if ($tab[0] eq "library"){next;}
      $info_sample{$tab[0]}{"eco"}=$tab[3];
      $info_sample{$tab[1]}{"eco"}=$tab[3];
}
close $tsv;

my %to_subsample;
open my $tsv,"<",$in_alpha_base;
print "Reading $in_alpha_base\n";
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "array"){next;}
    my $set=$tab[1]."_".$tab[0];
    my $total_cover=$tab[3];
    my $max_cover=$tab[4];
    if ($total_cover>=$min_cumulated_cover && $max_cover>=$min_max_cover){
        $to_subsample{$info_sample{$tab[1]}{"eco"}}{$tab[1]."_".$tab[0]}=1;
    }
}
close $tsv;

## Make the selection
my %check;
my %seen;
foreach my $eco (sort keys %to_subsample){
    print "Selecting for $eco\n";
    my @tab=keys %{$to_subsample{$eco}};
    if (scalar(@tab)<$to_subsample){
        if ($eco eq "Other"){
            print "For Other, we take everything\n";
            foreach my $set (@tab){
                if ($set=~/^(.*)_(Ac_\d+)$/){
                    my $sample=$1;
                    my $array=$2;
                    $check{$array}{$sample}=1;
                }
            }
        }
        else{
            die("Pblm, for $eco, we don't have enough data\n");
        }
    }
    print "We have ".scalar(@tab)." to pick from for $eco\n";
    &fisher_yates_shuffle(\@tab);
    my $n=$to_subsample;
    %seen=();
    my $i=0;
    while($n>0){
        if (!defined($tab[$i])){
            if ($eco eq "Other"){
                ## Ok, we take for the whole thing, but we won't show separately
                print "Stopped while still wanting $n for eco $eco\n";
                $n=0;
            }
            else{
                print "Stopped while still wanting $n for eco $eco\n";
                print "Came to the end of tab, and did not get 300 unique (probably the same array multiple times, enough to get us under the 300)\n";
                die("\n");
            }
        }
        else{
            if ($tab[$i]=~/^(.*)_(Ac_\d+)$/){
                my $sample=$1;
                my $array=$2;
                if (!defined($seen{$array})){
                    print "Taking $array - $sample\n";
                    $check{$array}{$sample}=1;
                    $seen{$array}=1;
                    $n--;
                }
                else{
                    print "We already saw array $array\n";
                    if ($eco eq "Human-associated_Respiratory-system"){
                        print "\t but for respiratory system we still take because otherwise we don't have enough\n";
                        die("we should not be there anymore\n");
                    }
                }
                $i++
            }
            else{
                die("pblm with code $tab[$i]\n");
            }
        }
    }
}

print "Reading $in_file_sets and $in_file_sets_bis\n";
my %store;
my $i=0;
open my $tsv,"cat $in_file_sets $in_file_sets_bis |";
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "cluster_id"){next;}
    if ($check{$tab[2]}{$tab[3]}==1){
        $store{$tab[2]}{$tab[3]}{$tab[0]}=$tab[1];
    }
    $i++;
    if ($i % 1000000 == 0){print ".. $i ..\n"};
}
close $tsv;

open my $s1,">",$out_file;
print "Preparing $out_file";
print $s1 "cluster_id\ttotal_coverage\trelative_coverage\trank\trelative_rank\tcrispr_array\ttaxo_lvl\tsra_run\tecosystem\n";
foreach my $array (keys %store){
    foreach my $sample (keys %{$store{$array}}){
        my $eco=$info_sample{$sample}{"eco"};
        print $array."\t".$sample."\t".$eco."\n";
        my %subs=();
        %subs=subsample_set($store{$array}{$sample},$to_subsample_spacers);
        my @t=sort {$subs{$b} <=> $subs{$a} or $a cmp $b} keys %subs;
        my $max=$subs{$t[0]};
        my $rank_max=scalar(@t);
        print $eco."\t".$array."\t".$sample."\t".$max."\t".$rank_max."\n";
        my $rank=0;
        foreach my $spacer (@t){
            $rank++;
            print $s1 $spacer."\t".$subs{$spacer}."\t".($subs{$spacer}/$max)."\t".$rank."\t".($rank/$rank_max)."\t".$array."\t".$info_array{$array}{"taxo_lvl"}."\t".$sample."\t".$eco."\n";
        }
        my $last=$rank+1;
        for (my $rank=$last;$rank<=$to_subsample_spacers;$rank++){
            print $s1 "Fake\t0\t0\t".$rank."\t".($rank/$rank_max)."\t".$array."\t".$info_array{$array}{"taxo_lvl"}."\t".$sample."\t".$eco."\n";
        }
    }
}
close $s1;

sub subsample_set(){
    my $hash=$_[0];
    my $total=$_[1];
    my @tmp=();
    foreach my $spacer (keys %{$hash}){
        for (my $i=0;$i<=$$hash{$spacer};$i++){
            push(@tmp,$spacer);
        }
    }
    &fisher_yates_shuffle(\@tmp);
    my %count;
    for (my $i=0;$i<$total;$i++){
        $count{$tmp[$i]}++;
    }
    return %count;
}

sub fisher_yates_shuffle(){
    my $deck = shift;  # $deck is a reference to an array
    return unless @$deck; # must not be empty!
    my $i = @$deck;
    while (--$i) {
            my $j = int rand ($i+1);
            @$deck[$i,$j] = @$deck[$j,$i];
    }
}

sub guess_level(){
    my $taxo=$_[0];
    my @t=split(";",$taxo);
    my $max_lvl="NA";
    foreach my $cell (@t){
        my @t2=split("__",$cell);
        if ($t2[1] ne "unclassified"){
            if (!defined($translate_level{$t2[0]})){die("Pblm with cell $t2[0] // $t2[1] // $cell // $taxo\n");}
            $max_lvl=$translate_level{$t2[0]};
        }
    }
    return $max_lvl;
}