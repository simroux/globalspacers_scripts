#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
my %translate_level=("d"=>"domain","p"=>"phylum","c"=>"class","o"=>"order","f"=>"family","g"=>"genus","s"=>"species");

GetOptions ('help' => \$h, 'h' => \$h); 
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
	print "# Script to prepare a file with all the info of beta diversity, and prepare a sub
# Arguments :
# run
";
	die "\n";
}

## For the plot we will take 1000 spacer sets (all from different arrays) with at least 100 spacers and a max cover of at least 10 across each ecosystem
my $min_spacer=10; ## min spacer per set to consider
my $min_maxcover=20; ## min max cover spacer to consider per set

# ## Min 10 samples, and we sub-sample exactly 10
my $min_samples=10; ## min sample per array to consider - this time to 10, for the subsample, to get a nice curve
my $n_sample_to_subsample=10; ## we will take 10 to be systematic
my $n_array_to_subsample=100;

my $out_file="spacer_sets_betadiv_forplot_subsamples.tsv";
my $out_file_stats="spacer_sets_betadiv_forplot_subsamples-summarized.tsv";

my $file_array_info="Data/Spacer_db/Array_info_filtered_for_db-Oct24-25.tsv";
my $file_sample_info="Data/Spacer_db/Runs_to_ecosystem_and_sequencing_and_study_for_db-Jul28-24.tsv";

my $in_spacer_base="Analyses/Spacer_database/spacer_clusters_and_metadata_for_alphadiv.tsv";
my $in_file_set_size="Analyses/Spacer_database/spacer_sets_alphadiv.tsv";
my $in_file_sets="spacer_clusters_and_metadata.tsv";

my %info_array;
open my $tsv,"<",$file_array_info;
print ".. reading $file_array_info ..\n";
while(<$tsv>){
      chomp($_);
      my @tab=split("\t",$_);
      if ($tab[0] eq "repeat_cluster"){next;}
      $info_array{$tab[0]}{"type"}=$tab[1];
      $info_array{$tab[0]}{"taxo_source"}=$tab[2];
      $info_array{$tab[0]}{"taxo_genus"}=$tab[6];
      if ($tab[3] eq "NA"){$info_array{$tab[0]}{"taxo_lvl"}="None";}
      else{$info_array{$tab[0]}{"taxo_lvl"}=&guess_level($tab[3]);}
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

my %to_subsample;
my %info_set;
open my $tsv,"<",$in_file_set_size;
print ".. reading $in_file_set_size ..\n";
my $i=0;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "array"){next;}
    if ($tab[2]>=$min_spacer && $tab[4]>=$min_maxcover){
        $to_subsample{"all"}{$tab[0]}{$tab[1]}=$tab[2];
        $info_set{$tab[0]}{$tab[1]}{"max_cover"}=$tab[4];
        if (defined($info_sample{$tab[1]}{"eco"}) && $info_sample{$tab[1]}{"eco"} ne "Unknown"){
            $to_subsample{$info_sample{$tab[1]}{"eco"}}{$tab[0]}{$tab[1]}=$tab[2];
        }
    }
    $i++;
    if ($i % 1000000 == 0){print " .. $i\n";}
}
close $tsv;

## Make the selection of combinations of array and samples for each array
my %check;
my %seen;
my @tab=();
foreach my $eco (sort keys %to_subsample){
    print "Selecting for $eco\n";
    @tab=();
    foreach my $array (keys %{$to_subsample{$eco}}){
        my $test=scalar(keys %{$to_subsample{$eco}{$array}});
        if ($test>=$min_samples){
            push(@tab,$array);
        }
    }
    if (scalar(@tab)<$n_array_to_subsample){
        print "Pblm, for $eco, we don't have enough data -- ".scalar(@tab)." vs $n_array_to_subsample\n";
        if ($eco eq "Other"){
            print "That's fine for Other\n";
        }
        else{
            die("Pblm, for $eco, we don't have enough data -- ".scalar(@tab)." vs $n_array_to_subsample\n");
        }
    }
    &fisher_yates_shuffle(\@tab);
    my $n=$n_array_to_subsample;
    %seen=();
    my $i=0;
    while($n>0){
        if (!defined($tab[$i])){
            if ($eco eq "Other"){
                ## Ok, we take for the whole thing, but we won't show separately
                print "Stopped while still wanting $n\n";
                $n=0;
            }
            else{
                print "Stopped while still wanting $n\n";
                die("Came to the end of tab, and did not get 1000 unique (probably the same array multiple times, enough to get us under the 1000)\n");
            }
        }
        else{
            my $array=$tab[$i];
            print "selected $array with eco $eco\n";
            print "now picking $n_sample_to_subsample samples from here\n";
            my @pot=keys %{$to_subsample{$eco}{$array}};
            if (scalar(@pot)>$n_sample_to_subsample){
                &fisher_yates_shuffle(\@pot);
                for (my $i=0;$i<$n_sample_to_subsample;$i++){
                    if (defined($pot[$i])){
                        $check{$array}{$eco}{$pot[$i]}=1;
                        if (!defined($info_set{$array}{$pot[$i]}{"max_cover"})){print "No max cover for $array - $pot[$i] ? ($eco)\n"; <STDIN>;}
                    }
                }
            }
            else{
                foreach my $sample (@pot){
                    $check{$array}{$eco}{$sample}=1;
                }
            }
            $n--;
        }
        $i++;
    }
}

print "Reading $in_spacer_base to get which spacers are common vs rare\n";
my %info_spacer;
open my $tsv,"<",$in_spacer_base;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "array"){next;}
    if (defined($check{$tab[2]})){
        foreach my $eco (sort keys %{$check{$tab[2]}}){
            if (defined($check{$tab[2]}{$eco}{$tab[3]}) && $check{$tab[2]}{$eco}{$tab[3]}==1){
                # print $tab[2]."\t".$tab[3]."\t".$info_set{$tab[2]}{$tab[3]}{"max_cover"}."\n";
                my $ratio=$tab[1]/$info_set{$tab[2]}{$tab[3]}{"max_cover"};
                my $type="unknown";
                if ($ratio>=0.5){$type="common";}
                elsif($ratio>=0.05){$type="intermediary";}
                elsif($tab[1]==1){$type="singleton";}
                else{$type="rare";}
                $info_spacer{$tab[0]}{$type}++;
            }
        }
    }
    $i++;
    if ($i % 1000000 == 0){print $i." ... \n";}
}
close $tsv;



my %tmp;
my %stats;
print ".. reading $in_file_sets ..\n";
open my $s1,">",$out_file;
print $s1 "array\tecosystem\tspacer\tspacer_type\tn_sample\tratio_sample\tn_common\tn_singleton\n";
open my $tsv,"<",$in_file_sets;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "cluster_id"){next;}
    if (defined($check{$tab[1]})){
        my $array=$tab[1];
        $tab[4]=~s/\[//g;
        $tab[4]=~s/\]//g;
        $tab[4]=~s/\s//g;
        my @t=split(",",$tab[4]);
        foreach my $eco (keys %{$check{$array}}){
            my $n=0;
            my $max=scalar(keys %{$check{$array}{$eco}}); ## should always be 10
            foreach my $sample (@t){
                if (defined($check{$array}{$eco}{$sample})){
                    $n++;
                }
            }
            if ($n>0){
                my $spacer_type="NA";
                if (!defined($info_spacer{$tab[0]}{"common"})){$info_spacer{$tab[0]}{"common"}=0;}
                if (!defined($info_spacer{$tab[0]}{"singleton"})){$info_spacer{$tab[0]}{"singleton"}=0;}
                if ($info_spacer{$tab[0]}{"common"}>0){$spacer_type="common";}
                elsif ($info_spacer{$tab[0]}{"singleton"}>0){$spacer_type="singleton";}
                else{$spacer_type="unknown";}
                print $s1 $array."\t".$eco."\t".$tab[0]."\t".$spacer_type."\t".$n."\t".($n/$max)."\t".$info_spacer{$tab[0]}{"common"}."\t".$info_spacer{$tab[0]}{"singleton"}."\n";
                my $ratio=sprintf("%.02f",($n/$max));
                $stats{$array}{$eco}{"total"}++;
                $stats{$array}{$eco}{"ratios"}{$ratio}++;
                $stats{$array}{$eco}{"n_sample_all"}{$n}++;
                if ($spacer_type ne "NA"){
                    $stats{$array}{$eco}{"n_sample_".$spacer_type}{$n}++;
                    $stats{$array}{$eco}{"total_".$spacer_type}++;
                }
                # print $array."\t".$eco."\t".$tab[0]."\t".$n."\t".($n/$max)."\n";
            }
        }
    }
}
close $s1;

open my $s2,">",$out_file_stats;
print $s2 "array\tecosystem\tn_samples\tcount_all\tpercentage_all\tcount_common\tpercentage_common\tcount_sgton\tpercentage_sgton\n";
foreach my $array (sort keys %stats){
    foreach my $eco (sort keys %{$stats{$array}}){
        foreach my $n_sample (sort {$a <=> $b} keys %{$stats{$array}{$eco}{"n_sample_all"}}){
            print $s2 $array."\t".$eco."\t".$n_sample."\t".
            $stats{$array}{$eco}{"n_sample_all"}{$n_sample}."\t".($stats{$array}{$eco}{"n_sample_all"}{$n_sample}/$stats{$array}{$eco}{"total"})."\t".
            $stats{$array}{$eco}{"n_sample_common"}{$n_sample}."\t".($stats{$array}{$eco}{"n_sample_common"}{$n_sample}/$stats{$array}{$eco}{"total_common"})."\t".
            $stats{$array}{$eco}{"n_sample_singleton"}{$n_sample}."\t".($stats{$array}{$eco}{"n_sample_singleton"}{$n_sample}/$stats{$array}{$eco}{"total_singleton"})."\n";
        }
    }
}
close $s2;

sub fisher_yates_shuffle {
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
            $max_lvl=$t2[0];
        }
    }
    return $max_lvl;
}