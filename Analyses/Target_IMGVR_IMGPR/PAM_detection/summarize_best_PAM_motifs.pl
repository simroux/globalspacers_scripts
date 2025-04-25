#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use File::Basename;
use Cwd;
my $h=0;
my $in_file_motifs="";
my $in_file_arrays="";
my $out_file="";
my $out_stat_file="";
my $plasmid;
GetOptions ('help' => \$h, 'h' => \$h, 'i=s'=>\$in_file_motifs, 'd=s'=>\$in_file_arrays, 'o=s'=>\$out_file, 's=s'=>\$out_stat_file);
if ($h==1 || $in_file_motifs eq "" || $in_file_arrays eq ""){ # If asked for help or did not set up any argument
	print "# Script to summarize the best PAM motifs for each repeat (if possible)
# Arguments :
# -i: input file of the de novo detection of putative motifs
# -d: input file with information about each repeat
# -o: output file - best motif by repeat
# -s: output file - statistics about motif detection
";
	die "\n";
}

if (-e $out_file || -e $out_stat_file){
    die("$out_file or $out_stat_file already exists, I refuse to do anything\n");
}


my %check_low_targets;
my %store;
open my $tsv,"<",$in_file_motifs;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "Array"){next;}
    my $array=$tab[0];
    if ($tab[2] eq "less_than_10_observations"){
        $check_low_targets{$array}=1;
        next;
    }
    $check_low_targets{$array}=2;
    my $location=$tab[1];
    my $motif=$tab[4];
    my $n_obs=$tab[5]/($tab[6]/100); ## total observations, i.e. number of sites used to define the motif divided by percentage used
    my $p_obs=$tab[6];
    my $count_n = $motif =~ tr/N//;
    if ($count_n<=4){ ## At least 1 bp motif
        if ($p_obs >= 50){ ## Built from at least 50% of the sequences
            $store{$array}{$location}{"n_obs"}=$n_obs;
            $store{$array}{$location}{"p_obs"}=$p_obs;
            $store{$array}{$location}{"ns"}=$count_n;
            $store{$array}{$location}{"motif"}=$motif;
            if ($count_n==4){$store{$array}{$location}{"mono"}=1;} ## Some custom flag to sort the mononucleotide motifs better
            else{$store{$array}{$location}{"mono"}=2;}
        }
    }
}
close $tsv;

my %stats;
open my $s1,">",$out_file;
open my $tsv,"<",$in_file_arrays;
print $s1 "Array\tType\tMotif\tLocation\tPercentage\tTotal_observations\n";
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "repeat_cluster"){next;}
    my $array=$tab[0];
    my $type=$tab[1];
    if (defined($store{$array})){
        my @t=sort {$store{$array}{$b}{"mono"} <=> $store{$array}{$a}{"mono"} or $store{$array}{$b}{"p_obs"} <=> $store{$array}{$a}{"p_obs"} or $store{$array}{$a}{"ns"} <=> $store{$array}{$b}{"ns"} or $a <=> $b} keys %{$store{$array}};
        ## Sorting by single first should always prioritize the non-mononucleotide, i.e. di-, tri-, etc, then we take the most common, and if tied then we take the longest
        my $best=$t[0];
        print $s1 $array."\t".$type."\t".$store{$array}{$best}{"motif"}."\t".$best."\t".$store{$array}{$best}{"p_obs"}."\t".$store{$array}{$best}{"n_obs"}."\n";
        print $array."\t".$type."\t".$store{$array}{$best}{"motif"}."\t".$best."\t".$store{$array}{$best}{"p_obs"}."\t".$store{$array}{$best}{"n_obs"}."\n";
        my $code=$best."_".$store{$array}{$best}{"motif"};
        $stats{$type}{$code}++;
    }
    elsif($check_low_targets{$array}==1){
        print $s1 $array."\t".$type."\tLow_target_n\tNA\tNA\tNA\n";
    }
    elsif($check_low_targets{$array}==2){
        print $s1 $array."\t".$type."\tNA\tNA\tNA\tNA\n";
        $stats{$type}{"NA"}++; ## It was in the list that we tried to get a motif for, and not flagged as low_target_n, so genuinely we did not find anything
    }
    else{
        print $s1 $array."\t".$type."\tNo_target_observed\tNA\tNA\tNA\n";
    }
}
close $tsv;
close $s1;


print "##### Summary by type ####\n";
my $top=0;
open my $s2,">",$out_stat_file;
print $s2 "Array_type\tMotif\tCount\n";
foreach my $type (sort keys %stats){
    $top=0;
    foreach my $motif (sort {$stats{$type}{$b} <=> $stats{$type}{$a} or $a cmp $b} keys %{$stats{$type}}){
        print $s2 $type."\t".$motif."\t".$stats{$type}{$motif}."\n";
        $top++;
        if ($top<=10){print $type."\t".$motif."\t".$stats{$type}{$motif}."\n";}
    }
}
close $s2;