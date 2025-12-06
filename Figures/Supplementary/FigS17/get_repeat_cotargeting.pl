#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
use Custom::Utils;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h); # , 'i=s'=>\$in_file, 'r=s'=>\$out_dir_root);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to get pairs of repeat co-targeting the same virus
# Arguments :
# run
";
    die "\n";
}


my $in_file_interest="../../Main/Fig_5/List_of_interest.tsv";
my $target_to_repeat_file="../../Main/Fig_5/Target_to_repeat.tsv";
my $repeat_to_lib_file="taxon_to_ecosystem.tsv";

my $out_file="cotargeting_vs_codetection.tsv";
my $out_file_2="cotargeting_vs_codetection_min100.tsv";
my $out_file_3="cotargeting_vs_codetection_min500.tsv";

my %check_virus;
print "Reading $in_file_interest ..\n";
open my $tsv,"<",$in_file_interest;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[1] eq "multi_class_with_high_pos" || $tab[1] eq "one_class_high_other_not_high" || $tab[1] eq "one_class_highpos_other_class_high"){
        $check_virus{$tab[0]}=1;
    }
}
close $tsv;

print "Reading $target_to_repeat_file ..\n";
my %info_repeat;
open my $tsv,"<",$target_to_repeat_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "target_id"){next;}
    if ($check_virus{$tab[0]}==1){
        if ($tab[2] eq "Genome_medium-confidence" || $tab[2] eq "Genome_high-confidence"){
            my @t=split(";",$tab[3]);
            my $class=$t[0].";".$t[1].";".$t[2];
            if ($t[2] ne "c__unclassified"){
                $info_repeat{$tab[1]}{"virus"}{$tab[0]}=1;
                $info_repeat{$tab[1]}{"taxo"}=$class;
            }
        }
    }
}
close $tsv;

print "Loading repeat to ecosystem from $repeat_to_lib_file ..\n";
open my $tsv,"<",$repeat_to_lib_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "crispr_array"){next;}
    if (defined($info_repeat{$tab[0]})){
        $info_repeat{$tab[0]}{"lib"}{$tab[2]}=1;
    }
}
close $tsv;

print "Looking at pairs of repeats from different classes, and checking if they co-target and how much they co-occur\n";
open my $s1,">",$out_file;
print $s1 "repeat_1\trepeat_1_taxo\trepeat_2\trepeat_2_taxo\tcotargeting\ttotal_obs\tcommon_obs\trepeat_1_obs\trepeat_2_obs\n";
open my $s2,">",$out_file_2;
print $s2 "repeat_1\trepeat_1_taxo\trepeat_2\trepeat_2_taxo\tcotargeting\ttotal_obs\tcommon_obs\trepeat_1_obs\trepeat_2_obs\n";
open my $s3,">",$out_file_3;
print $s3 "repeat_1\trepeat_1_taxo\trepeat_2\trepeat_2_taxo\tcotargeting\ttotal_obs\tcommon_obs\trepeat_1_obs\trepeat_2_obs\n";
my @list_repeats=sort keys %info_repeat;
for (my $i=0;$i<$#list_repeats;$i++){
    my $repeat_1=$list_repeats[$i];
    my @libs_1=keys %{$info_repeat{$repeat_1}{"lib"}};
    my $n_sample_1=scalar(@libs_1);
    my @virus_1=keys %{$info_repeat{$repeat_1}{"virus"}};
    for (my $k=$i+1;$k<=$#list_repeats;$k++){
        my $repeat_2=$list_repeats[$k];
        if ($info_repeat{$repeat_1}{"taxo"} ne $info_repeat{$repeat_2}{"taxo"}){
            my $n_sample_2=scalar(keys %{$info_repeat{$repeat_2}{"lib"}});
            my $n_common=0;
            foreach my $lib (@libs_1){
                if (defined($info_repeat{$repeat_2}{"lib"}{$lib})){
                    $n_common++;
                }
            }
            my $n_total=$n_sample_1+$n_sample_2;
            my $cotargeting="no";
            foreach my $virus (@virus_1){
                if (defined($info_repeat{$repeat_2}{"virus"}{$virus})){
                    $cotargeting="yes";
                    last;
                }
            }
            my $line=$repeat_1."\t".$info_repeat{$repeat_1}{"taxo"}."\t".$repeat_2."\t".$info_repeat{$repeat_2}{"taxo"}."\t".
            $cotargeting."\t".$n_total."\t".$n_common."\t".$n_sample_1."\t".$n_sample_2;
            print $s1 $line."\n";
            if ($n_sample_1>=100 && $n_sample_2>=100){
                print $s2 $line."\n";
                if ($n_sample_1>=500 && $n_sample_2>=500){
                    print $s3 $line."\n";
                }
            }
        }
    }
}
close $s1;
close $s2;
close $s3;

