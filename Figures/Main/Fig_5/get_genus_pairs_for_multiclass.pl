#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to get pairs of co-targeting genera
# Arguments :
# run
";
    die "\n";
}

my $in_list="Virus_to_targeting-range.tsv";
my $in_genus="Target_to_taxon.tsv";

my $genus_list="out_genus_list.txt";

my %check;
open my $tsv,"<",$in_list;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[6] eq "multi_class_balanced"){
        $check{$tab[0]}=1;
    }
}
close $tsv;

my %store;
open my $tsv,"<",$in_genus;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "target_id"){next;}
    if ($check{$tab[0]}==1){
        my @t=split(";",$tab[1]);
        if ($t[5] ne "g__unclassified"){
            $store{$tab[0]}{$tab[1]}=1;
        }
    }
}
close $tsv;

my %count;
foreach my $uvig (sort keys %store){
    my @t=sort {$a cmp $b} keys %{$store{$uvig}};
    for (my $i=0;$i<$#t;$i++){
        my @t1=split(";",$t[$i]);
        my $class_1=$t1[0].";".$t1[1].";".$t1[2];
        for (my $k=$i+1;$k<=$#t;$k++){
            my @t2=split(";",$t[$k]);
            my $class_2=$t2[0].";".$t2[1].";".$t2[2];
            if ($class_1 ne $class_2){
                $count{$t[$i]."|".$t[$k]}++;
            }
        }
    }
}
open my $s1,">",$genus_list;
foreach my $pair (sort {$count{$a} <=> $count{$b} or $a cmp $b} keys %count){
    print $pair."\t".$count{$pair}."\n";
    print $s1 $pair."\t".$count{$pair}."\n";
}
close $s1;