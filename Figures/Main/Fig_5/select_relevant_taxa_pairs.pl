#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to select pairs we think may be relevant to report for the 'multitaxa targeting' search
# Arguments :
# run
";
    die "\n";
}

## What we want as information in the final table is just a categorization of viruses 
# single class
# multi class - one dominant
# multi class - balanced
my $virus_info="../../../Data/Spacer_db/IMGVR_sequence_information_Oct17.tsv";

my $in_file="Target_to_taxon.tsv";
my $out_file="Virus_to_targeting-range.tsv";

print "#######################\n";
print "Load info about viruses from $virus_info to get only phages / archaeal viruses\n";
my %check_virus;
my %info_virus;
open my $tsv,"<",$virus_info;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "uvig"){next;}
    if ($tab[10] eq "phage" || $tab[10] eq "archaea"){
        $check_virus{$tab[0]}=1;
        $info_virus{$tab[0]}{"quality"}=$tab[5];
        $info_virus{$tab[0]}{"length"}=$tab[3];
        $info_virus{$tab[0]}{"taxo"}=$tab[6];
        $info_virus{$tab[0]}{"lifestyle"}=$tab[12];
        if ($tab[14] eq ""){$tab[14]="no";}
        $info_virus{$tab[0]}{"acr"}=$tab[14];
        if ($tab[15] eq ""){$tab[15]="no";}
        $info_virus{$tab[0]}{"dgr"}=$tab[15];
    }
}
close $tsv;

my %store;
print "Reading $in_file\n";
open my $tsv,"<",$in_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "target_id"){next;}
    if ($check_virus{$tab[0]}==1){
        my @t=split(";",$tab[1]);
        if ($t[2] eq "c__unclassified"){} ## We ignore, no class
        else{
            my $class=$t[0].";".$t[1].";".$t[2];
            $store{$tab[0]}{$class}{"n_clusters"}+=$tab[4];
        }
    }
}
close $tsv;

print "Processing viruses\n";
open my $s1,">",$out_file;
print $s1 "virus\tlist_class\ttotal_hits\tn_class\tbest_class\tbest_ratio\tcategory\tquality\tlength\ttaxo\tlifestyle\tacr\tdgr\n";
foreach my $virus (sort keys %store){
    if ($check_virus{$virus}==1){
        my @t=sort {$store{$virus}{$b}{"n_clusters"} <=> $store{$virus}{$a}{"n_clusters"} or $a cmp $b} keys %{$store{$virus}};
        my $n_class=scalar(@t);
        my $total=0;
        foreach my $class (@t){$total+=$store{$virus}{$class}{"n_clusters"};}
        my $best=$store{$virus}{$t[0]}{"n_clusters"}/$total;
        my $cat="unknown";
        if ($n_class==1){$cat="single_class";}
        elsif($best>=0.8){$cat="dominant_class";}
        else{$cat="multi_class_balanced";}
        print $s1 $virus."\t".join("|",@t)."\t".$total."\t".$n_class."\t".$t[0]."\t".$best."\t".$cat."\t".
        $info_virus{$virus}{"quality"}."\t".$info_virus{$virus}{"length"}."\t".$info_virus{$virus}{"taxo"}."\t".
        $info_virus{$virus}{"lifestyle"}."\t".$info_virus{$virus}{"acr"}."\t".$info_virus{$virus}{"dgr"}."\n";
    }
}
close $s1;
