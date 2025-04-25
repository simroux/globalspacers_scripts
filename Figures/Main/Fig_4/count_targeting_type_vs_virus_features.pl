#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h); 
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to get counts for the panel showing different targeting levels and types
# Arguments :
# run
";
    die "\n";
}

my $info_virus_file="../../../Data/Spacer_db/IMGVR_sequence_information_Oct17.tsv";
my $mismatch_profile_file="../../../Analyses/Target_IMGVR_IMGPR/Beyond_near_exact/Virus_to_array_hits_profile-medium_only_with_cc_with_virusinfo.tsv";
my $cover_file_base="../../../Analyses/Target_IMGVR_IMGPR/Target_coverage/uvig_hq_coverage_by_spacer-with_sp_info.tsv";

my $out_file="Final_counts_for_type_vs_level-and-dgr.tsv";

my %check_virus;
print "Reading $info_virus_file\n";
open my $tsv,"<",$info_virus_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "uvig"){next;}
    if ($tab[5] eq "High-quality" || $tab[5] eq "Reference"){
        if ($tab[10] eq "other" || $tab[10] eq "unknown"){next;}
        $check_virus{$tab[0]}=1;
        if ($tab[15] eq ""){$tab[15]="no";}
        $info_virus{$tab[0]}{"dgr"}=$tab[15];
    }
}
close $tsv;

my %store;

print "Reading $mismatch_profile_file\n";
open my $tsv,"<",$mismatch_profile_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "uvig"){next;}
    if (!defined($check_virus{$tab[0]})){next;}
    if ($tab[9] eq "NA"){$tab[9]="intermediary";}
    $store{$tab[0]}{$tab[1]}{"mismatch"}=$tab[9];
}
close $tsv;


my %count;
my %check;
print "Reading $cover_file_base\n";
open my $tsv,"<",$cover_file_base;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "uvig"){next;}
    if (!defined($check_virus{$tab[0]})){next;}
    my $type=$tab[13];
    if (!defined($check{"type"}{$type})){$check{"type"}{$type}=1;}
    my $race="NA";
    if (defined($store{$tab[0]}{$tab[1]}{"mismatch"})){$race=$store{$tab[0]}{$tab[1]}{"mismatch"};}
    if (!defined($check{"race"}{$race})){$check{"race"}{$race}=1;}
    ## Count the different category
    $count{"all"}{"all"}{"type"}{$type}++;
    $count{"all"}{"all"}{"race"}{$race}++;
    ## also add counts by dgr yes/no
    $count{"dgr"}{$info_virus{$tab[0]}{"dgr"}}{"type"}{$type}++;
    $count{"dgr"}{$info_virus{$tab[0]}{"dgr"}}{"race"}{$race}++;
}
close $tsv;

my @list_types=keys %{$check{"type"}};
my @list_races=keys %{$check{"race"}};


open my $s1,">",$out_file;
print $s1 "category\tvalue\tmetric\ttype\tcount\n";
foreach my $cat (sort keys %count){
    foreach my $value (sort keys %{$count{$cat}}){
        foreach my $type (@list_types){
            if (!defined($count{$cat}{$value}{"type"}{$type})){$count{$cat}{$value}{"type"}{$type}=0;}
            print $s1 $cat."\t".$value."\ttype\t".$type."\t".$count{$cat}{$value}{"type"}{$type}."\n";
        }
        foreach my $race (@list_races){
            if (!defined($count{$cat}{$value}{"race"}{$race})){$count{$cat}{$value}{"race"}{$race}=0;}
            print $s1 $cat."\t".$value."\trace\t".$race."\t".$count{$cat}{$value}{"race"}{$race}."\n";
        }
    }
}
close $s1;

