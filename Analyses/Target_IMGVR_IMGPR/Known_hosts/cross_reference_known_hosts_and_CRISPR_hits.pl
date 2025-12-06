#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h); 
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to check the targeting of phages for which we know the host
# Arguments :
# run
";
    die "\n";
}

my $info_virus_file="../../../Data/Spacer_db/IMGVR_sequence_information_Oct17.tsv";
my $known_host_file="../../../Data/Additional_data/Phages_with_known_hosts.tsv";
my $hit_file="Phages_with_known_hosts-all_hit_counts_to_repeats.tsv";
my $profile_file="../Beyond_near_exact/Virus_to_repeat_hits_profile-medium_with_cat.tsv";
my $out_file="Consistency_known_host_with_hits.tsv";

my %info_virus;
my %check;
print "Reading $info_virus_file\n";
open my $tsv,"<",$info_virus_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "uvig"){next;}
    $info_virus{$tab[0]}{"quality"}=$tab[5];
    if ($tab[10] eq "other" || $tab[10] eq "unknown"){next;}
    $check{$tab[0]}=1;
    $info_virus{$tab[0]}{"ecosystem"}=$tab[1];
    if ($tab[12] eq ""){$tab[12]="virulent_unknown";}
    $info_virus{$tab[0]}{"lifestyle"}=$tab[12];
    if ($tab[13] eq ""){$tab[13]="no";}
    $info_virus{$tab[0]}{"fibers"}=$tab[13];
    if ($tab[14] eq ""){$tab[14]="no";}
    $info_virus{$tab[0]}{"acr"}=$tab[14];
    if ($tab[15] eq ""){$tab[15]="no";}
    $info_virus{$tab[0]}{"dgr"}=$tab[15];
}
close $tsv;


my %known_host;
print "reading $known_host_file\n";
open my $tsv,"<",$known_host_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "uvig_id"){next;}
    $known_host{$tab[0]}{"host"}=$tab[1];
    $known_host{$tab[0]}{"source"}=$tab[2];
}
close $tsv;


print "reading $profile_file\n";
my %store_profiles;
open my $tsv,"<",$profile_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "uvig"){next;}
    if (defined($known_host{$tab[0]})){
        $store_profiles{$tab[0]}{$tab[1]}{"profile"}=$tab[8];
        $store_profiles{$tab[0]}{$tab[1]}{"n_hit_0"}=$tab[2];
        $store_profiles{$tab[0]}{$tab[1]}{"n_hit_1"}=$tab[3];
        $store_profiles{$tab[0]}{$tab[1]}{"n_hit_2"}=$tab[4];
        $store_profiles{$tab[0]}{$tab[1]}{"n_hit_3"}=$tab[5];
    }
}
close $tsv;

open my $s1,">",$out_file;
print $s1 "uvig\trepeat\tn_hits\tprofile\tn_hit_0\tn_hit_1\tn_hit_2\tn_hit_3\tresult\tref_host_source\tref_host\thit_host_source\thit_host\tquality\tecosystem\tlifestyle\tfibers\tacr\tdgr\n";
print "reading $hit_file\n";
open my $tsv,"<",$hit_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "target_id"){next;}
    if ($tab[4] eq "No-information_long" || $tab[4] eq "No-information_short"){next;}
    if ($tab[4] eq "Genome_high-confidence" || $tab[4] eq "Genome_medium-confidence"){
        my $taxo_hit=$tab[5];
        my $result=&compare($taxo_hit,$known_host{$tab[0]}{"host"});
        if ($result ne "NA"){ ## NA result only happens if ref or hit is unclassified at the domain level, i.e. we have 0 information at all..
            if (!defined($store_profiles{$tab[0]}{$tab[2]})){
                $store_profiles{$tab[0]}{$tab[2]}{"profile"}="NA";
                $store_profiles{$tab[0]}{$tab[2]}{"n_hit_0"}="NA";
                $store_profiles{$tab[0]}{$tab[2]}{"n_hit_1"}="NA";
                $store_profiles{$tab[0]}{$tab[2]}{"n_hit_2"}="NA";
                $store_profiles{$tab[0]}{$tab[2]}{"n_hit_3"}="NA";
            }
            print $s1 $tab[0]."\t".$tab[2]."\t".$tab[1]."\t".
            $store_profiles{$tab[0]}{$tab[2]}{"profile"}."\t".$store_profiles{$tab[0]}{$tab[2]}{"n_hit_0"}."\t".$store_profiles{$tab[0]}{$tab[2]}{"n_hit_1"}."\t".$store_profiles{$tab[0]}{$tab[2]}{"n_hit_2"}."\t".$store_profiles{$tab[0]}{$tab[2]}{"n_hit_3"}."\t".$result."\t".
            $known_host{$tab[0]}{"source"}."\t".$known_host{$tab[0]}{"host"}."\t".$tab[4]."\t".$taxo_hit."\t".
            $info_virus{$tab[0]}{"quality"}."\t".$info_virus{$tab[0]}{"ecosystem"}."\t".$info_virus{$tab[0]}{"lifestyle"}."\t".
            $info_virus{$tab[0]}{"fibers"}."\t".$info_virus{$tab[0]}{"acr"}."\t".$info_virus{$tab[0]}{"dgr"}."\n";
        }
    }
}
close $tsv;
close $s1;
print "Wrote $out_file\n";


sub compare(){
    my @t=split(";",$_[0]);
    my @ref=split(";",$_[1]);
    my $res="NA";
    for (my $i=0;$i<=5;$i++){
        if ($t[$i]=~/__unclassified$/ || $ref[$i]=~/__unclassified$/){
            $i=10;
            $res.="_before_genus";
        }
        elsif ($t[$i] ne $ref[$i]){
            $i=10;
            $res="inconsistent";
        }
        else{
            $res="consistent";
        }
    }
    if ($t[0] eq "d__unclassified" || $ref[0] eq "d__unclassified"){$res="NA";}
    return $res;
}