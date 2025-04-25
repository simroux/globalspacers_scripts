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

my $known_host_file="../../../Data/Additional_data/Phages_with_known_hosts.tsv";
my $hit_file="Phages_with_known_hosts-all_hit_counts_to_repeats.tsv";
my $profile_file="../Beyond_near_exact/Virus_to_array_hits_profile-medium_only_with_cc_with_virusinfo.tsv";
my $out_file="Consistency_known_host_with_hits.tsv";

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
        $store_profiles{$tab[0]}{$tab[1]}=$tab[9];
    }
}
close $tsv;

open my $s1,">",$out_file;
print $s1 "uvig\trepeat\tn_hits\tprofile\tresult\tref_host_source\tref_host\thit_host_source\thit_host\n";
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
        if ($result ne "NA"){
            if (!defined($store_profiles{$tab[0]}{$tab[2]})){$store_profiles{$tab[0]}{$tab[2]}="NA";}
            print $s1 $tab[0]."\t".$tab[2]."\t".$tab[1]."\t".$store_profiles{$tab[0]}{$tab[2]}."\t".$result."\t".
            $known_host{$tab[0]}{"source"}."\t".$known_host{$tab[0]}{"host"}."\t".$tab[4]."\t".$taxo_hit."\n";
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