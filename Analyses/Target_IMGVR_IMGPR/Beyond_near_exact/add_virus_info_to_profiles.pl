#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h); 
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to add virus information to the profile table. We need to do this after everything, because it needs information from the Known_hosts folder, which itself needs information from Virus_to_repeat_hits_profile-medium_with_cat.tsv
# Arguments :
# run
";
    die "\n";
}


my $profile_file="Virus_to_repeat_hits_profile-medium_with_cat.tsv";
my $info_virus_file="../../../Data/Spacer_db/IMGVR_sequence_information_Oct17.tsv";
my $known_hosts="../Known_hosts/Consistency_known_host_with_hits.tsv";
my $out_file_step3="Virus_to_repeat_hits_profile-medium_only_with_cat_with_virusinfo.tsv";


my %info_virus;
print "Reading $info_virus_file\n";
open my $tsv,"<",$info_virus_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "uvig"){next;}
    if ($tab[5] eq "High-quality" || $tab[5] eq "Reference"){
        if ($tab[10] eq "other" || $tab[10] eq "unknown"){next;}
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
}
close $tsv;

my %store_known;
print "Reading $known_hosts\n";
open my $tsv,"<",$known_hosts;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "uvig"){next;}
    if ($tab[4] ne ""){
        $store_known{$tab[0]}{$tab[1]}=$tab[4];
    }
}
close $tsv;



open my $tsv,"<",$profile_file;
open my $s1,">",$out_file_step3;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "uvig"){
        print $s1 $_."\tecosystem\tlifestyle\tfibers\tacr\tdgr\tknown_host\n";
    }
    my $known="NA";
    if (defined($store_known{$tab[0]}{$tab[1]})){
        $known=$store_known{$tab[0]}{$tab[1]};
    }
    my $line=$_."\t".$info_virus{$tab[0]}{"ecosystem"}."\t".$info_virus{$tab[0]}{"lifestyle"}."\t".$info_virus{$tab[0]}{"fibers"}."\t".$info_virus{$tab[0]}{"acr"}."\t".$info_virus{$tab[0]}{"dgr"}."\t".$known;
    print $s1 $line."\n";
}
close $s1;
close $tsv;

