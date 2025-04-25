#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h); 
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to get the different features of viruses targeted by multiple taxa
# Arguments :
# run
";
    die "\n";
}

my $info_virus_file="../../../Data/Spacer_db/IMGVR_sequence_information_Oct17.tsv";
my $in_file_multiclass="List_of_interest.tsv";
my $out_file="Multitaxa_features.tsv";

my %info_virus;
print "Reading $info_virus_file ..\n";
open my $tsv,"<",$info_virus_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "uvig"){next;}
    if ($tab[10] eq "phage" || $tab[10] eq "archaea"){
        $info_virus{$tab[0]}{"eco"}=$tab[1];
        $info_virus{$tab[0]}{"votu"}=$tab[2];
        $info_virus{$tab[0]}{"length"}=$tab[3];
        $info_virus{$tab[0]}{"quality"}=$tab[5];
        $info_virus{$tab[0]}{"taxo"}=$tab[6];
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

my %count;
print "Reading $in_file_multiclass ..\n";
open my $s1,">",$out_file;
print $s1 "uvig\tcategory\tquality\tlifestyle\tfibers\tacr\tdgr\tecosystem\tvotu\n";
open my $tsv,"<",$in_file_multiclass;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if (defined($info_virus{$tab[0]})){
        print $s1 $tab[0]."\t".$tab[1]."\t".$info_virus{$tab[0]}{"quality"}."\t".$info_virus{$tab[0]}{"lifestyle"}."\t".$info_virus{$tab[0]}{"fibers"}."\t".
        $info_virus{$tab[0]}{"acr"}."\t".$info_virus{$tab[0]}{"dgr"}."\t".$info_virus{$tab[0]}{"eco"}."\t".$info_virus{$tab[0]}{"votu"}."\n";
    }
}
close $tsv;
close $s1;