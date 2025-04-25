#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h); 
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to get the type of high targeting along with basic virus info
# Arguments :
# run
";
    die "\n";
}
my $n_bins=30;

my $info_virus_file="../../../Data/Spacer_db/IMGVR_sequence_information_Oct17.tsv";
my $targeting_file="../../../Analyses/Target_IMGVR_IMGPR/Target_coverage/uvig_hq_coverage_by_spacer-with_sp_info.tsv";
my $out_file="spacer_hit_distribution.tsv";

my %check;
print "Reading $info_virus_file\n";
open my $tsv,"<",$info_virus_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "uvig"){next;}
    if ($tab[5] eq "High-quality" || $tab[5] eq "Reference"){
        if ($tab[10] eq "other" || $tab[10] eq "unknown"){next;} ## Ignore eukaryotic and unknown viruses
        $check{$tab[0]}=1;
    }
}
close $tsv;

my %stat;
my $max=0;
open my $tsv,"<",$targeting_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "uvig"){next;}
    if (!defined($check{$tab[0]})){next;}
    my $summarize=0;
    if ($tab[5]<10){$summarize=$tab[5];}
    elsif($tab[5]<100){
        my $n=int($tab[5]/10)*10;
        my $min=$n;
        my $max=$n+9;
        $summarize=$min."-".$max;
    }
    else{
        $summarize=">= 100"
    }
    $stat{$summarize}++;
}
close $tsv;


open my $s1,">",$out_file;
print $s1 "n_spacers\tn_obs\n";
foreach my $n (sort {$a <=> $b} keys %stat){
    print $s1 $n."\t".$stat{$n}."\n";
}
close $s1;