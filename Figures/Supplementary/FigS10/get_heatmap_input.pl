#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h); 
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
	print "# Script to get the input for the heatmap matching target and spacer ecosystem
# Arguments :
# run
";
	die "\n";
}

my $in_vr="link_by_ecosystem.tsv";
my $out_file="vr_ecosystem_heatmap.tsv";

my $in_pr="link_by_ecosystem-pr.tsv";
my $out_file_pr="pr_ecosystem_heatmap.tsv";

my %tmp;
my $total=0;
open my $tsv,"<",$in_vr;
while(<$tsv>){
     chomp($_);
     my @tab=split("\t",$_);
     if ($tab[0] eq "virus_ecosystem"){next;}
     $tmp{"virus_total"}{$tab[0]}+=$tab[2];
     $tmp{"spacer_total"}{$tab[1]}+=$tab[2];
     $total+=$tab[2];
}
close $tsv;

open my $s1,">",$out_file;
print $s1 "virus_ecosystem\tspacer_ecosystem\tn_observation\tp_observation\n";
open my $tsv,"<",$in_vr;
while(<$tsv>){
     chomp($_);
     my @tab=split("\t",$_);
     if ($tab[0] eq "virus_ecosystem"){next;}
     if ($tab[0] eq "Other" || $tab[1] eq "Other"){next;}
     my $p_obs=$tab[2]/$tmp{"virus_total"}{$tab[0]}*100;
     print $s1 $_."\t".$p_obs."\n";
}
close $tsv;
close $s1;

%tmp=();
$total=0;
open my $tsv,"<",$in_pr;
while(<$tsv>){
     chomp($_);
     my @tab=split("\t",$_);
     if ($tab[0] eq "plasmid_ecosystem"){next;}
     $tmp{"plasmid_total"}{$tab[0]}+=$tab[2];
     $tmp{"spacer_total"}{$tab[1]}+=$tab[2];
     $total+=$tab[2];
}
close $tsv;

open my $s1,">",$out_file_pr;
print $s1 "plasmid_ecosystem\tspacer_ecosystem\tn_observation\tp_observation\n";
open my $tsv,"<",$in_pr;
while(<$tsv>){
     chomp($_);
     my @tab=split("\t",$_);
     if ($tab[0] eq "plasmid_ecosystem"){next;}
     if ($tab[0] eq "Other" || $tab[1] eq "Other"){next;}
     my $p_obs=$tab[2]/$tmp{"plasmid_total"}{$tab[0]}*100;
     print $s1 $_."\t".$p_obs."\n";
}
close $tsv;
close $s1;
