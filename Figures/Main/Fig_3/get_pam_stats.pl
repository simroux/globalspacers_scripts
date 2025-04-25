#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use File::Basename;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
	print "# Script to get some numbers about PAM detection across VR and PR 
# Arguments :
# toto
";
	die "\n";
}

my %stats;

my $array_clean="Array_to_PAM_assignment_clean.tsv";

print "Reading $array_clean .. \n";
open my $tsv,"<",$array_clean;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "array"){next;}
    $stats{"array"}{$tab[0]}=&translate($tab[3]);
}
close $tsv;

foreach my $array (sort keys %{$stats{"array"}}){
    if ($stats{"array"}{$array} eq "found"){
        $stats{"array_found"}++;
        $stats{"array_total"}++;
    }
    elsif($stats{"array"}{$array} eq "not_found"){
        $stats{"array_total"}++;
    }
}

print "Arrays: \n";
print "Total: ".$stats{"array_found"}."\n";
print "Wpam: ".$stats{"array_total"}."\n";
print "Ratio: ".sprintf("%.02f",($stats{"array_found"}/$stats{"array_total"}*100))."%\n";


sub translate(){
    my $code=$_[0];
    my $tr="unknown";
    if ($code eq "de_novo" || $code eq "expected" || $code eq "known_type_I" || $code eq "known_type_IV" || $code eq "known_types_I_III_V" || $code eq "known_types_I_II_V_VI" || $code eq "known_types_I_IV" || $code eq "known_types_I_IV_V" || $code eq "known_types_I_V" || $code eq "known_type_V" || $code eq "known_type_VI"){
        $tr="found";
    }
    elsif($code eq "no_motif_detected"){
        $tr="not_found";
    }
    elsif($code eq "no_motif_less_than_50" || $code eq "unexpected_less_than_50" || $code eq "unknown_less_than_50"){
        $tr="too_few";
    }
    return $tr;
}