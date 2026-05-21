#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
GetOptions ('help' => \$h, 'h' => \$h);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to get the input needed for Fig. S1 B - Cas gene detection summary
# Arguments :
# none
";
    die "\n";
}

my $in_file="../../../Analyses/Repeat_database/Array_info_filtered_for_db-updated_type_and_cas.tsv";
my $out_file="Cas_gene_overview_sourcedata.R";


my %store;
my %stats;
open my $tsv,"<",$in_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "repeat_cluster"){next;}
    $store{"all"}.="\"".$tab[0]."\",";
    $stats{"t"}++;
    if ($tab[1] ne "Unknown"){
        if ($tab[7] ne "Unknown"){
            $store{"from_repeat"}.="\"".$tab[0]."\",";
            $stats{"r"}++;
        }
        if ($tab[8] ne "Unknown" && $tab[8] ne "NA"){
            $store{"from_cas"}.="\"".$tab[0]."\",";
            $stats{"c"}++;
        }
        if ($tab[8] ne "Unknown" && $tab[8] ne "NA" && $tab[7] ne "Unknown"){
            $stats{"r_c"}++;
        }
    }
}
close $tsv;

## Note: for eulerr, the counts are exclusive, i.e. counts of "A" are count of elements that are A and only A. So counts for repeat and cas are 0, and for totals, we need to do some substractions
open my $s1,">",$out_file;
print $s1 "
cas_type_source <- c(
  \"A\" = ".($stats{"t"}-$stats{"r"}-$stats{"c"}+$stats{"r_c"}).",
  \"B\" = 0,
  \"C\" = 0,
  \"A&B\" = ".($stats{"r"}-$stats{"r_c"}).",
  \"A&C\" = ".($stats{"c"}-$stats{"r_c"}).",
  \"B&C\" = 0,
  \"A&B&C\" = ".$stats{"r_c"}."
)\n";
close $s1;
