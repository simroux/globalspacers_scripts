#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
	print "# Script to refine the link between spacer and ecosystem for complicated cases
# Arguments :
# run
";
	die "\n";
}

my $in_file_from_duckdb="spacer_clusters_to_ecosystem.tsv";
my $out_file_spacer_eco="spacer_to_eco-clean.tsv";

my %store;
my $i=0;
print "Reading $in_file_from_duckdb\n";
open my $tsv,"<",$in_file_from_duckdb;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "cluster_id"){next;}
    $store{$tab[0]}{$tab[1]}=$tab[2];
    $i++;
    if ($i % 1000000 == 0){
        print ".. $i \n";
    }
}
close $tsv;

print "Data loaded, now preparing output file $out_file_spacer_eco\n";
$i=0;
open my $s1,">",$out_file_spacer_eco;
print $s1 "cluster_id\tmaj_eco\n";
foreach my $spacer (keys %store){
    my @t=keys %{$store{$spacer}};
    if (scalar(@t)==1){
        print $s1 $spacer."\t".$t[0]."\n";
    }
    else{
        ## Majority rule
        @t=sort {$store{$spacer}{$b} <=> $store{$spacer}{$a} or $a cmp $b} keys %{$store{$spacer}};
        my $total=0;
        foreach my $eco (@t){$total+=$store{$spacer}{$eco};}
        my $test=$store{$spacer}{$t[0]}/$total;
        if ($test>=0.5){
            print $s1 $spacer."\t".$t[0]."\n";
        }
        else{
            print $s1 $spacer."\tmixed\n";
        }
    }
    if ($i % 1000000 == 0){
        print ".. $i \n";
    }
}
close $s1;
