#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use File::Basename;
use Cwd;
my $h=0;
my $in_minced="";
my $in_mapping="";
my $out_file="";
GetOptions ('help' => \$h, 'h' => \$h, 'i=s'=>\$in_minced, 'm=s'=>\$in_mapping, 'o=s'=>\$out_file);
if ($h==1 || $in_minced eq "" || $in_mapping eq "" || $out_file eq ""){ # If asked for help or did not set up any argument
	print "# Script to get the non-redundant regions to mask because they look like potential CRISPR arrays, so should not be considered for protospacers
# Arguments :
# Arguments :
# -i: input file from minced
# -m: input file from the mapping of repeats 
# -o: output file
";
	die "\n";
}

my %store;
my $i=0;
print "Reading $in_minced  ..\n";
open my $tsv,"<",$in_minced;
while(<$tsv>){
    chomp($_);
    if ($_=~/^#/){next;}
    my @tab=split("\t",$_);
    $i++;
    $store{$tab[0]}{$i}{"start"}=$tab[3];
    $store{$tab[0]}{$i}{"end"}=$tab[4];
}
close $tsv;
print "Reading $in_mapping  ..\n";
open my $tsv,"<",$in_mapping;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    $i++;
    $store{$tab[0]}{$i}{"start"}=$tab[2];
    $store{$tab[0]}{$i}{"end"}=$tab[3];
}
close $tsv;

print "Preparing the output file $out_file\n";
my $start_c=-1;
my $end_c=-1;
open my $s1,">",$out_file;
foreach my $seq (sort keys %store){
    my @list_regions=sort {$store{$seq}{$a}{"start"} <=> $store{$seq}{$b}{"start"} or $store{$seq}{$a}{"end"} <=> $store{$seq}{$b}{"end"} or $a cmp $b} keys %{$store{$seq}};
    print "Processing ".scalar(@list_regions)." for $seq\n";
    for (my $i=0;$i<=$#list_regions;$i++){
        $start_c=$store{$seq}{$list_regions[$i]}{"start"};
        $end_c=$store{$seq}{$list_regions[$i]}{"end"};
        print "\tRegion in $start_c -- $end_c\n";
        if ($i<$#list_regions){
            while($store{$seq}{$list_regions[$i+1]}{"start"} < $end_c){
                $i++;
                if ($store{$seq}{$list_regions[$i]}{"end"}>$end_c){
                    $end_c=$store{$seq}{$list_regions[$i]}{"end"};
                    print "\t\tExtended to $end_c\n";
                }
                if ($i==$#list_regions){last;}
            }
        }
        print $s1 $seq."\t".$start_c."\t".$end_c."\n";
    }
}
close $s1;