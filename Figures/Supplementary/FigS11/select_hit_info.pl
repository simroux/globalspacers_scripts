#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h);
if ($h==1 || $ARGV[0] eq ""){
	print "# Script to filter the big hit table to get a more reasonable input file for supp figure only focusing on atypical phages. Could be done in R as well, doesn't really matter.
# Arguments :
# run
";
	die "\n";
}
my $selec_file="selected_viruses.tsv";
my $in_file="../../Main/Fig_3/fig_hits_input_vr.tsv";
my $out_file="fig_hits_input_vr-atypical_phages.tsv";

my %check;
open my $tsv,"<",$selec_file;
while(<$tsv>){
     chomp($_);
     my @tab=split("\t",$_);
     if ($tab[0] eq "uvig"){next;}
     $check{$tab[0]}=1;
}
close $tsv;


open my $s1,">",$out_file;
open my $tsv,"<",$in_file;
while(<$tsv>){
     chomp($_);
     my @tab=split("\t",$_);
     if ($tab[0] eq "uvig"){
          print $s1 $_."\tcustom_tax\n";
     }
     else{
          my $full_taxo=$tab[6];
          my $class=$tab[7];
          my $type=$tab[10];
          if ($check{$tab[0]}==1){
               my $custom_tax=$class;
               if ($custom_tax eq ";;;"){
                    my @t=split(";",$full_taxo);
                    $custom_tax=";;;;;".$t[5];
                    print $full_taxo."\t".$class."\t".$type."\t".$custom_tax."\n";
               }
               print $s1 $_."\t".$custom_tax."\n";
          }
     }
}
close $tsv;
close $s1;
