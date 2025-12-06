#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
my %translate_level=("d"=>"domain","p"=>"phylum","c"=>"class","o"=>"order","f"=>"family","g"=>"genus","s"=>"species");
GetOptions ('help' => \$h, 'h' => \$h); 
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
	print "# Script to prepare a file for some distribution of spacer coverage for a subsample of spacer sets
# Arguments :
# run
";
	die "\n";
}

my $file_array_info="Data/Spacer_db/Array_info_filtered_for_db-Oct24-25.tsv";
my $file_sample_info="Data/Spacer_db/Runs_to_ecosystem_and_sequencing_and_study_for_db-Jul28-24.tsv";

my $in_alpha_base="Analyses/Spacer_database/spacer_sets_alphadiv.tsv";
my $out_file="Spacer_set_size_info.tsv";

my %info_array;
open my $tsv,"<",$file_array_info;
print "Reading $file_array_info\n";
while(<$tsv>){
      chomp($_);
      my @tab=split("\t",$_);
      if ($tab[0] eq "repeat_cluster"){next;}
      my @t=split("-",$tab[1]);
      $info_array{$tab[0]}{"type"}=$t[0];
      $info_array{$tab[0]}{"taxo_source"}=$tab[2];
      $info_array{$tab[0]}{"taxo_genus"}=$tab[6];
      my @t=split(";",$tab[6]);
      $info_array{$tab[0]}{"simple_genus"}=$t[0].";".$t[1].";".$t[5];
      $info_array{$tab[0]}{"simple_phylum"}=$t[0].";".$t[1];
      if ($tab[6] eq "NA"){
        $info_array{$tab[0]}{"simple_genus"}="NA";
        $info_array{$tab[0]}{"simple_phylum"}="NA";      
      }
}
close $tsv;

my %info_sample;
open my $tsv,"<",$file_sample_info;
print "Reading $file_sample_info\n";
while(<$tsv>){
      chomp($_);
      my @tab=split("\t",$_);
      if ($tab[0] eq "library"){next;}
      $info_sample{$tab[0]}{"eco"}=$tab[3];
      $info_sample{$tab[1]}{"eco"}=$tab[3];
}
close $tsv;

my $i=0;

open my $tsv,"<",$in_alpha_base;
print "Reading $in_alpha_base\n";
open my $s1,">",$out_file;
print $s1 "array\tsample\ttype\tecosystem\tgenus\tspacer_set_size\tmax_cover\tn_singleton\n";
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "array"){next;}
    print $s1 $tab[0]."\t".$tab[1]."\t".$info_array{$tab[0]}{"type"}."\t".$info_sample{$tab[1]}{"eco"}."\t".$info_array{$tab[0]}{"taxo_genus"}.
    "\t".$tab[2]."\t".$tab[4]."\t".$tab[6]."\n";
}
close $tsv;
close $s1;