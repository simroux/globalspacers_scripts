#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
use Custom::Utils;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h); # , 'i=s'=>\$in_file, 'r=s'=>\$out_dir_root);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to get the type of high targeting by sample
# Arguments :
# run
";
    die "\n";
}

my $info_virus_file="../../../Data/Spacer_db/IMGVR_sequence_information_Oct17.tsv";
my $targeting_file_global="../../../Analyses/Target_IMGVR_IMGPR/Target_coverage/uvig_hq_coverage_by_spacer-with_sp_info.tsv";
my $targeting_file="../../../Analyses/Target_IMGVR_IMGPR/Target_coverage/uvig_hq_coverage_by_spacer-bysample.tsv";
my $out_file="spacer_hit_distribution.tsv";

my %check;
print "Reading $info_virus_file\n";
open my $tsv,"<",$info_virus_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "uvig"){next;}
    if ($tab[5] eq "High-quality" || $tab[5] eq "Reference"){
        if ($tab[10] eq "other" || $tab[10] eq "unknown"){next;} ## Ignore eukaryotic viruses
        $check{$tab[0]}=1;
    }
}
close $tsv;

my %store;
my $i=0;
print "Reading $targeting_file_global ... \n";
open my $tsv,"<",$targeting_file_global;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "uvig"){next;}
    if (!defined($check{$tab[0]})){next;}
    my $type=$tab[13];
    if ($type eq "high"){
        $store{$tab[0]}{$tab[1]}{"global"}{"len"}=$tab[2];
        $store{$tab[0]}{$tab[1]}{"global"}{"type"}=$type;
        $store{$tab[0]}{$tab[1]}{"global"}{"coverage"}=$tab[4];
    }
    $i++;
    if ($i % 1000000 == 0){print ".. $i .. \n";}
}
close $tsv;

print "Reading $targeting_file ... \n";
$i=0;
open my $tsv,"<",$targeting_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "uvig"){next;}
    if (!defined($check{$tab[0]})){next;} ## Ignore the references we do not want
    if (!defined($store{$tab[0]}{$tab[1]})){
        if ($tab[4]<200 || $tab[6]<10){
            ## Expected, we ignore because it's not good in global coverage, and it's not good in local
            next;
        }
        else{
            die("That's weird, how did we end up with $tab[0] and $tab[1] ??\n");
        }
    }
    my $type="NA";
    if ($tab[4]<200 || $tab[6]<10){ ## Covering less than 20bp of the virus or less than 10 spacers
        $type="low";
    }
    elsif($tab[5]>=10 || $tab[4]>=1000){ ## Covering more than 10% (but no more than 50%) or more than 1kb (but less than 5kb)
        $type="high";
    }
    else{
        $type="medium";
    }
    if (!defined($store{$tab[0]}{$tab[1]}{"local"})){
        $store{$tab[0]}{$tab[1]}{"local"}{"type"}=$type;
        $store{$tab[0]}{$tab[1]}{"local"}{"coverage"}=$tab[5];
        $store{$tab[0]}{$tab[1]}{"local"}{"sample"}=$tab[2];
    }
    elsif($tab[5]>=$store{$tab[0]}{$tab[1]}{"local"}{"coverage"}){
        $store{$tab[0]}{$tab[1]}{"local"}{"type"}=$type;
        $store{$tab[0]}{$tab[1]}{"local"}{"coverage"}=$tab[5];
        $store{$tab[0]}{$tab[1]}{"local"}{"sample"}=$tab[2];
    }
    $i++;
    if ($i % 1000000 == 0){
        print ".. $i ..\n";
    }
}
close $tsv;

print "Writing $out_file\n";
open my $s1,">",$out_file;
print $s1 "virus\tlength\tarray\tglobal_type\tglobal_cover\tlocal_type\tlocal_cover\tlocal_sample\n";
foreach my $virus (sort keys %store){
    foreach my $array (sort keys %{$store{$virus}}){
        print $s1 $virus."\t".$store{$virus}{$array}{"global"}{"len"}."\t".$array."\t".
        $store{$virus}{$array}{"global"}{"type"}."\t".$store{$virus}{$array}{"global"}{"coverage"}."\t".
        $store{$virus}{$array}{"local"}{"type"}."\t".$store{$virus}{$array}{"local"}{"coverage"}."\t".$store{$virus}{$array}{"local"}{"sample"}."\n";
    }
}
close $s1;