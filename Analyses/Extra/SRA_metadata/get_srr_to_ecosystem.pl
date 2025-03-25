#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
use Text::Unidecode;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
	print "# Script to get a complete list of all metaG processed and add ecosystem classification next to it
# Arguments :
# none
";
	die "\n";
}

my $full_list="SRA_filtering_complete_Jul24.tsv";
my $info_SRA_run="SraRunInfo_MetaG_not_Amplicon_Public-All-Combined-Nov_18_2023.csv";
my $file_eco_img="Info_eco_IMG.tsv";
my $file_eco_mgn="MGnify_runs_to_ecosystem-summarized.tsv";
my $file_eco_sra="SRA_Eco_clean.tsv";
my $file_title="SRA_title.tsv";
my $file_amplified="list_runs_pred_amplified_from_spacer_array_data.txt";
my $file_sorting="list_cellsorting_runs.txt";
my $file_fastq_ex="excluded_list_fastq.txt";

my $out_file="Runs_to_ecosystem_and_sequencing_and_study.tsv";
my $out_file_ex="Excluded_runs_with_rationale.tsv";

my %remove_amp;
open my $tsv,"<",$file_amplified;
while(<$tsv>){
    chomp($_);
    $remove_amp{$_}=1;
}
close $tsv;
open my $tsv,"<",$file_sorting;
while(<$tsv>){
    chomp($_);
    $remove_amp{$_}=1;
}
close $tsv;

my %remove_lowdiv;
open my $tsv,"<",$file_fastq_ex;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "Run"){next;}
    if ($tab[1] eq "full_excluded"){
        if ($tab[2]==$tab[3]){
            $remove_lowdiv{$tab[0]}="low_diversity";
        }
        else{
            $remove_lowdiv{$tab[0]}="fastq_issue";
        }
    }
}
close $tsv;



my %check;
my %info;
print "Read list $full_list\n";
open my $s2,">",$out_file_ex;
open my $tsv,"<",$full_list;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if (defined($remove_amp{$tab[0]})){
        print $s2 $tab[0]."\tpredicted_amplified\n";
    }
    elsif(defined($remove_lowdiv{$tab[0]})){
        print $s2 $tab[0]."\t".$remove_lowdiv{$tab[0]}."\n";
    }
    elsif ($tab[1] eq "accepted" || $tab[1] eq "already_done"){
        $check{$tab[0]}=1;
        $info{$tab[0]}{"strategy"}=$tab[2];
        $info{$tab[0]}{"selection"}=$tab[3];
        $info{$tab[0]}{"sequencing"}=$tab[4];
        $info{$tab[0]}{"nspot"}=$tab[5];
        $info{$tab[0]}{"nbase"}=int($tab[5]*$tab[6]);
    }
    else{
        print $s2 $tab[0]."\t".$tab[1]."\n";
    }
}
close $tsv;
close $s2;

print "Read file $info_SRA_run\n";
open my $csv,"<",$info_SRA_run;
while(<$csv>){
    chomp($_);
    my $line=$_;
    while($line=~/(.*)\"([^\"]+)\"(.*)/){
        my $pref=$1;
        my $content=$2;
        my $suf=$3;
        $content=~s/,//g;
        $line=$pref.$content.$suf;
    }
    my @tab=split(",",$line);
    if ($check{$tab[0]}==1){
        $info{$tab[0]}{"nspot"}=$tab[3];
        $info{$tab[0]}{"nbase"}=$tab[4];
        $info{$tab[0]}{"sequencing"}=$tab[18];
        $info{$tab[0]}{"study"}=$tab[20];
        $info{$tab[0]}{"bioproject"}=$tab[21];
        $info{$tab[0]}{"srasample"}=$tab[24];
        $info{$tab[0]}{"biosample"}=$tab[25];
    }
}
close $csv;

my %check_exists;
print "Read SRA file $file_eco_sra\n";
open my $tsv,"<",$file_eco_sra;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "SRA Run"){next;}
    if ($check{$tab[0]}==1){
        my @t=split(";",$tab[1]);
        my $eco=$t[0].";".$t[1];
        $info{$tab[0]}{"eco"}=$eco;
        $info{$tab[0]}{"biofilm"}=$t[2];
        if (!defined($check_exists{$eco})){
            $check_exists{$eco}=1;
        }
    }
}
close $tsv;
my $in_key="extra_img_summarized.tsv";
my %translate_img;
open my $tsv,"<",$in_key;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if (defined($translate_img{$tab[0]})){
        die("Already a thing for $tab[0]?\n");
    }
    $translate_img{$tab[0]}=$tab[1];
    print $tab[0]."|||".$tab[1]."\n";
    if (!defined($check_exists{$tab[1]})){
        $check_exists{$tab[1]}=1;
    }
}
close $tsv;

print "Read IMG file $file_eco_img\n";
open my $tsv,"<",$file_eco_img;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "SRA Run"){next;}
    if ($check{$tab[0]}==1){
        if (!defined($info{$tab[0]}{"eco"})){
            print "$tab[0] had no eco, but img tells us $tab[10]\n";
            if (defined($translate_img{$tab[10]})){
                $tab[10]=$translate_img{$tab[10]};
                print "\t actually we count it as $tab[10]\n";
            }
            if ($check_exists{$tab[10]}==1){
                $info{$tab[0]}{"eco"}=$tab[10]; 
            }
            else{
                print "\t which means we should translate $tab[10] for img\n";
                <STDIN>;
            }
        }
        elsif(!($tab[10]=~/Unknown/) && ($info{$tab[0]}{"eco"}=~/Unknown/ || $info{$tab[0]}{"eco"}=~/XX/)){
            print "For $tab[0] maybe we could replace ".$info{$tab[0]}{"eco"}." by img ".$tab[10]."\n";
            if (defined($translate_img{$tab[10]})){
                $tab[10]=$translate_img{$tab[10]};
                print "\t actually we count it as $tab[10]\n";
            }
            if ($check_exists{$tab[10]}==1){
                $info{$tab[0]}{"eco"}=$tab[10];
            }
            else{
                print "\t which means we should translate $tab[10] for img\n";
                <STDIN>;
            }
        }
    }
    else{
        print "We excluded $tab[0] ??\n";
        # <STDIN>;
    }
}
close $tsv;
print "Read Mgn file $file_eco_mgn\n";
open my $tsv,"<",$file_eco_mgn;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "SRA Run"){next;}
    if ($check{$tab[0]}==1){
        if (!defined($info{$tab[0]}{"eco"})){
            print "$tab[0] had no eco, but MGNify tells us $tab[1]\n";
            if (defined($translate_img{$tab[1]})){
                $tab[1]=$translate_img{$tab[1]};
                print "\t actually we count it as $tab[1]\n";
            }
            if ($check_exists{$tab[1]}==1){
                $info{$tab[0]}{"eco"}=$tab[1]; 
            }
            else{
                print "\t which means we should translate $tab[1] for MGNify\n";
                <STDIN>;
            }
        }
        elsif(!($tab[1]=~/Unknown/) && ($info{$tab[0]}{"eco"}=~/Unknown/ || $info{$tab[0]}{"eco"}=~/XX/)){
            print "For $tab[0] maybe we could replace ".$info{$tab[0]}{"eco"}." by MGNify ".$tab[1]."\n";
            if (defined($translate_img{$tab[1]})){
                $tab[1]=$translate_img{$tab[1]};
                print "\t actually we count it as $tab[1]\n";
            }
            if ($check_exists{$tab[1]}==1){
                $info{$tab[0]}{"eco"}=$tab[1];
            }
            else{
                print "\t which means we should translate $tab[1] for MGNify\n";
                <STDIN>;
            }
        }
        # $info{$tab[0]}{"eco"}=$tab[1];
    }
}
close $tsv;

my %biosample_to_title;
print "Read file $file_title\n";
open my $tsv,"<",$file_title;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    $tab[1]=unidecode($tab[1]);
    $biosample_to_title{$tab[0]}=$tab[1];
}
close $tsv;


my @list_columns=("eco","biofilm","sequencing","bioproject","biosample","nspot","nbase");

open my $s1,">",$out_file;
print $s1 "SRA_run\tEcosystem\tBiofilm\tSequencing_platform\tBioProject\tBioSample\tN_spot\tN_bases\tTitle\n";
foreach my $run (sort keys %info){
    my $line=$run;
    foreach my $col (@list_columns){
        if (!defined($info{$run}{$col})){
            if ($col eq "eco"){
                $info{$run}{$col}="Unknown;Unknown";
            }
            elsif($col eq "biofilm"){
                $info{$run}{$col}="Unknown";
            }
            else{
                $info{$run}{$col}="NA";
            }
        }
        $line.="\t".$info{$run}{$col};
    }
    $line.="\t".$biosample_to_title{$info{$run}{"biosample"}};
    print $s1 $line."\n";
}
close $s1;

