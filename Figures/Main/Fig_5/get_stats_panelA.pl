#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
use Custom::Utils;
# use Parallel::ForkManager;
my $h=0;
# my $in_file;
# my $out_dir_root="Results/";
GetOptions ('help' => \$h, 'h' => \$h); # , 'i=s'=>\$in_file, 'r=s'=>\$out_dir_root);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to get the numbers for an overview of virus targeting by multiple taxa
# Arguments :
# toto
";
    die "\n";
}

my $info_virus_file="../../Data/IMGVR_sequence_information_Oct17.tsv";

my $in_file="Summary_multitaxa.tsv";
my $in_file_multiclass="Summary_multiclass_wprofile.tsv";

my $out_file="Counts_for_multitaxa_frequency.tsv";
my $out_file_interest="List_of_interest.tsv";

my %info_virus;
my %check_virus;
print "Reading $info_virus_file ..\n";
open my $tsv,"<",$info_virus_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "uvig"){next;}
    if ($tab[10] eq "other" || $tab[10] eq "unknown"){next;}
    if ($tab[10] eq "phage" || $tab[10] eq "archaea"){
        $info_virus{$tab[0]}{"quality"}=$tab[5];
        $info_virus{$tab[0]}{"length"}=$tab[3];
        $info_virus{$tab[0]}{"taxo"}=$tab[6];
        if ($tab[5] eq "Reference"){$tab[5]="High-quality";}
        $check_virus{$tab[0]}{"quality"}=$tab[5];
    }
    else{
        die("Unexpected type - $tab[10]\n");
    }
}
close $tsv;

my %seen;
my %count;
print "Reading $in_file_multiclass ..\n";
open my $s2,">",$out_file_interest;
open my $tsv,"<",$in_file_multiclass;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if (defined($check_virus{$tab[0]})){
        my $type_multiclass="";
        my $code_high="";
        my $code_high_mis="";
        if ($tab[10]>0){
            if ($tab[10]==1){
                $code_high="single_taxon";
                if ($tab[9] eq "phylum" || $tab[9] eq "domain"){$type_multiclass="no_class_with_high";}
                else{$type_multiclass="one_class_high_other_not_high";}
            }
            else{
                $code_high="multi_".$tab[9];
                if ($tab[10]>1 && ($tab[9] eq "class" || $tab[9] eq "phylum" || $tab[9] eq "domain")){
                    if ($tab[12]>0){
                        if ($tab[12]==1){
                            $code_high_mis="single_taxon";
                            if ($tab[11] eq "phylum" || $tab[11] eq "domain"){$type_multiclass="no_class_with_high_pos";}
                            else{$type_multiclass="one_class_highpos_other_class_high";}
                        }
                        else{
                            $code_high_mis="multi_".$tab[11];
                            if ($tab[11] eq "class" || $tab[11] eq "phylum" || $tab[11] eq "domain"){$type_multiclass="multi_class_with_high_pos";}
                            else{$type_multiclass="one_class_highpos_other_class_high";}
                        }
                    }
                    elsif($tab[12]==0){
                        $code_high_mis="no_positive";
                        $type_multiclass="no_class_with_high_pos";
                    }
                }
                else{
                    $type_multiclass="one_class_high_other_not_high"; ## We are in multi, but multi order / family / genus, so single class
                }
            }
        }
        else{
            $code_high="no_high";
            $type_multiclass="no_class_with_high";
        }
        $count{"hc_high"}{$code_high}++;
        if ($code_high_mis ne ""){
            $count{"hc_high_mis"}{$code_high_mis}++;
        }
        print $s2 $tab[0]."\t".$type_multiclass."\n";
        $seen{$tab[0]}=1;
    }
}
close $tsv;
print "Reading $in_file ..\n";
open my $tsv,"<",$in_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if (defined($check_virus{$tab[0]})){
        my $code="";
        if ($tab[8]>0){
            if ($tab[8]==1){$code="single_taxon";}
            else{$code="multi_".$tab[7];}
            $count{"hc"}{$code}++;
            if (!defined($seen{$tab[0]})){
                print $s2 $tab[0]."\tsingle_class\n";
            }
        }
    }
}
close $tsv;
close $s2;


open my $s1,">",$out_file;
print $s1 "type\tcat\tcount\n";
foreach my $type (keys %count){
    foreach my $cat (keys %{$count{$type}}){
        print $s1 $type."\t".$cat."\t".$count{$type}{$cat}."\n";
    }
}
close $s1;

