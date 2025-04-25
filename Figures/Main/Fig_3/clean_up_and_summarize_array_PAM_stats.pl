#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use File::Basename;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
	print "# Script to get the summarized stats on repeat <-> PAM connections
# Arguments :
# none
";
	die "\n";
}

## We use a relatively stringent cutoff:
# At least 50% of the neighborhood for 2 and more
# At least 80% for mononucleotide
# Also a minimum number of 50 neighborhood observed, except if PAM is expected based on type, or type unknown and canonical PAM

my $array_sample_file="../Fig_1/array_to_sample_counts.tsv";
my $curated_file="../../../Data/Additional_data/manually_curated_PAMs_per_type.tsv";

my $pam_file="../../../Analyses/Target_IMGVR_IMGPR/PAM_detection/Repeat_to_PAM_assignment.tsv";
my $out_file_clean="Repeat_to_PAM_assignment_clean.tsv";
my $out_file_clean_sup_tab="Repeat_to_PAM_assignment_for_sup_table.tsv";
my $out_file_sum_by_spacer="Summary_PAM_detection_by_spacer.tsv";
my $out_file_sum="Summary_PAM_detection_by_type.tsv";

if (-e $out_file_clean){
    die("$out_file_clean already exists, I refuse to do anything\n");
}

### First clean up the assignment per array
my %info_array;
open my $tsv,"<",$array_sample_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "repeat_cluster"){next;}
    $info_array{$tab[0]}{"n_sample"}=$tab[1];
    $info_array{$tab[0]}{"type"}=$tab[3];
    $tab[5]=~s/\"//g;
    $info_array{$tab[0]}{"lca"}=$tab[5];
}
close $tsv;


my %type_to_motif;
my %motif_to_type;
open my $tsv,"<",$curated_file;
while(<$tsv>){
	chomp($_);
	my @tab=split("\t",$_);
	if ($tab[0] eq "Type"){next;}
	## We store the link between the CRISPR type and the motif
	$type_to_motif{$tab[0]}{$tab[1]}=1;
    my @t=split("-",$tab[0]);
    $motif_to_type{$tab[1]}{"all"}{$t[0]}=1;
}
close $tsv;
foreach my $motif (sort keys %motif_to_type){
    my @t=sort keys %{$motif_to_type{$motif}{"all"}};
    if ($#t==0){
        $motif_to_type{$motif}{"summary"}="known_type_".$t[0];
    }
    else{
        $motif_to_type{$motif}{"summary"}="known_types_".$t[0];
        for (my $i=1;$i<=$#t;$i++){
            $motif_to_type{$motif}{"summary"}.="_".$t[$i];
        }
    }
}

my %store_clean_by_array;
open my $s1,">",$out_file_clean;
print $s1 "array\ttype\tpredicted_pam\tprediction_type\tpam_freq\tweighted_pam_freq\tn_obs\tarray_lca\n";
open my $sup,">",$out_file_clean_sup_tab;
print $sup "array\ttype\tn_obs\tpredicted_pam\tprediction_type\tpam_freq\tweighted_pam_freq\n";
open my $tsv,"<",$pam_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "Array"){next;}
    my $type=$info_array{$tab[0]}{"type"};
    my $ratio_count=$tab[10];
    my $ratio_spacers=$tab[8];
    my $motif=$tab[5];
    my $motif_key=$tab[6];
    my $denovo=$tab[2];
    my $n_n=$motif=~tr/N//;
    my $n_non_n=length($motif)-$n_n;
    if ($motif eq "NA"){$n_n=0; $n_non_n=0;}
    my $n_neigh=$tab[4];
    my $detect_type="";
    if ($n_neigh>=50){ ## At least 50 neighborhoods
        if ($motif eq "NA"){
            $detect_type="no_motif_detected";
        }
        else{
            if ($n_non_n>=2){ ## Motif of at least 2
                if ($ratio_count>=50){ ## Motif found in at least 50% of the neighborhoods
                    if ($motif eq $denovo){ ## this is a de novo motif, we clean up and report
                        $motif=&clean_denovomotif($motif);
                        $detect_type="de_novo";
                    }
                    else{ ## this is a known motif, we check if it's expected
                        if ($type eq "Unknown"){
                            ## we can't have expected, so we simply report where the motif came from
                            $detect_type=$motif_to_type{$motif}{"summary"};
                        }
                        else{
                            if ($type_to_motif{$type}{$motif}==1){
                                $detect_type="expected";
                            }
                            else{
                                $detect_type=$motif_to_type{$motif}{"summary"};
                            }
                        }
                    }
                }
                else{ ## Motif found in less than 50% of the neighborhoods
                    $motif="NA";
                    $detect_type="no_motif_detected";
                }
            }
            else{ ## Should be 1
                if ($n_non_n!=1){print "Pblm - n non n should be 1\n"; die("\n");}
                if ($ratio_count>=80){ ## Motif found in at least 80% of the neighborhoods (more stringent for 1-nucl motifs)
                    $motif=&clean_denovomotif($motif);
                    $detect_type="de_novo"; ## All the 1-nucleotide motifs are de novo
                }
                else{
                    $motif="NA";
                    $detect_type="no_motif_detected";
                }
            }
        }
    }
    else{ ## less than 50 neighborhood, we only consider if we have an expected canonical motif (or a canonical motif if no expectation)
        if ($motif eq "NA"){
            $detect_type="no_motif_less_than_50";
        }
        elsif ($type eq "Unknown"){
            $motif="NA";
            $detect_type="unknown_less_than_50";
        }
        else{
            if ($n_non_n>=2){ ## Motif of at least 2
                if ($ratio_count<50){ ## If less than 50%, we ignore
                    $motif="NA";
                    $detect_type="no_motif_less_than_50";
                }
                elsif ($type_to_motif{$type}{$motif}==1){ ## If more than 50% and expected, this is the only case we take
                    $detect_type="expected";
                }
                elsif(defined($motif_to_type{$motif})){
                    $motif="NA";
                    $detect_type="unexpected_less_than_50";
                }
                else{
                    $motif="NA";
                    $detect_type="no_motif_less_than_50";
                }
            }
            else{
                $motif="NA";
                $detect_type="no_motif_less_than_50";
            }
        }
    }
    if ($motif eq "NA"){
        $ratio_count="NA"; ## If no motif, we can't count
        $ratio_spacers="NA";
    }
    else{
        $store_clean_by_array{$tab[0]}=$motif_key;
    }
    print $s1 $tab[0]."\t".$tab[1]."\t".$motif."\t".$detect_type."\t".$ratio_count."\t".$ratio_spacers."\t".$n_neigh."\t".$info_array{$tab[0]}{"lca"}."\n";
    if ($n_neigh>=50 || $detect_type eq "expected"){
        print $sup $tab[0]."\t".$tab[1]."\t".$n_neigh."\t".$motif."\t".$detect_type."\t".$ratio_count."\t".$ratio_spacers."\n";
    }
}
close $tsv;
close $s1;
close $sup;

## Now we can read this cleaned-up file and prepare the summary
my %summary;
open my $s3,">",$out_file_sum_by_spacer;
print $s3 "array\ttype\tlvl_1_type\tcount_neigh\trefined_motif\tper_neigh\n";
open my $tsv,"<",$out_file_clean;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "array"){next};
    ## 
    my $type=$tab[1];
    if ($type eq "I-F_T"){$type="I-F";} ## For summary, we gather I-F_T with I-F
    my @t=split("-",$type);
    my $motif=$tab[2];
    my $det_type=$tab[3];
    my $freq_pam=$tab[4];
    my $n_neigh=$tab[6];
    ## 
    print $s3 $tab[0]."\t".$type."\t".$t[0]."\t".$n_neigh."\t".$motif."\t".$freq_pam."\n";
    ## We ignore cases without a motif and less than 50
    if ($det_type=~/_less_than_50/){next;}
    ## 
    $summary{$type}{"motifs"}{$motif}++;
    $summary{$type}{"total"}++;
    if ($n_neigh>=10){
        $summary{$type}{"motifs_min_10"}{$motif}++;
        if ($n_neigh>=100){ ## NOTE - ARBITRARY FOR NOW -> At least 100 distinct neighborhoods, 
            $summary{$type}{"motifs_min_100"}{$motif}++;
        }
    }
}
close $tsv;
close $s3;

open my $s1,">",$out_file_sum;
print $s1 "Type\tLvl_1_type\tMotif\tExpected\tCount\tCount_min_10\tCount_min_100\n";
my %tmp;
foreach my $type (sort keys %summary){
    %tmp=();
    my @t=split("-",$type);
    foreach my $motif (sort keys %{$summary{$type}{"motifs"}}){
        my $exp="other";
        if ($type_to_motif{$type}{$motif}==1){$exp="primary";}
        if ($motif eq "NA"){$exp="no_motif";}
        if (!defined($summary{$type}{"motifs_min_10"}{$motif})){$summary{$type}{"motifs_min_10"}{$motif}=0;}
        if (!defined($summary{$type}{"motifs_min_100"}{$motif})){$summary{$type}{"motifs_min_100"}{$motif}=0;}
        print $s1 $type."\t".$t[0]."\t".$motif."\t".$exp."\t".$summary{$type}{"motifs"}{$motif}."\t".
        $summary{$type}{"motifs_min_10"}{$motif}."\t".$summary{$type}{"motifs_min_100"}{$motif}."\n";
    }
}
close $s1;

print "Cleaning up neigh file VR\n";
open my $tsv,"<",$neigh_file_vr;
open my $s1,">",$out_file_neigh_vr;
print $s1 "array\tupstream\tdownstream\tpam_presence\tpam_motif\n";
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "Array"){next;}
    my @t=split(",",$tab[6]);
    foreach my $motif (@t){
        if ($store_clean_by_array{$tab[0]} eq $motif){
            print $s1 $tab[0]."\t".$tab[2]."\t".$tab[3]."\tyes\t".$motif."\n";
            last;
        }
    }
}
close $tsv;
close $s1;

print "Cleaning up neigh file PR\n";
open my $tsv,"<",$neigh_file_pr;
open my $s1,">",$out_file_neigh_pr;
print $s1 "array\tupstream\tdownstream\tpam_presence\tpam_motif\n";
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "Array"){next;}
    my @t=split(",",$tab[6]);
    foreach my $motif (@t){
        if ($store_clean_by_array{$tab[0]} eq $motif){
            print $s1 $tab[0]."\t".$tab[2]."\t".$tab[3]."\tyes\t".$motif."\n";
            last;
        }
    }
}
close $tsv;
close $s1;

sub clean_denovomotif(){
    my $motif=$_[0];
    while($motif=~/^N(.*)/){$motif=$1;}
    while(length($motif)<3){$motif="N".$motif;}
    return $motif;
}