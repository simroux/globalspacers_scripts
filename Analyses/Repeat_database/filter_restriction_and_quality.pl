#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
use Custom::Utils;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
	print "# Script to get a complete list of all repeats we can use
# Arguments :
# toto
";
	die "\n";
}

my $filtered_list="all_repeats_clstr-filtered-clean.tab";
my $info_and_complex="Repeat_QC/Array_info_with_complexity_for_db.tsv";
my $taxo_info="Phylodist/all_repeats_clstr-filtered-clean_wlca_phylodist.tab";
my $updated_confidence="Updated_LCA_confidence.tsv";
## Note - for taxo, we will expand the lca status to an lca origin -> goes from confidence (Low/high) to Genome_high, Genome_low, Contig
## We will also exclude everything that is Euk and Viruses outside of the Caudoviricetes


my $out_file="Array_info_filtered_for_db.tsv";

my %check;
my %store;
open my $tsv,"<",$filtered_list;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    $check{$tab[1]}=1;
    $store{$tab[1]}{"type"}=$tab[3];
    $store{$tab[1]}{"lca_confidence"}=$tab[8];
    $store{$tab[1]}{"lca"}=$tab[7];
}
close $tsv;

open my $tsv,"<",$updated_confidence;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    $store{$tab[0]}{"lca_confidence"}=$tab[1];
}
close $tsv;


# my %info_taxo;
open my $tsv,"<",$taxo_info;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "Array"){next;}
    if (defined($check{$tab[1]})){
        if ($tab[11]=~/^d__unclassified/){
            $store{$tab[1]}{"lca_meta"}="NA_long"; ## Necessarily long since we had at least some affi until we reached this useless LCA
        }
        else{
            $store{$tab[1]}{"lca_meta"}=$tab[11];
            if ($tab[11]=~/^d__Euk/ || ($tab[11]=~/^d__Viruses/ && !($tab[11]=~/Caudoviricetes/))){
                print "We will exclude $tab[1] -- $tab[11]\n";
                $store{$tab[1]}{"exclude_meta"}=1;
            }   
        }
    }
}
close $tsv;

open my $s1,">",$out_file;
open my $tsv,"<",$info_and_complex;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "repeat_cluster"){
        print $s1 "repeat_cluster\ttype\tlca_origin\tlca_full\tlca_class\tlca_family\tlca_genus\n";
    }
    else{
        if ($check{$tab[0]}==1){
            my $flag=$tab[7];
            if ($flag eq "NA"){
                ## Replacing type flag because it's better calculated now (including stuff removed because of restrictions)
                if ($tab[1] ne $store{$tab[0]}{"type"}){
                    print "replacing $tab[1] by ".$store{$tab[0]}{"type"}." for $tab[0] -- Type\n";
                    # <STDIN>;
                    $tab[1]=$store{$tab[0]}{"type"};
                }
                ## Updating taxo and confidence
                if ($tab[2] ne $store{$tab[0]}{"lca_confidence"}){
                    print "replacing $tab[2] by ".$store{$tab[0]}{"lca_confidence"}." for $tab[0] -- LCA Confidence\n";
                    if ($tab[2] eq "Medium-confidence" && $store{$tab[0]}{"lca_confidence"} eq "Low-confidence"){}
                    elsif ($tab[2] eq "High-confidence" && $store{$tab[0]}{"lca_confidence"} eq "Medium-confidence"){}
                    elsif ($tab[2] eq "High-confidence" && $store{$tab[0]}{"lca_confidence"} eq "Low-confidence"){}
                    elsif ($tab[2] eq "Low-confidence" && $store{$tab[0]}{"lca_confidence"} eq "Medium-confidence"){}
                    else{<STDIN>;}
                    $tab[2]=$store{$tab[0]}{"lca_confidence"};
                }
                if ($tab[3] ne $store{$tab[0]}{"lca"}){
                    print "replacing $tab[3] by ".$store{$tab[0]}{"lca"}." for $tab[0] -- LCA\n";
                    <STDIN>;
                    $tab[3]=$store{$tab[0]}{"lca"};
                }
                ## Adding meta taxo if needed
                if ($tab[3]=~/^d__unclassified/){ ## Not really an affiliation, we remove
                    $tab[2]="No-information_long"; ## These are long by definition because they were in a genome
                    $tab[3]="NA";
                    $tab[4]="NA";
                    $tab[5]="NA";
                    $tab[6]="NA";
                }
                elsif ($tab[2] eq "Low-confidence"){$tab[2]="Genome_low-confidence";}
                elsif ($tab[2] eq "Medium-confidence"){$tab[2]="Genome_medium-confidence";}
                elsif ($tab[2] eq "High-confidence"){$tab[2]="Genome_high-confidence";}
                elsif ($tab[2] eq "No-information"){
                    if (defined($store{$tab[0]}{"lca_meta"})){
                        if ($store{$tab[0]}{"exclude_meta"}==1){
                            print "Removing $tab[0] from final table because exclude meta\n";
                            next;
                        }
                        if ($store{$tab[0]}{"lca_meta"} eq "NA_short"){
                            $tab[2]="No-information_short";
                        }
                        elsif ($store{$tab[0]}{"lca_meta"} eq "NA_long"){
                            $tab[2]="No-information_long";
                        }
                        elsif ($store{$tab[0]}{"lca_meta"} eq "NA"){
                            print "THIS SHOULD NOT HAPPEN \n";
                            die("\n");
                        }
                        else{
                            $tab[2]="Contig";
                            my @t=split(";",$store{$tab[0]}{"lca_meta"});
                            $tab[3]=$store{$tab[0]}{"lca_meta"};
                            $tab[4]=join(";",@t[0..2]);
                            $tab[5]=join(";",@t[0..4]);
                            $tab[6]=join(";",@t[0..5]);
                        }
                    }
                }
                else{
                    print "Did not anticipate $tab[2]\n";
                    die("\n");
                }
                print $s1 join("\t",@tab[0..6])."\n";
            }
            else{
                print "we remove $tab[0] because of low complexity\n";
            }
        }
        else{
            print "we remove $tab[0] because of restrictions\n";
        }
        
    }
}
close $tsv;
close $s1;


