#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
my $h='';
my $cmd='';
my $out='';
GetOptions ('help' => \$h, 'h' => \$h); # , 'i=s'=>\$in_file, 'r=s'=>\$out_dir_root);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to create the new final array table for the database with updated types based on Cas detection
# Arguments :
# toto
";
    die "\n";
}

my $in_cas="cas_by_repeat_cluster.tsv";
my $in_db="Array_info_filtered_for_db-updated_type.tsv";
my $out_db="Array_info_filtered_for_db-updated_type_and_cas.tsv";


my %store;
open my $tsv,"<",$in_cas;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "repeat_cluster"){next;}
    if ($tab[2] eq "unknown"){ ## Ignoring cases where we could not identify a clear CRISPR type, to avoid reporting inconsistent sets of Cas genes
        $tab[2]="Unknown";
        $tab[5]="NA";
    } 
    if ($tab[2] eq "II-C2"){
        ## Have to change the II-C2 because the nomenclature changed a bit in 10.1038/s41564-025-02180-8
        if ($tab[5]=~/Cas4/){$tab[2]="II-D";} ## Presence of Cas4 is typical of II-D
        else{$tab[2]="II-C";} ## Without Cas4, it seems like this would be more a II-C (potentially)    
    }
    $store{$tab[0]}{"cas_n"}=$tab[1];
    $store{$tab[0]}{"cas_type"}=$tab[2];
    $store{$tab[0]}{"cas_genes"}=$tab[5];
}
close $tsv;



open my $s1,">",$out_db;
open my $tsv,"<",$in_db;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "repeat_cluster"){
        print $s1 $tab[0]."\t".$tab[1]."\t".$tab[2]."\t".$tab[3]."\t".$tab[4]."\t".$tab[5]."\t".$tab[6]."\ttype_from_repeat\ttype_from_Cas_genes\tlist_Cas_genes\tn_obs_Cas_loci\n";
        next;
    }
    if ($store{$tab[0]}{"cas_n"} eq "No_cas_gene_predicted"){
        print $s1 $tab[0]."\t".$tab[1]."\t".$tab[2]."\t".$tab[3]."\t".$tab[4]."\t".$tab[5]."\t".$tab[6]."\t".$tab[1]."\tNA\tNA\t0\n";
    }
    elsif($store{$tab[0]}{"cas_type"}=~/^Hybrid/ || $store{$tab[0]}{"cas_type"} eq "Unknown"){ ## For hybrid or unknown cas genes, we revert back to the prediction from repeats (small number anyway)
        print $s1 $tab[0]."\t".$tab[1]."\t".$tab[2]."\t".$tab[3]."\t".$tab[4]."\t".$tab[5]."\t".$tab[6]."\t".$tab[1]."\t".$store{$tab[0]}{"cas_type"}."\t".$store{$tab[0]}{"cas_genes"}."\t".$store{$tab[0]}{"cas_n"}."\n";
    }
    else{
        print $s1 $tab[0]."\t".$store{$tab[0]}{"cas_type"}."\t".$tab[2]."\t".$tab[3]."\t".$tab[4]."\t".$tab[5]."\t".$tab[6]."\t".$tab[1]."\t".$store{$tab[0]}{"cas_type"}."\t".$store{$tab[0]}{"cas_genes"}."\t".$store{$tab[0]}{"cas_n"}."\n";
    }
}
close $tsv;
close $s1;