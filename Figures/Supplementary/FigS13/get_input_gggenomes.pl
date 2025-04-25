#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
use Custom::Utils;
my $h=0;
my $uvig_id="";
GetOptions ('help' => \$h, 'h' => \$h, 'i=s'=>\$uvig_id); # , 'i=s'=>\$in_file, 'r=s'=>\$out_dir_root);
if ($h==1 || $ARGV[0] eq "" || $uvig_id eq ""){ # If asked for help or did not set up any argument
    print "# Script to get the input file for gggenomes
# Arguments :
# -i: uvig_id (e.g. IMGVR_UViG_3300045988_178767)
";
    die "\n";
}

my $info_file="../../../Data/Spacer_db/IMGVR_sequence_information_Oct17.tsv";

my $out_file_1=$uvig_id."_genome.tsv";
my $out_file_2=$uvig_id."_genes.tsv";

my $metag="";
my $scaffold_id="";
my $length=0;

print "Reading $info_file\n";
open my $tsv,"<",$info_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    my @t=split(/\|/,$tab[0]);
    if ($t[0] eq $uvig_id){
        $metag=$t[1];
        $scaffold_id=$t[2];
        $length=$tab[3];
    }
}
close $tsv;

my $gff_file="../../../Data/Additional_data/".$metag."_____".$scaffold_id.".gff";

open my $s1,">",$out_file_1;
print $s1 "seq_id\tlength\n";
print $s1 $uvig_id."\t".$length."\n";
close $s1;

print "Reading $gff_file\n";
open my $s2,">",$out_file_2;
print $s2 "seq_id\tstart\tend\ttype\tstrand\tfeat_id\tproduct\tsum_product\n";
open my $gff,"<",$gff_file;
while(<$gff>){
    chomp($_);
    if ($_=~/^#/){next;}
    my @tab=split("\t",$_);
    if ($tab[2] eq "CDS"){
        $tab[8]=~/ID=([^;]+)/;
        my $id=$1;
        $tab[8]=~/Product=([^;]+)/;
        my $prod=$1;
        my $sum_prod="NA";
        if ($prod eq "hypothetical protein"){$sum_prod="hypothetical protein";}
        elsif($prod=~/^pfam\S+ (\S+)/){
            $sum_prod=$1;
        }
        print $s2 $uvig_id."\t".$tab[3]."\t".$tab[4]."\t".$tab[2]."\t".$tab[6]."\t".$id."\t".$prod."\t".$sum_prod."\n";
    }
}
close $s2;