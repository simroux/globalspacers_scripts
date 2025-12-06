#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
use Custom::Utils;
my $h=0;
my $data_folder="";
GetOptions ('help' => \$h, 'h' => \$h); 
if ($h==1 || $ARGV[0] eq "" || 'd=s'=>\$data_folder){
    print "# Script to detect conserved positions (minimum 75% conservation per residue) in potential PAM regions, when at least 5 distinct regions
# Arguments :
# -d: Data folder (i.e. /.../Data/)
";
    die "\n";
}


my $in_file_arrays=$data_folder."/Spacer_db/Array_info_filtered_for_db-Oct24-25.tsv";
my $file_nei="Neighborhood_for_pam_motif_all_imgvr.tsv";
my $file_nei_pr="Neighborhood_for_pam_motif_all_imgpr.tsv";
my $out_file="Conserved_positions_neighborhoods.tsv";


my %ref_array;
print "Loading array type information\n";
open my $tsv,"<",$in_file_arrays;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "repeat_cluster"){next;}
    my $array=$tab[0];
    my $type=$tab[1];
    $ref_array{$array}{"type"}=$type;
}
close $tsv;



my $i=0;
my $c_a="";
my %store_all_nei;
print "## Process neighborhoods to look for conserved positions\n";
open my $tsv," cat $file_nei $file_nei_pr |";
while(<$tsv>){
	chomp($_);
	my @tab=split("\t",$_);
	if ($tab[0] eq "crispr_array"){next;}
	$i++;
	$c_a=$tab[0];
	# print "Looking for potential PAMs in $tab[1]\t$tab[2] for repeat $c_a\n";
    $store_all_nei{$c_a}{"total"}++;
    my @t=split("",$tab[1]);
    my @t2=split("",$tab[2]);
    for (my $i=0;$i<=$#t;$i++){
        $store_all_nei{$c_a}{"upstream"}{$i}{$t[$i]}++;
        $store_all_nei{$c_a}{"downstream"}{$i}{$t2[$i]}++;
    }
    if ($i % 10000 == 0){
        print "$i ... \n";
    }
}

my $ratio=0;
my @tup=();
my @td=();
my %count;
open my $s1,">",$out_file;
print $s1 "repeat_cluster\tpredicted_type\tn_distinct_neighborhoods\tconserved_upstream\tconserved_downstream\tmotif_size_up\tmotif_size_down\n";
foreach my $repeat (sort keys %store_all_nei){
    my $hash=$store_all_nei{$repeat};
    %count=();
    if ($$hash{"total"}>=5){
        print "Repeat cluster $repeat\n";
        $ratio=$$hash{"total"}*0.75;
        @tup=();
        @td=();
        $count{"up"}=0;
        $count{"down"}=0;
        for (my $i=0;$i<=9;$i++){
            my @t=sort {$$hash{"upstream"}{$i}{$b} <=> $$hash{"upstream"}{$i}{$a} or $a cmp $b} keys %{$$hash{"upstream"}{$i}};
            if ($$hash{"upstream"}{$i}{$t[0]}>=$ratio){push(@tup,$t[0]); $count{"up"}++;}
            else{push(@tup,"N");}
            my @t2=sort {$$hash{"downstream"}{$i}{$b} <=> $$hash{"downstream"}{$i}{$a} or $a cmp $b} keys %{$$hash{"downstream"}{$i}};
            if ($$hash{"downstream"}{$i}{$t2[0]}>=$ratio){push(@td,$t2[0]); $count{"down"}++;}
            else{push(@td,"N");}
        }
        print "\t".join("",@tup)."\t".join("",@td)."\n";
        print $s1 $repeat."\t".$ref_array{$repeat}{"type"}."\t".$$hash{"total"}."\t".join("",@tup)."\t".join("",@td)."\t".$count{"up"}."\t".$count{"down"}."\n";
    }
}
close $s1;