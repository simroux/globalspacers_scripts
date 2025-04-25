#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h); 
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
	print "# Script to pre-compute a list of common spacers, that will be useful later
# Arguments :
# run
";
	die "\n";
}


my $in_file_1="spacer_clusters_and_metadata_for_alphadiv.tsv";
my $in_file_2="spacer_clusters_and_metadata_for_alphadiv_md10.tsv";

my $out_file="Common_spacer_list.tsv";

my $last;
my %seen;
my %tmp;
my %store_common;
my $i=0;
open my $tsv,"cat $in_file_1 $in_file_2 |";
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "cluster_id"){next;}
    my $code=$tab[2].",".$tab[3];
    if ($last ne $code){
        if (defined($seen{$code})){
            print "Pblm, we already saw $tab[1] ??\n";
            die("\n");
        }
        if ($last ne ""){
            &process(\%tmp,$last,\%store_common);
            $i++;
            if ($i % 10000 == 0){
                print ".. $i -- $last ..\n";
            }
        }
        %tmp=();
        $last=$code;
        $seen{$code}=1;
    }
    if (defined($tmp{$tab[0]})){
        print "pblm, already saw $tab[0] for $code ?\n";
        die("\n");
    }
    $store_common{$tab[0]}{"total_sets"}++;
    $tmp{"by_spacer"}{$tab[0]}=$tab[1];
    $tmp{"total_cover"}+=$tab[1];
}
if ($last ne ""){
    &process(\%tmp,$last,\%store_common);
}
close $tsv;

open my $s1,">",$out_file;
print $s1 "spacer\tn_set_total\tn_common_set\tlist_common_sets\n";
foreach my $spacer (sort keys %store_common){
    if (defined($store_common{$spacer}{"common_sets"})){
        my @list_c=sort keys %{$store_common{$spacer}{"common_sets"}};
        my $n_common=scalar(@list_c);
        my $list_common=join("|",@list_c);
        print $s1 $spacer."\t".$store_common{$spacer}{"total_sets"}."\t".$n_common."\t".$list_common."\n";
    }
}
close $s1;

sub process(){
    my $hash=$_[0];
    my $code=$_[1];
    my $hash_result=$_[2];
    ## Alpha diversity
    my %stats;
    my @ordered_list=sort {$$hash{"by_spacer"}{$b} <=> $$hash{"by_spacer"}{$a} or $a cmp $b} keys %{$$hash{"by_spacer"}};
    my $total=scalar(@ordered_list);
    my $max_cover=$$hash{"by_spacer"}{$ordered_list[0]};
    ## We only look at cases where we have max cover >=20
    if ($max_cover>=20){
        my $th_fifty=0.5*$max_cover;
        foreach my $spacer (@ordered_list){
            my $p=$$hash{"by_spacer"}{$spacer}/$$hash{"total_cover"};
            # print $spacer."\t".$$hash{"by_spacer"}{$spacer}."\t".$p."\n";
            if ($$hash{"by_spacer"}{$spacer}>=$th_fifty){
                $$hash_result{$spacer}{"common_sets"}{$code}=1;
            }
        }
    }
}
