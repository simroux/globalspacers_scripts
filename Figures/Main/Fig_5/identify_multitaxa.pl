#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use File::Basename;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
	print "# Script to look for viruses targeted by multiple taxa 
# Arguments :
# run
";
	die "\n";
}

my %levels=("0"=>"domain","1"=>"phylum","2"=>"class","3"=>"order","4"=>"family","5"=>"genus","6"=>"species");

my $in_file="Target_to_repeat.tsv";
my $uvig_file="../../../Data/Spacer_db/IMGVR_sequence_information_Oct17.tsv";

my $out_file="Detailed_multitaxa.tsv";
my $out_file_bis="Summary_multitaxa.tsv";
my $list_for_further="List_multiclass_uvigs.txt";
my $list_for_further_nohq="List_multiclass_uvigs_nohq.txt";
my $min_count=10; ## Minimum number of spacers to consider "not low", per virus-repeat pair

my %store;
print "#######################\n";
print "Reading $in_file ...\n";
my $i=0;
open my $tsv,"<",$in_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "target_id"){next;}
    $store{$tab[0]}{$tab[1]}{"lca_confidence"}=$tab[2];
    $store{$tab[0]}{$tab[1]}{"lca_genus"}=$tab[3];
    $store{$tab[0]}{$tab[1]}{"count"}=$tab[4];
    $i++;
    if ($i % 1000000 == 0){
        print ".. $i ..\n";
    }
}
close $tsv;

print "#######################\n";
my %check_multiclass;
my %tmp;
open my $s2,">",$out_file_bis;
print $s2 "uvig\tn_array\tmax_rank_all\tn_taxon_max_rank\tn_array_hc\tmax_rank_hc\tn_taxon_max_rank_hc\tmax_rank_hc_10p\tn_taxon_max_rank_hc_10p\n";
open my $s1,">",$out_file;
print $s1 "uvig\tarray_type\trank\tn_array\tn_taxa\tlist_taxa\n";
foreach my $uvig (sort keys %store){
    %tmp=();
    print "## $uvig ##\n";
    my $n_array_total=0;
    my $n_array_hc=0;
    foreach my $array (sort keys %{$store{$uvig}}){
        $n_array_total++;
        $tmp{"all"}{$store{$uvig}{$array}{"lca_genus"}}+=$store{$uvig}{$array}{"count"};
        if ($store{$uvig}{$array}{"lca_confidence"} eq "Genome_high-confidence" || $store{$uvig}{$array}{"lca_confidence"} eq "Genome_medium-confidence"){
            $n_array_hc++;
            $tmp{"hc_only"}{$store{$uvig}{$array}{"lca_genus"}}+=$store{$uvig}{$array}{"count"};
            if ($store{$uvig}{$array}{"count"}>=$min_count){
                $tmp{"hc_and_10p"}{$store{$uvig}{$array}{"lca_genus"}}+=$store{$uvig}{$array}{"count"};
            }
        }
    }
    print "With all arrays\n";
    my ($max_rank,$n_taxon)=&get_count($tmp{"all"},$uvig,"all",$s1);
    print "\t".$max_rank."\t".$n_taxon."\n";
    print "Only hc arrays\n";
    my ($max_rank_hc,$n_taxon_hc)=&get_count($tmp{"hc_only"},$uvig,"high-confidence_only",$s1);
    print "\t".$max_rank_hc."\t".$n_taxon_hc."\n";
    print "Only hc arrays and at least 10 hits\n";
    my ($max_rank_hc_10p,$n_taxon_hc_10p)=&get_count($tmp{"hc_and_10p"},$uvig,"high-confidence_only_10p",$s1);
    print "\t".$max_rank_hc_10p."\t".$n_taxon_hc_10p."\n";
    print $s2 $uvig."\t".$n_array_total."\t".$max_rank."\t".$n_taxon."\t".$n_array_hc."\t".$max_rank_hc."\t".$n_taxon_hc."\t".$max_rank_hc_10p."\t".$n_taxon_hc_10p."\n";
    if ($n_taxon_hc_10p>1){
        if ($max_rank_hc_10p eq "class" || $max_rank_hc_10p eq "phylum" || $max_rank_hc_10p eq "domain"){
            $check_multiclass{$uvig}=1;
        }
    }
}
close $s1;
close $s2;

my %check_uvig;
print "Reading $uvig_file\n";
open my $tsv,"<",$uvig_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "uvig"){next;}
    if ($tab[5] eq "High-quality" || $tab[5] eq "Reference"){
    }
    elsif($tab[10] eq "phage" || $tab[10] eq "archaea"){
        $check_uvig{$tab[0]}=1;
    }
}
close $tsv;

open my $s3,">",$list_for_further;
open my $s4,">",$list_for_further_nohq;
print $s4 "uvig\n";
foreach my $uvig (keys %check_multiclass){
    print $s3 $uvig."\n";
    if ($check_uvig{$uvig}==1){
        print $s4 $uvig."\n";
    }
}
close $s3;
close $s4;


sub get_count(){
    my $hash=$_[0];
    my $uvig=$_[1];
    my $type=$_[2];
    my $sout=$_[3];
    my %count;
    foreach my $lca (keys %{$hash}){
        my @t=split(";",$lca);
        for (my $i=0;$i<=5;$i++){
            if ($t[$i]=~/__unclassified/){}
            else{
                my $taxo=join(";",@t[0..$i]);
                $count{$i}{$taxo}+=$$hash{$lca};
                # print $i."\t".$taxo."\t".$count{$i}{$taxo}."\n";
            }
        }
    }
    # print "=== summary ===\n";
    my $best_rank="None";
    my $best_rank_ntaxon=0;
    my $tag=0;
    my @test=keys %{$count{0}};
    if ($#test<0){} ## We end there, nothing that passed our min count cutoff
    else{
        for (my $i=0;$i<=5;$i++){
            my @list=sort {$count{$i}{$b} <=> $count{$i}{$a}} keys %{$count{$i}};
            if ($#list<0){ 
                ## We stop there
                $i=10;
            }
            else{
                my $line="";
                foreach my $taxo (@list){
                    $line.=$taxo." (".$count{$i}{$taxo}.") ";
                }
                chop($line);
                my $n_taxa=scalar(@list);
                if ($tag==0){
                    $best_rank=$levels{$i};
                    $best_rank_ntaxon=$n_taxa;
                    if ($n_taxa>1){
                        $tag=1;
                    }
                }
                print $s1 $uvig."\t".$type."\t".$levels{$i}."\t".$n_taxa."\t".$line."\n";
            }
        }
    }
    return ($best_rank,$best_rank_ntaxon);
}
