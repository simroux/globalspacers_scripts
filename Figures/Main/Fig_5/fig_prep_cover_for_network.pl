#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
my $panel="";
GetOptions ('help' => \$h, 'h' => \$h, 'p=s'=>\$panel); 
if ($h==1 || $panel eq ""){ # If asked for help or did not set up any argument
    print "# Script to import the sequences and info we want for the panels A or B in figure 5
# Arguments :
# -p: a or b depending on the panel we want
";
    die "\n";
}

my $w_dir="";
my $list_select;
if ($panel eq "a"){
    $w_dir="panel_a/";
    $list_select="list_panel_a";
}
elsif ($panel eq "b"){
    $w_dir="panel_b/";
    $list_select="list_panel_b";
}
else{die("I don't know this panel\n");}

if (!(-d $w_dir)){
    &run_cmd("mkdir $w_dir");
}


my $cl_to_repeat="../../../Analyses/Target_Custom_datasets/All_spacers_vs_selected_viruses/nr_spacers_hq_vs_selected_uvigs_db_all_hits_spacer_info.tsv";
my $hit_file="../../../Analyses/Target_Custom_datasets/All_spacers_vs_selected_viruses/nr_spacers_hq_vs_selected_uvigs_db_all_hits.tsv";
my $in_fasta="../../../Data/Additional_data/selected_uvigs_for_network.fna";

my $out_file=$w_dir."/selected_uvigs_coverage.tsv";
my $out_file_net_edges=$w_dir."/selected_uvigs_coverage_net-edges.tsv";
my $out_file_net_nodes=$w_dir."/selected_uvigs_coverage_net-nodes.tsv";

my %check;
open my $tsv,"<",$list_select;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    $check{$tab[0]}=1;
    print "selecting $tab[0]\n";
}
close $tsv;

print "Loading cluster to repeat and taxo link\n";
my %cl_to_sp;
my %repeat_to_taxon;
my $i=0;
open my $tsv,"<",$cl_to_repeat;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    ### Note: we ignore repeats that are without taxonomy or with only contig taxonomy
    if ($tab[25]=~/^Genome/){
        $cl_to_sp{$tab[0]}{$tab[20]}=1;
        if (!defined($repeat_to_taxon{$tab[20]})){
            $repeat_to_taxon{$tab[20]}=&clean($tab[26]);
        }
    }
    $i++;
    if ($i % 10000 == 0){
        print "... $i ...\n";
    }
}
close $tsv;


$i=0;
my %store_cover;
my %store_n;
my $last=0;
print "Loading cover info\n";
open my $tsv,"<",$hit_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "cluster_id"){next;}
    if ($check{$tab[1]}==1){
        if (!defined($store_cover{$tab[1]})){$store_cover{$tab[1]}={};}
        # print "Add to store for $tab[1]: $tab[2]/$tab[3]\n";
        foreach my $repeat (keys %{$cl_to_sp{$tab[0]}}){
            my $tax=$repeat_to_taxon{$repeat};
            if (!defined($store_cover{$tab[1]}{$tax})){$store_cover{$tab[1]}{$tax}={};}
            &add_to_store($store_cover{$tab[1]}{$tax},$tab[2],$tab[3]);
            $store_n{$tab[1]}{$tax}{"n_unique"}{$tab[0]}++;
            $store_n{$tab[1]}{$tax}{"n_repeat"}{$repeat}++;
            $store_n{$tab[1]}{$tax}{"number"}++;
            $store_n{$tab[1]}{$tax}{"mismatch"}{$tab[5]}++;
            $store_n{$tab[1]}{$tax}{"n_base"}+=($tab[3]-$tab[2])+1;
            $last=0;
            foreach my $key (sort {$a <=> $b} keys %{$store_cover{$tab[1]}{$tax}}){
                $last=$store_cover{$tab[1]}{$tax}{$key};
            }
        }
    }
    $i++;
    if ($i % 10000 == 0){
        print "... $i ...\n";
    }
}
close $tsv;

my %info_uvig;
my $c_c="";
open my $fa,"<",$in_fasta;
while(<$fa>){
    chomp($_);
    if ($_=~/^>(\S+)/){
        $c_c=$1;
    }
    else{
        $info_uvig{$c_c}{"length"}+=length($_);
    }
}
close $fa;

my %check_uvig;
my %check_tax;
print "Loading uvig info and calculating coverage\n";
open my $sedge,">",$out_file_net_edges;
print $sedge "node_1\tnode_2\tn_unique\tratio_covered\tmis_cat\n";
open my $s1,">",$out_file;
print $s1 "uvig\ttaxon\tlength\tcovered\tpercent\tn_repeat\tn_unique_spacer\tmedian_hit_per_spacer\tmax_hit\thit_profile\tprofile_cat\n";
foreach my $uvig (sort keys %info_uvig){
    if (defined($store_cover{$uvig})){
        print "Preparing output for $uvig ..\n";
        foreach my $tax (keys %{$store_cover{$uvig}}){
            my $covered=&get_cover_length($store_cover{$uvig}{$tax});
            my $len=$info_uvig{$uvig}{"length"};
            my $ratio=sprintf("%.02f",$covered/$len*100);
            my ($n_unique,$median_hit,$max_hit)=&get_stats($store_n{$uvig}{$tax}{"n_unique"});
            my $n_repeat=scalar(keys %{$store_n{$uvig}{$tax}{"n_repeat"}});
            if (!defined($store_n{$uvig}{$tax}{"mismatch"}{"0"})){$store_n{$uvig}{$tax}{"mismatch"}{"0"}=0;}
            if (!defined($store_n{$uvig}{$tax}{"mismatch"}{"1"})){$store_n{$uvig}{$tax}{"mismatch"}{"1"}=0;}
            if (!defined($store_n{$uvig}{$tax}{"mismatch"}{"2"})){$store_n{$uvig}{$tax}{"mismatch"}{"2"}=0;}
            if (!defined($store_n{$uvig}{$tax}{"mismatch"}{"3"})){$store_n{$uvig}{$tax}{"mismatch"}{"3"}=0;}
            my $profile=$store_n{$uvig}{$tax}{"mismatch"}{"0"}.";".$store_n{$uvig}{$tax}{"mismatch"}{"1"}.";".$store_n{$uvig}{$tax}{"mismatch"}{"2"}.";".$store_n{$uvig}{$tax}{"mismatch"}{"3"};
            my $n_good=$store_n{$uvig}{$tax}{"mismatch"}{"0"}+$store_n{$uvig}{$tax}{"mismatch"}{"1"};
            my $sum=$store_n{$uvig}{$tax}{"mismatch"}{"0"}+$store_n{$uvig}{$tax}{"mismatch"}{"1"}+$store_n{$uvig}{$tax}{"mismatch"}{"2"}+$store_n{$uvig}{$tax}{"mismatch"}{"3"};
            my $ratio_emp=$n_good/$sum;
            my $emp_cat="unknown";
            if ($ratio_emp<=0.5){$emp_cat="mostly_inexact";}
            elsif($ratio_emp<=0.8){$emp_cat="mixed";}
            elsif($ratio_emp>0.8){$emp_cat="mostly_exact";}
            print $uvig."\t".$tax."\t".$len."\t".$covered."\t".$ratio."\t".$n_repeat."\t".$n_unique."\t".$median_hit."\t".$max_hit."\t".$profile."\t".$emp_cat."\n";
            print $s1 $uvig."\t".$tax."\t".$len."\t".$covered."\t".$ratio."\t".$n_repeat."\t".$n_unique."\t".$median_hit."\t".$max_hit."\t".$profile."\t".$emp_cat."\n";
            if ($panel eq "a"){
                ## includes everything with 10 hits or more
                if ($n_unique>=10 && $n_good>=1){
                    my $clean_uvig=$uvig;
                    if ($uvig=~/IMGVR/){my @t=split(/\|/,$uvig); $clean_uvig=$t[0];}
                    my @t=split(";",$tax);
                    if ($t[1] eq "p__unclassified"){}
                    else{
                        $t[1]=~s/p__//;
                        $t[2]=~s/g__//;
                        my $clean_tax=$t[1].";".$t[2];
                        print $sedge $clean_uvig."\t".$clean_tax."\t".$n_unique."\t".$ratio."\t".$emp_cat."\n";
                        $check_uvig{$uvig}=$clean_uvig;
                        $check_tax{$tax}=$clean_tax;
                    }
                }
            }
            elsif($panel eq "b"){
                ## includes everything with 10 hits or more and covering at least 1%
                ## Also: no unclassified genus
                if ($n_unique>=10 && $n_good>=1 && $ratio>=1){
                    my $clean_uvig=$uvig;
                    if ($uvig=~/IMGVR/){my @t=split(/\|/,$uvig); $clean_uvig=$t[0];}
                    my @t=split(";",$tax);
                    if ($t[1] eq "p__unclassified" || $t[2] eq "g__unclassified"){}
                    else{
                        $t[1]=~s/p__//;
                        $t[2]=~s/g__//;
                        my $clean_tax=$t[1].";".$t[2];
                        print $sedge $clean_uvig."\t".$clean_tax."\t".$n_unique."\t".$ratio."\t".$emp_cat."\n";
                        $check_uvig{$uvig}=$clean_uvig;
                        $check_tax{$tax}=$clean_tax;
                    }
                }
            }
        }
    }
}
close $s1;
close $sedge;

open my $snode,">",$out_file_net_nodes;
print $snode "id\tfull_id\ttype\ttax\n";
foreach my $uvig (keys %check_uvig){
    print $snode $check_uvig{$uvig}."\t".$uvig."\tvirus\tvirus\n";
}
foreach my $tax (keys %check_tax){
    my @t=split(";",$tax);
    print $snode $check_tax{$tax}."\t".$tax."\thost\t".$t[0].";".$t[1]."\n";
}
close $snode;

sub clean(){
    my @t=split(";",$_[0]);
    my $clean=$t[0].";".$t[1].";".$t[5];
    return $clean;
}

sub run_cmd(){
    my $cmd=$_[0];
    if ($_[1] ne "veryquiet"){print "$cmd\n";}
    if ($_[1] ne "dummy"){
            my $out=`$cmd`;
            if ($_[1] ne "quiet" && $_[1] ne "veryquiet"){
                    if ($_[1] eq "stderr"){print STDERR "$out\n";}
                    else{print "$out\n";}
            }
            return($out);
    }
    else{
            print " ### dummy run\n";
    }
}


sub get_cover_length(){
    my $hash=$_[0];
    my $total=0;
    foreach my $key (sort {$a <=> $b} keys %{$hash}){$total+=$$hash{$key}-$key+1;}
    return $total;
}

sub get_stats(){
    my $hash=$_[0];
    my $n_unique=scalar(keys %{$hash});
    # print "\t$n_unique unique spacers\n";
    my @tab_hits=();
    my $max=0;
    foreach my $key (keys %{$hash}){
        push(@tab_hits,$$hash{$key});
        # print "\t\t$key matched the target on $$hash{$key} loci\n";
        if ($$hash{$key}>$max){$max=$$hash{$key};}
    }
    my $median=&median(\@tab_hits);
    return ($n_unique,$median,$max);
}

sub median(){
        my($data) = @_;
        if(@$data == 1){
            return $$data[0];
        }
        my $med="NA";
        @$data=sort {$a <=> $b} (@$data);
        if (scalar(@$data) % 2 ==0 ){
                $med=(@{$data}[scalar(@$data)/2]+@{$data}[scalar(@$data)/2-1])/2;
        }
        else{
                $med=@{$data}[int(scalar(@$data)/2)];
        }
        return $med;
}


sub add_to_store(){
    my $hash=$_[0];
    my $start=$_[1];
    my $end=$_[2];
    my $tag=0;
    # print "Adding $start - $end to the existing coverage\n";
    foreach my $key (sort {$a <=> $b} keys %{$hash}){ ## Need to do that in order, otherwise we may miss new spacers covering gaps between two existing windows
        if ($start>=$key && $start <=$$hash{$key}){
            ## Starting between $key and $$hash{$key}, maybe expand in 3'
            if ($end > $$hash{$key}){$$hash{$key}=$end;}
            $tag=1;
            $start=$key; ## New start for this one
            # print "New start - $start\n";
        }
        elsif($start<$key && $end>=$key){
            ## May update both start and end
            if ($end > $$hash{$key}){$$hash{$key}=$end;}
            if ($start < $key){
                $$hash{$start}=$$hash{$key};
                delete($$hash{$key});
            }
            $tag=1;
        }
    }
    if ($tag==0){
        $$hash{$start}=$end;
    }
}
