#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use File::Basename;
use Cwd;
my $h=0;
my $in_sp_to_target="";
my $n_cpu=6;
GetOptions ('help' => \$h, 'h' => \$h, 'i=s'=>\$in_sp_to_target, 'c=s'=>\$n_cpu);
if ($h==1 || $in_sp_to_target eq ""){ # If asked for help or did not set up any argument
	print "# Script to get the maximum distance (ANI) between targets for individual spacers
# Arguments :
# -i: sp_to_target file
# -c: number of threads (default = 6)
";
	die "\n";
}

my $in_imgvr="../../Data/Spacer_db/IMGVR_sequence_information_Oct17.tsv";
my $in_imgvr_fasta="IMGVR_all_nucleotides-high_confidence.fna"; ## Can be downloaded from IMG/VR v4 website
my $out_file=$in_sp_to_target;
$out_file=~s/\.tsv/_distance.tsv/;
if ($in_sp_to_target eq $out_file){
    die("Pblm, I could not guess the output file from $in_sp_to_target\n");
}

my $path_lzani="/tools/lz-ani-1.0.1_x64-linux/lz-ani";

$in_sp_to_target=~/\/([^\/]+)\.tsv/;
my $id_tmp=$1;
my $tmp_dir="/tmp/tmp_for_".$id_tmp;
&run_cmd("mkdir $tmp_dir");

my $i=0;
# vOTU information already included
print "Loading info about IMG/VR sequences\n";
my %check_hq;
open my $tsv,"<",$in_imgvr;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "UVIG"){next;}
    if ($tab[5] eq "High-quality" || $tab[5] eq "Reference"){$check_hq{$tab[0]}=1;}
    $i++;
    if ($i % 100000 == 0){print " .. ".$i."\n";}
}
close $tsv;

print "Loading relevant IMG/VR sequences\n";
my %store_seq;
my $c_c="";
my $tag=0;
$i=0;
open my $fna,"<",$in_imgvr_fasta;
while(<$fna>){
    chomp($_);
    if ($_=~/^>(\S+)/){
        $c_c=$1;
        $tag=0;
        my @t=split(/\|/,$c_c);
        if ($check_hq{$t[0]}==1){
            $tag=1;
        }
        $i++;
        if ($i % 10000 == 0){print " .. ".$i."\n";}
    }
    elsif($tag==1){
        $store_seq{$c_c}.=$_;
    }
}
close $fna;

print "Now looking at individual spacer clusters\n";
my %tmp;
$i=0;
my $c_clust="";
my $slog="";
open my $s1,">",$out_file;
print $s1 "sp_cluster_id\tn_uvig_target\tn_votu_target\tn_hq_uvig_target\tn_hq_votu_target\tmax_distance_between_hq_target\tn_missing_hq\n";
open my $tsv,"<",$in_sp_to_target;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "sp_cluster_id"){next;}
    if ($c_clust ne $tab[0]){
        if ($c_clust ne ""){
            print "Processing data acquired for $c_clust\n";
            &process(\%tmp,$slog,$tmp_dir,$s1);
        }
        $c_clust=$tab[0];
        %tmp=();
    }
    $tmp{"by_uvig"}{$tab[1]}=1;
    $tmp{"by_votu"}{$tab[2]}=1;
    my @t=split(/\|/,$tab[1]);
    if ($check_hq{$t[0]}==1){
        # print $t[0]." is hq\n";
        $tmp{"by_hq_uvig"}{$tab[1]}=1;
        $tmp{"by_hq_votu"}{$tab[2]}=1;
    }
    else{
        # print $t[0]." is not hq\n";
    }
    $i++;
    if ($i % 1000000 == 0){print " .. ".$i."\n";}
}
close $tsv;
if ($c_clust ne ""){
    print "Processing data acquired for $c_clust\n";
    &process(\%tmp,$slog,$tmp_dir,$s1);
}
close $s1;

sub process(){
    my $hash=$_[0];
    my $sout=$_[1];
    my $tmp_dir=$_[2]."/";
    my $stab=$_[3];
    my $n_uvig=scalar(keys %{$$hash{"by_uvig"}});
    my $n_votu=scalar(keys %{$$hash{"by_votu"}});
    my @list_hq_uvigs=keys %{$$hash{"by_hq_uvig"}};
    my $n_hq_uvig=scalar(@list_hq_uvigs);
    my $n_hq_votu=scalar(keys %{$$hash{"by_hq_votu"}});
    my $max_dist=0;
    my $n_missing=-1;
    if ($n_hq_uvig<=1){
        print "skipping this one, 0 or 1 hq_uvig\n";
        $max_dist="NA";
        $n_missing="NA";
    }
    else{
        ## Preparing the fasta file for which we want lz-ani distance
        my $tmp_fna=$tmp_dir."input.fna";
        open my $s1,">",$tmp_fna;
        foreach my $hq_uvig (@list_hq_uvigs){
            print $s1 ">".$hq_uvig."\n";
            if (!defined($store_seq{$hq_uvig})){
                die("Missing sequence for $hq_uvig\n");
            }
            print $s1 $store_seq{$hq_uvig}."\n";
        }
        close $s1;
        my $tmp_dist=$tmp_dir."output.dist";
        &run_cmd("$path_lzani all2all --in-fasta $tmp_fna --out $tmp_dist --out-format lite -t $n_cpu -V 5");
        my $test="";
        open my $tsv,"<",$tmp_dist;
        while(<$tsv>){
            chomp($_);
            my @tab=split("\t",$_);
            if ($tab[0] eq "idx1"){next;}
            # print $_."\n";
            my $dist=1-$tab[3];
            my $test=$tab[0]."_".$tab[1];
            my @tab_2=split("\t",<$tsv>); ## Get the next line, as this should be the reverse, i.e. seq 2 vs sseq 1, and we take the lowest distance between the two
            my $test_2=$tab_2[1]."_".$tab_2[0];
            if ($test ne $test_2){
                print join("\t",@tab)."\n";
                print join("\t",@tab_2)."\n";
                print $test."\n";
                print $test_2."\n";
                die("Pblm, we had two lines that should be the same pair but were not\n");
            }
            my $dist_2=1-$tab_2[3];
            if ($dist_2 < $dist){
                $dist=$dist_2; ## We take the lowest dist, i.e. ANI with respect to the shortest sequence
            }
            if ($dist>$max_dist){
                $max_dist=$dist;
            }
        }
        close $tsv;
        print "Final max dist: $max_dist\n";
    }
    print "\tn uvigs target: ".$n_uvig."\n";
    print "\tn votu target: ".$n_votu."\n";
    print "\tn uvigs target high-quality: ".$n_hq_uvig."\n";
    print "\tn votu target high-quality: ".$n_hq_votu."\n";
    print "\tmax distance between targets (hq only): ".$max_dist."\n";
    print "\tnumber of pairs without an ANI (very low ?): ".$n_missing."\n";
    print $stab $c_clust."\t".$n_uvig."\t".$n_votu."\t".$n_hq_uvig."\t".$n_hq_votu."\t".$max_dist."\t".$n_missing."\n";
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
