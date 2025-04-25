#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
use Custom::Utils;
my $h=0;
my $tag_multitaxa=0;
GetOptions ('help' => \$h, 'h' => \$h); # , 'i=s'=>\$in_file, 'r=s'=>\$out_dir_root);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to get uvig to repeat numbers of different match types, for multitaxa uvigs
# Arguments :
# run
";
    die "\n";
}

my $info_virus_file="../../../Data/Spacer_db/IMGVR_sequence_information_Oct17.tsv";

my $no_mis_sp="../Target_coverage/All_spacers_vr_additional_multitaxa.tsv";
my $mis_sp="All_spacers_vr_additional_multitaxa_2-3mis.tsv";

my $nomis_hits="../Target_coverage/All_hits_vr_additional_multitaxa.tsv";
my $mis_hits="All_hits_vr_additional_multitaxa_2-3mis.tsv";

my $out_file="Multitaxa_virus_to_array_hits_profile.tsv";
my $out_file_step2="Multitaxa_virus_to_array_hits_profile.with_cc.tsv";

my %sp_to_array;
my %count;
my $i=0;
print "Reading $mis_sp\n";
open my $tsv,"<",$mis_sp;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "cluster_id"){next;}
    $sp_to_array{$tab[0]}{$tab[6]}=1;
    $i++;
    if ($i % 1000000 == 0){print ".. $i ..\n";}
}
close $tsv;
print "Reading $mis_hits\n";
$i=0;
open my $tsv,"<",$mis_hits;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "cluster_id"){next;}
    foreach my $array (keys %{$sp_to_array{$tab[0]}}){
        $count{$tab[1]}{$array}{$tab[6]}++;
    }
    $i++;
    if ($i % 1000000 == 0){print ".. $i ..\n";}
}
close $tsv;

%sp_to_array=();
## Now for 0 and 1 mismatches
print "Reading $no_mis_sp\n";
$i=0;
open my $tsv,"<",$no_mis_sp;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "cluster_id"){next;}
    $sp_to_array{$tab[0]}{$tab[6]}=1;
    $i++;
    if ($i % 1000000 == 0){print ".. $i ..\n";}
}
close $tsv;
$i=0;
print "Reading $nomis_hits\n";
open my $tsv,"<",$nomis_hits;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "cluster_id"){next;}
    foreach my $array (keys %{$sp_to_array{$tab[0]}}){
        $count{$tab[1]}{$array}{$tab[6]}++;
    }
    $i++;
    if ($i % 1000000 == 0){print ".. $i ..\n";}
}
close $tsv;

## Need this for correlation computations
my $xm=1.5;
my $xs=5;
my $xs_sq=sqrt(5);
my @xt=(-1.5,-0.5,0.5,1.5);

print "######### Preparing output\n";
open my $s1,">",$out_file;
print $s1 "uvig\tarray\tn_hit_0\tn_hit_1\tn_hit_2\tn_hit_3\n";
open my $s2,">",$out_file_step2;
print $s2 "uvig\tarray\tn_hit_0\tn_hit_1\tn_hit_2\tn_hit_3\ttotal\tcc\tslope\tcategory\n";
foreach my $virus (keys %count){
    foreach my $array (keys %{$count{$virus}}){
        if (!defined($count{$virus}{$array}{0})){$count{$virus}{$array}{0}=0;}
        if (!defined($count{$virus}{$array}{1})){$count{$virus}{$array}{1}=0;}
        if (!defined($count{$virus}{$array}{2})){$count{$virus}{$array}{2}=0;}
        if (!defined($count{$virus}{$array}{3})){$count{$virus}{$array}{3}=0;}
        my $line=$virus."\t".$array."\t".$count{$virus}{$array}{0}."\t".$count{$virus}{$array}{1}."\t".$count{$virus}{$array}{2}."\t".$count{$virus}{$array}{3};
        print $s1 $line."\n";
        my @t=($count{$virus}{$array}{0},$count{$virus}{$array}{1},$count{$virus}{$array}{2},$count{$virus}{$array}{3});
        my $sum=$t[0]+$t[1]+$t[2]+$t[3];
        my $avg=$sum/4;
        if ($sum>=10){
            my $num=0;
            my $den=0;
            for (my $i=0;$i<=3;$i++){
                $num+=($t[$i]-$avg)*$xt[$i];
                $den+=($t[$i]-$avg)*($t[$i]-$avg);
            }
            $den=sqrt($den)*$xs_sq;
            my $cc=0;
            if ($den>0){$cc=$num/$den;}
            my $slope=$num/$xs;
            my $cat="NA";
            if ($cc>=0.7){$cat="positive";}
            elsif($cc<=-0.7){$cat="negative";}
            # print join(" ",@t)."\t".$cc."\t".$slope."\n";
            $line=$line."\t".$sum."\t".$cc."\t".$slope."\t".$cat;
            print $s2 $line."\n";
        }
    }
}
close $s1;
close $s2;