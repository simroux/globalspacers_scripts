#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h); 
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to get plasmid to repeat numbers of different match types
# Arguments :
# run
";
    die "\n";
}

my $no_mis_sp="../Target_coverage/All_spacers_pr_complete.tsv";
my $mis_sp="All_spacers_pr_hq_2-3mis.tsv";

my $nomis_hits="../Target_coverage/All_hits_pr_complete.tsv";
my $mis_hits="All_hits_pr_hq_2-3mis.tsv";

my $targeting_file="../Target_coverage/plasmid_complete_coverage_by_spacer.tsv";
my $info_plasmid_file="../../../Data/Spacer_db/IMGVR_sequence_information_Oct17.tsv";

my $out_file="Plasmid_to_repeat_hits_profile.tsv";
my $out_file_step2="Plasmid_to_repeat_hits_profile-medium_with_cc.tsv";

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

my %check;
my %store_type;
print "Reading $targeting_file\n";
open my $tsv,"<",$targeting_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    ### Try to look at everything with 10 spacers or more, 200bp or more, i.e. nothing "low"
    if ($tab[5]>=10 && $tab[3]>=200){
        $check{$tab[0]}{$tab[1]}=1;
    }
}
close $tsv;


## Need this for correlation computations
my $xm=1.5;
my $xs=5;
my $xs_sq=sqrt(5);
my @xt=(-1.5,-0.5,0.5,1.5);

print "######### Preparing output\n";
open my $s1,">",$out_file;
print $s1 "plasmid\tarray\tn_hit_0\tn_hit_1\tn_hit_2\tn_hit_3\n";
open my $s2,">",$out_file_step2;
print $s2 "plasmid\tarray\tn_hit_0\tn_hit_1\tn_hit_2\tn_hit_3\ttotal\tcc\tslope\tcategory\n";
foreach my $plasmid (keys %count){
    foreach my $array (keys %{$count{$plasmid}}){
        if (!defined($count{$plasmid}{$array}{0})){$count{$plasmid}{$array}{0}=0;}
        if (!defined($count{$plasmid}{$array}{1})){$count{$plasmid}{$array}{1}=0;}
        if (!defined($count{$plasmid}{$array}{2})){$count{$plasmid}{$array}{2}=0;}
        if (!defined($count{$plasmid}{$array}{3})){$count{$plasmid}{$array}{3}=0;}
        my $line=$plasmid."\t".$array."\t".$count{$plasmid}{$array}{0}."\t".$count{$plasmid}{$array}{1}."\t".$count{$plasmid}{$array}{2}."\t".$count{$plasmid}{$array}{3};
        print $s1 $line."\n";
        if ($check{$plasmid}{$array}==1){
            my @t=($count{$plasmid}{$array}{0},$count{$plasmid}{$array}{1},$count{$plasmid}{$array}{2},$count{$plasmid}{$array}{3});
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
}
close $s1;
close $s2;