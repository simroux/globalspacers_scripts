#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h); 
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to summarize sets for main figure
# Arguments :
# run
";
    die "\n";
}

my $in_file="Analyses/Extra/SRA_metadata/Overview_time_series.txt";
my $out_file="spacer_sets_betadiv_forplot_within-individuals.tsv";
my $n_sub=10;
my $min_spacer=10; ## min spacer per set to consider
my $min_maxcover=20; ## min max cover spacer to consider per set

my %info_series;
my %check;
open my $tsv,"<",$in_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "series_name"){next;}
    if ($tab[4]>=10){
        $check{$tab[0]}=1;
        $info_series{$tab[0]}{"eco"}=$tab[3];
    }
}
close $tsv;

my %info_lib;
my %store;
my %stats;
my %tmp;
my $tag=0;
open my $s1,">",$out_file;
print $s1 "series\tecosystem\tarray\tn_samples\tcount_all\tpercentage_all\tcount_common\tpercentage_common\tcount_sgton\tpercentage_sgton\n";
foreach my $series_set (keys %check){
    my $in_file="Analyses/Spacer_database/Individual_time_series/".$series_set."/".$series_set."_filt_spacers.tsv";
    if (!(-e $in_file)){next;}
    print "We look into $in_file\n";
    # 
    %info_lib=();
    %store=();
    open my $tsv,"<",$in_file;
    while(<$tsv>){
        chomp($_);
        my @tab=split("\t",$_);
        if ($tab[0] eq "spacer_id"){next;}
        my $spacer_id=$tab[0];
        my $array=$tab[5];
        my $library=$tab[7];
        if (!defined($info_lib{$library})){
            $info_lib{$library}{"date"}=$tab[10];
            my @t=split("-",$tab[10]);
            $info_lib{$library}{"datepseudodays"}=$t[0]*365+$t[1]*30+$t[2]; ## it's fake, but it should do the work for sorting and getting a rough sense of time between timepoints, i.e. 0 or >0 days
        }
        my $cover=$tab[3];
        my $sp_cluster=$tab[11];
        $store{$array}{$library}{$sp_cluster}{"coverage"}=$cover;
    }
    close $tsv;
    ##### 
    ## Get a list of all series_set for this subject
    my @lib_list=sort {$info_lib{$a}{"datepseudodays"} <=> $info_lib{$b}{"datepseudodays"} or $a cmp $b} keys %info_lib;
    my $n_lib=scalar(@lib_list);
    ## Prep some stats to determine if array is relevant
    %stats=();
    %tmp=();
    $tag=0;
    foreach my $array (sort keys %store){
        print "## checking array $array\n";
        my $line="";
        $tag=0;
        %stats=();
        %tmp=();
        for (my $i=0;$i<=$#lib_list;$i++){
            my $lib=$lib_list[$i];
            $line.="\t".$lib;
            if (defined($store{$array}{$lib})){
                $line.="\t".scalar(keys %{$store{$array}{$lib}})."\t";
                # my $total_cover=0;
                # my $total_sp=0;
                my @list_spacers=sort {$store{$array}{$lib}{$b}{"coverage"} <=> $store{$array}{$lib}{$a}{"coverage"} or $a cmp $b} keys %{$store{$array}{$lib}};
                my $th_cover=$store{$array}{$lib}{$list_spacers[0]}{"coverage"}*0.5;
                if ($store{$array}{$lib}{$list_spacers[0]}{"coverage"}<$min_maxcover || scalar(@list_spacers)<$min_spacer){
                    foreach my $spacer (@list_spacers){
                        # print "\t".$array."\t".$lib."\t".$spacer."\n";
                        $line.=$spacer." (".$store{$array}{$lib}{$spacer}{"coverage"}.") ";
                    }
                    chop($line);
                    ## Max cover is low, we only prepare for the output  
                } 
                else{
                    foreach my $spacer (@list_spacers){
                        # print "\t".$array."\t".$lib."\t".$spacer."\n";
                        $line.=$spacer." (".$store{$array}{$lib}{$spacer}{"coverage"}.") ";
                        # $total_sp++;
                        # $total_cover+=$store{$array}{$lib}{$spacer}{"coverage"};
                        $stats{$spacer}{"seen"}{$lib}=1;
                        if ($store{$array}{$lib}{$spacer}{"coverage"}>=$th_cover){
                            $stats{$spacer}{"common"}{$lib}=1;
                        }
                        elsif($store{$array}{$lib}{$spacer}{"coverage"}>1){
                            $stats{$spacer}{"not_common"}{$lib}=1;
                        }
                        elsif($store{$array}{$lib}{$spacer}{"coverage"}==1){
                            $stats{$spacer}{"singleton"}{$lib}=1;
                        }
                        $tmp{$info_lib{$lib}{"date"}}{$lib}=1;
                    }
                    chop($line);
                }
            }
            else{
                $line.="\t0\tNA";
                # $tag=1;
            }
            $line.="\n";
        }
        print $line."\n\n";
        ## Select the libs potentially interesting
        my @select_dates=keys %tmp;
        # print $series_set."\t".$array."\t".scalar(@select_lib)."\n";
        if (scalar(@select_dates)>=$n_sub){
            &fisher_yates_shuffle(\@select_dates);
            @select_dates=@select_dates[0..$n_sub-1];
            my @select_lib;
            foreach my $date (@select_dates){
                my @test=keys %{$tmp{$date}};
                &fisher_yates_shuffle(\@test);
                push(@select_lib,$test[0]);
            }
            my $n_lib=$n_sub;
            foreach my $spacer (sort keys %stats){
                $stats{$spacer}{"n_seen"}=0;
                $stats{$spacer}{"n_common"}=0;
                $stats{$spacer}{"n_not_common"}=0;
                $stats{$spacer}{"n_singleton"}=0;
                foreach my $lib (@select_lib){
                    if(defined($stats{$spacer}{"seen"}{$lib})){$stats{$spacer}{"n_seen"}++;}
                    if(defined($stats{$spacer}{"common"}{$lib})){$stats{$spacer}{"n_common"}++;}
                    if(defined($stats{$spacer}{"not_common"}{$lib})){$stats{$spacer}{"n_not_common"}++;}
                    if(defined($stats{$spacer}{"singleton"}{$lib})){$stats{$spacer}{"n_singleton"}++;}
                }
                ## 
                my $n_missing=$n_lib-$stats{$spacer}{"n_seen"};
                my $type="NA";
                if ($stats{$spacer}{"n_common"}>=1 || $stats{$spacer}{"n_not_common"}>=1){
                    $type="not_singleton";
                }
                elsif($stats{$spacer}{"singleton"}>=1){$type="singleton";}
                print $series_set."\t".$spacer."\t".$type."\t".$stats{$spacer}{"n_seen"}."\t".$n_missing."\t".$stats{$spacer}{"n_common"}."\t".$stats{$spacer}{"n_singleton"}."\n";
                if ($stats{$spacer}{"n_seen"}==0){} ## Ignore, spacer not seen in this sample subset
                else{
                    $stats{"summary"}{"total"}++;
                    $stats{"summary"}{$stats{$spacer}{"n_seen"}}{"all"}++;
                    if ($stats{$spacer}{"n_common"}>=1){
                        $stats{"summary"}{"total_common"}++;
                        $stats{"summary"}{$stats{$spacer}{"n_seen"}}{"common"}++;
                    }
                    if ($stats{$spacer}{"n_singleton"}>=1){
                        $stats{"summary"}{"total_singleton"}++;
                        $stats{"summary"}{$stats{$spacer}{"n_seen"}}{"singleton"}++;
                    }
                }
            }
            # <STDIN>;
            my $test=0;
            for (my $i=1;$i<=$n_sub;$i++){
                if (!defined($stats{"summary"}{$i}{"all"})){$stats{"summary"}{$i}{"all"}=0;}
                if (!defined($stats{"summary"}{$i}{"common"})){$stats{"summary"}{$i}{"common"}=0;}
                if (!defined($stats{"summary"}{$i}{"singleton"})){$stats{"summary"}{$i}{"singleton"}=0;}
                print $s1 $series_set."\t".$info_series{$series_set}{"eco"}."\t".$array."\t".$i."\t".$stats{"summary"}{$i}{"all"}."\t".$stats{"summary"}{$i}{"all"}/$stats{"summary"}{"total"}."\t".
                $stats{"summary"}{$i}{"common"}."\t".$stats{"summary"}{$i}{"common"}/$stats{"summary"}{"total_common"}."\t".
                $stats{"summary"}{$i}{"singleton"}."\t".$stats{"summary"}{$i}{"singleton"}/$stats{"summary"}{"total_singleton"}."\n";
                $test+=$stats{"summary"}{$i}{"all"}/$stats{"summary"}{"total"};
            }
            if ($test<0.9999999999999){
                print $stats{"summary"}{"total"}."\n";
                for (my $i=1;$i<=$n_sub;$i++){
                    print $stats{"summary"}{$i}{"all"}."\n";
                }
                print "PBLM -> WE ARE BELOW 100% ?\n";
                print $test."\n";
                die("\n");
            }
        }        
    }
}
close $s1;


