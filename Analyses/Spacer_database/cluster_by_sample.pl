#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
my $batch_n=-1;
GetOptions ('help' => \$h, 'h' => \$h, 'b=s'=>\$batch_n); 
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to cluster high-quality spacers at 95% (or 80%) for each set
# Arguments :
# none
";
    die "\n";
}

## Toggle this line to do 95 or 80% clustering
# my $min_id="95";
my $min_id="80";


my $in_file="Data/Spacer_db/All_spacers_info_filtered-Jul19-24.tsv";
my $out_file="Spacer_clustering_at_".$min_id."_per_sample.tsv";
my $tmp_dir="/tmp/";

my %check_sample;
if ($batch_n!=-1){
    my $batch_file="batch_".$batch_n.".tsv";
    $out_file="batch_".$batch_n."_clustering_at_".$min_id."_per_sample.tsv";
    # $prev_out_file="Batches/batch_".$batch_n."_clustering_at_".$min_id."_per_sample.tsv";
    $tmp_dir="/tmp/batch_".$batch_n."/";
    if (!(-d $tmp_dir)){&run_cmd("mkdir $tmp_dir");}
    print "Reading $batch_file for $batch_n, results will now go into $out_file, and tmp dir is $tmp_dir\n";
    open my $tsv,"<",$batch_file;
    while(<$tsv>){
        chomp($_);
        $check_sample{$_}=1;
        # print "Checking sample $_\n";
    }
    close $tsv;
}

my %done;

my $c_c="";
my %tmp;
my %seen;
my $i=0;
open my $s1,">",$out_file;
print $s1 "sample\tarray\tn_hq_spacer\tmax_cover\tn_snr\tn_snr_s".$min_id."\tn_snr_c".$min_id."\tn_snr_c".$min_id."-sonly\tn_snr_c".$min_id."-nons\n";
open my $tsv,"<",$in_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "spacer_id"){next;}
    if ($batch_n==-1 || $check_sample{$tab[7]}==1){
        my $code=$tab[7]; ## Samples are ordered, but arrays are not. The code is only the sample
        if ($c_c eq ""){$c_c=$code;}
        if ($code ne $c_c){
            print "Processing for $code\n";
            foreach my $array (sort keys %tmp){
                if ($tmp{$array}{"n_hq"}>1){
                    print "\tarray $array\n";
                    if ($done{$c_c}{$array}==1){
                        print "$c_c and $array already processed\n";
                    }
                    else{
                        my ($sg_s,$sg_c,$sg_c_os,$sg_c_ns)=&process($tmp{$array}{"spacers"});
                        print $s1 $c_c."\t".$array."\t".$tmp{$array}{"n_hq"}."\t".$tmp{$array}{"max_cover"}."\t".$tmp{$array}{"n_singleton"}."\t".$sg_s."\t".$sg_c."\t".$sg_c_os."\t".$sg_c_ns."\n";
                    }
                }
            }
            %tmp=();
            $c_c=$code;
            if ($seen{$code}==1){
                die("we already saw this code $code - That's a problem\n");
            }
            $seen{$code}=1;
        }
        if ($tab[6]==1){
            $tmp{$tab[5]}{"n_hq"}++;
            if ($tab[3]==1){$tmp{$tab[5]}{"n_singleton"}++;}
            $tmp{$tab[5]}{"spacers"}{$tab[0]}{"seq"}=$tab[1];
            $tmp{$tab[5]}{"spacers"}{$tab[0]}{"cover"}=$tab[3];
            if (!defined($tmp{$tab[5]}{"max_cover"})){$tmp{$tab[5]}{"max_cover"}=$tab[3];}
            elsif ($tab[3] > $tmp{$tab[5]}{"max_cover"}){$tmp{$tab[5]}{"max_cover"}=$tab[3];}
        }
    }
    else{
        # print $tab[7]." is of no interest\n";
    }
}
close $tsv;
## Process the last
foreach my $array (sort keys %tmp){
    if ($tmp{$array}{"n_hq"}>1){
        if ($done{$c_c}{$array}==1){
            print "$c_c and $array already processed\n";
        }
        else{
            print "\tarray $array\n";
            my ($sg_s,$sg_c,$sg_c_os,$sg_c_ns)=&process($tmp{$array}{"spacers"});
            print $s1 $c_c."\t".$array."\t".$tmp{$array}{"n_hq"}."\t".$tmp{$array}{"max_cover"}."\t".$tmp{$array}{"n_singleton"}."\t".$sg_s."\t".$sg_c."\t".$sg_c_os."\t".$sg_c_ns."\n";
            # <STDIN>;
        }
    }
}
close $s1;

sub process(){
    my $hash=$_[0];
    my %stats;
    # print "\tcreating a fasta file\n";
    my $tmp_fasta=$tmp_dir."tmp.fasta";
    open my $s1,">",$tmp_fasta;
    foreach my $sp (keys %{$hash}){
        print $s1 ">".$sp."\n".$$hash{$sp}{"seq"}."\n";
        $stats{"singletons_total"}++;
    }
    close $s1;
    # print "\trunning cd-hit\n";
    my $tmp_nr=$tmp_dir."tmp.nr.fasta";
    my $tmp_clstr=$tmp_dir."tmp.nr.fasta.clstr";
    &run_cmd("cd-hit-est -i $tmp_fasta -o $tmp_nr -c 0.$min_id -T 4 -d 0","veryquiet");
    my $tmp_txt=$tmp_dir."tmp.nr.fasta.clstr.txt";
    &run_cmd("clstr2txt.pl $tmp_clstr > $tmp_txt","veryquiet");
    # print "\tnow parsing\n";
    open my $tsv,"<",$tmp_txt;
    while(<$tsv>){
        chomp($_);
        my @tab=split("\t",$_);
        if ($tab[0] eq "id"){next;}
        if ($tab[2]==1){
            ## Singleton even after clustering
            if ($$hash{$tab[0]}{"cover"}==1){$stats{"singletons_not_clustered"}++;}
        }
        else{
            if ($$hash{$tab[0]}{"cover"}==1){$stats{"singletons_clustered"}++;} ## Here, just singleton clustered vs unclustered, we don't care if this is a group of all singletons, or if it clusters with a more highly covered spacer
            if ($$hash{$tab[0]}{"cover"}>1){
                $stats{"clusters"}{$tab[1]}{"non_singleton"}++;
            }
            else{
                $stats{"clusters"}{$tab[1]}{"singletons"}++;
            }
        }
    }
    close $tsv;
    foreach my $cl (keys %{$stats{"clusters"}}){
        if (defined($stats{"clusters"}{$cl}{"non_singleton"}) && ($stats{"clusters"}{$cl}{"non_singleton"}>0)){
            $stats{"singletons_clustered-nonsingleton"}+=$stats{"clusters"}{$cl}{"singletons"}; ## All singletons in this cluster were clustered with a non-singleton
        }
        else{
            $stats{"singletons_clustered-singleton_only"}+=$stats{"clusters"}{$cl}{"singletons"};
        }
    }
    if (!defined($stats{"singletons_not_clustered"})){$stats{"singletons_not_clustered"}=0;}
    if (!defined($stats{"singletons_clustered"})){$stats{"singletons_clustered"}=0;}
    if (!defined($stats{"singletons_clustered-singleton_only"})){$stats{"singletons_clustered-singleton_only"}=0;}
    if (!defined($stats{"singletons_clustered-nonsingleton"})){$stats{"singletons_clustered-nonsingleton"}=0;}
    return($stats{"singletons_not_clustered"},$stats{"singletons_clustered"},$stats{"singletons_clustered-singleton_only"},$stats{"singletons_clustered-nonsingleton"});
}


sub run_cmd {
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

