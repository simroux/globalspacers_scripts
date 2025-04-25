#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use File::Basename;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
	print "# Script to calculate how many positions are covered by spacers for each UViG vs array and sample
# Arguments :
# run
";
	die "\n";
}

my $info_file="../../../Data/Spacer_db/IMGVR_sequence_information_Aug26.tsv";
my $cl_to_array="All_spacers_vr_hq.tsv";
my $hit_file="All_hits_vr_hq.tsv";
my $global_cover_file="uvig_hq_coverage_by_spacer.tsv";

my $out_file="uvig_hq_coverage_by_spacer-bysample.tsv";
my $log_file="uvig_hq_coverage-bysample_list_pblm.txt";


print "Making the list of arrays and viruses we care about\n";
my %check_virus;
my %check_array;
open my $tsv,"<",$global_cover_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "uvig"){next;}
    if ($tab[4]>=10 || $tab[3]>=1000){ # 10% coverage or 1kb
        $check_virus{$tab[0]}{$tab[1]}=1;
        $check_array{$tab[1]}=1;
    }
}
close $tsv;


print "Loading cluster to array and sample link\n";
my %cl_to_sp;
my $i=0;
open my $tsv,"<",$cl_to_array;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if (defined($check_array{$tab[6]})){
        $cl_to_sp{$tab[0]}{$tab[6]}{$tab[8]}=1;
        $i++;
        if ($i % 1000000 == 0){
            print "... $i ...\n";
        }
    }
}
close $tsv;

$i=0;
my %store_cover;
my %store_n;
my $last=0;
print "Loading cover info\n";
open my $slog,">",$log_file;
open my $tsv,"<",$hit_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "cluster_id"){next;}
    if (defined($check_virus{$tab[1]}) && defined($cl_to_sp{$tab[0]})){ ## If we care about this virus and this spacer
        foreach my $array (keys %{$cl_to_sp{$tab[0]}}){
            if (defined($check_virus{$tab[1]}) && ($check_virus{$tab[1]}{$array}==1)){ ## And if we care specifically about this array for this virus
                foreach my $sample (keys %{$cl_to_sp{$tab[0]}{$array}}){
                    if (!defined($store_cover{$tab[1]}{$array}{$sample})){$store_cover{$tab[1]}{$array}{$sample}={};}
                    &add_to_store($store_cover{$tab[1]}{$array}{$sample},$tab[3],$tab[4]);
                    $store_n{$tab[1]}{$array}{$sample}{"n_unique"}{$tab[0]}++;
                    $store_n{$tab[1]}{$array}{$sample}{"number"}++;
                    $store_n{$tab[1]}{$array}{$sample}{"n_base"}+=($tab[4]-$tab[3])+1;
                    # print "Current coverage for $tab[1] is\n";
                    $last=0;
                    foreach my $key (sort {$a <=> $b} keys %{$store_cover{$tab[1]}{$array}{$sample}}){
                        # print "\t".$key."\t".$store_cover{$tab[1]}{$key}."\n";
                        if ($key<$last){
                            print $slog "## pblm after adding $tab[0] - $tab[1] - $tab[3] - $tab[4]\n";
                            print $slog "We have $key but the last stop was $last\n";
                        }
                        $last=$store_cover{$tab[1]}{$array}{$sample}{$key};
                    }
                }
            }
        }
        $i++;
        if ($i % 100000 == 0){
            print "... $i ...\n";
        }
        # if ($i==10000000){last;}
    }
}
close $tsv;
close $slog;

my %info_uvig;
print "Loading uvig info and calculating coverage\n";
open my $s1,">",$out_file;
print $s1 "uvig\tarray\tsample\tlength\tcovered\tpercent\tn_unique_spacer\tmedian_hit_per_spacer\tmax_hit\n";
open my $tsv,"<",$info_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "uvig"){next;}
    # $info_uvig{$tab[0]}{"length"}=$tab[3];
    if (defined($store_cover{$tab[0]})){
        foreach my $array (keys %{$store_cover{$tab[0]}}){
            foreach my $sample (keys %{$store_cover{$tab[0]}{$array}}){
                my $covered=&get_cover_length($store_cover{$tab[0]}{$array}{$sample});
                my $ratio=sprintf("%.02f",$covered/$tab[3]*100);
                my ($n_unique,$median_hit,$max_hit)=&get_stats($store_n{$tab[0]}{$array}{$sample}{"n_unique"});
                # print "Preparing stats for $tab[0] // $array\n";
                if ($covered > 0){
                    print $tab[0]."\t".$array."\t".$sample."\t".$tab[3]."\t".$covered."\t".$ratio."\t".$n_unique."\t".$median_hit."\t".$max_hit."\n";
                    print $s1 $tab[0]."\t".$array."\t".$sample."\t".$tab[3]."\t".$covered."\t".$ratio."\t".$n_unique."\t".$median_hit."\t".$max_hit."\n";
                }
            }
        }
    }
}
close $tsv;
close $s1;

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
	if (scalar(@$data) % 2 == 0 ){
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
