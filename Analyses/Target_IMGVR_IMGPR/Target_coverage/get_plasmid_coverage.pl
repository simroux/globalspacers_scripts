#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use File::Basename;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
	print "# Script to calculate how many positions are covered by spacers for each (predicted complete) plasmid
# Arguments :
# run
";
	die "\n";
}

my $info_file="../../../Data/Spacer_db/IMGPR_sequence_information_Aug26.tsv";
my $cl_to_array="All_spacers_pr_complete.tsv";
my $hit_file="All_hits_pr_complete.tsv";

my $out_file="plasmid_complete_coverage_by_spacer.tsv";
my $log_file="plasmid_complete_coverage_list_pblm.txt";

print "Loading cluster to array link\n";
my %cl_to_sp;
my $i=0;
open my $tsv,"<",$cl_to_array;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    $cl_to_sp{$tab[0]}{$tab[6]}=1;
    $i++;
    if ($i % 1000000 == 0){
        print "... $i ...\n";
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
    if (!defined($store_cover{$tab[1]})){$store_cover{$tab[1]}={};}
    # print "Add to store: $tab[2]/$tab[3]\n";
    foreach my $array (keys %{$cl_to_sp{$tab[0]}}){
        if (!defined($store_cover{$tab[1]}{$array})){$store_cover{$tab[1]}{$array}={};}
        &add_to_store($store_cover{$tab[1]}{$array},$tab[3],$tab[4]);
        $store_n{$tab[1]}{$array}{"n_unique"}{$tab[0]}++;
        $store_n{$tab[1]}{$array}{"number"}++;
        $store_n{$tab[1]}{$array}{"n_base"}+=($tab[4]-$tab[3])+1;
        # print "Current coverage for $tab[1] is\n";
        $last=0;
        foreach my $key (sort {$a <=> $b} keys %{$store_cover{$tab[1]}{$array}}){
            # print "\t".$key."\t".$store_cover{$tab[1]}{$key}."\n";
            if ($key<$last){
                print $slog "## pblm after adding $tab[0] - $tab[1] - $tab[3] - $tab[4]\n";
                print $slog "We have $key but the last stop was $last\n";
                
            }
            $last=$store_cover{$tab[1]}{$array}{$key};
        }
    }
    $i++;
    if ($i % 100000 == 0){
        print "... $i ...\n";
    }
}
close $tsv;
close $slog;

my %info_uvig;
print "Loading uvig info and calculating coverage\n";
open my $s1,">",$out_file;
print $s1 "plasmid\tarray\tlength\tcovered\tpercent\tn_unique_spacer\tmedian_hit_per_spacer\tmax_hit\n";
open my $tsv,"<",$info_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "full_plasmid_id"){next;}
    if (defined($store_cover{$tab[0]})){
        foreach my $array (keys %{$store_cover{$tab[0]}}){
            my $covered=&get_cover_length($store_cover{$tab[0]}{$array});
            my $ratio=sprintf("%.02f",$covered/$tab[8]*100);
            my ($n_unique,$median_hit,$max_hit)=&get_stats($store_n{$tab[0]}{$array}{"n_unique"});
            # print "Preparing stats for $tab[0] // $array\n";
            if ($covered > 0){
                print $tab[0]."\t".$array."\t".$tab[8]."\t".$covered."\t".$ratio."\t".$n_unique."\t".$median_hit."\t".$max_hit."\n";
                print $s1 $tab[0]."\t".$array."\t".$tab[8]."\t".$covered."\t".$ratio."\t".$n_unique."\t".$median_hit."\t".$max_hit."\n";
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
    my @tab_hits=();
    my $max=0;
    foreach my $key (keys %{$hash}){
        push(@tab_hits,$$hash{$key});
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
    print "Adding $start - $end to the existing coverage\n";
    foreach my $key (sort {$a <=> $b} keys %{$hash}){ ## Need to do that in order, otherwise we may miss new spacers covering gaps between two existing windows
        if ($start>=$key && $start <=$$hash{$key}){
            ## Starting between $key and $$hash{$key}, maybe expand in 3'
            if ($end > $$hash{$key}){$$hash{$key}=$end;}
            $tag=1;
            $start=$key; ## New start for this one
            print "New start - $start\n";
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
