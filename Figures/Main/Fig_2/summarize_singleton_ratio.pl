#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h); # , 'i=s'=>\$in_file, 'r=s'=>\$out_dir_root);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to get a master table with all singleton info for all repeat/sample combinations, i.e. all sets
# Arguments :
# none
";
    die "\n";
}

my $out_file="Singleton_ratio_by_set.tsv";
my $out_file_2="input_figure_singletons.tsv";

my %store;
foreach my $in_file (<Analyses/Spacer_database/*per_sample.tsv>){
    print $in_file." .. \n";
    my $lvl="NA";
    if ($in_file=~/80_per_sample/){$lvl="n_sing_80";}
    elsif ($in_file=~/95_per_sample/){$lvl="n_sing_95";}
    ## Percentage of singletons at different levels
    open my $tsv,"<",$in_file;
    while(<$tsv>){
        chomp($_);
        my @tab=split("\t",$_);
        if ($tab[0] eq "sample"){next;}
        if (!defined($store{$tab[1]}{$tab[0]})){
            $store{$tab[0]}{$tab[1]}{"n_hq_spacer"}=$tab[2];
            $store{$tab[0]}{$tab[1]}{"max_cover"}=$tab[3];
            $store{$tab[0]}{$tab[1]}{"n_sing_100"}=$tab[4];
        }
        $store{$tab[0]}{$tab[1]}{$lvl}=$tab[5];
    }
    close $tsv;
}

open my $s1,">",$out_file;
print $s1 "repeat\tsample\ttotal_spacer\tmax_cover\tratio_s_100\tratio_s_95\tratio_s_80\n";
foreach my $repeat (sort keys %store){
    foreach my $sample (sort keys %{$store{$repeat}}){
        print $s1 $repeat."\t".$sample."\t".$store{$repeat}{$sample}{"n_hq_spacer"}."\t".$store{$repeat}{$sample}{"max_cover"}.
        "\t".$store{$repeat}{$sample}{"n_sing_100"}/$store{$repeat}{$sample}{"n_hq_spacer"}.
        "\t".$store{$repeat}{$sample}{"n_sing_95"}/$store{$repeat}{$sample}{"n_hq_spacer"}.
        "\t".$store{$repeat}{$sample}{"n_sing_80"}/$store{$repeat}{$sample}{"n_hq_spacer"}."\n";
    }
}
close $s1;

my @ids=("id_100","id_95","id_80");

my %store=();
my %stats=();
open my $tsv,"<",$out_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "repeat"){next;}
    $stats{$tab[2]}{"total"}++;
    $stats{$tab[2]}{"id_100"}{$tab[4]}++;
    $stats{$tab[2]}{"id_95"}{$tab[5]}++;
    $stats{$tab[2]}{"id_80"}{$tab[6]}++;
}
close $tsv;

open my $s1,">",$out_file_2;
print $s1 "total_spacer\tid\tmedian\tstderr_min\tstderr_max\tstdev_min\tstdev_max\n";
for (my $i=5;$i<=1000;$i++){
    if (!defined($stats{$i})){
        print "?!?! 0 sets with exactly $i spacers ?!\n";
    }
    foreach my $id (@ids){
        my $med=&weighted_median($stats{$i}{$id});
        my $stdev=&weighted_stdev($stats{$i}{$id});
        my $stderr=$stdev/sqrt($stats{$i}{"total"});
        print $s1 $i."\t".$id."\t".$med."\t".($med-$stderr)."\t".($med+$stderr)."\t".($med-$stdev)."\t".($med+$stdev)."\n";
    }
}
close $s1;




sub weighted_median {
        my $hash=$_[0];
        my $med="NA";
        my $verbose=0;
#         $verbose=1;
	my @t=sort {$a <=> $b} keys %{$hash};
	## Verify that all elements have a positive weight (or null);
	my $total_weight=0;
	foreach my $key (@t){
		if (defined($$hash{$key}) && $$hash{$key} ne "" && $$hash{$key}>=0){
			$total_weight+=$$hash{$key};
			if ($verbose==1){print $key."\t".$$hash{$key}."\n";}
		}
		else{
			print "######## THIS IS A PROBLEM - $key -> $$hash{$key}\n";
			die("\n");
		}
	}
	if ($verbose==1){print "Total weight: $total_weight\n";}
	my $th=$total_weight/2;
	if ($verbose==1){print "we're looking to go up to rank $th\n";}

	if ($#t==0){ ## Simple case with 1 element
		$med=$t[0];
	}
	elsif ($#t==1){ ## Simple case with 2 elements
		if ($$hash{$t[0]} == $$hash{$t[1]}){ ## If same weight, return the average
			$med=($t[0]+$t[1])/2;
		}
		elsif($$hash{$t[0]} > $$hash{$t[1]}){ ## Else take the larger weight
			$med=$t[0];
		}
		elsif($$hash{$t[1]} > $$hash{$t[0]}){
			$med=$t[1];
		}
		else{
			print "//// SHOULD NOT BE HERE\n";
		}
	}
	else{ ## We do the not super efficient but reliable "count up until you reach half of the total weight
		my $total;
		my $prev_key=0;
		for (my $i=0;$i<=$#t;$i++){
			$total+=$$hash{$t[$i]};
# 			print "$t[$i] goes to $total\n";
			if ($total>$th){
				if ($verbose==1){print "we are over $th, so we found our median, it's $t[$i]\n";}
				$med=$t[$i];
				last;
			}
			elsif($total==$th){
				if ($verbose==1){print "we are just on $th, so the median will be the average with the next key\n";}
				$med=($t[$i]+$t[$i+1])/2;
				last;
			}
		}
	}
	return($med)
}

sub weighted_stdev {
        my $data=$_[0];
        if (not $data) {
                die("Empty hash\n");
        }
        my $average = &weighted_average($data);
        my $sqtotal = 0;
        my $total_weight = 0;
        my $n_weight = 0;
        foreach my $key (keys %{$data}){
                $sqtotal += $$data{$key} * (($average-$key) ** 2);
                $total_weight += $$data{$key};
                $n_weight ++;
        }
        my $std = ($sqtotal / ($total_weight*(($n_weight-1)/$n_weight))) ** 0.5;
        return $std;
}



sub weighted_average {
        my $data=$_[0];
        if (not $data) {
                die("Empty hash\n");
        }
        my $total = 0;
        my $total_weight = 0;
        foreach my $key (keys %{$data}){
                $total += $key*$$data{$key};
                $total_weight += $$data{$key};
        }
        my $average = $total / $total_weight;
        return $average;
}

