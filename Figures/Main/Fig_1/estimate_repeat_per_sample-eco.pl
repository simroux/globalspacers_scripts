#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
	print "# Script to estimate the median number of CRISPR arrays per metagenome, with resampling to get variability
# Arguments :
# run
";
	die "\n";
}

my $in_file="sample_to_array_counts.tsv";
my $out_file="sample_to_array_counts-summarized_by_ecosystem.tsv";
my $out_file_bis="sample_to_array_counts-summarized_by_ecosystem-detailed.tsv";


my %store;
my %store_bis;
my %translate;
print "Reading $in_file\n";
open my $tsv,"<",$in_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "sra_run"){next;}
    if ($tab[1] eq "NA"){$tab[1]=0;}
    my $n_base=$tab[11];
    if ($n_base < 500000000){ 
        next;
    }
    my $array_per_gb=$tab[1]/$tab[11]*1000000000;
    $store{$tab[5]}{$tab[2]}=$array_per_gb;
    $store_bis{$tab[4]}{$tab[2]}=$array_per_gb;
    if (!defined($translate{$tab[4]})){
        $translate{$tab[4]}=$tab[5];
    }
}
close $tsv;

my $n_batch=50;
my $s_batch=1000;
my @tmp;
my %store_results;
foreach my $eco (sort keys %store){
    print "## $eco ##\n";
    my @t=keys %{$store{$eco}};
    my $max=scalar(@t);
    if ($max<$s_batch){die("Pblm - $max vs $s_batch\n");}
    my @tab=();
    $store_results{$eco}=\@tab;
    print "Max is $max\n";
    for (my $i=0;$i<$n_batch;$i++){
        @tmp=();
        print "\tbatch $i\n";
        &fisher_yates_shuffle(\@t);
        for (my $k=0;$k<$s_batch;$k++){
            push(@tmp,$store{$eco}{$t[$k]});
        }
        my $median=&median(\@tmp);
        push(@{$store_results{$eco}},$median);
        print "\t\t"."\t".$median."\n";
    }
}

open my $s1,">",$out_file;
print $s1 "ecosystem\taverage_median\tmin_median\tmax_median\n";
foreach my $eco (sort keys %store_results){
    my @sorted=sort {$a <=> $b} @{$store_results{$eco}};
    my $average=&average(\@sorted);
    print $eco."\t".$average."\t".$sorted[0]."\t".$sorted[$#sorted]."\n";
    print $s1 $eco."\t".$average."\t".$sorted[0]."\t".$sorted[$#sorted]."\n";
}
close $s1;

$n_batch=50;
$s_batch=1000;
@tmp=();
%store_results=();
my %info_batch;
foreach my $eco (sort keys %store_bis){
    print "## $eco ##\n";
    my $c_batch=$s_batch;
    my @t=keys %{$store_bis{$eco}};
    my $max=scalar(@t);
    if ($max<$c_batch){
        print "$eco - $max -> trying the next $c_batch\n";
        if ($max > 500){
            $c_batch=500;
            print "using 500\n";
        }
        elsif ($max > 50){
            $c_batch=50;
            print "using 50\n";
        }
        else{
            die("Pblm - $max vs $s_batch\n");
        }
    }
    $info_batch{$eco}{"total_samples"}=$max;
    $info_batch{$eco}{"subsample_size"}=$c_batch;
    my @tab=();
    $store_results{$eco}=\@tab;
    for (my $i=0;$i<$n_batch;$i++){
        @tmp=();
        print "\tbatch $i\n";
        &fisher_yates_shuffle(\@t);
        for (my $k=0;$k<$c_batch;$k++){
            push(@tmp,$store_bis{$eco}{$t[$k]});
        }
        my $median=&median(\@tmp);
        push(@{$store_results{$eco}},$median);
        print "\t\t"."\t".$median."\n";
    }
}

open my $s1,">",$out_file_bis;
print $s1 "ecosystem_sum\tecosystem_detailed\tsubsample_size\taverage_median\tmin_median\tmax_median\n";
foreach my $eco (sort keys %store_results){
    my @sorted=sort {$a <=> $b} @{$store_results{$eco}};
    my $average=&average(\@sorted);
    print $eco."\t".$average."\t".$sorted[0]."\t".$sorted[$#sorted]."\n";
    print $s1 $translate{$eco}."\t".$eco."\t".$info_batch{$eco}{"subsample_size"}."\t".$average."\t".$sorted[0]."\t".$sorted[$#sorted]."\n";
}
close $s1;



sub fisher_yates_shuffle {
	my $deck = shift;  # $deck is a reference to an array
	return unless @$deck; # must not be empty!
	my $i = @$deck;
	while (--$i) {
		my $j = int rand ($i+1);
		@$deck[$i,$j] = @$deck[$j,$i];
	}
}

sub average {
        my($data) = @_;
        if (not @$data) {
                die("Empty arrayn");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}

sub stdev {
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}

sub median {
        my($data) = @_;
        if(@$data == 1){
			return $$data[0];
            #     return 0;
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