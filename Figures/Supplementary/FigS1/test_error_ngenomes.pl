#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h); # , 'i=s'=>\$in_file, 'r=s'=>\$out_dir_root);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to get a nice estimation of the taxonomic assignment error rate based on subsampling
# Arguments :
# none
";
    die "\n";
}


my $array_file="rData/Additional_data/epeats_to_genome_contigs.tsv";
my $repeat_file="Data/Spacer_db/Array_info_filtered_for_db-Nov1-24.tsv";
my $genome_file="Data/Additional_data/all_genome_info.nometag.tsv";
my $out_file_aggregated="aggregated_stats_subsampling-all.tsv";

my %levels=("0"=>"domain","1"=>"phylum","2"=>"class","3"=>"order","4"=>"family","5"=>"genus","6"=>"species");
my %prefixes=("0"=>"d__","1"=>"p__","2"=>"c__","3"=>"o__","4"=>"f__","5"=>"g__","6"=>"s__");


my @n_to_sample=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15);

my %info_genome;
print "Reading genome file $genome_file\n";
open my $tsv,"<",$genome_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "Contig"){next;}
    $info_genome{$tab[1]}{"type"}=$tab[3];
    $info_genome{$tab[1]}{"taxo"}=$tab[4];
}
close $tsv;

my %store;
print "Reading repeat file $repeat_file\n";
open my $tsv,"<",$repeat_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "repeat_cluster"){next;}
    $store{$tab[0]}{"taxo"}=$tab[3];
}
close $tsv;
print "Reading array file $array_file\n";
open my $tsv,"<",$array_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "Cluster"){next;}
    if ($tab[8] eq "NA"){next;} ## Unclassified contig
    my $array=$tab[0];
    if (!defined($store{$tab[0]})){next;} ## These are repeats we removed later in the QC process, so we ignore
    my @list_genomes=split(";",$tab[2]);
    foreach my $genome (@list_genomes){
        if ($info_genome{$genome}{"taxo"} ne "NA"){
            if ($info_genome{$genome}{"type"} eq "Bin" || $info_genome{$genome}{"type"} eq "Isolate"){
                $store{$array}{"genome_list"}{$genome}=1;
            }
        }
    }
}
close $tsv;

my $max_n=500;

my %stats;
my %store_tmp;
print "Getting ready to subsamples\n";
for (my $iter=0;$iter<100;$iter++){
    print "## ITERATION $iter ##\n";
    %store_tmp=();
    foreach my $array (sort keys %store){
        my $true_lca=$store{$array}{"taxo"};
        my @t_true_lca=split(";",$true_lca);
        my @t_genomes=keys %{$store{$array}{"genome_list"}};
        my $n_genome=scalar(@t_genomes);
        if ($n_genome > 1){
            # print "############### array $array -- $n_genome ###########\n";
            foreach my $to_sample (@n_to_sample){
                ## At least 5 more than the number to sample
                if ($n_genome>($to_sample+5)){
                    &fisher_yates_shuffle(\@t_genomes);
                    my @tmp_list=@t_genomes[0..($to_sample-1)];
                    # print "List of genomes sampled:".join(" ",@tmp_list)."\n";
                    my $sample_lca=&get_lca(\@tmp_list);
                    my @t_sample_lca=split(";",$sample_lca);
                    # print "LCA of this array is: ".$true_lca."\n";
                    # print "By sampling $to_sample genomes, we get: ".$sample_lca."\n";
                    my $tag=0;
                    for (my $i=0;$i<=6;$i++){
                        my $ori=join(";",@t_true_lca[0..$i]);
                        my $sample=join(";",@t_sample_lca[0..$i]);
                        if ($ori eq $sample){
                            # print "Match at $i - ".$levels{$i}."\n";
                            $store_tmp{$to_sample}{$array}{$i}="yes";
                        }
                        else{
                            # print "Disagreement at $i - ".$levels{$i}."\n";
                            $store_tmp{$to_sample}{$array}{$i}="no";
                        }
                    }
                    # <STDIN>;
                }
            }
        }
    }

    my %tmp;
    foreach my $to_sample (sort {$a <=> $b} keys %store_tmp){
        my @t=keys %{$store_tmp{$to_sample}};
        &fisher_yates_shuffle(\@t);
        %tmp=();
        for (my $i=0;$i<=$max_n;$i++){
            if (!defined($t[$i])){
                print "DID NOT HAVE $i examples of $to_sample ??\n";
                <STDIN>;
            }
            for (my $rank=0;$rank<=6;$rank++){
                $tmp{$rank}{"total"}++;
                # print $t[$i]."\t".$to_sample."\t".$rank."\t".$store_tmp{$to_sample}{$t[$i]}{$rank}."\n";
                if ($store_tmp{$to_sample}{$t[$i]}{$rank} eq "yes"){
                    $tmp{$rank}{"agreed"}++;
                }
            }
        }
        for (my $rank=0;$rank<=6;$rank++){
            my $ratio=$tmp{$rank}{"agreed"}/$tmp{$rank}{"total"};
            my $error=1-$ratio;
            $stats{$to_sample}{$rank}{$iter}=$error;
            # print $to_sample."\t".$levels{$rank}."\t".$ratio."\t".$error."\n";
        }
    }
}

my @tmp=();
open my $s2,">",$out_file_aggregated;
print $s2 "n_sample\trank\tavg_error\tmed_error\tsd_error\n";
foreach my $n_sample (sort {$a <=> $b} keys %stats){
    for (my $rank=0;$rank<=6;$rank++){
        @tmp=();
        foreach my $iter (keys %{$stats{$n_sample}{$rank}}){
            push(@tmp,$stats{$n_sample}{$rank}{$iter});
        }
        my $avg=&average(\@tmp);
        my $med=&median(\@tmp);
        my $stdev=&stdev(\@tmp);
        print $s2 $n_sample."\t".$levels{$rank}."\t".$avg."\t".$med."\t".$stdev."\n";
    }
}
close $s2;

    
sub get_lca(){
    my @list=@{$_[0]};
    # print "List of genomes sampled:".join(" ",@list)."\n";
    my %tmp;
    foreach my $genome (@list){
        # print "-- genome $genome\n";
        my @t=split(";",$info_genome{$genome}{"taxo"});
        for (my $i=0;$i<=6;$i++){
            $tmp{$i}{$t[$i]}++;
        }
    }
    my $tag=0;
    my @lca=();
    for (my $i=0;$i<=6;$i++){
        my @t=keys %{$tmp{$i}};
        # print $i."\t".scalar(@t)."\t".$t[0]."\t".$tmp{$i}{$t[0]}."\n";
        # <STDIN>;
        if ($tag==1){
            $lca[$i]=$prefixes{$i}."unclassified";
        }
        elsif (scalar(@t)==1){
            $lca[$i]=$t[0];
        }
        else{
            $tag=1;
            $lca[$i]=$prefixes{$i}."unclassified";
        }
    }
    return join(";",@lca);
}

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