#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to get a nice estimation of the taxonomic assignment error rate based on subsampling
# Arguments :
# run
";
    die "\n";
}

my $repeat_file="Data/Spacer_db/Array_info_filtered_for_db-Nov1-24.tsv";
my $array_file="Data/Additional_data/repeats_to_genome_contigs.tsv";
my $genome_file="Data/Additional_data/all_genome_info.nometag.tsv";
my $out_file_lca="lca_per_type_and_len.tsv";
my $log_file="lca_per_type_and_len.log";

my %levels=("-1"=>"none","0"=>"domain","1"=>"phylum","2"=>"class","3"=>"order","4"=>"family","5"=>"genus","6"=>"species");
my %prefixes=("0"=>"d__","1"=>"p__","2"=>"c__","3"=>"o__","4"=>"f__","5"=>"g__","6"=>"s__");

my %info_contig;
print "Reading genome file $genome_file\n";
open my $tsv,"<",$genome_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "Contig"){next;}
    $info_contig{$tab[0]}{"genome"}=$tab[1];
    $info_contig{$tab[0]}{"length"}=$tab[2];
    $info_contig{$tab[0]}{"type"}=$tab[3];
    $info_contig{$tab[0]}{"taxo"}=$tab[4];
}
close $tsv;

my %check;
print "Reading repeat file $repeat_file\n";
open my $tsv,"<",$repeat_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "repeat_cluster"){next;}
    $check{$tab[0]}=1;
}
close $tsv;


my %stats;
print "Reading array file $array_file\n";
open my $slog,">",$log_file;
open my $tsv,"<",$array_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "Cluster"){next;}
    if ($tab[1] eq "NA"){next;} ## No information, all repeats were found in metagenome contigs
    my $array=$tab[0];
    if (!defined($check{$tab[0]})){next;} ## These are repeats we removed later in the QC process, so we ignore
    my @list_contigs;
    foreach my $contig (split(";",$tab[1])){
        if ($info_contig{$contig}{"taxo"} ne "NA"){
            push(@list_contigs,$contig);
        }
    }
    my $n_contigs=scalar(@list_contigs);
    if ($n_contigs>1){
        print $array." ... \n";
        if ($n_contigs>50){ ## cap it at 50 contigs
            &fisher_yates_shuffle(\@list_contigs);
            @list_contigs=@list_contigs[0..49];
            $n_contigs=scalar(@list_contigs);
        }
        my $last=$n_contigs-1;
        for (my $i=0;$i<$last;$i++){
            for (my $k=$i+1;$k<=$last;$k++){
                my ($type,$len,$lca_rank)=&get_info($list_contigs[$i],$list_contigs[$k]);
                my $print=$list_contigs[$i]."\t".$info_contig{$list_contigs[$i]}{"taxo"}."\t".$info_contig{$list_contigs[$i]}{"genome"}."\t".$info_contig{$list_contigs[$i]}{"type"}."\t".$info_contig{$list_contigs[$i]}{"length"}."\n";
                $print.=$list_contigs[$k]."\t".$info_contig{$list_contigs[$k]}{"taxo"}."\t".$info_contig{$list_contigs[$k]}{"genome"}."\t".$info_contig{$list_contigs[$k]}{"type"}."\t".$info_contig{$list_contigs[$k]}{"length"}."\n";
                $print.=$type."\t".$len."\t".$lca_rank."\n";
                print $slog $print;
                # <STDIN>;
                if ($lca_rank<-1 || $lca_rank>6){
                    print "Pblm - $lca_rank\n";
                    print $print;
                    <STDIN>;
                }
                $stats{$type}{$len}{$lca_rank}++;
            }
        }
    }
}
close $tsv;
close $slog;

open my $s1,">",$out_file_lca;
print $s1 "type,length,lca_rank,count\n";
foreach my $type (sort keys %stats){
    foreach my $len (sort keys %{$stats{$type}}){
        foreach my $rank (sort keys %{$stats{$type}{$len}}){
            print $s1 $type.",".$len.",".$levels{$rank}.",".$stats{$type}{$len}{$rank}."\n";
        }
    }
}
close $s1;

sub get_info(){
    my $go=$_[0];
    my $gt=$_[1];
    my $lo=$info_contig{$go}{"length"};
    my $lt=$info_contig{$gt}{"length"};
    my $to=$info_contig{$go}{"type"};
    my $tt=$info_contig{$gt}{"type"};
    my @tax_o=split(";",$info_contig{$go}{"taxo"});
    my @tax_t=split(";",$info_contig{$gt}{"taxo"});
    my $type="NA";
    if ($to eq "Bin" || $tt eq "Bin"){$type="Bin";}
    elsif($to eq "Isolate" && $tt eq "Isolate"){$type="Isolate";}
    my $length=0;
    if ($lo<$lt){$length=$lo;}
    else{$length=$lt;}
    $length=&get_cat($length);
    my $max_rank="NA";
    for (my $i=0;$i<=6;$i++){
        my $test_o=join(";",@tax_o[0..$i]);
        my $test_t=join(";",@tax_t[0..$i]);
        # print $test_o."\t..vs..\n".$test_t."\n";
        if ($test_o ne $test_t){
            $max_rank=$i-1; ## Previous rank is the max one
            $i=10;
        }
    }
    if ($max_rank eq "NA"){$max_rank=6;} ## We went through the whole thing without any issue
    return ($type,$length,$max_rank);
}

sub get_cat(){
    my $length=$_[0];
    my $cat="NA";
    if ($length>=0 && $length<10000){$cat="0-10kb";}
    elsif ($length>=10000 && $length<100000){$cat="10-100kb";}
    elsif ($length>=100000 && $length<500000){$cat="100-500kb";}
    elsif ($length>=500000){$cat="500kb+";}
    return $cat;
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
