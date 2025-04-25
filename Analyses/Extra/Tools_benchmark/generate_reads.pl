#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
use Parallel::ForkManager;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h); 
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to get the selected genomes and generate the reads, using a local directory with the corresponding IMG genomes
# Arguments :
# run
";
    die "\n";
}

my $in_tsv="List_genomes_sub200.tsv";
my $out_dir="Reads_sub200";


my $in_dir="IMG.taxon.fna/";

my %check;
my $i=0;
open my $tsv,"<",$in_tsv;
while(<$tsv>){
     chomp($_);
     my @tab=split("\t",$_);
     $check{$tab[0]}=1;
}
close $tsv;

my $n_processes = 24;
my $pm = Parallel::ForkManager->new( $n_processes );
foreach my $genome (keys %check){
     print "Taking reads from $genome\n";
     my $in_file=$in_dir.$genome.".fna";
     $pm->start and next;
     &run_cmd("art_illumina -ss HS25 -p -l 150 -f 0.8 -m 200 -s 10 -i $in_file -o $out_dir/Cover_0.8_HS25_$genome");
     &run_cmd("art_illumina -ss HS25 -p -l 150 -f 10 -m 200 -s 10 -i $in_file -o $out_dir/Cover_10_HS25_$genome");
     &run_cmd("art_illumina -ss HS25 -p -l 150 -f 100 -m 200 -s 10 -i $in_file -o $out_dir/Cover_100_HS25_$genome");
     &run_cmd("art_illumina -ss HS25 -qs 3 -qs2 3 -p -l 150 -f 0.8 -m 200 -s 10 -i $in_file -o $out_dir/Cover_0.8_HS25qs3_$genome");
     &run_cmd("art_illumina -ss HS25 -qs 3 -qs2 3 -p -l 150 -f 10 -m 200 -s 10 -i $in_file -o $out_dir/Cover_10_HS25qs3_$genome");
     &run_cmd("art_illumina -ss HS25 -qs 3 -qs2 3 -p -l 150 -f 100 -m 200 -s 10 -i $in_file -o $out_dir/Cover_100_HS25qs3_$genome");
     $pm->finish;
}
$pm->wait_all_children;

sub run_cmd{
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
