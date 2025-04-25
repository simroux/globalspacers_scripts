#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use File::Basename;
use Cwd;
my $h=0;
my $infile="";
my $infile_cluster="";
my $out_file="";
GetOptions ('help' => \$h, 'h' => \$h, 'i=s'=>\$infile, 'c=s'=>\$infile_cluster, 'o=s'=>\$out_file);
if ($h==1 || $infile eq "" || $infile_cluster eq "" || $out_file eq ""){ # If asked for help or did not set up any argument
	print "# Script to prepare input for duckdb (change column name, aggregate tsv files if you split spacers and process it by batch)
# Arguments :
# -i: input file
# -c: table file listing spacers currently in the duckdb database
# -o: output file
";
	die "\n";
}


my %check_cluster;
my $i=0;
print "Reading $infile_cluster\n";
open my $tsv,"<",$infile_cluster;
while(<$tsv>){
	chomp($_);
	my @tab=split("\t",$_);
	if ($tab[0] eq "cluster_id"){next;}
	$check_cluster{$tab[0]}=1;
	$i++;
	if ($i % 1000000 == 0){print "\t $i ... \n";}
}
close $tsv;

open my $s1,">",$out_file;
$i=0;
print $s1 "cluster_id\ttarget_id\thit_start\thit_end\thit_strand\tn_mismatches\tCIGAR\tMD\tspacer\tprotospacer\tupstream\tdownstream\n";
print "Processing $infile ...\n";
open my $tsv,"gunzip -c $infile |";
$i=0;
while(<$tsv>){
	chomp($_);
	my @tab=split("\t",$_);
	if ($tab[0] eq "Spacer id"){next;}
	if ($tab[12] eq ""){
		if ($check_cluster{$tab[0]}==1){
			print $s1 join("\t",@tab[0..11])."\n";
			$i++;
			if ($i % 500000 == 0){print "\t $i ... \n";}
		}
		else{
			print "$tab[0] is not interesting anymore because not in the final db\n";
		}
	}
}
close $tsv;
close $s1;