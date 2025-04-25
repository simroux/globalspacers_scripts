#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use File::Basename;
use Cwd;
my $h=0;
my $dbdir="";
my $result_file="";
my $n_cpu=8;
my $quiet=0;
GetOptions ('help' => \$h, 'h' => \$h, 'd=s'=>\$dbdir, 'o=s' => \$result_file, 'q'=>\$quiet);
if ($h==1 || $dbdir eq "" || $result_file eq ""){ # If asked for help or did not set up any argument
	print "# Script to post-process the results of spacermap
# Arguments :
# -d: database directory (used in spacermap computation)
# -o: output file from spacermap (xxx.matches.tsv)
# -q: quiet mode
";
	die "\n";
}

## Guess the ids and check all the required files are here
print "## Loading information about the database and checking that all required files are present\n";
## Guess name of dustmasker file and mapping file to arrays
my $input_dust="";
my $input_arrays="";
my @list_fna=<$dbdir/*.fna.bgz>;
if ($#list_fna==0){
    my $input_fasta=$list_fna[0];
	my $test_name=basename($input_fasta);
	if ($test_name=~/(.*\.f[^\/]+)\.bgz$/){
		my $id=$1;
		my $test_1=$dbdir."/".$id.".dustmasker";
		if (-e $test_1){
			$input_dust=$test_1;
		}
		else{
			print "\t\tMissing $test_1\n";
			die("\n");
		}
		my $test_2=$dbdir."/".$id.".predicted_arrays.tsv";
		if (-e $test_2){
			$input_arrays=$test_2;
		}
		else{
			print "\t\tMissing $test_2\n";
			die("\n");
		}
	}
	else{print "I was not able to extract an id from $input_fasta\n"; die("\n");}
    
}
else{
    print "\tFound following files when looking for fna:\n";
    foreach my $file (@list_fna){
        print "\t\t$file\n";
    }
    die("\tPblm, I have multiple files ending with fna.bgz, I can't know which one is the relevant one here or 0 files which is also not good\n");
}
## Prep the output file path
## Final output file - add "_clean" before tsv, warning if same file (do not want to overwrite)
my $out_file=$result_file;
$out_file=~s/\.tsv\.gz/_clean.tsv\.gz/;
if ($out_file eq $result_file){
	print "Pblm - I was not able to guess the name I should use for the output file: \n";
	print $result_file."\n";
	print $out_file."\n";
	die("\n");
}


## Load dustmasker results per contig
print "## Loading dustmasker info from $input_dust ..\n";
my $count=0;
my %info_dust;
my $c_id="";
my $n_block=0;
open my $fts,"<",$input_dust;
while(<$fts>){
	chomp($_);
	if ($_=~/^>(\S+)/){
		$c_id=$1;
		$n_block=0;
		$count++;
		if ($count % 1000000 == 0){print "$count ... \n";}
	}
	else{
		$_=~s/\s//g;
		my @t=split("-",$_);
		$info_dust{$c_id}{$n_block}{"start"}=$t[0];
		$info_dust{$c_id}{$n_block}{"end"}=$t[1];
		# print $c_id."\t".$n_block."\t".$info_dust{$c_id}{$n_block}{"start"}.":".$info_dust{$c_id}{$n_block}{"end"}."\n";
		$n_block++;
	}
}
close $fts;

## Load dustmasker results per contig
print "## Loading info about MGE-encoded arrays from $input_arrays ..\n";
$count=0;
my %info_arrays;
$c_id="";
$n_block=0;
open my $fts,"<",$input_arrays;
while(<$fts>){
	chomp($_);
	my @tab=split("\t",$_);
	$count++;
	if ($count % 1000000 == 0){print "$count ... \n";}
	$info_arrays{$tab[0]}{$n_block}{"start"}=$tab[1];
	$info_arrays{$tab[0]}{$n_block}{"end"}=$tab[2];
	# print $tab[0]."\t".$n_block."\t".$info_arrays{$tab[0]}{$n_block}{"start"}.":".$info_arrays{$tab[0]}{$n_block}{"end"}."\n";
	$n_block++;
}
close $fts;


## Then go through the results for each contig and flag cases that are potentially problematic
print "Reading $result_file and writing results to $out_file\n";
$count=0;
my %tmp;
my $i=0;
$c_id="";
## Changing to compressing on the fly
open my $s1,"| gzip -c > $out_file";
my $tsv;
if ($result_file=~/\.gz$/){open $tsv,"gunzip -c $result_file |";}
else{open $tsv,"<",$result_file;}
while(<$tsv>){
	chomp($_);
	$count++;
	if ($count % 100000 == 0){print "$count ... \n";}
	# print $_."\n";
	my @tab=split("\t",$_);
	if ($tab[0] eq "Spacer id"){
		print $s1 $_."\tFlags\n";
		next;
	}
	if ($tab[1] ne $c_id){
		if ($c_id ne ""){
			&process($c_id,\%tmp,$s1);
		}
		%tmp=();
		$c_id=$tab[1];
		$i=0;
	}
	$tmp{$i}{"line"}=$_;
	$i++;
}
close $tsv;
## PROCESS LAST ONE
if ($c_id ne ""){
	&process($c_id,\%tmp,$s1);
}
close $s1;

sub process(){
	my $id=$_[0];
	my $hash=$_[1];
	my $sout=$_[2];
	# my %startlist;
	# print "Processing $id\n";
	## Check which hits are overlapping with dust -> Flag
	## Also flag high mismatches
	foreach my $hit (sort keys %{$hash}){
		my @t=split("\t",$$hash{$hit}{"line"});
		my $start=$t[2];
		my $end=$t[3];
		my $mis=$t[5];
		if ($end<$start){
			print "Pblm with start / end == $start / $end\n";
			die("\n");
		}
		# $startlist{$start}{"count"}++;
		if (defined($info_dust{$id})){
			foreach my $block (keys %{$info_dust{$id}}){
				# print $id."\t".$block."\n";
				if ($info_dust{$id}{$block}{"end"}<$info_dust{$id}{$block}{"start"}){
					print "Pblm with start / end on dust == ".$info_dust{$id}{$block}{"start"}." / ".$info_dust{$id}{$block}{"end"}."\n";
					die("\n");
				}
				if ($end<$info_dust{$id}{$block}{"start"}){}
				elsif ($start>$info_dust{$id}{$block}{"end"}){}
				else{
					if ($quiet!=1){
						print "Found an overlap with dust region: \n";
						print $t[0]."\t".$t[1]."\t".$start."\t".$end."\n";
						print $id."\t".$info_dust{$id}{$block}{"start"}."\t".$info_dust{$id}{$block}{"end"}."\n";
					}
					# <STDIN>;
					$$hash{$hit}{"flags"}{"low_complexity"}=1;
				}
			}
		}
		if (defined($info_arrays{$id})){
			foreach my $block (keys %{$info_arrays{$id}}){
				# print $id."\t".$block."\n";
				if ($info_arrays{$id}{$block}{"end"}<$info_arrays{$id}{$block}{"start"}){
					print "Pblm with start / end on arrays == ".$info_arrays{$id}{$block}{"start"}." / ".$info_arrays{$id}{$block}{"end"}."\n";
					die("\n");
				}
				if ($end<$info_arrays{$id}{$block}{"start"}){}
				elsif ($start>$info_arrays{$id}{$block}{"end"}){}
				else{
					if ($quiet!=1){
						print "Found an overlap with array region: \n";
						print $t[0]."\t".$t[1]."\t".$start."\t".$end."\n";
						print $id."\t".$info_arrays{$id}{$block}{"start"}."\t".$info_arrays{$id}{$block}{"end"}."\n";
					}
					# <STDIN>;
					$$hash{$hit}{"flags"}{"potential_array"}=1;
				}
			}
		}
		if ($mis>1){
			$$hash{$hit}{"flags"}{"mismatches"}=1;
		}
		## This is only for when we don't do additional filtering based on multiple hits, etc, 
		my $combined_flags=join(";",sort keys %{$$hash{$hit}{"flags"}});
		print $sout $$hash{$hit}{"line"}."\t".$combined_flags."\n";
		# print $$hash{$hit}{"line"}."\t".$combined_flags."\n";
	}
}