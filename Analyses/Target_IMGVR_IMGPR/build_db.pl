#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
my $in_file="";
my $out_dir_root="";
my $n_cpu=8;
GetOptions ('help' => \$h, 'h' => \$h, 'i=s'=>\$in_file, 'd=s'=>\$out_dir_root);
if ($h==1 || $in_file eq "" || $out_dir_root eq ""){ # If asked for help or did not set up any argument
	print "# Script to build a database for spacermap
# Arguments :
# -i: input fasta file of virus sequences (fasta or fasta.gz)
# -d: database directory (must not exists)
";
	die "\n";
}

if ($out_dir_root=~/\/$/){}
else{$out_dir_root.="/";}

if (-d $out_dir_root){
    die("The folder $out_dir_root exists already");
}
&run_cmd("mkdir $out_dir_root");

my $out_file_fasta=$out_dir_root."Input_sequences.fna.bgz";

print "Importing the fasta file as a bgz compressed\n";
my $n=0;
open my $s1,"| bgzip > $out_file_fasta";
open my $fna,"<",$in_file;
while(<$fna>){
    chomp($_);
    if ($_=~/^>/){
        $n++;
        print $s1 $_."\n";
    }
    else{
        print $s1 $_."\n";
    }
}
close $fna;
close $s1;

print "Running dustmasker (will be used later, in post-processing)\n";
my $dust_out=$out_dir_root."Input_sequences.fna.dustmasker";
&run_cmd("dustmasker -in $in_file -out $dust_out");

print "Building the index\n";
&run_cmd("samtools faidx $out_file_fasta");

print "Building the bowtie-1 database\n";
if ($n>=200000){
    print "We have more than 200,000 sequences, we will break this down by batches of 200k sequences\n";
	my $tmp_folder=$out_dir_root."tmp_batches/";
	&run_cmd("mkdir $tmp_folder");
	my $i=0;
	my $c_count=0;
	my $s1;
	my $out_file=$tmp_folder."batch_".sprintf("%05d",$i).".fna";
	open $s1,">",$out_file;
	open my $fna,"<",$in_file;
	while(<$fna>){
		chomp($_);
		if ($_=~/^>/){
			if ($c_count==200000){
				close $s1;
				$i++;
				print "Switching to batch ".sprintf("%05d",$i)."\n";
				my $out_file=$tmp_folder."batch_".sprintf("%05d",$i).".fna";
				open $s1,">",$out_file;
				$c_count=0;
			}
			$c_count++;
		}
		print $s1 $_."\n";
	}
	close $fna;
	close $s1;
    for (my $k=0;$k<=$i;$k++){
		my $in_file=$tmp_folder."batch_".sprintf("%05d",$k).".fna";
		&run_cmd("bowtie-build -o 2 --threads $n_cpu $in_file $out_dir_root/batch_".sprintf("%05d",$k)." ");	
	}
	print "Should be able to remove $tmp_folder now if everything went right\n";
}
else{
    &run_cmd("bowtie-build -o 2 --threads $n_cpu $in_file $out_dir_root/batch_00000");
}

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
