#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
use Parallel::ForkManager;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to get reads as interleaved and compressed
# Arguments :
# run
";
    die "\n";
}


my $in_tsv="List_genomes_sub200.tsv";
my $in_dir="Reads_sub200";

my $out_dir="Selected_Reads_sub200";



my %check;
my $i=0;
open my $tsv,"<",$in_tsv;
while(<$tsv>){
     chomp($_);
     my @tab=split("\t",$_);
     $check{$tab[0]}=1;
}
close $tsv;

my @cover=("0.8","10","100");
my @profiles=("HS25","HS25qs3");

my $n_processes = 4;
my $pm = Parallel::ForkManager->new( $n_processes );
foreach my $genome (keys %check){
     $pm->start and next;
     print "Taking reads from $genome\n";
     foreach my $cov (@cover){
          foreach my $prof (@profiles){
               my $file_r1=$in_dir."/Cover_".$cov."_".$prof."_".$genome."1.fq";
               my $file_r2=$in_dir."/Cover_".$cov."_".$prof."_".$genome."2.fq";
               my $out_file=$out_dir."/Cover_".$cov."_".$prof."_".$genome."_il.fastq.gz";
               &run_cmd("shifter -i bryce911/bbtools reformat.sh in=$file_r1 in2=$file_r2 out=$out_file");
          }
     }
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
