#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h); 
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to run SE on the simulations
# Arguments :
# run
";
    die "\n";
}

my $in_list="List_genomes_sub200.tsv";
my $read_dir="Selected_Reads_sub200/";
my $out_dir="Results_SE_sub200/";

my $path_db="./SE_Db_r214GBIMG_0.9";

my %check;
open my $tsv,"<",$in_list;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "taxon"){next;}
    $check{$tab[0]}=1;
}
close $tsv;

my @cover=("0.8","10","100");
my @profiles=("HS25","HS25qs3");

foreach my $genome (keys %check){
     print "Processing $genome\n";
     foreach my $cov (@cover){
          foreach my $prof (@profiles){
             my $out_dir_run=$out_dir."/".$genome."_Cover_".$cov."_".$prof."/";
             my $log_file=$out_dir."/".$genome."_Cover_".$cov."_".$prof.".log.txt";
             my $read_file=$read_dir."/Cover_".$cov."_".$prof."_".$genome."_il.fastq.gz";
             if (-d $out_dir_run){
                 print $out_dir_run." already exists\n";
             }
             else{
                 &run_cmd("spacerextractor extract_spacers --repeat_db_dir $path_db --input_fastq $read_file --n_threads 6 -o $out_dir_run > $log_file 2> $log_file");
                 &run_cmd("spacerextractor filter_hq_spacers -d $path_db --wdir $out_dir_run -t 6 -m 0 >> $log_file 2>> $log_file");
             }
        }
   }
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
