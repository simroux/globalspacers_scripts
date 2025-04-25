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
my $out_dir="Results_Crass_sub200/";

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
                 print $out_dir." already exists\n";
             }
             else{
                  my $timestamp = getLoggingTime();
                  open my $slog,">",$log_file;
                  print $slog "Start:\t".$timestamp."\n";
                  close $slog;
                  &run_cmd("crass -o $out_dir_run -l 4 $read_file >> $log_file 2>> $log_file");
                  my $timestamp = getLoggingTime();
                  open my $slog,">>",$log_file;
                  print $slog "End:\t".$timestamp."\n";
                  close $slog;
             }
             my $out_file_final=$out_dir."/".$genome."_Cover_".$cov."_".$prof.".crass_results.tsv";
             my $in_file_crass_result=$out_dir_run."/crass.crispr";
             if (-e $in_file_crass_result){
                  my $c_c="";
                  open my $s1,">",$out_file_final;
                  open my $tsv,"<",$in_file_crass_result;
                  while(<$tsv>){
                       chomp($_);
                       if ($_=~/drseq=\"(\S+)\"/){$c_c=$1;}
                       elsif($_=~/spacer cov=\"(\d+)\" seq=\"(\S+)\"/){
                            print $s1 $c_c."\t".$2."\t".$1."\n";
                       }
                  }
                  close $tsv;
                  close $s1;
             }
             # last;
        }
   }
}

sub getLoggingTime {

    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
    my $nice_timestamp = sprintf ( "%04d%02d%02d %02d:%02d:%02d",
                                   $year+1900,$mon+1,$mday,$hour,$min,$sec);
    return $nice_timestamp;
}

sub run_cmd {
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
