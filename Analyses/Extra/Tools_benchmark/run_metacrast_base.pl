#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h); # , 'i=s'=>\$in_file, 'r=s'=>\$out_dir_root);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to run SE on the simulations
# Arguments :
# none
";
    die "\n";
}
my $in_list="List_genomes_sub200.tsv";
my $read_dir="Selected_Reads_sub200/";
my $spacer_file="Spacer_table_sub200.tsv";
my $out_dir="Results_MetaCrast_sub200/";


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
     my $repeat_file=$out_dir."/".$genome."_repeats.fna";
     if (!(-e $repeat_file)){
          my %seen;
          open my $s1,">",$repeat_file;
          open my $tsv,"<",$spacer_file;
          while(<$tsv>){
               chomp($_);
               my @tab=split("\t",$_);
               if ($tab[0] eq $genome){
                    if (!defined($seen{$tab[2]})){
                         print $s1 ">".$tab[1]."\n";
                         print $s1 $tab[2]."\n";
                         $seen{$tab[2]}=1;
                    }
               }
          }
          close $tsv;
          close $s1;
     }
     foreach my $cov (@cover){
          foreach my $prof (@profiles){
             my $out_dir_run=$out_dir."/".$genome."_Cover_".$cov."_".$prof."/";
             my $tmp_dir_run=$out_dir."/".$genome."_Cover_".$cov."_".$prof."-tmp/";
             my $log_file=$out_dir."/".$genome."_Cover_".$cov."_".$prof.".log.txt";
             my $read_file=$read_dir."/Cover_".$cov."_".$prof."_".$genome."_il.fastq.gz";
             if (-d $out_dir_run){
                 print $out_dir." already exists\n";
             }
             else{
                  my $tmp_fastq=$genome."_Cover_".$cov."_".$prof.".tmp.fastq";
                  &run_cmd("gunzip -c $read_file > $tmp_fastq");
                  my $timestamp = getLoggingTime();
                  open my $slog,">",$log_file;
                  print $slog "Start:\t".$timestamp."\n";
                  close $slog;
                  &run_cmd("MetaCRAST -p $repeat_file -i $tmp_fastq -o $out_dir_run -d 0 -n 6 -q -t $tmp_dir_run");
                  my $timestamp = getLoggingTime();
                  open my $slog,">>",$log_file;
                  print $slog "End:\t".$timestamp."\n";
                  close $slog;
                  &run_cmd("rm $tmp_fastq");
             }
             my $out_file_final=$out_dir."/".$genome."_Cover_".$cov."_".$prof.".metacrast_results.tsv";
             open my $s1,">",$out_file_final;
             foreach my $file (<$out_dir_run/Spacers-*.fa>){
                  $file=~/Spacers-\d+-(\S+)\.fa/;
                  my $dr=$1;
                  my $c_seq="";
                  open my $fa,"<",$file;
                  while(<$fa>){
                       chomp($_);
                       if ($_=~/^>/){
                            if ($c_seq ne ""){
                                 print $s1 $dr."\t".$c_seq."\n";
                            }
                            $c_seq="";
                       }
                       else{
                            $c_seq.=$_;
                       }
                  }
                  close $fa;
                  if ($c_seq ne ""){
                       print $s1 $dr."\t".$c_seq."\n";
                  }
             }
             close $s1;
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
