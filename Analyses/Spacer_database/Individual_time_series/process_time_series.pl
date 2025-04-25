#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h); # , 'i=s'=>\$in_file, 'r=s'=>\$out_dir_root);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to process all time series
# Arguments :
# none
";
    die "\n";
}

my $in_file="Analyses/Extra/SRA_metadata/Overview_time_series.txt";

my %check;
open my $tsv,"<",$in_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "series_name"){next;}
    $check{$tab[0]}=1;
}
close $tsv;

foreach my $series_name (sort keys %check){
    my $spacer_file="Individual_time_series/".$series_name."/".$series_name."_filt_spacers.tsv";
    my $hit_file="Individual_time_series/".$series_name."/".$series_name."_filt_hits.tsv";
    my $cover_file="Individual_time_series/".$series_name."/".$series_name."_filt_hits_coverage.tsv";
    if (!(-e $spacer_file) || !(-e $hit_file)){
        &run_cmd("python extract_time_series_data.py -id $series_name");
    }
    else{print "Spacer and hit file already there for $series_name\n";}
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
