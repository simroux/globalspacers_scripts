#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h); # , 'i=s'=>\$in_file, 'r=s'=>\$out_dir_root);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to get a subsample of uvig_hq_coverage_by_spacer-with_sp_info.tsv to plot
# Arguments :
# run
";
    die "\n";
}

my $in_file="../../../Analyses/Target_IMGVR_IMGPR/Target_coverage/uvig_hq_coverage_by_spacer-with_sp_info.tsv";
my $out_file="uvig_hq_cover-vs-n_spacer_subsampled.tsv";

my @t=();
my $f_l="";
print "reading $in_file\n";
open my $tsv,"<",$in_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "uvig"){
        $f_l=$_;
        next;
    }
    my $code=$tab[0].";".$tab[1];
    push(@t,$code);
}
close $tsv;


&fisher_yates_shuffle(\@t);
print "subsampling\n";
my %check;
for (my $i=0;$i<100000;$i++){
    my @tab=split(";",$t[$i]);
    $check{$tab[0]}{$tab[1]}=1;
}

print "Writing the new subsampled table\n";
open my $s1,">",$out_file;
print $s1 $f_l."\n";
open my $tsv,"<",$in_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "uvig"){
        next;
    }
    if ($check{$tab[0]}{$tab[1]}==1){
        print $s1 $_."\n";
    }
}
close $tsv;
close $s1;


sub fisher_yates_shuffle {
	my $deck = shift;  # $deck is a reference to an array
	return unless @$deck; # must not be empty!
	my $i = @$deck;
	while (--$i) {
		my $j = int rand ($i+1);
		@$deck[$i,$j] = @$deck[$j,$i];
	}
}
