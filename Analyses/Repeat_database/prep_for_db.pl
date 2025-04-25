#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
	print "# Script to prepare the table of runs for the database
# Arguments :
# run
";
	die "\n";
}


my $in_file="Array_info_with_complexity.tsv";
my $out_file="Array_info_with_complexity_for_db.tsv";

open my $s1,">",$out_file;
open my $tsv,"<",$in_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "Cluster"){
        print $s1 "repeat_cluster\ttype\tlca_status\tlca_full\tlca_class\tlca_family\tlca_genus\tflag\n";
    }
    else{
        print $s1 $tab[0]."\t".$tab[3]."\t".$tab[6]."\t".$tab[5]."\t".&get_class($tab[5])."\t".&get_family($tab[5])."\t".&get_genus($tab[5])."\t".$tab[12]."\n";
    }
}
close $tsv;
close $s1;


sub get_class(){
    my $taxo=$_[0];
    if ($taxo eq "NA"){return "NA";}
    else{
        my @t=split(";",$taxo);
        if ($t[2]=~/^c__/){
            return join(";",@t[0..2]);
        }
        else{
            print "Pblm with $taxo\n";
            die("\n");
        }
    }
}


sub get_family(){
    my $taxo=$_[0];
    if ($taxo eq "NA"){return "NA";}
    else{
        my @t=split(";",$taxo);
        if ($t[4]=~/^f__/){
            return join(";",@t[0..4]);
        }
        else{
            print "Pblm with $taxo\n";
            die("\n");
        }
    }
}


sub get_genus(){
    my $taxo=$_[0];
    if ($taxo eq "NA"){return "NA";}
    else{
        my @t=split(";",$taxo);
        if ($t[5]=~/^g__/){
            return join(";",@t[0..5]);
        }
        else{
            print "Pblm with $taxo\n";
            die("\n");
        }
    }
}