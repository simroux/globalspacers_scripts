#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h); # , 'i=s'=>\$in_file, 'r=s'=>\$out_dir_root);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to re-run repeattyper on repeat sequences to get updated Type predictions
# Arguments :
# none
";
    die "\n";
}


my $in_file="all_repeats_clstr.tab";
my $in_file_tab="Array_info_filtered_for_db.tsv";
my $out_file_tab="Array_info_filtered_for_db-updated_type.tsv";

my $out_file="all_repeats_tmp.txt";
my $out_file_tmp="all_repeats_tmp_pred.tsv";

my %seq_to_id;
print "Reading $in_file\n";
open my $s1,">",$out_file;
open my $tsv,"<",$in_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "Array"){next;}
    if ($tab[2]=~/N/ || $tab[2]=~/Y/ || $tab[2]=~/R/  || $tab[2]=~/S/){next;} ## Ignore repeats with Ns or other ambiguities
    if ($tab[2] eq "TTTCCGTCTCCTTGTGGAGTATTATTCATTCTTAt"){ ## Fix a lower-case sequence
        $tab[2]="TTTCCGTCTCCTTGTGGAGTATTATTCATTCTTAT";
    }
    print $s1 $tab[2]."\n";
    if (defined($seq_to_id{$tab[2]})){
        die("Pblm, already an id for $tab[2]\n");
    }
    $seq_to_id{$tab[2]}=$tab[1];
}
close $tsv;
close $s1;

&run_cmd("repeatType $out_file > $out_file_tmp");


my %store;
print "Reading $out_file_tmp\n";
open my $tsv,"<",$out_file_tmp;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if (defined($seq_to_id{$tab[0]})){
        if ($tab[2]>=0.75){
            $store{$seq_to_id{$tab[0]}}=$tab[1];
        }
    }
    else{
        die("Pblm, no id for $tab[0] ?\n");
    }
}
close $tsv;

print "Preparing output\n";
open my $tsv,"<",$in_file_tab;
open my $s1,">",$out_file_tab;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "repeat_cluster"){
        print $s1 $_."\n";
        next;
    }
    if (defined($store{$tab[0]})){
        if ($tab[1] ne $store{$tab[0]}){
            print $tab[1]." becomes ".$store{$tab[0]}."\n";
        }
        $tab[1]=$store{$tab[0]};
    }
    else{
        if ($tab[1] ne "Unknown"){
            print $tab[1]." becomes Unknown\n";
        }
        $tab[1]="Unknown";
    }
    print $s1 $tab[0]."\t".$tab[1]."\t".$tab[2]."\t".$tab[3]."\t".$tab[4]."\t".$tab[5]."\t".$tab[6]."\n";
}
close $tsv;
close $s1;

sub run_cmd(){
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

