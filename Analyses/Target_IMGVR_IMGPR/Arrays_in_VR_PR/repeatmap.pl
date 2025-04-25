#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use File::Basename;
use Cwd;
my $h=0;
my $in_file="";
my $dbdir="";
my $outdir="";
my $tmp_dir="tmp_searches/";
my $n_cpu=12;
my $quiet="";
GetOptions ('help' => \$h, 'h' => \$h, 'i=s'=>\$in_file, 'd=s'=>\$dbdir, 'o=s'=>\$outdir, 't=s'=>\$tmp_dir, 'c=s'=>\$n_cpu, 'q'=>\$quiet);
if ($h==1 || $in_file eq "" || $dbdir eq "" || $outdir eq ""){ # If asked for help or did not set up any argument
	print "# Script to compare a set of sequences (fasta file input) to IMG/VR or PR
# Arguments :
# -i: input file (repeats - fasta)
# -d: database folder (folder path)
# -o: output folder
## Optional
# -t: tmp dir (default: tmp_searches/)
# -c: number of cpu/threads (default: 12)
# -q: quiet mode (will be less verbose)
";
	die "\n";
}

if (!($dbdir=~/\/$/)){$dbdir.="/";}

if ($quiet==1){
    print "This will be a quiet run of repeatmap\n";
}


## 
print "## Loading information about the database and checking that all required files are present\n";
## Get number of batches
my @list_batches=<$dbdir/*.1.ebwt*>;
my $max_batch=-1;
foreach my $file (@list_batches){
    if ($file=~/rev/){next;}
    else{
        if ($quiet!=1){print $file."\n";}
        $max_batch++;
    }
}
## Guess name
my $db_name="Dummy_db";
my $test_name=basename($dbdir);
if ($test_name ne ""){$db_name=$test_name;}
else{print "\tI could not guess the name of the db, so it will be called $db_name\n";}
print "\tMax batch number is predicted to be --> $max_batch\n";
## Check that we have the bgz fasta indexed, and store path for later
my $input_fasta="";
my @list_fna=<$dbdir/*.fna.bgz>;
if ($#list_fna==0){
    $input_fasta=$list_fna[0];
    my $test_1=$input_fasta.".fai";
    my $test_2=$input_fasta.".gzi";
    if (-e $test_1 && -e $test_2){}
    else{
        print "\t\tMissing $test_1 or $test_2\n";
        die("\n");
    }
}
else{
    print "\tFound following files when looking for fna:\n";
    foreach my $file (@list_fna){
        print "\t\t$file\n";
    }
    die("\tPblm, I have multiple files ending with fna.bgz, I can't know which one is the relevant one here or 0 files which is also not good\n");
}

## 
print "## Getting information about input file and preparing output paths\n";
## Extract id from file name
my $id="Dummy";
my $test_name=basename($in_file);
if ($test_name=~/(.*)\.f[^\/]+$/){$id=$1;}
else{print "I was not able to extract an id from $in_file so I'll use the id $id\n";}
## Prep the different output files
my $out_file=$outdir."/".$id."_vs_".$db_name."_all.sam"; ## Combined sam (combined because sometimes you have multiple batches)
my $out_file_final=$outdir."/".$id."_vs_".$db_name."_all.regions.tsv";

## Get a single sam file against the whole database
if (-e $out_file){
    print "## $out_file is already here, we don't re-generate\n";
}
else{
    print "## Running the mapping step\n";
    ## Run mapping
    my @list_files;
    for (my $i=0;$i<=$max_batch;$i++){
        print "\t Mapping spacers against batch $i\n";
        my $db_path=$dbdir."batch_".sprintf("%05d",$i);
        my $batch_file=$tmp_dir."/".$id."_vs_batch_".sprintf("%05d",$i).".sam";
        if (-e $batch_file){
            print "$batch_file is already there\n";
        }
        else{
            &run_cmd("bowtie -x $db_path -f $in_file -a -v 3 --threads $n_cpu -t --sam | samtools view -F 4 > $batch_file");
        }
        push(@list_files,$batch_file);
    }
    ## Concatenate all sam files and remove tmp files
    print "\n\nConcatenating all sam files\n";
    &run_cmd("cat ".join(" ",@list_files)." > $out_file");
}

## Load info about each it
my $n=0;
print "Parsing sam file - $out_file\n";
my %store;
open my $sam,"<",$out_file;
while(<$sam>){
    chomp($_);
    $n++;
    if ($n % 1000000 == 0){
        print "... $n ...\n";
    }
    my @tab=split(" ",$_);
    my $repeat=$tab[0];
    my $target=$tab[2];
    my $start=$tab[3];
    $store{$target}{$repeat}{$start}=1;
}
close $sam;

open my $s1,">",$out_file_final;
foreach my $target (keys %store){
    print "Processing hits for target $target\n";
    &process($target,$store{$target},$s1);
    # <STDIN>;
}
close $s1;

sub process(){
    my $target=$_[0];
    my $hash=$_[1];
    my $sout=$_[2];
    my $max=0;
    my $start_c="";
    my $end_c="";
    my $tag=0;
    foreach my $repeat (sort keys %{$hash}){
        print "Looking at repeat $repeat\n";
        my @ordered_list=sort {$a <=> $b} keys %{$$hash{$repeat}};
        $start_c=$ordered_list[0];
        $end_c=$ordered_list[0];
        $max=$start_c+200;
        print "\tWe start from $start_c // $end_c\n";
        for (my $i=1;$i<=$#ordered_list;$i++){
            my $start=$ordered_list[$i];
            while($start<$max && $i<=$#ordered_list){
                $max=$start+200;
                $end_c=$start;
                $start=$ordered_list[$i];
                print "\t\tWe try to extend to $start\n";
                $i++;
            }
            if ($start_c != $end_c){ ## We found a window
                print "\tWe found a window from $start_c to $end_c, we make it + 50 because this is a bit fuzzy just in case\n";
                print $target."\t".($start_c-50)."\t".($end_c+50)."\n";
                print $sout $target."\t".$repeat."\t".($start_c-50)."\t".($end_c+50)."\n";
            }
            $start_c=$ordered_list[$i];
            $end_c=$ordered_list[$i]; ## Starting back from there
            $max=$start_c+200;
            print "\tWe start back from $start_c // $end_c\n";
        }
    }
}


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
