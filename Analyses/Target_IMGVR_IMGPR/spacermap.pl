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
	print "# Script to compare a set of sequences (fasta file input) to IMG/VR
# Arguments :
# -i: input file (spacer - fasta)
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

## Defaults parameters
my $n_upstream=10;
my $n_downstream=10;

if ($quiet==1){
    print "This will be a quiet run of spacermap\n";
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
my $out_file_final=$outdir."/".$id."_vs_".$db_name."_all.matches.tsv";
my $max_per_region=10000000;
my $tmp_region_file=$tmp_dir."/regions.txt";
my @tab_region_files=($tmp_region_file);
## Checking that the tmp dir and out dir exists
if (!(-d $tmp_dir)){
    die("Pblm, $tmp_dir does not seem to be accessible\n");
}
if (!(-d $outdir)){
    die("Pblm, $outdir does not seem to be accessible\n");
}

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
my %store;
my $n=0;
print "Parsing sam file - $out_file\n";
open my $sam,"<",$out_file;
while(<$sam>){
    chomp($_);
    $n++;
    if ($n % 100000 == 0){
        print "... $n ...\n";
    }
    ## Sam format:
    # Query (spacer)
    # Flag - expects 0 or 16 for plus and minus strand
    # Target (target)
    # Hit start (on target)
    # Quality (ignored)
    # CIGAR
    # Rnext / Pnext / Tlen - all 3 ignored (paired-end related)
    # Sequence
    # Sequence quality (we ignore)
    # From there on, various fflags, including number of mismatches
    my @tab=split(" ",$_);
    my $spacer=$tab[0];
    my $target=$tab[2];
    my $flag=$tab[1];
    my $strand="NA";
    if ($flag==0){$strand="plus";}
    elsif($flag==16){$strand="minus";}
    else{
        die("Unexpected flag $flag\n$_\n");
    }
    my $start=$tab[3];
    if (defined($store{$target}{$spacer}{$start}{$strand})){
        print "!!! Already a hit for $target - $spacer - $start - $strand ?!?! \n";
        die("\n");
        next;
    }
    my $cigar=$tab[5];
    my $md_string="NA";
    my $n_mis="NA";
    for (my $i=11;$i<=$#tab;$i++){
        if ($tab[$i]=~/NM:i:(\d+)/){$n_mis=$1;}
        if ($tab[$i]=~/MD:Z:(.+)/){$md_string=$1;}
    }
    $store{$target}{$spacer}{$start}{$strand}{"seq"}=$tab[9];
    $store{$target}{$spacer}{$start}{$strand}{"CIGAR"}=$cigar;
    $store{$target}{$spacer}{$start}{$strand}{"n_mis"}=$n_mis;
    $store{$target}{$spacer}{$start}{$strand}{"md_string"}=$md_string;
    ## Store relevant info including
    # Target id
    # Spacer id
    # Start
    # Strand
    # Sequence
    # N mismatches
    # CIGAR
    # MD string
}
close $sam;

print "Extracting regions upstream and downstream of potential protospacer\n";
my %seen;
my $n=0;
my $n_batch=1;
open my $s1,">",$tmp_region_file;
foreach my $target (sort keys %store){
    foreach my $spacer (sort keys %{$store{$target}}){
        foreach my $start (sort keys %{$store{$target}{$spacer}}){
            foreach my $strand (sort keys %{$store{$target}{$spacer}{$start}}){
                my $real_start=$start-10; 
                if ($real_start<1){$real_start=1;} ## We need -1 because of the way sequences are extracted later
                my $real_end=$start+length($store{$target}{$spacer}{$start}{$strand}{"seq"})+9; ## We only add 9 because length already add +1 (because of the way coordinates work)
                $store{$target}{$spacer}{$start}{$strand}{"real_start"}=$real_start;
                $store{$target}{$spacer}{$start}{$strand}{"real_end"}=$real_end;
                my $test=$target."_".$real_start."_".$real_end;
                if (defined($seen{$test})){
                    # print "$test was already here\n";
                }
                else{
                    $seen{$test}=1;
                    print $s1 $target.":".$real_start."-".$real_end."\n";
                    $n++;
                    if ($n % $max_per_region == 0){
                        close $s1;
                        my $new_tmp_region_file=$tmp_region_file.$n_batch.".txt";
                        push(@tab_region_files,$new_tmp_region_file);
                        print "adding regions to $new_tmp_region_file ...\n";
                        $n_batch++;
                        $n=1;
                        open $s1,">",$new_tmp_region_file;
                    }
                }
            }
        }
    }
}
close $s1;

# chr:from-to
###
## Use samtools faidx as a stream, then process sequences on the fly assuming only 1 per line (which should always be the case since we use -n 0), that way you don't store too many additional things (already stored all the hits, that's enough..)
## If processing on the fly seems to cause memory issues, we could instead generate a temp fasta file and get the same result
$n=0;
open my $s2,">",$out_file_final;
## header
print $s2 "Spacer id\tTarget id\tStart\tEnd\tStrand\tN mismatches\tCIGAR\tMD\tspacer\tprotospacer\tupstream\tdownstream\n";
foreach my $r_file (@tab_region_files){
    print "Processing $r_file ...\n";
    open my $fna,"samtools faidx $input_fasta -r $r_file -n 0 -c |";
    while(<$fna>){
        chomp($_);
        if ($_=~/^>(\S+)/){
            my $seq_id=$1;
            my $seq=<$fna>;
            chomp($seq);
            $seq=~tr/[a-z]/[A-Z]/; ## Bowtie1 uses only upper case, and if we don't do that we may see mismatches where there are none
            if ($quiet!=1){
                print $seq_id."\n";
                print $seq."\n";
            }
            $n++;
            if ($n % 100000 == 0){
                print "... $n ...\n";
            }
            if ($seq_id=~/^(.*):(\d+)-(\d+)$/){
                my $target=$1;
                my $real_start=$2;
                my $real_end=$3;
                &process_proto($real_start,$real_end,$seq,$store{$target},$target,$s2);
            }
            else{
                die("Pblm with format of id $seq_id\n");
            }
            # <STDIN>;
            
        }
    }
    close $fna;
}
close $s2;

sub process_proto(){
    my $r_start=$_[0];
    my $r_end=$_[1];
    my $seq=$_[2];
    my $hash=$_[3];
    my $target=$_[4];
    my $sout=$_[5];
    foreach my $spacer (sort keys %{$hash}){
        foreach my $start (sort keys %{$$hash{$spacer}}){
            foreach my $strand (sort keys %{$$hash{$spacer}{$start}}){
                if ($r_start == $$hash{$spacer}{$start}{$strand}{"real_start"} && $r_end == $$hash{$spacer}{$start}{$strand}{"real_end"}){
                    if ($quiet!=1){print "We found one that matches $r_start and $r_end ==> $start / --$strand-- \n";}
                    if ($$hash{$spacer}{$start}{$strand}{"seen"}==1){
                        print "We found one that matches $r_start and $r_end ==> $start / --$strand-- \n";
                        print "?!?! But we found it before ? THAT IS NOT EXPECTED !\n";
                        die("\n");
                    }
                    my $spacer_seq=$$hash{$spacer}{$start}{$strand}{"seq"};
                    my $proto=$seq;
                    my $exp_len=length($spacer_seq);
                    my $upstream="NA";
                    my $protospacer="NA";
                    my $downstream="NA";
                    if ($proto=~/^(\w{10})(\w{$exp_len})(\w{10})/){
                        $upstream=$1;
                        $protospacer=$2;
                        $downstream=$3;
                    }
                    else{
                        print "We have something unexpected:\n";
                        print $proto."\n";
                        print $spacer_seq."\n";
                        ## Looking if this is truncated in start;
                        if ($r_start==1){
                            ## It is, we can back-calculate by how much by comparing $start to $r_start
                            print "Probably truncated in start -> $start / $r_start\n";
                            my $short_upstream=$start-$r_start;
                            if ($proto=~/^(\w{$short_upstream})(\w{$exp_len})(\w{10})/){
                                $upstream=$1;
                                $protospacer=$2;
                                $downstream=$3;
                            }
                        }
                        elsif ($proto=~/^(\w{10})(\w{$exp_len})(\w*)/){
                            ## It is not truncated at start, so we assume it was truncated in the end, because it was the end of the contig
                            $upstream=$1;
                            $protospacer=$2;
                            $downstream=$3;
                        }
                    }
                    if ($upstream ne "NA" && $downstream ne "NA" && $protospacer ne "NA"){
                        if ($strand eq "minus"){ ## If match on the minus strand, everything must be rev-comped, because by default bowtie rev-comp the input (spacer) but to get the correct upstream/downstream you need to rev-comp the protospacer
                            my $tmp=$downstream;
                            $downstream=&revcomp($upstream);
                            $upstream=&revcomp($tmp);
                            $protospacer=&revcomp($protospacer);
                            $proto=&revcomp($proto);
                            $spacer_seq=&revcomp($spacer_seq); 
                            ## Note -> You may also reverse-completement MD String if you want, but here we opted to keep the original MD from bowtie1
                        }
                        $protospacer=&get_pretty($protospacer,$spacer_seq,$$hash{$spacer}{$start}{$strand}{"n_mis"});
                        if ($quiet!=1){
                            print "\tInfo: ".$$hash{$spacer}{$start}{$strand}{"CIGAR"}." ".$$hash{$spacer}{$start}{$strand}{"n_mis"}." ".$$hash{$spacer}{$start}{$strand}{"md_string"}."\n";
                            print "\tComparing spacers and protospacer:\n";
                            print "\t\t".$proto."\n";
                            print "\t\t".(" " x length($upstream)).$protospacer.(" " x length($downstream))."\n";
                            print "\t\t".(" " x length($upstream)).$spacer_seq.(" " x length($downstream))."\n";
                            print "Upstream: ".$upstream."\n";
                            print "Downstream: ".$downstream."\n";
                        }
                        ## back-calculating end coordinate of hit
                        my $end=$start+length($spacer_seq)-1;
                        ## also padding upstream / downstream if incomplete
                        while(length($upstream)<10){$upstream="N".$upstream;}
                        while(length($downstream)<10){$downstream.="N";}
                        print $sout $spacer."\t".$target."\t".$start."\t".$end."\t".$strand."\t".$$hash{$spacer}{$start}{$strand}{"n_mis"}."\t".
                        $$hash{$spacer}{$start}{$strand}{"CIGAR"}."\t".$$hash{$spacer}{$start}{$strand}{"md_string"}."\t".
                        $spacer_seq."\t".$protospacer."\t".$upstream."\t".$downstream."\n";
                    }
                    else{
                        die("I do no understand protospacer $proto\n");
                    }
                    $$hash{$spacer}{$start}{$strand}{"seen"}=1; 
                }
            }
        }
    }
}

sub get_pretty(){
    my @t_proto=split("",$_[0]);
    my @t_seq=split("",$_[1]);
    my $exp_mis=$_[2];
    my $ret="";
    my $flag=0;
    if ($#t_proto != $#t_seq){$flag=4;}
    for (my $i=0;$i<=$#t_proto;$i++){
        if ($t_proto[$i] eq $t_seq[$i]){$ret.=".";}
        else{$ret.=$t_proto[$i]; $flag++;}
    }
    if ($flag>3){
        print "!!!! THIS SHOULD NOT HAPPEN -> MORE THAN 3 MISMATCHES BETWEEN SPACER AND PROTOSPACER ??\n";
        print $_[0]."\n";
        print $_[1]."\n";
        die("\n");
    }
    elsif ($flag != $exp_mis){
        print "!!!! THIS SHOULD NOT HAPPEN -> WE EXPECT $exp_mis mismatches but we found $flag ??\n";
        print $_[0]."\n";
        print $_[1]."\n";
        die("\n");
    }
    return $ret;
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

sub revcomp(){
	my $nuc=$_[0];
	$nuc=~tr/atcguryswkmbdhvn/tagcayrswmkvhdbn/;
	$nuc=~tr/ATCGURYSWKMBDHVN/TAGCAYRSWMKVHDBN/;
	$nuc=reverse($nuc);
	return $nuc;
}