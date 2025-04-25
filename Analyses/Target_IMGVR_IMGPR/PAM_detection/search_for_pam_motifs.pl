#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use File::Basename;
use Cwd;
my $h=0;
my $in_file="";
my $out_file="";
GetOptions ('help' => \$h, 'h' => \$h, 'i=s'=>\$in_file, 'o=s'=>\$out_file);
if ($h==1 || $in_file eq "" || $out_file eq ""){ # If asked for help or did not set up any argument
	print "# Script to look for PAM motifs from the export of DuckDB
# Arguments :
# -i: input file (default: Neighborhood_for_pam_motif_all.test)
# -o: out_file (default :Motifs_from_neighborhood_all.test)
## Run -> to run it
";
	die "\n";
}

my $tmp_file="tmp_meme/tmp.fna";
my $tmp_meme="tmp_meme/tmp_output";

my %store;
my $c_id="";
my $i=0;
my %tmp;
open my $s1,">",$out_file;
print $s1 "Array\tSide\tMotif_name\tMotif\tMotif_clean\tN_obs\tP_obs\n";
open my $tsv,"<",$in_file;
while(<$tsv>){
	chomp($_);
	my @tab=split("\t",$_);
	if ($tab[0] eq "crispr_array"){next;}
	if ($c_id ne $tab[0]){
		if ($c_id ne ""){
			print "## $c_id ##\n";
			print "up\n";
			my $motif=&find_motif($tmp{"up"});
            print $s1 $c_id."\tupstream\t".$motif."\n";
			print "down\n";
			my $motif=&find_motif($tmp{"do"});
            print $s1 $c_id."\tdownstream\t".$motif."\n";
		}
		$c_id=$tab[0];
		%tmp=();
		$i=0;
	}
	$i++;

	$tmp{"up"}{$i}=substr($tab[1],-5); ## last 5bp upstream
	$tmp{"do"}{$i}=substr($tab[2],0,5); ## first 5bp downstream
}
close $tsv;
if ($c_id ne ""){
	print "## $c_id ##\n";
	print "up\n";
    my $motif=&find_motif($tmp{"up"});
    print $s1 $c_id."\tupstream\t".$motif."\n";
    print "down\n";
    my $motif=&find_motif($tmp{"do"});
    print $s1 $c_id."\tdownstream\t".$motif."\n";
}
close $s1;


sub find_motif(){
	my $hash=$_[0];
    my $return="NA";
	## Make the fasta
    my $n=0;
	open my $s1,">",$tmp_file;
	foreach my $seq (sort keys %{$hash}){
		print $s1 ">Seq_".$seq."\n";
		print $s1 $$hash{$seq}."\n";
        $n++;
	}
	close $s1;
    if ($n<10){
        return "less_than_10_observations";
    }
	## Run meme
    print "meme $tmp_file -oc $tmp_meme -dna -p 12 -minw 5 -maxw 5\n";
	&run_cmd("meme $tmp_file -oc $tmp_meme -dna -p 12 -minw 5 -maxw 5 2>&1","veryquiet"); ## minw5 and maxw5, the motif should always be the full sequences we provide as input, essentially
	## Parse xml
	my %tmp;
	open my $xml,"<",$tmp_meme."/meme.xml";
	while(<$xml>){
		chomp($_);
	      if ($_=~/<training_set .* primary_count=\"(\d+)\"/){
	            $tmp{"primary_count"}=$1;
	      }
	      if($_=~/<motif id="([^\"]+)" name="([^\"]+)" .* sites="(\d+)" /){
	            my $id=$1;
	            my $motif=$2;
	            my $count=$3;
	            $tmp{"motifs"}{$id}{"motif"}=$motif;
	            $tmp{"motifs"}{$id}{"count"}=$count;
	      }
	}
	close $xml;
	foreach my $id (sort keys %{$tmp{"motifs"}}){
	      $tmp{"motifs"}{$id}{"clean_motif"}=$tmp{"motifs"}{$id}{"motif"};
	      $tmp{"motifs"}{$id}{"clean_motif"}=~s/[^ACGT]/N/g;
	      print $id."\t".$tmp{"motifs"}{$id}{"motif"}."\t".$tmp{"motifs"}{$id}{"clean_motif"}."\t".$tmp{"motifs"}{$id}{"count"}."\t".($tmp{"motifs"}{$id}{"count"}/$tmp{"primary_count"}*100)."%\n";
          if ($return ne "NA"){
            die("Multiple motifs ???\n");
          }
          $return=$id."\t".$tmp{"motifs"}{$id}{"motif"}."\t".$tmp{"motifs"}{$id}{"clean_motif"}."\t".$tmp{"motifs"}{$id}{"count"}."\t".($tmp{"motifs"}{$id}{"count"}/$tmp{"primary_count"}*100);
	}
	## Remove directory
	&run_cmd("rm -r $tmp_meme");
    ## return
    return $return;
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
