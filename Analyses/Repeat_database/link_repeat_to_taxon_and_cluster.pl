#!/usr/bin/env perl
use strict;
use autodie;
use Custom::Utils;
use Getopt::Long;
my $h='';
my $cmd='';
my $out='';
my $batch_n="";
GetOptions ('help' => \$h, 'h' => \$h);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
	print "# Script to link repeat to cluster and to taxon (as much as possible)
# Arguments :
# none
";
	die "\n";
}



my %info_contig;
my %info_genome;
my %info_repeat;
my %confidence_tier=("3"=>"High-confidence","2"=>"Medium-confidence","1"=>"Low-confidence","0"=>"No-information");

my $info_repeat_file="all_repeats.tsv";
my $ori_fna="all_repeats.fna";
my $repeat_dedup="all_repeats_dedup.fna.clstr";
my $info_genome_file="all_genome_info.tsv";

my $out_dedup_repeats_tab="all_repeats_dedup.tab";
my $out_dedup_fna="all_repeats_dedup-custom-reps.fna";
my $out_dedup_log="all_repeats_dedup-custom-reps.log";
my $out_repeat_clstr_fna="all_repeats_dedup_nr.fna";
my $out_repeat_clstr="all_repeats_dedup_nr.fna.clstr";
my $out_final_reps="all_repeats_clstr.tab";
my $out_final_reps_log="all_repeats_clstr.log";

## Load repeat info
my %check_c;
print "Reading $info_repeat_file\n";
open my $tsv,"<",$info_repeat_file;
while(<$tsv>){
      chomp($_);
      my @tab=split("\t",$_);
      if ($tab[0] eq "Contig"){next;}
	my $repeat=$tab[1];
	my $contig="";
	if($repeat=~/^(.*)-\d+_\d+/){
		$contig=$1;
	}
	elsif($repeat=~/^(.*)_\d+/){
		$contig=$1;
	}
	else{
		die("Pblm, I don't understand $contig\n$_\n");
	}
	$check_c{$contig}=1;
	$info_repeat{$repeat}{"contig"}=$contig;
      $info_repeat{$repeat}{"consensus"}=$tab[4];
      $info_repeat{$repeat}{"trusted"}=$tab[11];
	$info_repeat{$repeat}{"n_repeat"}=$tab[5];
      if ($tab[14]>=0.8){ ### NOTE - 0.8 is totally arbitrary for now, but it's to assign a type which is "fine" for now (i.e. it does not impact the repeats themselves), default seems to be 0.75 but unclear why some with 0.75 get assigned and other don't
            $info_repeat{$tab[1]}{"type"}=$tab[13];
      }
	else{
		$info_repeat{$tab[1]}{"type"}="Unknown";
	}
}
close $tsv;

## Load repeat seq
my $c_c="";
my %store_seq;
print "Loading repeat sequences from $ori_fna\n";
open my $fna,"<",$ori_fna;
while(<$fna>){
	chomp($_);
	if ($_=~/^>(\S+)/){$c_c=$1;}
	else{$store_seq{$c_c}.=$_;}
}
close $fna;

my %info_contigs;
print "Reading $info_genome_file for information on contigs and genomes\n";
open my $tsv,"<",$info_genome_file;
while(<$tsv>){
	chomp($_);
	my @tab=split("\t",$_);
	if ($tab[0] eq "Contig"){next;}
	$info_contig{$tab[0]}{"genome"}=$tab[1];
	$info_contig{$tab[0]}{"length"}=$tab[2];
	$info_contig{$tab[0]}{"type"}=$tab[3];
	$info_contig{$tab[0]}{"taxo"}=$tab[4];
}
close $tsv;

## Add taxo to each repeat
print "Transfer genome and taxonomy information to repeats, and calculating the quality tier\n";
foreach my $repeat (keys %info_repeat){
      my $contig=$info_repeat{$repeat}{"contig"};
	$info_repeat{$repeat}{"genome"}=$info_contig{$contig}{"genome"};
	$info_repeat{$repeat}{"contig_length"}=$info_contig{$contig}{"length"};
	$info_repeat{$repeat}{"taxo"}=$info_contig{$contig}{"taxo"};
	$info_repeat{$repeat}{"confidence_tier"}=0;
	if ($info_contig{$contig}{"type"} eq "Bin"){
		if ($info_repeat{$repeat}{"contig_length"}>=500000){ $info_repeat{$repeat}{"confidence_tier"}=2;}
		else{ $info_repeat{$repeat}{"confidence_tier"}=1;}
	}
	elsif($info_contig{$contig}{"type"} eq "Meta"){
		$info_repeat{$repeat}{"confidence_tier"}=0;
	}
	else{ $info_repeat{$repeat}{"confidence_tier"}=3;}
}


## Deduplicate - store dedup in separate hash
print "Processing the first (dedupe) clustering to select the right representative\n";
open my $s1,">",$out_dedup_repeats_tab;
print $s1 "Rep\tConfidence_tier\tN_arrays\tN_genome\tGenome_list\tTaxo_LCA\tTotal_repeats_all_arrays\tList_arrays\n";
open my $s2,">",$out_dedup_fna;
open my $s3,">",$out_dedup_log;
my %tmp;
my %info_dedup_repeat;
my %check_rep;
print "Reading $repeat_dedup\n";
open my $fts,"<",$repeat_dedup;
while(<$fts>){
	chomp($_);
      if($_=~/^>/){
		&process_first_step(\%tmp,$s1,$s2,$s3);
		%tmp=();
	}
      else{
		my @tab=split("\t",$_);
		if ($tab[1]=~/(\d+)nt, >(\S+).\.\./){
			$tmp{$2}{"len"}=$1;
		}
	}
}
close $fts;
&process_first_step(\%tmp,$s1,$s2,$s3);
close $s1;
close $s2;
close $s3;

&run_cmd("cd-hit-est -M 10000 -d 0 -c 1 -i $out_dedup_fna -o $out_repeat_clstr_fna");

### Identify which clusters must be refined
my %selected_rep;
# open my $s5,">",$out_list_to_refine;
print "Parsing the second clustering to select the best rep\n";
open my $s4,">",$out_final_reps_log;
open my $fts,"<",$out_repeat_clstr;
my $ntmp=0;
while(<$fts>){
      chomp($_);
	if($_=~/^>/){
		$ntmp++;
		&process_second_step(\%tmp,$s4,$ntmp);
		%tmp=();
	}
	else{
		my @tab=split("\t",$_);
		if ($tab[1]=~/(\d+)nt, >(\S+).\.\./){
			$tmp{$2}{"len"}=$1;
		}
	}
}
&process_second_step(\%tmp,$s4,$ntmp);
close $fts;
close $s4;
# close $s5;

## Now preparing the output -> full table of repr with cluster names
open my $sf,">",$out_final_reps;
print $sf "Array\tCluster id\tConsensus sequence\tType\tTotal repeats\tTotal arrays\tTotal genomes\tGTDB taxonomy\tTaxonomy confidence\tList genomes\tList arrays\n";
my $n_cluster=0;
foreach my $rep (sort keys %selected_rep){
	$n_cluster++;
	my $name_cluster="Ac_".sprintf("%05d",$n_cluster);
	print $sf $rep."\t".$name_cluster."\t".$store_seq{$rep}."\t".$info_repeat{$rep}{"type"}."\t".$info_dedup_repeat{$rep}{"total_repeats"}."\t".$info_dedup_repeat{$rep}{"total_array"}."\t".$info_dedup_repeat{$rep}{"total_genome"}."\t".
	$info_dedup_repeat{$rep}{"taxo"}."\t".$confidence_tier{$info_repeat{$rep}{"confidence_tier"}}."\t".$info_dedup_repeat{$rep}{"genome_list"}."\t".$info_dedup_repeat{$rep}{"array_list"}."\n";
}
close $sf;


sub process_first_step{
	my %tmp=%{$_[0]};
	my $s1=$_[1];
	my $s2=$_[2];
	my $s3=$_[3];
	my $test=scalar(keys %tmp);
	if ($test<1){}
	elsif ($test==1){
		my $rep=(keys %tmp)[0];
		$info_dedup_repeat{$rep}=$info_repeat{$rep};
		print $s1 $rep."\t".$confidence_tier{$info_repeat{$rep}{"confidence_tier"}}."\t1\t1\t".$info_repeat{$rep}{"genome"}."\t".$info_repeat{$rep}{"taxo"}."\t".$info_repeat{$rep}{"n_repeat"}."\t".$rep."\n";
		print $s2 ">".$rep."\n".$store_seq{$rep}."\n";
		$info_dedup_repeat{$rep}{"total_array"}=1;
		$info_dedup_repeat{$rep}{"array_list"}=$rep;
		$info_dedup_repeat{$rep}{"total_genome"}=1;
		$info_dedup_repeat{$rep}{"genome_list"}=$info_repeat{$rep}{"genome"};
		$info_dedup_repeat{$rep}{"total_repeats"}=$info_repeat{$rep}{"n_repeat"};
		$info_dedup_repeat{$rep}{"taxo"}=$info_repeat{$rep}{"taxo"};
	}
	elsif($test>1){
		print "Checking a cluster\n";
		print $s3 "##### New cluster\n";
		my ($rep,$n_genome)=&get_rep(\%tmp,$s3);
		$info_dedup_repeat{$rep}=$info_repeat{$rep};
		my @list_arrays=split(";",$info_repeat{$rep}{"arrays"});
		my $tmp_total=0;
		foreach my $repeat (@list_arrays){
			print $repeat."\t".$info_repeat{$repeat}{"n_repeat"}."\n";
			$tmp_total+=$info_repeat{$repeat}{"n_repeat"};
			my $genome=$info_repeat{$repeat}{"genome"};
			# print $genome."\t".$info_repeat{$repeat}{"taxo"}."\t".$info_genome{$genome}{"type"}."\t".$info_repeat{$repeat}{"contig_length"}."\n";
		}
		if ($tmp_total==0){
			print $rep."\t".$confidence_tier{$info_repeat{$rep}{"confidence_tier"}}."\t".scalar(@list_arrays)."\t".$n_genome."\t".$info_repeat{$rep}{"genome"}."\t".
			$info_repeat{$rep}{"taxo"}."\t".$tmp_total."\t".$info_repeat{$rep}{"arrays"}."\n";
			print "Pblm, there should never be any case with 0 repeat\n";
			die("\n");
		}
		print $s1 $rep."\t".$confidence_tier{$info_repeat{$rep}{"confidence_tier"}}."\t".scalar(@list_arrays)."\t".$n_genome."\t".$info_repeat{$rep}{"genome"}."\t".
		$info_repeat{$rep}{"taxo"}."\t".$tmp_total."\t".$info_repeat{$rep}{"arrays"}."\n";
		print $s2 ">".$rep."\n".$store_seq{$rep}."\n";
		$info_dedup_repeat{$rep}{"total_array"}=scalar(@list_arrays);
		$info_dedup_repeat{$rep}{"array_list"}=$info_repeat{$rep}{"arrays"};
		$info_dedup_repeat{$rep}{"total_genome"}=$n_genome;
		$info_dedup_repeat{$rep}{"genome_list"}=$info_repeat{$rep}{"genome"};
		$info_dedup_repeat{$rep}{"total_repeats"}=$tmp_total;
		$info_dedup_repeat{$rep}{"taxo"}=$info_repeat{$rep}{"taxo"};
		# <STDIN>;
	}
}


sub process_second_step(){
	my %tmp=%{$_[0]};
	my $s4=$_[1];
	my $ntmp=$_[3];
	my $test=scalar(keys %tmp);
	if ($test<1){}
	elsif ($test==1){
		my $rep=(keys %tmp)[0];
		print "This is a single rep, all good\n";
		$selected_rep{$rep}=1;
	}
	elsif($test>1){
		print "Checking cluster $ntmp\n";
		print $s4 "#### New cluster $ntmp\n";
		my %tmp2;
		## We filter the NAs out if there are any
		my %filtered_tmp;
		print "########## STEP 2 ##########\n";
		foreach my $repeat (keys %tmp){
			print $repeat."\t".$info_repeat{$repeat}{"confidence_tier"}."\t".$info_repeat{$repeat}{"n_repeat"}."\t".$info_repeat{$repeat}{"taxo"}."\n";
			if ($info_repeat{$repeat}{"taxo"} ne "NA"){$filtered_tmp{$repeat}=$tmp{$repeat};}
		}
		if (scalar(keys %filtered_tmp)>=1){
			%tmp=%filtered_tmp;
			print "## We filtered the list to: ".join(";",keys %filtered_tmp)."\n";
			# <STDIN>;
		}
		else{
			print "## This will be a full NA, so we keep all\n";
		}
		## Sort these by genus, because we otherwise keep them separated if across genus
		foreach my $rep (keys %tmp){
			my @t=split(";",$info_dedup_repeat{$rep}{"taxo"});
			my $genus=join(";",@t[0..5]);
			$tmp2{$genus}{$rep}=1;
			print "\t".$rep."\t".$genus."\n";
			print $s4 "\t".$rep."\t".$genus."\n";
		}
		foreach my $genus (sort keys %tmp2){
			my @sorted_rep=sort { $info_dedup_repeat{$b}{"total_repeats"} <=> $info_dedup_repeat{$a}{"total_repeats"} or $info_dedup_repeat{$b}{"total_genome"} <=> $info_dedup_repeat{$a}{"total_genome"} or $info_dedup_repeat{$b}{"total_array"} <=> $info_dedup_repeat{$a}{"total_array"} or $a cmp $b} keys %{$tmp2{$genus}};
			if (scalar(@sorted_rep)==1){
				print "\tThis is a cluster but we have only one repeat for this genus - $genus, so we keep this one\n";
				print $s4 "\tThis is a cluster but we have only one repeat for this genus - $genus, so we keep this one\n";
				$selected_rep{$sorted_rep[0]}=1;
			}
			else{
				my $code=join(";",@sorted_rep);
				print "\tThis is for genus $genus -> Here are the different repeats\n";
				print $s4 "\tThis is for genus $genus -> Here are the different repeats\n";
				my $i=0;
				my %tmp_consensus;
				foreach my $rep (@sorted_rep){
					$i++;
					print "\t\t$i   ".$rep."\t".$store_seq{$rep}."\t".$info_dedup_repeat{$rep}{"total_repeats"}."\t".$info_dedup_repeat{$rep}{"total_array"}."\t".$info_dedup_repeat{$rep}{"genome_list"}."\n";
					$tmp_consensus{"total_array"}+=$info_dedup_repeat{$rep}{"total_array"};
					$tmp_consensus{"total_repeats"}+=$info_dedup_repeat{$rep}{"total_repeats"};
					foreach my $genome (split(";",$info_dedup_repeat{$rep}{"genome_list"})){
						$tmp_consensus{"genome"}{$genome}++;
					}
					foreach my $array (split(";",$info_dedup_repeat{$rep}{"array_list"})){
						$tmp_consensus{"array"}{$array}++;
					}
					$tmp_consensus{"taxo"}{$info_dedup_repeat{$rep}{"taxo"}}=1;
				}

				my $n_genome=scalar(keys %{$tmp_consensus{"genome"}});
				my $lca_taxo=&get_lca($tmp_consensus{"taxo"});
				my $genome_list=join(";",sort keys %{$tmp_consensus{"genome"}});
				my $array_list=join(";",sort keys %{$tmp_consensus{"array"}});
				print "\t\t with the automated system, we would pick ".$sorted_rep[0]."\n";
				print $s4 "\t\t with the automated system, we would pick ".$sorted_rep[0]."\n";
				my $ratio=$info_dedup_repeat{$sorted_rep[0]}{"total_repeats"}/$tmp_consensus{"total_repeats"};
				if ($ratio>=0.75){
					print "\t\t And this would represent > 75% of the repeats, so we take automatically, we don't review further\n";
					print $s4 "\t\t And this would represent > 75% of the repeats, so we take automatically, we don't review further\n";
					$selected_rep{$sorted_rep[0]}=1;
				}
				else{
					print "\t\t And this would not represent > 75% of the repeats, so we take all of them because we're not sure\n";
					print $s4 "We are not sure so we take all of them\n";
					foreach my $rep (@sorted_rep){
						$selected_rep{$rep}=1;
					}
					# print $s5 $ntmp."\t".$genus."\t".$code."\t".$lca_taxo."\t".$genome_list."\t".$n_genome."\t".$array_list."\t".$tmp_consensus{"total_array"}."\t".$tmp_consensus{"total_repeats"}."\n";
				}
			}
		}
	}
}




sub get_rep(){
	my $hash=$_[0];
	my $sout=$_[1];
	## First, check if we have some taxo. If so, we only consider for rep the ones with taxo
	my @list_repeats=sort {$info_repeat{$b}{"confidence_tier"} <=> $info_repeat{$a}{"confidence_tier"} or $info_repeat{$b}{"n_repeat"} <=> $info_repeat{$a}{"n_repeat"} or $a <=> $b} keys %{$hash};
	my @filtered_repeats;
	print "########## GET REP ##########\n";
	foreach my $repeat (@list_repeats){
		print $repeat."\t".$info_repeat{$repeat}{"confidence_tier"}."\t".$info_repeat{$repeat}{"n_repeat"}."\t".$info_repeat{$repeat}{"taxo"}."\n";
		if ($info_repeat{$repeat}{"taxo"} ne "NA"){push(@filtered_repeats,$repeat);}
	}
	if (scalar(@filtered_repeats)>=1){
		@list_repeats=@filtered_repeats;
		print "## We filtered the list to: ".join(";",@list_repeats)."\n";
		# <STDIN>;
	}
	else{
		print "## This will be a full NA, so we keep all\n";
	}
	# my @list_repeats=sort {$info_repeat{$b}{"confidence_tier"} <=> $info_repeat{$a}{"confidence_tier"} or $info_repeat{$b}{"n_repeat"} <=> $info_repeat{$a}{"n_repeat"} or $a <=> $b} keys %{$hash};
	my $rep=$list_repeats[0];
	my $best_conf=$info_repeat{$rep}{"confidence_tier"};
	my $rep_len=$$hash{$rep}{"len"};
	# print "Ref: $rep / $rep_len\n";
	print "############################ LOOKING FOR A NEW REP #####################################\n";
	print $sout "Ref: $rep / $rep_len\n";
	my %list_genomes;
	my %list_arrays;
	my %list_taxo;
	foreach my $repeat (@list_repeats){
		# print "\t".$repeat."\t".$$hash{$repeat}{"len"}."\t".$info_repeat{$repeat}{"n_repeat"}."\t".$info_repeat{$repeat}{"genome"}."\t".$info_repeat{$repeat}{"taxo"}."\n";
		print $sout "\t".$repeat."\t".$$hash{$repeat}{"len"}."\t".$info_repeat{$repeat}{"n_repeat"}."\t".$info_repeat{$repeat}{"genome"}."\t".$info_repeat{$repeat}{"taxo"}."\t".
		$info_repeat{$repeat}{"contig_length"}."\t".$info_genome{$info_repeat{$repeat}{"genome"}}{"type_raw"}."\t".$info_repeat{$repeat}{"confidence_tier"}."\n";
		if ($$hash{$repeat}{"len"}!=$rep_len){
			die("PBLM  ---- $repeat vs $rep\n");
		}
		## In a first step, only consider CRISPR arrays of the same quality tier as the one you picked as rep (i.e. the "best" one)
		if ($info_repeat{$repeat}{"confidence_tier"} == $best_conf){
			$list_genomes{$info_repeat{$repeat}{"genome"}}=1;
			$list_taxo{$info_repeat{$repeat}{"taxo"}}=1;
			$list_arrays{$repeat}=1;
		}
	}
	print "First LCA with best repeats\n";
	foreach my $tax (keys %list_taxo){print $tax."\n";}
	my $first_lca=&get_lca(\%list_taxo);
	print "\t\t ==> $first_lca\n";
	print "Now adding other repeats if class is consistent\n";
	my @t=split(";",$first_lca);
	my $lca_class=$t[2];
	foreach my $repeat (@list_repeats){
		## Only for cases where the repeat is not as good as the first one
		if ($info_repeat{$repeat}{"confidence_tier"} < $best_conf){
			my @t=split(";",$info_repeat{$repeat}{"taxo"});
			if ($t[2] eq $lca_class){
				print "\tSame class, we also consider $repeat / ".$info_repeat{$repeat}{"taxo"}."\n";
				$list_genomes{$info_repeat{$repeat}{"genome"}}=1;
				$list_taxo{$info_repeat{$repeat}{"taxo"}}=1;
				$list_arrays{$repeat}=1;
			}
			else{
				print "\tDifferent class, we ignore $repeat / ".$info_repeat{$repeat}{"taxo"}."\n";
			}
		}
	}
	print "Now we have a final list of genomes and we put together our refined LCA\n";
	my $list_genomes_cat=join(";",keys %list_genomes);
	my $list_arrays_cat=join(";",keys %list_arrays);
	my $n_genome=scalar(keys %list_genomes);
	$info_repeat{$rep}{"genome"}=$list_genomes_cat;
	$info_repeat{$rep}{"arrays"}=$list_arrays_cat;
	if ($list_taxo{"NA"}==1){
		## Note, this should be only when we have 0 taxo for the whole thing
		my $test=scalar(keys %list_taxo);
		if ($test == 1){
			$info_repeat{$rep}{"taxo"}="NA";
		}
		else{
			print "## Pblm, we have this list of taxo, which includes NAs, should not exist\n";
			print join(";",keys %list_taxo)."\n";
			die("dying now\n");
		}
	}
	else{
		$info_repeat{$rep}{"taxo"}=&get_lca(\%list_taxo);
	}
	$info_repeat{$rep}{"confidence"}=$best_conf;
	print $info_repeat{$rep}{"taxo"}."\n";
	# print $n_genome." genome(s)\tLCA => ".$info_repeat{$rep}{"taxo"}."\n";
	print $sout $n_genome." genome(s)\tLCA => ".$info_repeat{$rep}{"taxo"}."\n";
	# if ($info_repeat{$rep}{"taxo"}=~/s__unclassified/){<STDIN>;}
	# if ($info_repeat{$rep}{"taxo"}=~/g__unclassified/){<STDIN>;}
	return ($rep,$n_genome);
}




sub get_lca(){
	my @list_ranks=("d","p","c","o","f","g","s");
	my $hash=$_[0];
	my %tmp;
	foreach my $taxo (sort keys %{$hash}){
		foreach my $taxon (split(";",$taxo)){
			$taxon=~/^(\w)__/;
			my $lvl=$1;
			$tmp{$lvl}{$taxon}++;
			# print $lvl."\t".$taxon."\n";
		}
	}
	my $tag=0;
	my $lca="";
	foreach my $rank (@list_ranks){
		my @test=keys %{$tmp{$rank}};
		if ($#test==0){
			$lca.=$test[0].";";
		}
		else{
			$lca.=$rank."__unclassified;";
		}
	}
	chop($lca);
	return $lca;
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
