#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use File::Basename;
use Cwd;
my $h=0;
my $data_folder="";
GetOptions ('help' => \$h, 'h' => \$h, 'd=s'=>\$data_folder);
if ($h==1 || $data_folder eq ""){ # If asked for help or did not set up any argument
	print "# Script to associate a predicted PAM motif for each array
# Arguments :
# -d: Data folder (i.e. /.../Data/)
";
	die "\n";
}

my $file_ref_pam=$data_folder."/Additional_data/manually_curated_PAMs_per_type.tsv";
my $in_file_arrays=$data_folder."/Spacer_db/Array_info_filtered_for_db-Apr12-24.tsv";
my $in_file_denovo="Potential_additional_denovo_by_array_combined.tsv";
my $file_nei="Neighborhood_for_pam_motif_all_imgvr.tsv";
my $file_nei_pr="Neighborhood_for_pam_motif_all_imgpr.tsv";

my $out_file="Repeat_to_PAM_assignment.tsv";
my $out_file_ter="Neighborhood_to_PAM_assignment_imgvr.tsv";
my $out_file_ter_pr="Neighborhood_to_PAM_assignment_imgpr.tsv";

if (-e $out_file || -e $out_file_ter || -e $out_file_ter_pr){
    die("$out_file or $out_file_ter or $out_file_ter_pr already exists, I refuse to do anything\n");
}

my %type_to_motif;
my %canonical_motif;
my %info_motif;
print "Loading canonical motifs\n";
open my $tsv,"<",$file_ref_pam;
while(<$tsv>){
	chomp($_);
	my @tab=split("\t",$_);
	if ($tab[0] eq "Type"){next;}
	## We store the link between the CRISPR type and the motif
	$type_to_motif{$tab[0]}{$tab[1]}=1;
	## If we don't already have the info on this motif, we load/compute it
	if (!defined($info_motif{$tab[1]})){
		print $tab[1]."\tcanonical\n";
		my ($regexp_up,$regexp_down)=&motif_to_regexp($tab[1]);
		$info_motif{$tab[1]}{"regexp_up"}=$regexp_up;
		$info_motif{$tab[1]}{"regexp_down"}=$regexp_down;
		$info_motif{$tab[1]}{"size"}=$tab[3];
		foreach my $motif (sort keys %info_motif){
			if ($motif ne $tab[1]){
				my $reg_test=$info_motif{$motif}{"regexp_up"};
				# print $motif."\t".$tab[1]."\t".$regexp_up."\t".$reg_test."\n";
				if ($motif=~/$regexp_up$/){
					print $tab[1]." is a sub-motif of ".$motif."\n";
					$info_motif{$motif}{"submotifs"}{$tab[1]}=1;
				}
				if ($tab[1]=~/$reg_test$/){
					print $motif." is a sub-motif of ".$tab[1]."\n";
					$info_motif{$tab[1]}{"submotifs"}{$motif}=1;
				}
			}
		}
		$canonical_motif{$tab[1]}=1;
	}
}
close $tsv;

my %ref_array;
print "Loading array type information\n";
open my $tsv,"<",$in_file_arrays;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "repeat_cluster"){next;}
    my $array=$tab[0];
    my $type=$tab[1];
    $ref_array{$array}{"type"}=$type;
}
close $tsv;

my %denovo_motif;
print "Loading motif from denovo prediction\n";
open my $tsv,"<",$in_file_denovo;
while(<$tsv>){
	chomp($_);
	my @tab=split("\t",$_);
	$denovo_motif{$tab[0]}{"motif"}=$tab[1];
	if (!defined($info_motif{$tab[1]})){
		$info_motif{$tab[1]}{"regexp_up"}=$tab[3];
		$info_motif{$tab[1]}{"regexp_down"}=$tab[4];
		$info_motif{$tab[1]}{"size"}=$tab[5];
		foreach my $sub (split(/\|/,$tab[6])){
			$info_motif{$tab[1]}{"submotifs"}{$sub}=1;
		}
	}
}
close $tsv;

my $i=0;
my $c_t="";
my $c_a="";
my %store_all_nei;
print "## Based on these motifs, process neighborhoods of VR\n";
open my $s3,">",$out_file_ter;
print $s3 "Array\tType\tUpstream\tDownstream\tMotif\tScore\n";
open my $tsv,"<",$file_nei;
while(<$tsv>){
	chomp($_);
	my @tab=split("\t",$_);
	if ($tab[0] eq "crispr_array"){next;}
	$i++;
	$c_a=$tab[0];
	if (defined($ref_array{$c_a}{"type"})){$c_t=$ref_array{$c_a}{"type"};}
	else{$c_t="Unknown";}
	print "Looking for potential PAMs in $tab[1]\t$tab[2] for type $c_t and potential de novo motif ".$denovo_motif{$c_a}{"motif"}."\n";
	my ($best_motif,$best_score,@list)=&check_pam($tab[1],$tab[2],$c_t,$denovo_motif{$c_a}{"motif"});
	print $tab[1]."\t".$tab[2]."\t".$best_motif."\t".$best_score."\t".join(",",@list)."\t".$denovo_motif{$c_a}{"motif"}."\n";
	$store_all_nei{$tab[0]}{"totals"}{"distinct"}++;
	$store_all_nei{$tab[0]}{"totals"}{"cover"}+=$tab[3];
	if ($best_score>0){
		$store_all_nei{$tab[0]}{"neighs"}{$i}{"n_motif"}=scalar(@list);
		$store_all_nei{$tab[0]}{"neighs"}{$i}{"best_motif"}=$best_motif;
		$store_all_nei{$tab[0]}{"neighs"}{$i}{"cover"}=$tab[3];
		foreach my $motif (@list){
			$store_all_nei{$tab[0]}{"neighs"}{$i}{"pot_motif"}{$motif}=1;
		}
		print $s3 $c_a."\t".$c_t."\t".$tab[1]."\t".$tab[2]."\t".$best_motif."\t".$best_score."\t".join(",",@list)."\n";
		print $c_a."\t".$c_t."\t".$tab[1]."\t".$tab[2]."\t".$best_motif."\t".$best_score."\t".join(",",@list)."\n";
	}
	else{
		print $s3 $c_a."\t".$c_t."\t".$tab[1]."\t".$tab[2]."\tNA\tNA\tNA\n";
	}
	# <STDIN>;
}
close $s3;

print "## Then the neighborhoods of PR\n";
open my $s3,">",$out_file_ter_pr;
print $s3 "Array\tType\tUpstream\tDownstream\tMotif\tScore\n";
open my $tsv,"<",$file_nei_pr;
while(<$tsv>){
	chomp($_);
	my @tab=split("\t",$_);
	if ($tab[0] eq "crispr_array"){next;}
	$i++;
	$c_a=$tab[0];
	if (defined($ref_array{$c_a}{"type"})){$c_t=$ref_array{$c_a}{"type"};}
	else{$c_t="Unknown";}
	print "Looking for potential PAMs in $tab[1]\t$tab[2] for type $c_t and potential de novo motif ".$denovo_motif{$c_a}{"motif"}."\n";
	my ($best_motif,$best_score,@list)=&check_pam($tab[1],$tab[2],$c_t,$denovo_motif{$c_a}{"motif"});
	print $tab[1]."\t".$tab[2]."\t".$best_motif."\t".$best_score."\t".join(",",@list)."\t".$denovo_motif{$c_a}{"motif"}."\n";
	$store_all_nei{$tab[0]}{"totals"}{"distinct"}++;
	$store_all_nei{$tab[0]}{"totals"}{"cover"}+=$tab[3];
	if ($best_score>0){
		$store_all_nei{$tab[0]}{"neighs"}{$i}{"n_motif"}=scalar(@list);
		$store_all_nei{$tab[0]}{"neighs"}{$i}{"best_motif"}=$best_motif;
		$store_all_nei{$tab[0]}{"neighs"}{$i}{"cover"}=$tab[3];
		foreach my $motif (@list){
			$store_all_nei{$tab[0]}{"neighs"}{$i}{"pot_motif"}{$motif}=1;
		}
		print $s3 $c_a."\t".$c_t."\t".$tab[1]."\t".$tab[2]."\t".$best_motif."\t".$best_score."\t".join(",",@list)."\n";
	}
	else{
		print $s3 $c_a."\t".$c_t."\t".$tab[1]."\t".$tab[2]."\tNA\tNA\tNA\n";
	}
}
close $s3;

open my $s1,">",$out_file;
print $s1 "Array\tType\tDe_novo_motif\tTotal_cover\tTotal_counts\tSelected_motif\tMotif\tMotif_cover\tRatio_cover\tMotif_count\tRatio_count\n";
foreach my $c_a (sort keys %store_all_nei){
	my $c_t="";
	if (defined($ref_array{$c_a}{"type"})){$c_t=$ref_array{$c_a}{"type"};}
	else{$c_t="Unknown";}
	my $tmp=$store_all_nei{$c_a};
	print "##########################\t".$c_a."\t".$c_t."\n";
	print "Totals\t".$$tmp{"totals"}{"distinct"}."\t".$$tmp{"totals"}{"cover"}."\n";
	## First, get counts for unique hits
	my $cand_motif=&get_best_cand($$tmp{"neighs"});
	if ($cand_motif ne "NA"){
		my ($total_hits,$total_cover)=&get_full_count($cand_motif,$$tmp{"neighs"});
		my @t=split("_",$cand_motif);
		my $motif=$t[0];
		print $s1 $c_a."\t".$c_t."\t".$denovo_motif{$c_a}{"motif"}."\t".$$tmp{"totals"}{"cover"}."\t".$$tmp{"totals"}{"distinct"}."\t".$motif."\t".$cand_motif."\t".$total_cover."\t".sprintf("%.2f",$total_cover/$$tmp{"totals"}{"cover"}*100)."\t".$total_hits."\t".sprintf("%.2f",$total_hits/$$tmp{"totals"}{"distinct"}*100)."\n";
		### 
		print $c_a."\t".$c_t."\t".$denovo_motif{$c_a}{"motif"}."\t".$$tmp{"totals"}{"cover"}."\t".$$tmp{"totals"}{"distinct"}."\t".$motif."\t".$cand_motif."\t".$total_cover."\t".sprintf("%.2f",$total_cover/$$tmp{"totals"}{"cover"}*100)."\t".$total_hits."\t".sprintf("%.2f",$total_hits/$$tmp{"totals"}{"distinct"}*100)."\n";
	}
	else{
		print "No candidate motif, we skip\n";
		print $s1 $c_a."\t".$c_t."\t".$denovo_motif{$c_a}{"motif"}."\t".$$tmp{"totals"}{"cover"}."\t".$$tmp{"totals"}{"distinct"}."\tNA\tNA\tNA\tNA\tNA\tNA\n";
		print $c_a."\t".$c_t."\t".$denovo_motif{$c_a}{"motif"}."\t".$$tmp{"totals"}{"cover"}."\t".$$tmp{"totals"}{"distinct"}."\tNA\tNA\tNA\tNA\tNA\tNA\n";
	}
}
close $s1;

sub get_full_count(){
	my $motif=$_[0];
	my $hash=$_[1];
	my $count=0;
	my $cover=0;
	foreach my $i (keys %{$hash}){
		if ($$hash{$i}{"n_motif"}>1){
			foreach my $motif_2 (keys %{$$hash{$i}{"pot_motif"}}){
				if ($motif eq $motif_2){
					$count++;
					$cover+=$$hash{$i}{"cover"};
				}
			}
		}
		elsif ($motif eq $$hash{$i}{"best_motif"}){
			$count++;
			$cover+=$$hash{$i}{"cover"};
		}
	}
	return ($count,$cover);
}

sub get_best_cand(){
	my $hash=$_[0];
	my $max=0;
	my %count;
	foreach my $i (keys %{$hash}){
		if ($$hash{$i}{"n_motif"}==1){
			$max++;
			$count{$$hash{$i}{"best_motif"}}++;
		}
	}
	my @t=sort {$count{$b} <=> $count{$a} or $a cmp $b} keys %count;
	if ($max==0){ ## 0 unique candidates, we can't say without manual curation
		return "NA";
	}
	elsif (($count{$t[0]}/$max)>=0.25 && ($count{$t[0]}>$count{$t[1]})){
		## We have a clear best candidate found in more than 25% of the unique hits, and a clear number 1 (i.e. no ties)
		return $t[0]; 
	}
	else{
		return "NA";
	}
}

sub check_pam(){
	# my ($best_motif,$best_score,@list)=&check_pam($tab[1],$tab[2],$c_t,$denovo_motif{$c_a});
	my %tmp;
	my $up=$_[0];
	my $down=$_[1];
	my $type=$_[2];
	my $denovo=$_[3];
	## Always test all canonicals
	foreach my $motif (keys %canonical_motif){
		# print "### $motif ###\n";
		my $score=&perfect_match($up,$info_motif{$motif}{"regexp_up"},"end");
		if ($score>0){
			$tmp{$motif."_up"}=$score;
			print "up\t".$up."\t".$motif."\t".$info_motif{$motif}{"regexp_up"}."\t"."end"."\t".$score."\n";
		}
		$score=&perfect_match($down,$info_motif{$motif}{"regexp_down"},"start");
		if ($score>0){
			$tmp{$motif."_down"}=$score;
			print "down\t".$down."\t".$motif."\t".$info_motif{$motif}{"regexp_down"}."\t"."start"."\t".$score."\n";
		}
	}
	if ($denovo ne ""){
		# print "### $denovo ###\n";
		my $score=&perfect_match($up,$info_motif{$denovo}{"regexp_up"},"end");
		if ($score>0){
			$tmp{$denovo."_up"}=$score;
			print "up\t".$up."\t".$denovo."\t".$info_motif{$denovo}{"regexp_up"}."\t"."end"."\t".$score."\n";
		}
		
		$score=&perfect_match($down,$info_motif{$denovo}{"regexp_down"},"start");
		if ($score>0){
			$tmp{$denovo."_down"}=$score;
			print "down\t".$down."\t".$denovo."\t".$info_motif{$denovo}{"regexp_down"}."\t"."start"."\t".$score."\n";
		}
	}
	# <STDIN>;
	my $n_hits=scalar(keys %tmp);
	if ($n_hits==0){
		print "No potential PAM here, we return nothing\n";
		return("NA",0,());
	}
	else{
		print "Some potential PAM here, we check if one or several, and if several, try to simplify\n";
		if ($n_hits>1){
			print "\tmultiple PAMs, we first simplify by removing all the ones that are not top score\n";
			my @t=sort {$tmp{$b} <=> $tmp{$a} or $a cmp $b} keys %tmp;
			foreach my $motif_found (@t){
				if ($tmp{$motif_found}<$tmp{$t[0]}){
					print "\t\tRemoving $motif_found because its score is $tmp{$motif_found} which is not as good at $tmp{$t[0]} for $t[0]\n";
					delete($tmp{$motif_found});
				}
			}
			print "\tmultiple PAMs, we check if we should simplify based on submotif\n";
			@t=sort {$tmp{$b} <=> $tmp{$a} or $info_motif{$b}{"size"} <=> $info_motif{$a}{"size"} or $a cmp $b} keys %tmp;
			foreach my $motif_found (@t){
				my @t1=split("_",$motif_found);
				my $motif=$t1[0];
				my $orient=$t1[1];
				print "We found motif $motif_found\n";
				if (defined($info_motif{$motif}{"submotifs"})){
					print "\t\t$motif has submotifs..\n";
					foreach my $motif_2 (sort keys %{$info_motif{$motif}{"submotifs"}}){
						my $test=$motif_2."_".$orient;
						print "\t\t\t$motif_2 -- $test\n";
						if (!defined($tmp{$test})){
							print "But the sub-motif was not actually found, probably because it's not canonical and only found in de novo\n";
							if ($canonical_motif{$test}==1){
								print "\t???? ACTUALLY NO, $test is FLAGGED AS CANONICAL ??\n";
								die("\n");
							}
						}
						else{
							print "We remove $test since it's covered by $motif_found\n";
							delete($tmp{$test});
						}
					}
				}
			}
		}
		## Now can take the best hits
		my @r=sort {$tmp{$b} <=> $tmp{$a} or $info_motif{$b}{"size"} <=> $info_motif{$a}{"size"} or $a cmp $b} keys %tmp; 
		my $best_motif=$r[0];
		my $best_score=$tmp{$r[0]};
		# if ($best_motif eq $denovo."_up" || $best_motif eq $denovo."_down"){
			# print "!! Yeah, the best motif is a de novo one \\o\/ \n";
			# <STDIN>;
		# }
		## Take best and list
		my @list;
		if ($best_score>0){
			foreach my $pot (@r){
				if ($tmp{$pot}==$best_score){
					push(@list,$pot);
				}
			}
		}
		print $up."\t".$down."\t".$type."\t".$best_motif."\t".$best_score."\n";
		return($best_motif,$best_score,@list);
	}
}

sub perfect_match(){
	my $seq=$_[0];
	my $motif=$_[1];
	my $location=$_[2];
	if ($location eq "start"){
		if ($seq=~/^$motif/){
			return 2;
		}
		else{
			return 0;
		}
	}
	elsif($location eq "end"){
		if ($seq=~/$motif$/){
			return 2;
		}
		else{
			return 0;
		}
	}
	else{
		return -1;
	}
}

sub motif_to_regexp(){
	my $motif=$_[0];
	my @t=split("",$motif);
	my @r=();
	foreach my $cell (@t){
		if ($cell eq "N"){push(@r,"[ATCGN]");}
		else{push(@r,$cell);}
	}
	my $up=join("",@r);
	@r=();
	for (my $i=$#t;$i>=0;$i--){
		my $cell=$t[$i];
		if ($cell eq "N"){push(@r,"[ATCGN]");}
		else{push(@r,&revcomp($cell));}
	}
	my $down=join("",@r);
	# print $motif."\t".$up."\t".$down."\n";
	# <STDIN>;
	return($up,$down);
}


sub revcomp(){
	my $nuc=$_[0];
	$nuc=~tr/atcguryswkmbdhvn/tagcayrswmkvhdbn/;
	$nuc=~tr/ATCGURYSWKMBDHVN/TAGCAYRSWMKVHDBN/;
	$nuc=reverse($nuc);
	return $nuc;
}





# print "Loading de novo motifs from VR\n";
# open my $tsv,"<",$file_denovo; 
# while(<$tsv>){
# 	chomp($_);
# 	my @tab=split("\t",$_);
# 	if ($tab[0] eq "Array"){next;}
# 	my $n_n=$tab[4]=~tr/N//;
# 	my $n_site=$tab[5];
# 	print "### $tab[0] - $tab[4] - $n_n - $n_site - $tab[6]\n";
# 	if ($n_n<=4){
# 		if ($tab[6] > 50){
# 			## potential motif
# 			print "Potential motif $tab[4] for $tab[0]\n";
# 			if (defined($denovo_motif{$tab[0]})){
# 				print "This is a potential replacement, we check if it should be a replacement\n";
# 				print $tab[6]." vs ".$denovo_motif{$tab[0]}{"perc"}."\n";
# 				if ($tab[6] < $denovo_motif{$tab[0]}{"perc"}){
# 					print "\twe actually skip this new one because it's percent $tab[6] is not as good as ".$denovo_motif{$tab[0]}{"perc"}."\n";
# 					# <STDIN>;
# 					next;
# 				}
# 				if ($n_n==4){
# 					if ($denovo_motif{$tab[0]}{"n_n"}<4){
# 						print "\twe actually skip this new one because it's a mononucleotide while the previous one is longer\n";
# 						# <STDIN>;
# 						next;
# 					}
# 				}
# 				print "We take as a replacement\n";
# 				# <STDIN>;
# 			}
# 			## If we get there, that means we want this motif
# 			my $motif="";
# 			if($tab[1] eq "upstream"){}
# 			else{
# 				$tab[4]=&revcomp($tab[4]);
# 			}
# 			my @t=split("",$tab[4]);
# 			# my $n_n=$motif=~tr/N//;
# 			my $motif=$tab[4]; ## We take the 5 bases
# 			my $check_motif=$t[2].$t[3].$t[4]; ## to check canonical
# 			my $check_motif_bis=$t[1].$t[2].$t[3].$t[4]; ## to check canonical
# 			if (defined($canonical_motif{$motif}) || defined($canonical_motif{$check_motif}) || defined($canonical_motif{$check_motif_bis})){
# 				print "motif $motif is already in the list of canonical, we skip\n";
# 				# <STDIN>;
# 			}
# 			else{
# 				$denovo_motif{$tab[0]}{"motif"}=$motif;
# 				$denovo_motif{$tab[0]}{"n_obs"}=$tab[5];
# 				$denovo_motif{$tab[0]}{"perc"}=$tab[6];
# 				$denovo_motif{$tab[0]}{"n_n"}=$n_n;
# 				$denovo_motif{$tab[0]}{"from_vr"}=1;
# 				if (defined($info_motif{$motif})){
# 					print "we already know about motif $motif, so we just link it to this array\n";
# 				}
# 				else{
# 					print "$motif is new, we take its characteristics\n";
# 					my ($regexp_up,$regexp_down)=&motif_to_regexp($motif);
# 					$info_motif{$motif}{"regexp_up"}=$regexp_up;
# 					$info_motif{$motif}{"regexp_down"}=$regexp_down;
# 					$info_motif{$motif}{"size"}=3-$n_n;
# 					foreach my $motif_2 (sort keys %info_motif){
# 						if ($motif ne $motif_2){
# 							my $reg_test=$info_motif{$motif_2}{"regexp_up"};
# 							# print $motif_2."\t".$motif."\t".$regexp_up."\t".$reg_test."\n";
# 							if ($motif_2=~/$regexp_up$/){
# 								# print $motif." is a sub-motif of ".$motif_2."\n";
# 								$info_motif{$motif_2}{"submotifs"}{$motif}=1;
# 								# <STDIN>;
# 							}
# 							if ($motif=~/$reg_test$/){
# 								# print $motif_2." is a sub-motif of ".$motif."\n";
# 								$info_motif{$motif}{"submotifs"}{$motif_2}=1;
# 								# <STDIN>;
# 							}
# 						}
# 					}
# 				}
# 			}
# 		}
# 	}
# }
# close $tsv;
#
#
# print "Loading additional de novo motifs from PR if we have some\n";
# open my $tsv,"<",$file_denovo_pr; 
# while(<$tsv>){
# 	chomp($_);
# 	my @tab=split("\t",$_);
# 	if ($tab[0] eq "Array"){next;}
# 	my $n_n=$tab[4]=~tr/N//;
# 	my $n_site=$tab[5];
# 	print "### $tab[0] - $tab[4] - $n_n - $n_site - $tab[6]\n";
# 	if ($n_n<=4){
# 		if ($tab[6] > 50){
# 			## potential motif
# 			print "Potential motif $tab[4] for $tab[0]\n";
# 			if (defined($denovo_motif{$tab[0]})){
# 				print "This is a potential replacement, we check if it should be a replacement\n";
# 				if ($denovo_motif{$tab[0]}{"from_vr"}==1){
# 					print "We had a VR motif for $tab[0]\n";
# 					print $tab[5]." vs ".$denovo_motif{$tab[0]}{"n_obs"}."\n";
# 					print $tab[6]." vs ".$denovo_motif{$tab[0]}{"perc"}."\n";
# 					if ($tab[5] < $denovo_motif{$tab[0]}{"n_obs"}){
# 						print "We had more observations in VR, so we trust VR more\n";
# 						next;
# 					}
# 					else{
# 						print "We have more observations in PR, so we look into this\n";
# 						# <STDIN>;
# 					}
# 				}
# 				if ($tab[6] < $denovo_motif{$tab[0]}{"perc"}){
# 					print "\twe actually skip this new one because it's percent $tab[6] is not as good as ".$denovo_motif{$tab[0]}{"perc"}."\n";
# 					# <STDIN>;
# 					next;
# 				}
# 				if ($n_n==4){
# 					if ($denovo_motif{$tab[0]}{"n_n"}<4){
# 						print "\twe actually skip this new one because it's a mononucleotide while the previous one is longer\n";
# 						# <STDIN>;
# 						next;
# 					}
# 				}
# 				print "We take as a replacement\n";
# 			}
# 			## If we get there, that means we want this motif
# 			my $motif="";
# 			if($tab[1] eq "upstream"){}
# 			else{
# 				$tab[4]=&revcomp($tab[4]);
# 			}
# 			my @t=split("",$tab[4]);
# 			# my $n_n=$motif=~tr/N//;
# 			my $motif=$tab[4]; ## We take the 5 bases
# 			my $check_motif=$t[2].$t[3].$t[4]; ## to check canonical
# 			my $check_motif_bis=$t[1].$t[2].$t[3].$t[4]; ## to check canonical
# 			if (defined($canonical_motif{$motif}) || defined($canonical_motif{$check_motif}) || defined($canonical_motif{$check_motif_bis})){
# 				print "motif $motif is already in the list of canonical, we skip\n";
# 				# <STDIN>;
# 			}
# 			else{
# 				$denovo_motif{$tab[0]}{"motif"}=$motif;
# 				$denovo_motif{$tab[0]}{"n_obs"}=$tab[5];
# 				$denovo_motif{$tab[0]}{"perc"}=$tab[6];
# 				$denovo_motif{$tab[0]}{"n_n"}=$n_n;
# 				$denovo_motif{$tab[0]}{"from_vr"}=2;
# 				if (defined($info_motif{$motif})){
# 					print "we already know about motif $motif, so we just link it to this array\n";
# 				}
# 				else{
# 					print "$motif is new, we take its characteristics\n";
# 					my ($regexp_up,$regexp_down)=&motif_to_regexp($motif);
# 					$info_motif{$motif}{"regexp_up"}=$regexp_up;
# 					$info_motif{$motif}{"regexp_down"}=$regexp_down;
# 					$info_motif{$motif}{"size"}=3-$n_n;
# 					foreach my $motif_2 (sort keys %info_motif){
# 						if ($motif ne $motif_2){
# 							my $reg_test=$info_motif{$motif_2}{"regexp_up"};
# 							# print $motif_2."\t".$motif."\t".$regexp_up."\t".$reg_test."\n";
# 							if ($motif_2=~/$regexp_up$/){
# 								# print $motif." is a sub-motif of ".$motif_2."\n";
# 								$info_motif{$motif_2}{"submotifs"}{$motif}=1;
# 								# <STDIN>;
# 							}
# 							if ($motif=~/$reg_test$/){
# 								# print $motif_2." is a sub-motif of ".$motif."\n";
# 								$info_motif{$motif}{"submotifs"}{$motif_2}=1;
# 								# <STDIN>;
# 							}
# 						}
# 					}
# 				}
# 			}
# 		}
# 	}
# }
# close $tsv;



# sub fuzzy_match(){
# 	my $seq=$_[0];
# 	my $motif=$_[1];
# 	my $location=$_[2];
# 	if ($location eq "start"){
# 		if ($seq=~/^$motif/){
# 			return 2;
# 		}
# 		elsif($seq=~/^[ATCGN]$motif/){
# 			return 1;
# 		}
# 		else{
# 			return 0;
# 		}
# 	}
# 	elsif($location eq "end"){
# 		if ($seq=~/$motif$/){
# 			return 2;
# 		}
# 		elsif($seq=~/$motif[ATCGN]$/){
# 			return 1;
# 		}
# 		else{
# 			return 0;
# 		}
# 	}
# 	else{
# 		return -1;
# 	}
# }


# sub simplify(){
# 	my @in_tab=@{$_[0]};
# 	my %tmp;
# 	foreach my $motif (@in_tab){

# 	}
# }

# sub check_pam(){
# 	my %tmp;
# 	my $up=$_[0];
# 	my $down=$_[1];
# 	my $type=$_[2];
# 	my $hash=$_[3];
# 	my $denovo=$_[4];
# 	if ($type ne ""){
# 		foreach my $motif (keys %{$$hash{$type}}){
# 			# print "### $motif ###\n";
# 			$tmp{$motif."_down"}=&fuzzy_match($down,$$hash{$type}{$motif}{"regexp_down"},"start");
# 			# print "down\t".$down."\t".$$hash{$type}{$motif}{"regexp_down"}."\t"."start"."\t".$tmp{$motif."_down"}."\n";
# 			$tmp{$motif."_up"}=&fuzzy_match($up,$$hash{$type}{$motif}{"regexp_up"},"end");
# 			# print "up\t".$up."\t".$$hash{$type}{$motif}{"regexp_up"}."\t"."end"."\t".$tmp{$motif."_up"}."\n";
# 		}
# 	}
# 	my @r=sort {$tmp{$b} <=> $tmp{$a} or $sort_motif{$b} <=> $sort_motif{$a} or $a cmp $b} keys %tmp;
# 	if ($tmp{$r[0]}==0){
# 		## No match to expected PAM, look for other pam, including the de novo one
# 		foreach my $type (sort keys %{$hash}){
# 			foreach my $motif (keys %{$$hash{$type}}){
# 				# print "### $motif ###\n";
# 				$tmp{$motif."_down"}=&fuzzy_match($down,$$hash{$type}{$motif}{"regexp_down"},"start");
# 				# print "down\t".$down."\t".$$hash{$type}{$motif}{"regexp_down"}."\t"."start"."\t".$tmp{$motif."_down"}."\n";
# 				$tmp{$motif."_up"}=&fuzzy_match($up,$$hash{$type}{$motif}{"regexp_up"},"end");
# 				# print "up\t".$up."\t".$$hash{$type}{$motif}{"regexp_up"}."\t"."end"."\t".$tmp{$motif."_up"}."\n";
# 			}
# 		}
# 		if(defined($$denovo{"motif_ori"})){
# 			my $motif=$$denovo{"motif_ori"};
# 			# print "### $motif ###\n";
# 			$tmp{$motif."_down"}=&fuzzy_match($down,$$denovo{"regexp_down"},"start");
# 			# print "down\t".$down."\t".$$denovo{"regexp_down"}."\t"."start"."\t".$tmp{$motif."_down"}."\n";
# 			$tmp{$motif."_up"}=&fuzzy_match($up,$$denovo{"regexp_up"},"end");
# 			# print "up\t".$up."\t".$$denovo{"regexp_up"}."\t"."end"."\t".$tmp{$motif."_up"}."\n";
# 		}
# 		## Include the de novo one
# 		@r=sort {$tmp{$b} <=> $tmp{$a} or $sort_motif{$b} <=> $sort_motif{$a} or $a cmp $b} keys %tmp;
# 	}
# 	my $best_motif=$r[0];
# 	my $best_score=$tmp{$r[0]};
# 	my @list;
# 	foreach my $pot (@r){
# 		if ($tmp{$pot}==$best_score){
# 			push(@list,$pot);
# 		}
# 	}
# 	# print $up."\t".$down."\t".$type."\t".$best_motif."\t".$best_score."\n";
# 	return($best_motif,$best_score,@list);
# }