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
	print "# Script to combine the motifs detected from IMG/VR and IMG/PR hits
# Arguments :
# -d: Data folder (i.e. /.../Data/)
";
	die "\n";
}

my $out_file="Potential_additional_denovo_by_array_combined.tsv";

my $file_ref_pam=$data_folder."/Additional_data/manually_curated_PAMs_per_type.tsv";
my $in_file_arrays=$data_folder."/Spacer_db/Array_info_filtered_for_db-Apr12-24.tsv";
my $file_denovo="Motifs_from_neighborhood_all_imgvr.tsv";
my $file_denovo_pr="Motifs_from_neighborhood_all_imgpr.tsv";

if (-e $out_file){
    die("$out_file already exists, I refuse to do anything\n");
}


my %canonical_motif;
my %info_motif;
print "Loading canonical motifs\n";
open my $tsv,"<",$file_ref_pam;
while(<$tsv>){
	chomp($_);
	my @tab=split("\t",$_);
	if ($tab[0] eq "Type"){next;}
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
print "Loading de novo motifs from VR\n";
open my $tsv,"<",$file_denovo; 
while(<$tsv>){
	chomp($_);
	my @tab=split("\t",$_);
	if ($tab[0] eq "Array"){next;}
	my $n_n=$tab[4]=~tr/N//;
	my $n_site=$tab[5];
	print "### $tab[0] - $tab[4] - $n_n - $n_site - $tab[6]\n";
	if ($n_n<=4){
		if ($tab[6] > 50){
			## potential motif
			print "Potential motif $tab[4] for $tab[0]\n";
			if (defined($denovo_motif{$tab[0]})){
				print "This is a potential replacement, we check if it should be a replacement\n";
				print $tab[6]." vs ".$denovo_motif{$tab[0]}{"perc"}."\n";
				if ($tab[6] < $denovo_motif{$tab[0]}{"perc"}){
					print "\twe actually skip this new one because it's percent $tab[6] is not as good as ".$denovo_motif{$tab[0]}{"perc"}."\n";
					next;
				}
				if ($n_n==4){
					if ($denovo_motif{$tab[0]}{"n_n"}<4){
						print "\twe actually skip this new one because it's a mononucleotide while the previous one is longer\n";
						next;
					}
				}
				print "We take as a replacement\n";
				# <STDIN>;
			}
			## If we get there, that means we want this motif
			my $motif="";
			if($tab[1] eq "upstream"){}
			else{
				$tab[4]=&revcomp($tab[4]);
			}
			my @t=split("",$tab[4]);
			# my $n_n=$motif=~tr/N//;
			my $motif=$tab[4]; ## We take the 5 bases
			my $check_motif=$t[2].$t[3].$t[4]; ## to check canonical
			my $check_motif_bis=$t[1].$t[2].$t[3].$t[4]; ## to check canonical
			if (defined($canonical_motif{$motif}) || defined($canonical_motif{$check_motif}) || defined($canonical_motif{$check_motif_bis})){
				print "motif $motif is already in the list of canonical, we skip\n";
				# <STDIN>;
			}
			else{
				$denovo_motif{$tab[0]}{"motif"}=$motif;
				$denovo_motif{$tab[0]}{"n_obs"}=$tab[5];
				$denovo_motif{$tab[0]}{"perc"}=$tab[6];
				$denovo_motif{$tab[0]}{"n_n"}=$n_n;
				$denovo_motif{$tab[0]}{"from_vr"}=1;
				if (defined($info_motif{$motif})){
					print "we already know about motif $motif, so we just link it to this array\n";
				}
				else{
					print "$motif is new, we take its characteristics\n";
					my ($regexp_up,$regexp_down)=&motif_to_regexp($motif);
					$info_motif{$motif}{"regexp_up"}=$regexp_up;
					$info_motif{$motif}{"regexp_down"}=$regexp_down;
					$info_motif{$motif}{"size"}=5-$n_n;
					foreach my $motif_2 (sort keys %info_motif){
						if ($motif ne $motif_2){
							my $reg_test=$info_motif{$motif_2}{"regexp_up"};
							# print $motif_2."\t".$motif."\t".$regexp_up."\t".$reg_test."\n";
							if ($motif_2=~/$regexp_up$/){
								# print $motif." is a sub-motif of ".$motif_2."\n";
								$info_motif{$motif_2}{"submotifs"}{$motif}=1;
								# <STDIN>;
							}
							if ($motif=~/$reg_test$/){
								# print $motif_2." is a sub-motif of ".$motif."\n";
								$info_motif{$motif}{"submotifs"}{$motif_2}=1;
								# <STDIN>;
							}
						}
					}
				}
			}
		}
	}
}
close $tsv;

print "Loading additional de novo motifs from PR if we have some\n";
open my $tsv,"<",$file_denovo_pr; 
while(<$tsv>){
	chomp($_);
	my @tab=split("\t",$_);
	if ($tab[0] eq "Array"){next;}
	my $n_n=$tab[4]=~tr/N//;
	my $n_site=$tab[5];
	print "### $tab[0] - $tab[4] - $n_n - $n_site - $tab[6]\n";
	if ($n_n<=4){
		if ($tab[6] > 50){
			## potential motif
			print "Potential motif $tab[4] for $tab[0]\n";
			if (defined($denovo_motif{$tab[0]})){
				print "This is a potential replacement, we check if it should be a replacement\n";
				if ($denovo_motif{$tab[0]}{"from_vr"}==1){
					print "We had a VR motif for $tab[0]\n";
					print $tab[5]." vs ".$denovo_motif{$tab[0]}{"n_obs"}."\n";
					print $tab[6]." vs ".$denovo_motif{$tab[0]}{"perc"}."\n";
					if ($tab[5] < $denovo_motif{$tab[0]}{"n_obs"}){
						print "We had more observations in VR, so we trust VR more\n";
						next;
					}
					else{
						print "We have more observations in PR, so we look into this\n";
						# <STDIN>;
					}
				}
				if ($tab[6] < $denovo_motif{$tab[0]}{"perc"}){
					print "\twe actually skip this new one because it's percent $tab[6] is not as good as ".$denovo_motif{$tab[0]}{"perc"}."\n";
					# <STDIN>;
					next;
				}
				if ($n_n==4){
					if ($denovo_motif{$tab[0]}{"n_n"}<4){
						print "\twe actually skip this new one because it's a mononucleotide while the previous one is longer\n";
						# <STDIN>;
						next;
					}
				}
				print "We take as a replacement\n";
			}
			## If we get there, that means we want this motif
			my $motif="";
			if($tab[1] eq "upstream"){}
			else{
				$tab[4]=&revcomp($tab[4]);
			}
			my @t=split("",$tab[4]);
			my $motif=$tab[4]; ## We take the 5 bases
			my $check_motif=$t[2].$t[3].$t[4]; ## to check canonical
			my $check_motif_bis=$t[1].$t[2].$t[3].$t[4]; ## to check canonical
			if (defined($canonical_motif{$motif}) || defined($canonical_motif{$check_motif}) || defined($canonical_motif{$check_motif_bis})){
				print "motif $motif is already in the list of canonical, we skip\n";
				# <STDIN>;
			}
			else{
				$denovo_motif{$tab[0]}{"motif"}=$motif;
				$denovo_motif{$tab[0]}{"n_obs"}=$tab[5];
				$denovo_motif{$tab[0]}{"perc"}=$tab[6];
				$denovo_motif{$tab[0]}{"n_n"}=$n_n;
				$denovo_motif{$tab[0]}{"from_vr"}=2;
				if (defined($info_motif{$motif})){
					print "we already know about motif $motif, so we just link it to this array\n";
				}
				else{
					print "$motif is new, we take its characteristics\n";
					my ($regexp_up,$regexp_down)=&motif_to_regexp($motif);
					$info_motif{$motif}{"regexp_up"}=$regexp_up;
					$info_motif{$motif}{"regexp_down"}=$regexp_down;
					$info_motif{$motif}{"size"}=5-$n_n;
					foreach my $motif_2 (sort keys %info_motif){
						if ($motif ne $motif_2){
							my $reg_test=$info_motif{$motif_2}{"regexp_up"};
							if ($motif_2=~/$regexp_up$/){
								$info_motif{$motif_2}{"submotifs"}{$motif}=1;
								# <STDIN>;
							}
							if ($motif=~/$reg_test$/){
								$info_motif{$motif}{"submotifs"}{$motif_2}=1;
								# <STDIN>;
							}
						}
					}
				}
			}
		}
	}
}
close $tsv;

print "Writing the output file\n";
open my $s1,">",$out_file;
print $s1 "array\tmotif\tp_obs\tregexp_up\tregexp_down\tsize\tsubmotifs\n";
foreach my $array (sort keys %denovo_motif){
    my $motif=$denovo_motif{$array}{"motif"};
    my $submotifs=join("|",sort keys %{$info_motif{$motif}{"submotifs"}});
    print $array."\t".$motif."\t".$denovo_motif{$array}{"perc"}."\t".
    $info_motif{$motif}{"regexp_up"}."\t".$info_motif{$motif}{"regexp_down"}."\t".$info_motif{$motif}{"size"}."\t".$submotifs."\n";
    print $s1 $array."\t".$motif."\t".$denovo_motif{$array}{"perc"}."\t".
    $info_motif{$motif}{"regexp_up"}."\t".$info_motif{$motif}{"regexp_down"}."\t".$info_motif{$motif}{"size"}."\t".$submotifs."\n";
}
close $s1;



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

