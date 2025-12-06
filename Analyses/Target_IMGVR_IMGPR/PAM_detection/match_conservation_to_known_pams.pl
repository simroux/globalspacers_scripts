#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
use Custom::Utils;
my $h=0;
my $data_folder="";
GetOptions ('help' => \$h, 'h' => \$h); 
if ($h==1 || $ARGV[0] eq "" || 'd=s'=>\$data_folder){
    print "# Script to check if conserved residue upstream or downstream of spacer hits match known PAMs (fuzzy matching if needed)
# Arguments :
# -d: Data folder (i.e. /.../Data/)
";
    die "\n";
}


my $file_ref_pam=$data_folder."/Additional_data/manually_curated_PAMs_per_type.tsv";
my $in_file_arrays=$data_folder."/Spacer_db/Array_info_filtered_for_db-Oct24-25.tsv";
my $in_file_cons="Conserved_positions_neighborhoods.tsv";

my $log_extra="pam_prediction.log";
my $out_file="Stat_motif_detection.tsv";
my $out_file_bis="Detailed_stats_motifs.tsv";

my %type_to_motif;
my %info_motif;
print "Loading canonical motifs\n";
open my $tsv,"<",$file_ref_pam;
while(<$tsv>){
	chomp($_);
	my @tab=split("\t",$_);
	if ($tab[0] eq "Type"){next;}
    if ($_=~/^#/){next;}
	## We store the link between the CRISPR type and the motif
	$type_to_motif{$tab[0]}{$tab[1]}=1;
	## If we don't already have the info on this motif, we load/compute it
	if (!defined($info_motif{$tab[1]})){
		print $tab[1]."\tcanonical\n";
		my $regexp_up=&motif_to_regexp($tab[1]);
        my $regexp_down=&motif_to_regexp($tab[2]);
		$info_motif{$tab[1]}{"regexp_up"}=$regexp_up;
		$info_motif{$tab[1]}{"regexp_down"}=$regexp_down;
		$info_motif{$tab[1]}{"size"}=$tab[3];
        $info_motif{$tab[1]}{"types"}{$tab[0]}=1;
	}
    else{
        $info_motif{$tab[1]}{"types"}{$tab[0]}=1;
    }
}
close $tsv;


## Summarize type by motif
foreach my $motif (sort keys %info_motif){
    my @t=sort keys %{$info_motif{$motif}{"types"}};
    $info_motif{$motif}{"summary_types"}=join(";",@t);
}

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

my %stats;
my %store_motifs;
open my $slog,">",$log_extra;
print $slog "Repeat\tPredicted type\tN distinct neighborhoods\tPredicted PAM\tPrediction type\n";
open my $tsv,"<",$in_file_cons;
my %tmp;
my $pick="";
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "repeat_cluster"){next;}
    %tmp=();
    $pick="";
    my $c_a=$tab[0];
    my $n_d_n=$tab[2];
    if ($n_d_n < 10){
        ## Ignoring if less than 10 distinct neighborhoods, we don't trust
        next;
    }
    my $up=$tab[3];
    my $down=$tab[4];
    my $cons_level=$tab[5]+$tab[6];
    my $c_t=$ref_array{$c_a}{"type"};
    if ($cons_level>=6){
        print $c_a." neighborhoods are too conserved, we don't trust\n";
        print $slog $c_a."\t".$c_t."\t".$n_d_n."\tNA\ttoo_conserved\n";
        # $stats{$c_t}{"extra"}++; ## We ignore these cases
        next;
    }
    elsif($cons_level==0){
        print $c_a." neighborhoods have 0 conservation, we can't do anything\n";
        print $slog $c_a."\t".$c_t."\t".$n_d_n."\tNA\tno_conservation\n";
        $stats{$c_t}{"no_conserved_position"}++;
        next;
    }
    print "Looking for potential PAMs for $c_a / $c_t, looking in \t$up\t$down\n";
    foreach my $motif (sort keys %info_motif){
        if ($up=~/$info_motif{$motif}{"regexp_up"}$/){
            print "\tfound $motif in upstream: ".$info_motif{$motif}{"regexp_up"}."\n";
            $tmp{$motif}{"upstream"}=1;
        }
        if ($down=~/^$info_motif{$motif}{"regexp_down"}/){
            print "\tfound $motif in downstream: ".$info_motif{$motif}{"regexp_down"}."\n";
            $tmp{$motif}{"downstream"}=1;
        }
    }
     ## We sort based on score, if ties, we take the longer motifs, or an arbitrary order just to be reproducible
    my @raw_list=sort {$info_motif{$b}{"size"} <=> $info_motif{$a}{"size"} or $a cmp $b} keys %tmp;
    ## We ignore motifs of size 1 unless we are in $c_t II-A, II-B, or II-C
    my @list=();
    foreach my $motif (@raw_list){
        if ($c_t eq "II-A" || $c_t eq "II-B" || $c_t eq "II-C"){push(@list,$motif);}
        else{ 
            if ($info_motif{$motif}{"size"}==1){
                print "We ignore $motif because of size ".$info_motif{$motif}{"size"}." and c-t $c_t\n";
            } ## We ignore motif of size 1, too much noise, only exception is for CRISPRs of Type II where we expect them
            else{push(@list,$motif);}
        }
    }
    ### 
    if (scalar(@list)>=1){
        if ($c_t eq "Unknown"){
            print "## No predicted type for this repeat cluster, we pick the highest score, unless it's a 1-base motif which we don't trust\n";
            $pick=$list[0];
            if ($pick eq ""){
                print "!!! LIST IS EMPTY ?\n";
                die("\n");
            }
        }
        else{ ## We look first into the expected
            print "This repeat is predicted to be of type $c_t, we'll favor motifs associated with this CRISPR type (if any)\n";
            foreach my $motif (@list){
                if ($pick ne ""){last;} ## We found a pick, no need to move further down the list
                # print "\tchecking for motif $motif\n";
                foreach my $type (keys %{$info_motif{$motif}{"types"}}){
                    if ($pick ne ""){last;} ## We found a pick, no need to move further down the list
                    # print "\t\tit is linked to $type\n";
                    if ($type eq $c_t){
                        print "## We found a motif expected for type $c_t, that is great !\n";
                        $pick=$motif;
                        print "We can stop there\n";
                        last;
                    }
                }
            }
            if ($pick eq ""){ ## No match for this type, we take the "best of the rest"
                $pick=$list[0];
            }
        }
        my @t2=keys %{$tmp{$pick}};
        my $orient="";
        if ($#t2==0){$orient=$t2[0];}
        else{
            print "The same motif is found upstream and downstream ?! That's interesting, but unexpected \n";
            $orient="both";
            # die("\n");
        }
        print $c_a."\t".$c_t."\t".$pick."\t".$info_motif{$pick}{"summary_types"}."\t".$orient."\t".$up."\t".$down."\n";
        my $match="known_maybe";
        foreach my $type (keys %{$info_motif{$pick}{"types"}}){
            if ($type eq $c_t){
                $match="known_expected";
            }
            elsif($match eq "known_maybe" && $c_t ne "Unknown"){
                $match="known_other";
            }
        }
        print $slog $c_a."\t".$c_t."\t".$n_d_n."\t".$pick."\t".$match."\n";
        $store_motifs{$c_t}{$pick}++;
        $stats{$c_t}{"conserved_position"}++;
        $stats{$c_t}{"matches"}{$match}++;
    }
    else{
        print "No motif here, listing the potential de novo ones\n";
        my $tag=0;
        my $pick="";
        my $up_short=substr($up,6,4);
        my $n = () = $up_short=~/N/g;
        if ($n<=2){
            print "\tpotential motif just upstream: ".$up_short."\n";
            print $slog $c_a."\t".$c_t."\t".$n_d_n."\t".$up_short."\tnew\n";
            $tag=1;
            $pick=$up_short;
        }
        elsif($n==3 && $c_t=~/^II-/){
            print "\tpotential motif just upstream of a single base, which is ok for types II: ".$up_short."\n";
            print $slog $c_a."\t".$c_t."\t".$n_d_n."\t".$up_short."\tnew\n";
            $pick=$up_short;
            $tag=1;
        }
        my $down_short=substr($down,0,4);
        my $n = () = $down_short=~/N/g;
        if ($n<=2){
            print "\tpotential motif just downstream: ".$down_short."\n";
            my $down_short_rc=&revcomp($down_short);
            print "\twe reverse-complement: ".$down_short_rc."\n";
            print $slog $c_a."\t".$c_t."\t".$n_d_n."\t".$down_short_rc."\tnew\n";
            $pick=$down_short_rc;
            $tag=1;
        }
        elsif($n==3 && $c_t=~/^II-/){
            print "\tpotential motif just downstream of a single base, which is ok for types II: ".$down_short."\n";
            my $down_short_rc=&revcomp($down_short);
            print "\twe reverse-complement: ".$down_short_rc."\n";
            print $slog $c_a."\t".$c_t."\t".$n_d_n."\t".$down_short_rc."\tnew\n";
            $pick=$down_short_rc;
            $tag=1;
        }
        if ($tag==0){
            print $slog $c_a."\t".$c_t."\t".$n_d_n."\tNA\tno_conservation\n";
            $stats{$c_t}{"no_conserved_position"}++;
        }
        else{
            $stats{$c_t}{"matches"}{"new_motif"}++;
            $stats{$c_t}{"conserved_position"}++;
            $store_motifs{$c_t}{$pick}++;
        }
    }
}
close $tsv;
close $slog;

open my $s1,">",$out_file;
print $s1 "type\ttotal\tconserved_pos\tno_conserved_pos\texpected\tother_known\tpotential_new\tpcent_conserved\tpcent_known_motif\tpcent_all_motif\n";
open my $s2,">",$out_file_bis;
print $s2 "type\ttotal\tmotif\tcount\tpcent\n";
foreach my $type (sort keys %stats){
    if (!defined($stats{$type}{"no_conserved_position"})){$stats{$type}{"no_conserved_position"}=0;}
    if (!defined($stats{$type}{"conserved_position"})){$stats{$type}{"conserved_position"}=0;}
    my $total=$stats{$type}{"no_conserved_position"}+$stats{$type}{"conserved_position"};
    my $line=$type."\t".$total."\t".$stats{$type}{"conserved_position"}."\t".$stats{$type}{"no_conserved_position"};
    if (!defined($stats{$type}{"matches"}{"known_expected"})){$stats{$type}{"matches"}{"known_expected"}=0;}
    if (!defined($stats{$type}{"matches"}{"known_other"})){$stats{$type}{"matches"}{"known_other"}=0;}
    if (!defined($stats{$type}{"matches"}{"new_motif"})){$stats{$type}{"matches"}{"new_motif"}=0;}
    $line.="\t".$stats{$type}{"matches"}{"known_expected"}."\t".$stats{$type}{"matches"}{"known_other"}."\t".$stats{$type}{"matches"}{"new_motif"};
    if ($total==0){
        print "NOTHING FOR TYPE $type ?\n";
        $line.="\tNA\tNA\tNA";
    }
    else{
        my $p_cons=sprintf("%0.2f",$stats{$type}{"conserved_position"}/$total*100);
        $line.="\t".$p_cons;
        if ($stats{$type}{"conserved_position"} == 0){
            $line.="\tNA\tNA";
        }
        else{
            my $p_known=sprintf("%0.2f",$stats{$type}{"matches"}{"known_expected"}/$stats{$type}{"conserved_position"}*100);
            my $p_motif=sprintf("%0.2f",($stats{$type}{"matches"}{"known_expected"}+$stats{$type}{"matches"}{"known_other"})/$stats{$type}{"conserved_position"}*100);
            $line.="\t".$p_known."\t".$p_motif;
        }
    }
    print $s1 $line."\n";
    ## Also put out the list of motifs with high-enough frequency to be relevant
    my @list_motifs=(sort {$store_motifs{$type}{$b} <=> $store_motifs{$type}{$a} or $a cmp $b} keys %{$store_motifs{$type}});
    ### Below we would take the percentage amongst neighborhood with conserved position, but we decided we preferred to take percentage of all neighborhoods
    # my $total;
    # foreach my $motif (@list_motifs){$total+=$store_motifs{$type}{$motif};}
    # my $check_total=$stats{$type}{"matches"}{"known_expected"}+$stats{$type}{"matches"}{"known_other"}+$stats{$type}{"matches"}{"new_motif"};
    # print "## CHECK - $total should be equal to $check_total\n";
    # if ($total != $check_total && $type ne "Unknown"){die("pblm");}
    foreach my $motif (@list_motifs){
        my $ratio=$store_motifs{$type}{$motif}/$total*100;
        if ($ratio>=1){
            print $type."\t".$motif."\t".$store_motifs{$type}{$motif}."\t".$ratio."\n";
            print $s2 $type."\t".$total."\t".$motif."\t".$store_motifs{$type}{$motif}."\t".$ratio."\n";
        }
    }
}
close $s1;
close $s2;

sub motif_to_regexp(){
	my $motif=$_[0];
	my @t=split("",$motif);
	my @r=();
	foreach my $cell (@t){
		if ($cell eq "N"){push(@r,"[ATCGN]");}
		else{push(@r,$cell);}
	}
	my $motif_regex=join("",@r);
	return($motif_regex);
}
