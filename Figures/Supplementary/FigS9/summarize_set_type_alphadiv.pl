#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
my %translate_level=("d"=>"domain","p"=>"phylum","c"=>"class","o"=>"order","f"=>"family","g"=>"genus","s"=>"species");
GetOptions ('help' => \$h, 'h' => \$h); 
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
	print "# Script to make some counts of each type of set (based on alphadiv)
# Arguments :
# run
";
	die "\n";
}

my $file_array_info="Data/Spacer_db/Array_info_filtered_for_db-Nov1-24.tsv";
my $file_sample_info="Data/Spacer_db/Runs_to_ecosystem_and_sequencing_and_study_for_db-Jul28-24.tsv";
my $in_alpha_base="Analyses/Spacer_database/spacer_sets_alphadiv.tsv";

my $out_file_frozen="spacer_sets_frozen.tsv";
my $out_file_stats_frozen="spacer_sets_frozen_summary_for_percentage.tsv";

my %info_array;
open my $tsv,"<",$file_array_info;
print "Reading $file_array_info\n";
while(<$tsv>){
      chomp($_);
      my @tab=split("\t",$_);
      if ($tab[0] eq "repeat_cluster"){next;}
      my @t=split("-",$tab[1]);
      $info_array{$tab[0]}{"type"}=$t[0];
      $info_array{$tab[0]}{"taxo_source"}=$tab[2];
      $info_array{$tab[0]}{"taxo_genus"}=$tab[6];
      my @t=split((";",$tab[6]));
      $info_array{$tab[0]}{"taxo_phylum"}=$t[0].";".$t[1];
      if ($tab[3] eq "NA"){$info_array{$tab[0]}{"taxo_lvl"}="none";}
      else{$info_array{$tab[0]}{"taxo_lvl"}=&guess_level($tab[3]);}
}
close $tsv;

my %info_sample;
open my $tsv,"<",$file_sample_info;
print "Reading $file_sample_info\n";
while(<$tsv>){
      chomp($_);
      my @tab=split("\t",$_);
      if ($tab[0] eq "library"){next;}
      $info_sample{$tab[0]}{"eco"}=$tab[3];
      $info_sample{$tab[1]}{"eco"}=$tab[3];
}
close $tsv;

my $i=0;
my %stats;
my %count_by_array;
my %count_by_sample;
my %array_to_eco;
open my $tsv,"<",$in_alpha_base;
print "Reading $in_alpha_base\n";
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "array"){
        # print join("\t",@tab[0..6])."\t".$tab[12]."\t".$tab[13]."\tP_min50\tP_max5\ttype\n";
        next;
    }
    my $p_min50=$tab[12]/$tab[2];
    my $p_max5=$tab[13]/$tab[2];
    my $n_spacer=$tab[2];
    my $max_cover=$tab[4];
    my $type="NA";
    if ($max_cover < 20){
        $type="low_coverage";
    }
    elsif($p_min50<=0.25){
        $type="standard";
    }
    elsif($p_min50>=0.50){
        $type="frozen";
    }
    else{
        $type="other";
    }
    # print join("\t",@tab[0..6])."\t".$tab[12]."\t".$tab[13]."\t".$p_min50."\t".$p_max5."\t".$type."\n";
    $stats{$type}++;
    $count_by_array{$tab[0]}{$type}++;
    $count_by_sample{$tab[1]}{$type}++;
    if (defined($info_sample{$tab[1]}{"eco"}) && $info_sample{$tab[1]}{"eco"} ne "Unknown"){
        $array_to_eco{$tab[0]}{$info_sample{$tab[1]}{"eco"}}++;
        $array_to_eco{$tab[0]}{"total"}++;
    }
    # <STDIN>;
    $i++;
    if ($i % 100000 == 0){print $i." ... \n";}
    # if ($i>=500000){last;}
}
close $tsv;
print "## General stats\n";
foreach my $type (sort keys %stats){
    print $type."\t".$stats{$type}."\n";
}
print "#####################\n";
my %stats_summary;
my %stats_taxo;
open my $s1,">",$out_file_frozen;
print $s1 "array\tlow_diversity\tlow_coverage\tother\tstandard\ttaxo\tarray_type\tmaj_ecosystem\n";
foreach my $array (sort keys %count_by_array){
    if (!defined($count_by_array{$array}{"frozen"})){$count_by_array{$array}{"frozen"}=0;}
    if (!defined($count_by_array{$array}{"low_coverage"})){$count_by_array{$array}{"low_coverage"}=0;}
    if (!defined($count_by_array{$array}{"standard"})){$count_by_array{$array}{"standard"}=0;}
    if (!defined($count_by_array{$array}{"other"})){$count_by_array{$array}{"other"}=0;}
    $stats_summary{"total"}++;
    my $total_qualifying=$count_by_array{$array}{"frozen"}+$count_by_array{$array}{"standard"};
    my $type_frozen="";
    if ($total_qualifying<10){
        $type_frozen="low_samples";
    }
    elsif($total_qualifying>=10){
        $type_frozen="high_samples";
    }
    if ($count_by_array{$array}{"frozen"}>0){
        my $maj_ecosystem=&get_maj_eco($array_to_eco{$array});
        print $s1 $array."\t".$count_by_array{$array}{"frozen"}."\t".$count_by_array{$array}{"low_coverage"}."\t".
        $count_by_array{$array}{"other"}."\t".$count_by_array{$array}{"standard"}."\t".$info_array{$array}{"taxo_genus"}."\t".
        $info_array{$array}{"type"}."\t".$maj_ecosystem."\n";
        my $ratio=$count_by_array{$array}{"frozen"}/$total_qualifying;
        if ($ratio>=0.25){
            $type_frozen.="_common_frozen";
        }
        else{
            $type_frozen.="_rare_frozen";
        }
    }
    else{
        $type_frozen.="_no_frozen";
    }
    $stats_summary{$type_frozen}++;
    ## Now for some stats:
    # ignoring low coverage
    $stats_taxo{"global"}{"standard"}+=$count_by_array{$array}{"standard"};
    $stats_taxo{"global"}{"frozen"}+=$count_by_array{$array}{"frozen"};
    if (defined($info_array{$array}{"taxo_genus"}) && ($info_array{$array}{"taxo_genus"} ne "NA") && ($info_array{$array}{"taxo_genus"} ne "")){
        $stats_taxo{"genus"}{$info_array{$array}{"taxo_genus"}}{"standard"}+=$count_by_array{$array}{"standard"};
        $stats_taxo{"genus"}{$info_array{$array}{"taxo_genus"}}{"frozen"}+=$count_by_array{$array}{"frozen"};
    }
    if (defined($info_array{$array}{"taxo_phylum"}) && ($info_array{$array}{"taxo_phylum"} ne "NA") && ($info_array{$array}{"taxo_phylum"} ne "")){
        $stats_taxo{"phylum"}{$info_array{$array}{"taxo_phylum"}}{"standard"}+=$count_by_array{$array}{"standard"};
        $stats_taxo{"phylum"}{$info_array{$array}{"taxo_phylum"}}{"frozen"}+=$count_by_array{$array}{"frozen"};
    }
}
close $s1;
open my $s2,">",$out_file_stats_frozen;
print $s2 "type\tcount\n";
foreach my $type (sort keys %stats_summary){
    print $s2 $type."\t".$stats_summary{$type}."\n";
}
close $s2;

$stats_taxo{"global"}{"ratio"}=$stats_taxo{"global"}{"frozen"}/($stats_taxo{"global"}{"frozen"}+$stats_taxo{"global"}{"standard"});
foreach my $genus (sort keys %{$stats_taxo{"genus"}}){
    if ($stats_taxo{"genus"}{$genus}{"frozen"}>0 || $stats_taxo{"genus"}{$genus}{"standard"}>0){
        $stats_taxo{"genus"}{$genus}{"ratio"}=$stats_taxo{"genus"}{$genus}{"frozen"}/($stats_taxo{"genus"}{$genus}{"frozen"}+$stats_taxo{"genus"}{$genus}{"standard"});
    }
}
foreach my $phylum (sort keys %{$stats_taxo{"phylum"}}){
    if ($stats_taxo{"phylum"}{$phylum}{"frozen"}>0 || $stats_taxo{"phylum"}{$phylum}{"standard"}>0){
        $stats_taxo{"phylum"}{$phylum}{"ratio"}=$stats_taxo{"phylum"}{$phylum}{"frozen"}/($stats_taxo{"phylum"}{$phylum}{"frozen"}+$stats_taxo{"phylum"}{$phylum}{"standard"});
    }
}
print "##############################################################################################################################\n";
foreach my $genus (sort {$stats_taxo{"genus"}{$b}{"ratio"} <=> $stats_taxo{"genus"}{$a}{"ratio"} or $stats_taxo{"genus"}{$b}{"frozen"} <=> $stats_taxo{"genus"}{$a}{"frozen"}} keys %{$stats_taxo{"genus"}}){
    print $genus."\t".$stats_taxo{"genus"}{$genus}{"frozen"}."\t".$stats_taxo{"genus"}{$genus}{"ratio"}."\t vs ".$stats_taxo{"global"}{"ratio"}."\n";
}
print "##############################################################################################################################\n";
foreach my $phylum (sort {$stats_taxo{"phylum"}{$b}{"ratio"} <=> $stats_taxo{"phylum"}{$a}{"ratio"} or $stats_taxo{"phylum"}{$b}{"frozen"} <=> $stats_taxo{"phylum"}{$a}{"frozen"}} keys %{$stats_taxo{"phylum"}}){
    print $phylum."\t".$stats_taxo{"phylum"}{$phylum}{"frozen"}."\t".$stats_taxo{"phylum"}{$phylum}{"ratio"}."\t vs ".$stats_taxo{"global"}{"ratio"}."\n";
}

### Now similar stats by ecosystem
my %stats_eco;
foreach my $sample (sort keys %count_by_sample){
    ## Now for some stats:
    # ignoring low coverage
    $stats_eco{"global"}{"standard"}+=$count_by_sample{$sample}{"standard"};
    $stats_eco{"global"}{"frozen"}+=$count_by_sample{$sample}{"frozen"};
    if (defined($info_sample{$sample}{"eco"}) && ($info_sample{$sample}{"eco"} ne "NA") && ($info_sample{$sample}{"eco"} ne "")){
        $stats_eco{$info_sample{$sample}{"eco"}}{"standard"}+=$count_by_sample{$sample}{"standard"};
        $stats_eco{$info_sample{$sample}{"eco"}}{"frozen"}+=$count_by_sample{$sample}{"frozen"};
    }
}

foreach my $eco (sort keys %stats_eco){
    $stats_eco{$eco}{"ratio"}=$stats_eco{$eco}{"frozen"}/($stats_eco{$eco}{"frozen"}+$stats_eco{$eco}{"standard"});
}
print "##############################################################################################################################\n";
foreach my $eco (sort {$stats_eco{$b}{"ratio"} <=> $stats_eco{$a}{"ratio"} or $stats_eco{$b}{"frozen"} <=> $stats_eco{$a}{"frozen"}} keys %stats_eco){
    print $eco."\t".$stats_eco{$eco}{"frozen"}."\t".$stats_eco{$eco}{"ratio"}."\t vs ".$stats_eco{"global"}{"ratio"}."\n";
}




sub guess_level(){
    my $taxo=$_[0];
    my @t=split(";",$taxo);
    my $max_lvl="NA";
    foreach my $cell (@t){
        my @t2=split("__",$cell);
        if ($t2[1] ne "unclassified"){
            if (!defined($translate_level{$t2[0]})){die("Pblm with cell $t2[0] // $t2[1] // $cell // $taxo\n");}
            $max_lvl=$translate_level{$t2[0]};
        }
    }
    return $max_lvl;
}

sub get_maj_eco(){
    my $hash=$_[0];
    my @t=sort {$$hash{$b} <=> $$hash{$a} or $b <=> $a} keys %{$hash};
    if ($t[0] eq "total"){shift(@t);}
    my $maj="unknown";
    if ($$hash{$t[0]}>=(0.5*$$hash{"total"})){
        print "Maj eco is $t[0] -- ".$$hash{$t[0]}." vs ".$$hash{"total"}."\n";
        $maj=$t[0];
    }
    else{
        print "Maj eco is mixed because $t[0] -- ".$$hash{$t[0]}." vs ".$$hash{"total"}."\n";
        $maj="mixed";
    }
    return $maj;
}