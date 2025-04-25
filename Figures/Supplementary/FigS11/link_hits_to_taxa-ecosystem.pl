#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h);
if ($h==1 || $ARGV[0] eq ""){
    print "# Script to cross-reference atypical phages and their targeting
# Arguments :
# run
";
    die "\n";
}


my $info_virus_file="selected_viruses.tsv";
my $spacer_info_file="selected_spacer_info.tsv";
my $spacer_hit_file="selected_spacer_hits.tsv";

my $int_file="Hit_taxon_and_eco_by_uvig.tsv";
my $out_file="Input_for_atypical_taxa_ecosystem.tsv";

my %info_virus;
print "Reading $info_virus_file..\n";
open my $tsv,"<",$info_virus_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "uvig"){next;}
    $info_virus{$tab[0]}{"taxo"}=$tab[6];
}
close $tsv;

my %info_spacer;
print "Reading $spacer_info_file..\n";
open my $tsv,"<",$spacer_info_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "spacer_id"){next;}
    $info_spacer{$tab[0]}{"lca_origin"}=$tab[11];
    $info_spacer{$tab[0]}{"lca"}=$tab[12];
    $info_spacer{$tab[0]}{"eco"}=$tab[19];
}
close $tsv;

my $i=0;
my %stats;
my %seen;
print "Reading $spacer_hit_file...\n";
open my $tsv,"<",$spacer_hit_file;
while(<$tsv>){
    chomp($_);
    $i++;
    my @tab=split("\t",$_);
    if ($tab[0] eq "uvig"){next;}
    my $uvig=$tab[0];
    my $cl=$tab[1];
    my $sp=$tab[7];
    if (!defined($seen{$uvig}{$cl})){
        $stats{$uvig}{"total_hits"}++;
        $stats{$uvig}{"lca"}{$info_spacer{$sp}{"lca"}}++;
        $stats{$uvig}{"eco"}{$info_spacer{$sp}{"eco"}}++;
        $seen{$uvig}{$cl}{$info_spacer{$sp}{"lca"}}=1;
        $seen{$uvig}{$cl}{$info_spacer{$sp}{"eco"}}=1;
    }
    elsif(!defined($seen{$uvig}{$cl}{$info_spacer{$sp}{"lca"}})){
        $stats{$uvig}{"lca"}{$info_spacer{$sp}{"lca"}}++;
        $seen{$uvig}{$cl}{$info_spacer{$sp}{"lca"}}=1;
    }
    elsif(!defined($seen{$uvig}{$cl}{$info_spacer{$sp}{"eco"}})){
        $stats{$uvig}{"eco"}{$info_spacer{$sp}{"eco"}}++;
        $seen{$uvig}{$cl}{$info_spacer{$sp}{"eco"}}=1;
    }
    if ($i % 100000 == 0){print ".. $i ..\n";}
}
close $tsv;

open my $s1,">",$int_file;
print $s1 "uvig\ttaxo\ttype\tcat\tcount\n";
my %store;
my %check;
foreach my $uvig (sort keys %stats){
    print $s1 $uvig."\t".$info_virus{$uvig}{"taxo"}."\ttotal\ttotal\t".$stats{$uvig}{"total_hits"}."\n";
    print $uvig."\t".$info_virus{$uvig}{"taxo"}."\ttotal\ttotal\t".$stats{$uvig}{"total_hits"}."\n";
    ## Preparing stats for next step as well
    my @t=split(";",$info_virus{$uvig}{"taxo"});
    my $custom_tax="NA";
    if ($info_virus{$uvig}{"taxo"}=~/^;;;/){
        $custom_tax=";;;;;".$t[5];
    }
    else{
        $custom_tax=$t[0].";".$t[1].";".$t[2].";".$t[3];
    }
    $check{$custom_tax}{$uvig}=1;
    foreach my $lca (keys %{$stats{$uvig}{"lca"}}){
        print $s1 $uvig."\t".$info_virus{$uvig}{"taxo"}."\tlca\t".$lca."\t".$stats{$uvig}{"lca"}{$lca}."\n";
        my @t2=split(";",$lca);
        if ($lca eq "NA" || $lca eq ""){}
        elsif($t2[1] eq "" || $t2[1] eq "p__unclassified"){}
        else{
            my $sum_tax=$t2[0].";".$t2[1];
            $store{$custom_tax}{$uvig}{"lca"}{$sum_tax}+=$stats{$uvig}{"lca"}{$lca};
        }
    }
    foreach my $eco (keys %{$stats{$uvig}{"eco"}}){
        print $s1 $uvig."\t".$info_virus{$uvig}{"taxo"}."\teco\t".$eco."\t".$stats{$uvig}{"eco"}{$eco}."\n";
        if ($eco eq ""){}
        elsif ($eco eq "Unknown"){}
        else{
            $store{$custom_tax}{$uvig}{"eco"}{$eco}+=$stats{$uvig}{"eco"}{$eco};
        }
    }
}
close $s1;

open my $s2,">",$out_file;
print $s2 "tax\ttype\tcat\tcount\n";
my %tmp;
foreach my $tax (sort keys %check){
    my $n_uvig=scalar(keys %{$check{$tax}});
    my $cutoff_taxa=0.05;
    if ($n_uvig>300){$cutoff_taxa=0.01;}
    %tmp=();
    foreach my $uvig (keys %{$check{$tax}}){
        my $tag_lca=0;
        foreach my $lca (keys %{$store{$tax}{$uvig}{"lca"}}){
            if ($store{$tax}{$uvig}{"lca"}{$lca}>=2){
                $tmp{"lca"}{$lca}++;
                $tag_lca=1;
            }
        }
        if ($tag_lca==0){
            $tmp{"none_lca"}++;
        }
        my $tag_eco=0;
        foreach my $eco (keys %{$store{$tax}{$uvig}{"eco"}}){
            if ($store{$tax}{$uvig}{"eco"}{$eco}>=2){
                $tmp{"eco"}{$eco}++;
                $tag_eco=1;
            }
        }
        if ($tag_eco==0){
            $tmp{"none_eco"}++;
        }
    }
    print $tax."\t".$n_uvig."\n";
    my $i=0; ## We show the top 3 systematically
    foreach my $lca (sort {$tmp{"lca"}{$b} <=> $tmp{"lca"}{$a} or $a cmp $b} keys %{$tmp{"lca"}}){
        print "\t".$lca."\t".$tmp{"lca"}{$lca}."\n";
        $tmp{"total_lca"}+=$tmp{"lca"}{$lca};
        my $ratio=$tmp{"lca"}{$lca}/$n_uvig;
        $i++;
        if ($ratio>=$cutoff_taxa || $i<=3){
            print $s2 $tax."\ttax\t".$lca."\t".$tmp{"lca"}{$lca}."\n";
        }
        else{
            $tmp{"other_lca"}+=$tmp{"lca"}{$lca};
        }
    }
    if($tmp{"other_lca"}>0){print $s2 $tax."\ttax\tother\t".$tmp{"other_lca"}."\n";}
    print $s2 $tax."\ttax\tnone\t".$tmp{"none_lca"}."\n";
    foreach my $eco (sort {$tmp{"eco"}{$b} <=> $tmp{"eco"}{$a} or $a cmp $b} keys %{$tmp{"eco"}}){
        print "\t\t".$eco."\t".$tmp{"eco"}{$eco}."\n";
        $tmp{"total_eco"}+=$tmp{"eco"}{$eco};
        my $ratio=$tmp{"eco"}{$eco}/$n_uvig;
        if ($ratio>=0){
            print $s2 $tax."\teco\t".$eco."\t".$tmp{"eco"}{$eco}."\n";
        }
        else{
            $tmp{"other_eco"}+=$tmp{"eco"}{$eco};
        }
    }
    if ($tmp{"other_eco"}>0){print $s2 $tax."\teco\tother\t".$tmp{"other_eco"}."\n";}
    print $s2 $tax."\teco\tnone\t".$tmp{"none_eco"}."\n";
}
close $s2;


