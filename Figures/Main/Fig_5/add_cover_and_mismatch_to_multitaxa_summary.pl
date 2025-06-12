#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use File::Basename;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
	print "# Script to look for viruses targeted by multiple phyla 
# Arguments :
# run
";
	die "\n";
}

my %levels=("0"=>"domain","1"=>"phylum","2"=>"class","3"=>"order","4"=>"family","5"=>"genus","6"=>"species");

my $virus_info="../../../Data/Spacer_db/IMGVR_sequence_information_Oct17.tsv";
my $in_file="Target_to_repeat.tsv";
my $summary_file="Summary_multitaxa.tsv";
my @profile_files=("../../../Analyses/Target_IMGVR_IMGPR/Beyond_near_exact/Virus_to_array_hits_profile-10_hits_with_cc.tsv","../../../Analyses/Target_IMGVR_IMGPR/Beyond_near_exact/Multitaxa_virus_to_array_hits_profile.with_cc.tsv");
my @cover_files=("../../../Analyses/Target_IMGVR_IMGPR/Target_coverage/uvig_hq_coverage_by_spacer-with_sp_info.tsv","../../../Analyses/Target_IMGVR_IMGPR/Target_coverage/uvigadditional_multitaxa_coverage_by_spacer-with_sp_info.tsv");

my $out_file="Summary_multiclass_wprofile.tsv";
my $out_missing_cover="List_missing_coverage.txt";

my $min_count=10; ## Minimum number of spacers to consider "not low", per virus-repeat pair

print "#######################\n";
print "Load info about viruses from $virus_info\n";
my %info_virus;
my %check_virus;
open my $tsv,"<",$virus_info;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "uvig"){next;}
    if ($tab[10] eq "phage" || $tab[10] eq "archaea"){
        $info_virus{$tab[0]}{"quality"}=$tab[5];
        $info_virus{$tab[0]}{"length"}=$tab[3];
        $info_virus{$tab[0]}{"taxo"}=$tab[6];
        if ($tab[5] eq "Reference"){$tab[5]="High-quality";}
        $check_virus{$tab[0]}{"quality"}=$tab[5];
    }

}
close $tsv;

print "#######################\n";
my %info_link;
print "Reading $summary_file\n";
open my $tsv,"<",$summary_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if (defined($check_virus{$tab[0]})){
        if ($tab[8]>1 && ($tab[7] eq "class" || $tab[7] eq "phylum" || $tab[7] eq "domain")){
            $check_virus{$tab[0]}{"multiclass"}=1;
        }
    }
}
close $tsv;


print "#######################\n";
print "Reading $in_file ...\n";
my $i=0;
open my $tsv,"<",$in_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if (defined($check_virus{$tab[0]}) && $check_virus{$tab[0]}{"multiclass"}==1){ ## We ignore uvigs not classified as phages or archaea
        $info_link{$tab[0]}{$tab[1]}{"lca_confidence"}=$tab[2];
        $info_link{$tab[0]}{$tab[1]}{"lca_genus"}=$tab[3];
        $info_link{$tab[0]}{$tab[1]}{"count"}=$tab[4];
        $i++;
        if ($i % 1000000 == 0){
            print ".. $i ..\n";
        }
    }
}
close $tsv;


print "#######################\n";
foreach my $file (@cover_files){
    print "Reading $file\n";
    open my $tsv,"<",$file;
    while(<$tsv>){
        chomp($_);
        my @tab=split("\t",$_);
        if ($tab[0] eq "uvig"){next;}
        if (!defined($info_link{$tab[0]})){next;} ## Ignore viruses we don't want to consider
        if (defined($info_link{$tab[0]}{$tab[1]})){
            $info_link{$tab[0]}{$tab[1]}{"coverage"}=$tab[4];
            $info_link{$tab[0]}{$tab[1]}{"covered"}=$tab[3];
            my $coverage_level="NA";
            ## Only consider coverage (10% + ) if the genome is high quality, i.e. if we have a complete genome of 5kb, we consider 500bp targeted as high targeting, but not if we have a fragment of 5kb
            if ($tab[3]>=1000 || ($tab[4]>=10 && $check_virus{$tab[0]}{"quality"} eq "High-quality")){$coverage_level="high";}
            elsif($tab[5]>=10){$coverage_level="medium";}
            else{$coverage_level="low";}
            $info_link{$tab[0]}{$tab[1]}{"coverage_level"}=$coverage_level;
            $info_link{$tab[0]}{$tab[1]}{"n_common"}=$tab[8];
            $info_link{$tab[0]}{$tab[1]}{"n_unique"}=$tab[11];
        }
    }
    close $tsv;
}

print "#######################\n";
foreach my $profile_file (@profile_files){
    print "Loading mismatch profiles from $profile_file\n";
    open my $tsv,"<",$profile_file;
    while(<$tsv>){
        chomp($_);
        my @tab=split("\t",$_);
        if ($tab[0] eq "uvig"){next;}
        if (defined($info_link{$tab[0]})){
            if (defined($info_link{$tab[0]}{$tab[1]})){
                $info_link{$tab[0]}{$tab[1]}{"mis_profile_cat"}=$tab[9];
                $info_link{$tab[0]}{$tab[1]}{"mis_profile"}=$tab[2].";".$tab[3].";".$tab[4].";". $tab[5];

            }
        }
    }
    close $tsv;
}


print "#######################\n";
open my $smis,">",$out_missing_cover;
my %tmp;
open my $s2,">",$out_file;
print $s2 "uvig\tn_array\tmax_rank_all\tn_taxon_max_rank\tn_array_hc\tmax_rank_hc\tn_taxon_max_rank_hc\tmax_rank_hc_10p\tn_taxon_max_rank_hc_10p\tmax_rank_hc_high\tn_taxon_max_rank_hc_high\tmax_rank_hc_high_positive\tn_taxon_max_rank_hc_high_positive\n";
foreach my $uvig (sort keys %info_link){
    %tmp=();
    print "## $uvig ##\n";
    my $n_array_total=0;
    my $n_array_hc=0;
    my $n_array_10p=0;
    my $n_array_ar=0;
    foreach my $array (sort keys %{$info_link{$uvig}}){
        $n_array_total++;
        # print "\t".$array."\t".$store{$uvig}{$array}{"lca_genus"}."\t".$store{$uvig}{$array}{"lca_confidence"}."\t".$store{$uvig}{$array}{"count"}."\n";
        $tmp{"all"}{$info_link{$uvig}{$array}{"lca_genus"}}+=$info_link{$uvig}{$array}{"count"};
        if ($info_link{$uvig}{$array}{"lca_confidence"} eq "Genome_high-confidence" || $info_link{$uvig}{$array}{"lca_confidence"} eq "Genome_medium-confidence"){
            $n_array_hc++;
            $tmp{"hc_only"}{$info_link{$uvig}{$array}{"lca_genus"}}+=$info_link{$uvig}{$array}{"count"};
            if ($info_link{$uvig}{$array}{"count"}>=$min_count){
                $tmp{"hc_and_10p"}{$info_link{$uvig}{$array}{"lca_genus"}}+=$info_link{$uvig}{$array}{"count"};
                if (!defined($info_link{$uvig}{$array}{"coverage_level"})){
                    print $smis $uvig."\t".$array."\tno coverage\n";
                }
                else{
                    if ($info_link{$uvig}{$array}{"coverage_level"} eq "high"){
                        $tmp{"hc_and_high"}{$info_link{$uvig}{$array}{"lca_genus"}}+=$info_link{$uvig}{$array}{"count"};
                        if (!defined($info_link{$uvig}{$array}{"mis_profile_cat"})){
                            print $smis $uvig."\t".$array."\tno profile\n";
                        }
                        if ($info_link{$uvig}{$array}{"mis_profile_cat"} eq "positive"){
                            $tmp{"hc_and_high_and_mis"}{$info_link{$uvig}{$array}{"lca_genus"}}+=$info_link{$uvig}{$array}{"count"};
                        }
                    }
                }
            }
        }
    }
    print "With all arrays\n";
    my ($max_rank,$n_taxon)=&get_count($tmp{"all"},$uvig,"all");
    print "\t".$max_rank."\t".$n_taxon."\n";
    print "Only hc arrays\n";
    my ($max_rank_hc,$n_taxon_hc)=&get_count($tmp{"hc_only"},$uvig,"high-confidence_only");
    print "\t".$max_rank_hc."\t".$n_taxon_hc."\n";
    print "Only hc arrays and at least 10 hits\n";
    my ($max_rank_hc_10p,$n_taxon_hc_10p)=&get_count($tmp{"hc_and_10p"},$uvig,"high-confidence_only_10p");
    print "\t".$max_rank_hc_10p."\t".$n_taxon_hc_10p."\n";
    print "Only hc arrays and at least 1kb (or 10% if hq)\n";
    my ($max_rank_hc_high,$n_taxon_hc_high)=&get_count($tmp{"hc_and_high"},$uvig,"high-confidence_only_high");
    print "\t".$max_rank_hc_high."\t".$n_taxon_hc_high."\n";
    print "Only hc arrays and at least 1kb (or 10% if hq) and only positive\n";
    my ($max_rank_hc_high_mis,$n_taxon_hc_high_mis)=&get_count($tmp{"hc_and_high_and_mis"},$uvig,"high-confidence_only_high");
    print "\t".$max_rank_hc_high_mis."\t".$n_taxon_hc_high_mis."\n";
    print $s2 $uvig."\t".$n_array_total."\t".$max_rank."\t".$n_taxon."\t".$n_array_hc."\t".$max_rank_hc."\t".$n_taxon_hc."\t".$max_rank_hc_10p."\t".$n_taxon_hc_10p."\t".$max_rank_hc_high."\t".$n_taxon_hc_high."\t".$max_rank_hc_high_mis."\t".$n_taxon_hc_high_mis."\n";
}
close $s2;
close $smis;

sub get_count(){
    my $hash=$_[0];
    my $uvig=$_[1];
    my $type=$_[2];
    my %count;
    foreach my $lca (keys %{$hash}){
        my @t=split(";",$lca);
        for (my $i=0;$i<=5;$i++){
            if ($t[$i]=~/__unclassified/){}
            else{
                my $taxo=join(";",@t[0..$i]);
                $count{$i}{$taxo}+=$$hash{$lca};
                # print $i."\t".$taxo."\t".$count{$i}{$taxo}."\n";
            }
        }
    }
    # print "=== summary ===\n";
    my $best_rank="None";
    my $best_rank_ntaxon=0;
    my $tag=0;
    my @test=keys %{$count{0}};
    if ($#test<0){} ## We end there, nothing that passed our min count cutoff
    else{
        for (my $i=0;$i<=5;$i++){
            my @list=sort {$count{$i}{$b} <=> $count{$i}{$a}} keys %{$count{$i}};
            if ($#list<0){ 
                ## We stop there
                $i=10;
            }
            else{
                my $line="";
                foreach my $taxo (@list){
                    $line.=$taxo." (".$count{$i}{$taxo}.") ";
                }
                chop($line);
                my $n_taxa=scalar(@list);
                if ($tag==0){
                    $best_rank=$levels{$i};
                    $best_rank_ntaxon=$n_taxa;
                    if ($n_taxa>1){
                        $tag=1;
                    }
                }
            }
        }
    }
    return ($best_rank,$best_rank_ntaxon);
}