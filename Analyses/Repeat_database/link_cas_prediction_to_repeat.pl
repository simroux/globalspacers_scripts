#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h); # , 'i=s'=>\$in_file, 'r=s'=>\$out_dir_root);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to parse the Cas prediction from CCtyper to get complementary Type predictions
# Arguments :
# none
";
    die "\n";
}

my $out_file="cas_by_repeat_cluster.tsv";

## First, linking individual repeats to their closest cas operon, if not already done
my @list_to_process;
foreach my $cct_folder (<all_cct/*/>){
    my $operon_file=$cct_folder."/cas_operons.tab";
    my $repeat_file=$cct_folder."/crisprs_all.tab";
    my $out_file=$cct_folder."/crisprs_all_to_closest_cas_operon.tsv";
    push(@list_to_process,$out_file);
    if (-e $out_file){
        print "$out_file already here\n";
        next;
    }
    if (!(-e $operon_file)){
        print "No operon detected in contigs, we stop\n";
        open my $tsv,"<",$repeat_file;
        open my $s1,">",$out_file;
        while(<$tsv>){
            chomp($_);
            my @tab=split("\t",$_);
            if ($tab[0] eq "Contig"){
                print $s1 $tab[0]."\t".$tab[1]."\t".$tab[2]."\t".$tab[3]."\tdistance\tCas\n";
            }
            else{
                print $s1 $tab[0]."\t".$tab[1]."\t".$tab[2]."\t".$tab[3]."\tNA\tNA\n";
            }
        }
        close $s1;
        close $tsv;
    }
    else{
        my %store;
        my $header="";
        open my $tsv,"<",$operon_file;
        while(<$tsv>){
            chomp($_);
            my @tab=split("\t",$_);
            if ($tab[0] eq "Contig"){$header=$_; next;}
            $store{$tab[0]}{$tab[1]}{"start"}=$tab[2];
            $store{$tab[0]}{$tab[1]}{"end"}=$tab[3];
            $store{$tab[0]}{$tab[1]}{"line"}=$_;
        }
        close $tsv;

        my %tmp;
        open my $tsv,"<",$repeat_file;
        open my $s1,">",$out_file;
        while(<$tsv>){
            chomp($_);
            my @tab=split("\t",$_);
            if ($tab[0] eq "Contig"){
                print $s1 $tab[0]."\t".$tab[1]."\t".$tab[2]."\t".$tab[3]."\tdistance\t".$header."\n";
            }
            else{
                %tmp=();
                foreach my $operon (sort keys %{$store{$tab[0]}}){
                    my $n1=abs($store{$tab[0]}{$operon}{"start"}-$tab[2]);
                    my $n2=abs($store{$tab[0]}{$operon}{"start"}-$tab[3]);
                    my $n3=abs($store{$tab[0]}{$operon}{"end"}-$tab[2]);
                    my $n4=abs($store{$tab[0]}{$operon}{"end"}-$tab[3]);
                    my @t=($n1,$n2,$n3,$n4);
                    print join("\t",@t)."\n";
                    @t=sort {$a <=> $b} (@t);
                    print join("\t",@t)."\n";
                    if ($t[0] <= 50000){
                        print $tab[1]."\t".$operon."\t".$tab[2].";".$tab[3]."\t".$store{$tab[0]}{$operon}{"start"}.";".$store{$tab[0]}{$operon}{"end"}."\t".$t[0]."\n";
                        $tmp{$operon}=$t[0];
                    }
                }
                my @list=sort {$tmp{$a} <=> $tmp{$b} or $a cmp $b} keys %tmp;
                if (scalar(@list)>0){
                    print "## Best one: ".$list[0]."\n";
                    print $s1 $tab[0]."\t".$tab[1]."\t".$tab[2]."\t".$tab[3]."\t".$tmp{$list[0]}."\t".$store{$tab[0]}{$list[0]}{"line"}."\n";
                    # <STDIN>;
                }
                else{
                    print $s1 $tab[0]."\t".$tab[1]."\t".$tab[2]."\t".$tab[3]."\tNA\tNA\n";
                }
            }
        }
        close $s1;
        close $tsv;
    }
    print "output written in $out_file\n";
}

## Now load all the repeat to cluster links
my $in_file="all_repeats_clstr.tab";
my %array_to_cluster;
print "Reading $in_file\n";
open my $tsv,"<",$in_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "Array"){next;}
    my @t=split(";",$tab[10]);
    foreach my $array (@t){
        $array_to_cluster{$array}=$tab[1];
    }
}
close $tsv;


## Now load all these Cas information

## 6 -> cas_locus    $store{$code}{"cas"}.
## 7 distance_bp    $store{$code}{"dist"}
## 8 interference   $store{$code}{"interference"}
## 9 adaptation $store{$code}{"adaptation"}
## 10 cas_gene_list list
## 11 array_size_n_spacer   $store{$code}{"size"}

## Contig	CRISPR	Start	End	distance	Contig	Operon	Start	End	Prediction	Complete_Interference	Complete_Adaptation	Best_type	Best_score	Genes	Positions	E-values	CoverageSeq	CoverageHMM	Strand_Interference	Strand_Adaptation
# 

## 
# $store{$code}{"dist"}=$tab[4];
#             $store{$code}{"interference"}=$tab[10];
#             $store{$code}{"adaptation"}=$tab[11];
#             my @t=split(",",$tab[14]);
#             foreach my $cas (@t){
#                 $cas=~s/^\[//;
#                 $cas=~s/\]$//;
#                 $cas=~s/^\s+//g;
#                 $cas=~s/^'//;
#                 $cas=~s/'$//;
#                 $store{$code}{"cas_list"}{$cas}++;
#             }

my %tmp;
my %store;
foreach my $file (@list_to_process){
    print ".. $file ..\n";
    ## 
    open my $tsv,"<",$file;
    while(<$tsv>){
        chomp($_);
        my @tab=split("\t",$_);
        if ($tab[0] eq "Contig"){next;}
        if (defined($array_to_cluster{$tab[1]})){
            my $code=$array_to_cluster{$tab[1]};
            # print $code."\n";
            # <STDIN>;
            $store{$code}{"total"}++;
            if ($tab[5] ne "NA"){
                $store{$code}{"typed"}++;
                $store{$code}{"type"}{$tab[9]}++;
                if ($tab[9]=~/^Hybrid/){
                    if ($tab[10]=~/\d',/){
                        $tab[10]=~s/',/\%',/g;
                        $tab[10]=~s/\%\%/\%/g;
                    }
                    $store{$code}{"interference_raw"}{$tab[8]}++;
                    $tab[10]=~s/^\[//;
                    $tab[10]=~s/\]$//;
                    $tab[10]=~s/\s//g;
                    %tmp=();
                    my @t=split(",",$tab[10]);
                    foreach my $check (@t){
                        $check=~s/\'//g;
                        my $test=&translate($check);
                        if ($test ne "NA"){$tmp{$test}++;}
                    }
                    my @t2=keys %tmp;
                    my $n=scalar(@t2);
                    if ($n==0){$store{$code}{"interference"}{"unknown"}++;}
                    elsif($n==1){$store{$code}{"interference"}{$t2[0]}++;}
                    else{$store{$code}{"interference"}{"mixed"}++;}
                    ## 
                    if ($tab[11]=~/\d',/){
                        $tab[11]=~s/',/\%',/g;
                        $tab[11]=~s/\%\%/\%/g;
                    }
                    $store{$code}{"adaptation_raw"}{$tab[9]}++;
                    $tab[11]=~s/^\[//;
                    $tab[11]=~s/\]$//;
                    $tab[11]=~s/\s//g;
                    %tmp=();
                    my @t=split(",",$tab[11]);
                    foreach my $check (@t){
                        $check=~s/\'//g;
                        my $test=&translate($check);
                        if ($test ne "NA"){$tmp{$test}++;}
                    }
                    my @t2=keys %tmp;
                    my $n=scalar(@t2);
                    if ($n==0){$store{$code}{"adaptation"}{"unknown"}++;}
                    elsif($n==1){$store{$code}{"adaptation"}{$t2[0]}++;}
                    else{$store{$code}{"adaptation"}{"mixed"}++;}
                }
                else{
                    if (!($tab[10]=~/\%$/) && $tab[10] ne "NA"){$tab[10].="%";}
                    my $interference=&translate($tab[10]);
                    $store{$code}{"interference"}{$interference}++;
                    $store{$code}{"interference_raw"}{$tab[10]}++;
                    if (!($tab[11]=~/\%$/) && $tab[11] ne "NA"){$tab[11].="%";}
                    my $adaptation=&translate($tab[11]);
                    $store{$code}{"adaptation"}{$adaptation}++;
                    $store{$code}{"adaptation_raw"}{$tab[11]}++;
                }
                my @t=split(",",$tab[14]);
                foreach my $cas (@t){
                    $cas=~s/^\[//;
                    $cas=~s/\]$//;
                    $cas=~s/^\s+//g;
                    $cas=~s/^'//;
                    $cas=~s/'$//;
                    $store{$code}{"cas_list"}{$cas}++;
                }
            }
        }
    }
    close $tsv;
}

open my $s1,">",$out_file;
print $s1 "repeat_cluster\tn_cas_obs\tc_type\tc_interference\tc_adaptation\tc_cas_genes\tall_type\tall_interference\tall_adaptation\tall_cas_genes\n";
foreach my $cluster (sort keys %store){
    print "## $cluster \n";
    if (!defined($store{$cluster}{"typed"})){
        print $store{$cluster}{"total"}."\t0\tNO CAS GENE\n";
        print $s1 $cluster."\tNo_cas_gene_predicted\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
    }
    else{
        print $store{$cluster}{"total"}."\t".$store{$cluster}{"typed"}."\n";
        my $th=$store{$cluster}{"typed"}*0.5;
        ## type
        print "\ttype\n";
        my $cons_type="unknown";
        my $list_type="";
        foreach my $type (sort {$store{$cluster}{"type"}{$b} <=> $store{$cluster}{"type"}{$a} or $a cmp $b} keys %{$store{$cluster}{"type"}}){
            print "\t\t".$type."\t".$store{$cluster}{"type"}{$type}."\n";
            if ($store{$cluster}{"type"}{$type}>=$th){$cons_type=$type;}
            $list_type.=$type." (".$store{$cluster}{"type"}{$type}.") ";
        }
        chop($list_type);
        print "\tconsensus: ".$cons_type."\n";
        print "\tlist: ".$list_type."\n";
        ## interference
        print "\tinterference\n";
        my $cons_interference="unknown";
        my $list_interference="";
        foreach my $interf (sort {$store{$cluster}{"interference_raw"}{$b} <=> $store{$cluster}{"interference_raw"}{$a} or $a cmp $b} keys %{$store{$cluster}{"interference_raw"}}){
            $list_interference.=$interf." (".$store{$cluster}{"interference_raw"}{$interf}.") ";
        }
        foreach my $interf (sort {$store{$cluster}{"interference"}{$b} <=> $store{$cluster}{"interference"}{$a} or $a cmp $b} keys %{$store{$cluster}{"interference"}}){
            if ($store{$cluster}{"interference"}{$interf}>=$th){$cons_interference=$interf;}
        }
        print "\tconsensus: ".$cons_interference."\n";
        print "\tlist: ".$list_interference."\n";
        ## adaptation
        print "\tadaptation\n";
        my $cons_adaptation="unknown";
        my $list_adaptation="";
        foreach my $adapt (sort {$store{$cluster}{"adaptation_raw"}{$b} <=> $store{$cluster}{"adaptation_raw"}{$a} or $a cmp $b} keys %{$store{$cluster}{"adaptation_raw"}}){
            $list_adaptation.=$adapt." (".$store{$cluster}{"adaptation_raw"}{$adapt}.") ";
        }
        foreach my $adapt (sort {$store{$cluster}{"adaptation"}{$b} <=> $store{$cluster}{"adaptation"}{$a} or $a cmp $b} keys %{$store{$cluster}{"adaptation"}}){
            if ($store{$cluster}{"adaptation"}{$adapt}>=$th){$cons_adaptation=$adapt;}
        }
        print "\tconsensus: ".$cons_adaptation."\n";
        print "\tlist: ".$list_adaptation."\n";
        ## cas list
        my $cas_list="";
        my $full_list="";
        foreach my $test (sort {$store{$cluster}{"cas_list"}{$b} <=> $store{$cluster}{"cas_list"}{$a} or $a cmp $b} keys %{$store{$cluster}{"cas_list"}}){
            $full_list.=$test." (".$store{$cluster}{"cas_list"}{$test}.") ";
            if ($store{$cluster}{"cas_list"}{$test}>=$th){
                $cas_list.=$test.";";
            }
        }
        chop($cas_list);
        chop($full_list);
        print $cas_list."\n";
        print $full_list."\n";
        # if ($store{$cluster}{"typed"}>=5){<STDIN>;}  
        print $s1 $cluster."\t".$store{$cluster}{"typed"}."\t".$cons_type."\t".$cons_interference."\t".$cons_adaptation."\t".$cas_list."\t".$list_type."\t".$list_interference."\t".$list_adaptation."\t".$full_list."\n";
    }
}
close $s1;

sub translate(){
    my $result;
    my $input=$_[0];
    if ($input=~/^\d+$/){
        $input.="%";
    }
    if ($_[0] eq "100%"){$result="complete";}
    elsif ($_[0] eq "0%"){$result="absent";}
    elsif($_[0] eq "NA"){$result="NA";}
    else{$result="partial";}
    return $result;
}