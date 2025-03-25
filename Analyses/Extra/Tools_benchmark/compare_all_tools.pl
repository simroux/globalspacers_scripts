#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h); # , 'i=s'=>\$in_file, 'r=s'=>\$out_dir_root);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to compare the results obtained on all simulations
# Arguments :
# none
";
    die "\n";
}

my $ref_file="Spacer_table_sub200.tsv";
my $out_file="Compare_spacer_recovery_sub200.tsv";
my $dir_crass="Results_Crass_sub200";
my $dir_se="Results_SE_sub200";
my $dir_mc="Results_MetaCrast_sub200";

my $array_info="./SE_Db_r214GBIMG_0.9/Arrays_info.tsv";
my %array_to_seq;
my %check_repeat;
open my $tsv,"<",$array_info;
while(<$tsv>){
     chomp($_);
     my @tab=split("\t",$_);
     if ($tab[0] eq "Cluster"){next;}
     $array_to_seq{$tab[0]}=$tab[2];
     $check_repeat{$tab[2]}=1;
     my $revcomp=&revcomp($tab[2]);
     $check_repeat{$revcomp}=1;
}
close $tsv;



my %expected;
my %info;
my $n=0;
open my $tsv,"<",$ref_file;
while(<$tsv>){
     chomp($_);
     $n++;
     my @tab=split("\t",$_);
     my $revcomp_repeat=&revcomp($tab[2]);
     if (defined($check_repeat{$tab[2]})){
          my $revcomp_spacer=&revcomp($tab[3]);
          ## Do we already know this repeat for this genome ?
          if (defined($info{$tab[0]})){
               if (defined($info{$tab[0]}{$tab[2]})){} ## We do, and it's this one, all good
               if(defined($info{$tab[0]}{$revcomp_repeat})){
                    ## We know the revcomp of this repeat, we should not know tab[2], and change it moving forward
                    if (defined($info{$tab[0]}{$tab[2]})){
                         print "!!! PBLM - WE SHOULD NOT STORE BOTH $tab[2] and $revcomp_spacer for $tab[0] !!\n";
                         die("\n");
                    }
                    $tab[2]=$revcomp_repeat;
               }
          }
          if (!defined($info{$tab[0]}{$tab[2]}{$tab[3]}) && !defined($info{$tab[0]}{$tab[2]}{$revcomp_spacer})){
               ## Only consider additional spacer if not identical (or RC) to a previous one
               $expected{$tab[0]}{$tab[2]}++;
               $info{$tab[0]}{$tab[2]}{$tab[3]}=1;
               $info{$tab[0]}{$tab[2]}{$revcomp_spacer}=1;
          }
     }
     else{
          print "Ignoring repeat $tab[2] -- not reliable apparently since not in our database\n";
     }
}
close $tsv;

my %check_code;
my %store;
my %ignored;
my %seen;
foreach my $log_file (<$dir_crass/*.log.txt>){
     if ($log_file=~/.*\/(\d+)_Cover_([^_]+)_([^_]+)\.log.txt/){
          my $genome=$1;
          my $cov=$2;
          my $prof=$3;
          my $code=$cov."_".$prof;
          $check_code{$code}{$genome}++;
          if (!defined($info{$genome})){
               print "0 repeats available for $genome\n";
               $ignored{$genome}=1;
          }
          my $result_file=$dir_crass."/".$genome."_Cover_".$cov."_".$prof.".crass_results.tsv";
          print "Reading $result_file -- $code // $genome \n";
          %seen=();
          if (-e $result_file){
               open my $tsv,"<",$result_file;
               while(<$tsv>){
                    chomp($_);
                    my @tab=split("\t",$_);
                    my $repeat=$tab[0];
                    my $revcomp=&revcomp($tab[0]);
                    if (defined($info{$genome}{$repeat})){} ## We know this repeat for this genome
                    elsif (defined($info{$genome}{$revcomp})){ $repeat=$revcomp;} ## Switching the repeat to its revcomp, which is the one we used
                    else{$repeat="";} ## Ignoring since we don't know this repeat or its revcomp for this genome
                    if ($repeat ne ""){
                         my $spacer=$tab[1];
                         if ($seen{$repeat}{$spacer}==1){} ## We ignore, we already saw this spacer (or its RC) for this repeat in this result file
                         else{
                              if (defined($info{$genome}{$repeat}{$spacer})){
                                   print "Found $tab[1] - expected for $tab[0]\n";
                                   $store{$genome}{$tab[0]}{$code}{"crass"}{"found"}++;
                              }
                              else{
                                   $store{$genome}{$tab[0]}{$code}{"crass"}{"extra"}++;
                              }
                              my $revcomp_spacer=&revcomp($spacer);
                              $seen{$repeat}{$spacer}=1;
                              $seen{$repeat}{$revcomp_spacer}=1;
                         }
                    }
               }
               close $tsv;
          }
          else{
               print "\t$code // $genome processed, but result file empty for CRASS \n";
          }
     }
     else{
          print "Pblm with $log_file\n";
          die("\n");
     }
}
## SE
%seen=();
foreach my $log_file (<$dir_se/*.log.txt>){
     if ($log_file=~/.*\/(\d+)_Cover_([^_]+)_([^_]+).log.txt/){
          my $genome=$1;
          my $cov=$2;
          my $prof=$3;
          my $code=$cov."_".$prof;
          $check_code{$code}{$genome}++;
          my $result_file=$dir_se."/".$genome."_Cover_".$cov."_".$prof."/Final_spacers.nr.denoised.hq.info.tsv";
          print "Reading $result_file -- $code // $genome \n";
          %seen=();
          if (-e $result_file){
               open my $tsv,"<",$result_file;
               while(<$tsv>){
                    chomp($_);
                    my @tab=split("\t",$_);
                    if ($tab[0] eq "Spacer id"){next;}
                    if (!defined($array_to_seq{$tab[5]})){
                         die("Pblm with $tab[5]\n");
                    }
                    my $tag_sing=0;
                    if ($tab[3]>1){$tag_sing=1;}
                    my $repeat=$array_to_seq{$tab[5]};
                    my $revcomp=&revcomp($repeat);

                    if (defined($info{$genome}{$repeat})){} ## We know this repeat for this genome
                    elsif (defined($info{$genome}{$revcomp})){ $repeat=$revcomp;} ## Switching the repeat to its revcomp, which is the one we used
                    else{$repeat="";} ## Ignoring since we don't know this repeat or its revcomp for this genome
                    if ($repeat ne ""){
                         my $spacer=$tab[1];
                         if ($seen{$repeat}{$spacer}==1){} ## We ignore, we already saw this spacer (or its RC) for this repeat in this result file
                         else{
                              if (defined($info{$genome}{$repeat}{$spacer})){
                                   print "Found $tab[1] - expected for $tab[0]\n";
                                   $store{$genome}{$repeat}{$code}{"spacer_extractor"}{"found"}++;
                                   if ($tag_sing==1){$store{$genome}{$repeat}{$code}{"spacer_extractor_nosing"}{"found"}++;}
                              }
                              else{
                                   $store{$genome}{$repeat}{$code}{"spacer_extractor"}{"extra"}++;
                                   if ($tag_sing==1){$store{$genome}{$repeat}{$code}{"spacer_extractor_nosing"}{"extra"}++;}
                              }
                              my $revcomp_spacer=&revcomp($spacer);
                              $seen{$repeat}{$spacer}=1;
                              $seen{$repeat}{$revcomp_spacer}=1;
                         }
                    }
               }
               close $tsv;
          }
          else{
               print "\t$code // $genome processed, but result file empty for SE \n";
          }
          ## Also look at multi repeats if it exists
          my $result_file_2=$dir_se."/".$genome."_Cover_".$cov."_".$prof."/Final_spacers.nr.denoised.multiarray.info.tsv";
          if (-e $result_file_2){
               print "Also reading $result_file_2 -- $code // $genome \n";
               open my $tsv,"<",$result_file_2;
               while(<$tsv>){
                    chomp($_);
                    my @tab=split("\t",$_);
                    if ($tab[0] eq "Spacer id"){next;}
                    my @tab_arrays=split(";",$tab[5]);
                    foreach my $array_code (@tab_arrays){
                         if (!defined($array_to_seq{$array_code})){
                              die("Pblm with $array_code\n");
                         }
                         my $tag_sing=0;
                         if ($tab[3]>1){$tag_sing=1;}
                         my $repeat=$array_to_seq{$array_code};
                         my $revcomp=&revcomp($repeat);

                         if (defined($info{$genome}{$repeat})){} ## We know this repeat for this genome
                         elsif (defined($info{$genome}{$revcomp})){ $repeat=$revcomp;} ## Switching the repeat to its revcomp, which is the one we used
                         else{$repeat="";} ## Ignoring since we don't know this repeat or its revcomp for this genome
                         if ($repeat ne ""){
                              my $spacer=$tab[1];
                              if ($seen{$repeat}{$spacer}==1){} ## We ignore, we already saw this spacer (or its RC) for this repeat in this result file
                              else{
                                   if (defined($info{$genome}{$repeat}{$spacer})){
                                        print "Found $tab[1] - expected for $tab[0]\n";
                                        $store{$genome}{$repeat}{$code}{"spacer_extractor"}{"found"}++;
                                        if ($tag_sing==1){$store{$genome}{$repeat}{$code}{"spacer_extractor_nosing"}{"found"}++;}
                                   }
                                   else{
                                        $store{$genome}{$repeat}{$code}{"spacer_extractor"}{"extra"}++;
                                        if ($tag_sing==1){$store{$genome}{$repeat}{$code}{"spacer_extractor_nosing"}{"extra"}++;}
                                   }
                                   my $revcomp_spacer=&revcomp($spacer);
                                   $seen{$repeat}{$spacer}=1;
                                   $seen{$repeat}{$revcomp_spacer}=1;
                              }
                         }
                    }
               }
               close $tsv;
          }
          else{
               print "\t$code // $genome processed, but result file empty for SE \n";
          }
     }
     else{
          print "Pblm with $log_file\n";
          die("\n");
     }
}
## MetaCRAST
%seen=();
foreach my $log_file (<$dir_mc/*.log.txt>){
     if ($log_file=~/.*\/(\d+)_Cover_([^_]+)_([^_]+)\.log.txt/){
          my $genome=$1;
          my $cov=$2;
          my $prof=$3;
          my $code=$cov."_".$prof;
          $check_code{$code}{$genome}++;
          # .metacrast_results.tsv
          my $result_file=$dir_mc."/".$genome."_Cover_".$cov."_".$prof.".metacrast_results.tsv";
          %seen=();
          print "Reading $result_file -- $code // $genome \n";
          if (-e $result_file){
               open my $tsv,"<",$result_file;
               while(<$tsv>){
                    chomp($_);
                    my @tab=split("\t",$_);
                    my $repeat=$tab[0];
                    my $revcomp=&revcomp($repeat);

                    if (defined($info{$genome}{$repeat})){} ## We know this repeat for this genome
                    elsif (defined($info{$genome}{$revcomp})){ $repeat=$revcomp;} ## Switching the repeat to its revcomp, which is the one we used
                    else{$repeat="";} ## Ignoring since we don't know this repeat or its revcomp for this genome
                    if ($repeat ne ""){
                         my $spacer=$tab[1];
                         if ($seen{$repeat}{$spacer}==1){} ## We ignore, we already saw this spacer (or its RC) for this repeat in this result file
                         else{
                              if (defined($info{$genome}{$repeat}{$spacer})){
                                   print "Found $tab[1] - expected for $tab[0]\n";
                                   $store{$genome}{$repeat}{$code}{"metacrast"}{"found"}++;
                              }
                              else{
                                   $store{$genome}{$repeat}{$code}{"metacrast"}{"extra"}++;
                              }
                              my $revcomp_spacer=&revcomp($spacer);
                              $seen{$repeat}{$spacer}=1;
                              $seen{$repeat}{$revcomp_spacer}=1;
                         }
                    }
               }
               close $tsv;
          }
          else{
               print "\t$code // $genome processed, but result file empty for MetaCRAST \n";
          }
     }
     else{
          print "Pblm with $log_file\n";
          die("\n");
     }
}

my %count;
my @method=("crass","spacer_extractor","spacer_extractor_nosing","metacrast");
open my $s1,">",$out_file;
print $s1 "code\tgenome\tarray\texpected\tcrass\tspacer_extractor\tspacer_extractor_nosing\tmetacrast\n";
foreach my $code (sort keys %check_code){
     foreach my $genome (sort keys %info){
          print $code."\t".$genome."\t".$check_code{$code}{$genome}."\n";
          if ($check_code{$code}{$genome}==""){$check_code{$code}{$genome}=0;}
          $count{$check_code{$code}{$genome}}++;
          if ($check_code{$code}{$genome}==3){ ## Processed and results from all 3
               foreach my $array (sort keys %{$info{$genome}}){
                    if (defined($store{$genome}{$array}{$code})){
                         # print "## We found some spacers for $genome - $array - $code at least once\n";
                         my $line=$code."\t".$genome."\t".$array."\t".$expected{$genome}{$array};
                         foreach my $met (@method){
                              if (!defined($store{$genome}{$array}{$code}{$met}{"found"})){$store{$genome}{$array}{$code}{$met}{"found"}=0;}
                              if (!defined($store{$genome}{$array}{$code}{$met}{"extra"})){$store{$genome}{$array}{$code}{$met}{"extra"}=0;}
                              $line.="\t".$store{$genome}{$array}{$code}{$met}{"found"}.";".$store{$genome}{$array}{$code}{$met}{"extra"};
                         }
                         print $line."\n";
                         print $s1 $line."\n";
                    }
                    else{
                         print "## We never found any spacer for $genome - $array - $code\n";
                    }
               }
          }
     }
}
close $s1;


print "## Simulations processed:\n";
foreach my $n (sort keys %count){
     print $n." methods: ".$count{$n}." simulations\n";
}


foreach my $genome (keys %ignored){
     print "Also we ignored $genome because no info matching our database\n";
}


sub revcomp(){
	my $nuc=$_[0];
	$nuc=~tr/atcguryswkmbdhvn/tagcayrswmkvhdbn/;
	$nuc=~tr/ATCGURYSWKMBDHVN/TAGCAYRSWMKVHDBN/;
	$nuc=reverse($nuc);
	return $nuc;
}
