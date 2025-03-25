#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
	print "# Script to look at the results of biosample API and classify biosample in simple ecosystem categories
# Arguments :
# none
";
	die "\n";
}

my $in_file="All_biosample_info_NCBIAPI.tsv";
my $wanted_file="SraRunInfo_MetaG_not_Amplicon_Public-All-Combined-Nov_18_2023.csv";
my $curated_keyword_file="Final_keyword_to_category.tsv";
my $file_info_category="Ordered_category_scores.txt";
my $combination_ok="Acceptable_combination.txt";
my $out_file="SRA_Eco_clean.tsv";
my $out_empty="Dump_missing_eco.txt";
my $out_missing="list_missing_biosample_new.txt";
my $out_contradictions="list_contradictions.txt";


print "Loading curated keyword to metadata link\n";
## Start with information from the 'well-formatted' fields like env_broad and env_local and env_medium
## Then look at Env package, including whether or not they used a specific env package (which gives a hint at what type of sample this is)
## Then move on to best guesses based on host , isolation source, description, or title
## NOTE - we try another approach below which is less strict, and is based on a score to prioritize different categories (and if ties, then we use the order of the keywords to solve)
## Now taken from a file
# my @ordered_keywords=("DATA_env_local_scale","DATA_env_medium","DATA_Env_Package","DATA_env_broad_scale",
# "DATA_air environmental package","DATA_built environment environmental package","DATA_humam oral environmental package","DATA_human associated environmental package","DATA_human-associated environmental package","DATA_human gut environmental package","DATA_human-gut environmental package","DATA_human oral environmental package","DATA_human skin environmental package","DATA_human vaginal environmental package","DATA_soil environmental package","DATA_wastewater/sludge environmental package","DATA_water environmental package","DATA_microbial mat/biofilm environmental package","DATA_plant-associated environmental package","DATA_sediment environmental package",
# "Package",
# "DATA_sample comment","DATA_host","DATA_isolation_source","DATA_Organism","DATA_Descr","Title");
my @ordered_keywords;
my %score_keywords;
print "loading info about category importance from $file_info_category\n";
open my $tsv,"<",$file_info_category;
while(<$tsv>){
      chomp($_);
      my @tab=split("\t",$_);
      if ($tab[0] eq ""){next;}
      if ($tab[1] eq ""){next;}
      if ($tab[0]=~/^#/){next;}
      push(@ordered_keywords,$tab[0]);
      $score_keywords{$tab[0]}=$tab[1];
}
close $tsv;

my %keyword_ref;
open my $tsv,"<",$curated_keyword_file;
while(<$tsv>){
      chomp($_);
      my @tab=split("\t",$_);
      if ($tab[0] eq "Category"){next;}
      if ($tab[0] ne "" && $tab[1] ne ""){
            if ($tab[2] ne ""){$keyword_ref{$tab[0]}{$tab[1]}{"lvl_1"}=$tab[2];}
            if ($tab[3] ne ""){$keyword_ref{$tab[0]}{$tab[1]}{"lvl_2"}=$tab[3];}
            if ($tab[4] ne ""){$keyword_ref{$tab[0]}{$tab[1]}{"biofilm"}=$tab[4];}
            if ($tab[0] eq "Title"){
                  my $alt="DATA_Title";
                  if ($tab[2] ne ""){$keyword_ref{$alt}{$tab[1]}{"lvl_1"}=$tab[2];}
                  if ($tab[3] ne ""){$keyword_ref{$alt}{$tab[1]}{"lvl_2"}=$tab[3];}
                  if ($tab[4] ne ""){$keyword_ref{$alt}{$tab[1]}{"biofilm"}=$tab[4];}
            }
            if ($tab[0] eq "DATA_Organism"){
                  my $alt="Organism";
                  if ($tab[2] ne ""){$keyword_ref{$alt}{$tab[1]}{"lvl_1"}=$tab[2];}
                  if ($tab[3] ne ""){$keyword_ref{$alt}{$tab[1]}{"lvl_2"}=$tab[3];}
                  if ($tab[4] ne ""){$keyword_ref{$alt}{$tab[1]}{"biofilm"}=$tab[4];}
            }
      }
}
close $tsv;

my %links_ok;
open my $tsv,"<",$combination_ok;
while(<$tsv>){
      chomp($_);
      my @tab=split(";",$_);
      $links_ok{$tab[1]}{$tab[0]}=1;
}
close $tsv;

my %check_biosample;
open my $tsv,"<",$wanted_file;
while(<$tsv>){
	chomp($_);
      my $line=$_;
      ## Need to clean-up because comma-separated but sometimes also comma in text
      while($line=~/(.*)\"([^\"]+)\"(.*)/){
            my $pref=$1;
            my $content=$2;
            my $suf=$3;
            $content=~s/,//g;
            $line=$pref.$content.$suf;
      }
	my @tab=split(",",$line);
      $check_biosample{$tab[25]}=2;
}
close $tsv;


print "Linking BioSamples to simplified ecosystems\n";
my $c_b="";
my %tmp;
open my $s3,">",$out_contradictions;
open my $s2,">",$out_empty;
open my $s1,">",$out_file;
open my $tsv,"<",$in_file;
while(<$tsv>){
      chomp($_);
      my @tab=split("\t",$_);
      if ($tab[0] eq "Query"){next;}
      if ($c_b ne $tab[0]){
            if ($c_b ne "" && defined($check_biosample{$c_b}) && $check_biosample{$c_b}>=2){
                  print "Processing $c_b\n";
                  my ($eco,$reason)=&guess_eco_byscore(\%tmp,\%keyword_ref,\@ordered_keywords,\%score_keywords,$s3);
                  print $s1 $c_b."\t".$eco."\t".$reason."\n";
                  print $c_b."\t".$eco."\t".$reason."\n";
                  my @t=split(";",$eco);
                  if ($#t>=1){}
                  else{
                        foreach my $key (sort keys %tmp){
                              print $s2 $c_b."\t".$key."\t".$tmp{$key}."\n";
                        }
                        # <STDIN>;
                  }
		  $check_biosample{$c_b}=3;
            }
            %tmp=();
            $c_b=$tab[0];
      }
      $tmp{$tab[1]}=$tab[2];
}
close $tsv;
if ($c_b ne "" &&  defined($check_biosample{$c_b}) && $check_biosample{$c_b}>=2){
      print "Processing $c_b\n";
      my ($eco,$reason)=&guess_eco_byscore(\%tmp,\%keyword_ref,\@ordered_keywords,\%score_keywords,$s3);
      print $c_b."\t".$eco."\t".$reason."\n";
      my @t=split(";",$eco);
      if ($#t>=1){}
      else{
            foreach my $key (sort keys %tmp){
                  print $s2 $c_b."\t".$key."\t".$tmp{$key}."\n";
            }
            # <STDIN>;
      }
      $check_biosample{$c_b}=3;
}
else{
      print "Not processing the last: --".$c_b."--\n";
}
close $s1;
close $s2;
close $s3;
open my $s2,">",$out_missing;
foreach my $c_b (sort keys %check_biosample){
	if ($check_biosample{$c_b}==2){
		print "--".$c_b."-- is missing\n";
		print $s2 $c_b."\n";
	}
}
close $s2;




sub guess_eco_byscore(){
      my $hash=$_[0];
      my $keyword_hash=$_[1];
      my $keyword_list=$_[2];
      my $keyword_score=$_[3];
      my $out_c=$_[4];
      my $lvl_1="";
      my $lvl_2="";
      my $biofilm="";
      my $reason="";
      my $tag=0;
      my %local_score;
      my $i=0;
      foreach my $category(@{$keyword_list}){
            $tag=0;
            $i++;
            if (defined($$hash{$category})){
                  # print "We found $category --> $$hash{$category} --> ".$$keyword_hash{$category}{$$hash{$category}}{"lvl_1"}.";".$$keyword_hash{$category}{$$hash{$category}}{"lvl_2"}.";".$$keyword_hash{$category}{$$hash{$category}}{"biofilm"}."\n";
                  my $field=$category.":".$$hash{$category};
                  my $value_1="";
                  my $value_2="";
                  my $value_biofilm="";
                  ### Loading value corresponding to this keyword in this category, or to any keyword in this category for the few metadata where the presence of the package itself is the metadata
                  if ($$keyword_hash{$category}{$$hash{$category}}{"lvl_1"} ne ""){$value_1=$$keyword_hash{$category}{$$hash{$category}}{"lvl_1"};}
                  elsif(defined($$keyword_hash{$category}{"Package_itself"}) && $$keyword_hash{$category}{"Package_itself"}{"lvl_1"} ne ""){$value_1=$$keyword_hash{$category}{"Package_itself"}{"lvl_1"};}
                  if ($$keyword_hash{$category}{$$hash{$category}}{"lvl_2"} ne ""){$value_2=$$keyword_hash{$category}{$$hash{$category}}{"lvl_2"};}
                  elsif(defined($$keyword_hash{$category}{"Package_itself"}) && $$keyword_hash{$category}{"Package_itself"}{"lvl_2"} ne ""){$value_2=$$keyword_hash{$category}{"Package_itself"}{"lvl_2"};}
                  if ($$keyword_hash{$category}{$$hash{$category}}{"biofilm"} ne ""){$value_biofilm=$$keyword_hash{$category}{$$hash{$category}}{"biofilm"};}
                  elsif(defined($$keyword_hash{$category}{"Package_itself"}) && $$keyword_hash{$category}{"Package_itself"}{"biofilm"} ne ""){$value_biofilm=$$keyword_hash{$category}{"Package_itself"}{"biofilm"};}
                  if ($value_1 ne ""){
                        $local_score{"lvl1"}{$value_1}{"total_score"}+=$$keyword_score{$category};
                        $local_score{"lvl1"}{$value_1}{"n_cat"}+=$$keyword_score{$category};
                        $local_score{"lvl1"}{$value_1}{"order"}=$i; ## Because the keywords are ordered, the first value to come in will be taken in case of tie
                        $local_score{"lvl1"}{$value_1}{"rationale"}.=$field.";";
                  }
                  if ($value_2 ne ""){
                        $local_score{"lvl2"}{$value_2}{"total_score"}+=$$keyword_score{$category};
                        $local_score{"lvl2"}{$value_2}{"n_cat"}+=$$keyword_score{$category};
                        $local_score{"lvl2"}{$value_2}{"order"}=$i; ## Because the keywords are ordered, the first value to come in will be taken in case of tie
                        $local_score{"lvl2"}{$value_2}{"rationale"}.=$field.";";
                  }
                  if ($value_biofilm eq "yes" && $biofilm eq ""){
                        $biofilm="yes";
                  }
            }
      }
      if (defined($local_score{"lvl1"})){
            ## Sorting is:
            # First on score (higher is better)
            # Second on the number of consistent categories (higher is better)
            # Third on the order which will be unique for each category (because it's the keyword order in the ordered list, lower is better)
            my @list_lvl_1=sort {$local_score{"lvl1"}{$b}{"total_score"} <=> $local_score{"lvl1"}{$a}{"total_score"} or $local_score{"lvl1"}{$b}{"n_cat"} <=> $local_score{"lvl1"}{$a}{"n_cat"} or $local_score{"lvl1"}{$a}{"order"} <=> $local_score{"lvl1"}{$b}{"order"}} keys %{$local_score{"lvl1"}};
            if (scalar(@list_lvl_1)==1){
                  print "All lvl 1 point to ".$list_lvl_1[0]." - that's easy\n";
            }
            else{
                  print "We have multiple lvl 1 possible - ".join(";",@list_lvl_1)."\n";
                  foreach my $value (@list_lvl_1){
                        print "\t".$value."\t".$local_score{"lvl1"}{$value}{"total_score"}."\t".$local_score{"lvl1"}{$value}{"n_cat"}."\t".$local_score{"lvl1"}{$value}{"order"}."\n";
                  }
                  if ($list_lvl_1[0] eq "Animal-associated"){ ## Check if we want to be more specific
                        if ($list_lvl_1[1] eq "Human-associated" || $list_lvl_1[1] eq "Birds-associated" || $list_lvl_1[1] eq "Fish-associated" || $list_lvl_1[1] eq "Mammals-associated"){
                              print "We authorized the replacement of $list_lvl_1[0] by $list_lvl_1[1]\n";
                              $list_lvl_1[0]=$list_lvl_1[1];
                        }
                  }
                  print "We take the best score: ".$list_lvl_1[0]."\n";
                  # <STDIN>;
            }
            $lvl_1=$list_lvl_1[0];
            chop($local_score{"lvl1"}{$lvl_1}{"rationale"});
            $reason=$local_score{"lvl1"}{$lvl_1}{"rationale"}."|";
      }
      else{
            print "No indication for lvl 1\n";
      }
      if (defined($local_score{"lvl2"})){
            my @list_lvl_2=sort {$local_score{"lvl2"}{$b}{"total_score"} <=> $local_score{"lvl2"}{$a}{"total_score"} or $local_score{"lvl2"}{$b}{"n_cat"} <=> $local_score{"lvl2"}{$a}{"n_cat"} or $local_score{"lvl2"}{$a}{"order"} <=> $local_score{"lvl2"}{$b}{"order"}} keys %{$local_score{"lvl2"}};
            if (scalar(@list_lvl_2)==1){
                  print "All lvl 2 point to ".$list_lvl_2[0]." - that's easy\n";
                  if ($lvl_1 ne ""){
                        if ($links_ok{$list_lvl_2[0]}{$lvl_1}==1){
                              print "\t and this is a known combination with lvl 1 $lvl_1, all good\n";
                              $lvl_2=$list_lvl_2[0];
                              chop($local_score{"lvl2"}{$lvl_2}{"rationale"});
                              $reason.=$local_score{"lvl2"}{$lvl_2}{"rationale"};
                        }
                        else{
                              print "Trying to add $list_lvl_2[0] to an existing lvl 1 $lvl_1, this is not allowed\n";
                              print $out_c $list_lvl_2[0]."\t".$lvl_1."\n";
                        }
                  }
                  else{
                        $lvl_2=$list_lvl_2[0];
                        chop($local_score{"lvl2"}{$lvl_2}{"rationale"});
                        $reason.=$local_score{"lvl2"}{$lvl_2}{"rationale"};
                  }
            }
            else{
                  print "We have multiple lvl 2 possible - ".join(";",@list_lvl_2)."\n";
                  foreach my $value (@list_lvl_2){
                        print $value."\t".$local_score{"lvl2"}{$value}{"total_score"}."\t".$local_score{"lvl2"}{$value}{"n_cat"}."\t".$local_score{"lvl2"}{$value}{"order"}."\n";
                  }
                  print "We would initially take the best score: ".$list_lvl_2[0]."\n";
                  if ($lvl_1 ne ""){
                        if ($links_ok{$list_lvl_2[0]}{$lvl_1}==1){
                              print "\t and this is a known combination with lvl 1 $lvl_1, all good\n";
                              $lvl_2=$list_lvl_2[0];
                        }
                        else{
                              print "Trying to add $list_lvl_2[0] to an existing lvl 1 $lvl_1, this is not allowed\n";
                              print $out_c $list_lvl_2[0]."\t".$lvl_1."\n";
                              print "Checking if we have another lvl 2 compatible with $lvl_1 maybe\n";
                              foreach my $value (@list_lvl_2){
                                    if ($links_ok{$value}{$lvl_1}==1){
                                          print "\t we found one - $value ! we take this one\n";
                                          $lvl_2=$value;
                                    }
                              }
                        }
                  }
                  else{
                        print "\tand we take it because we don't have lvl 1 anyway\n";
                        $lvl_2=$list_lvl_2[0];
                  }
                  # <STDIN>;
                  if ($lvl_2 ne ""){
                        chop($local_score{"lvl2"}{$lvl_2}{"rationale"});
                        $reason.=$local_score{"lvl2"}{$lvl_2}{"rationale"};
                  }
            }
      }
      else{
            print "No indication for lvl 2\n";
      }
      if ($reason=~/\|$/){chop($reason);}
      if ($lvl_1 eq ""){
            if ($lvl_2 eq ""){
                  $lvl_1="Unknown";
                  $lvl_2="Unknown";
            }
            else{
                  print "We have lvl 2 $lvl_2 without lvl 1 $lvl_1\n";
                  print $reason."\n";
                  print "We should be able to guess lvl 1\n";
                  if ($lvl_2=~/system/ || $lvl_2 eq "Skin"){
                        $lvl_1="Animal-associated";
                  }
                  if ($lvl_1 eq ""){die("\n");}
            }
      }
      elsif($lvl_2 eq ""){
            $lvl_2="Unknown";
      }
      elsif($lvl_2 eq "Other"){
            $lvl_2="Unknown"; ## We put together Other & Unknown
      }
      # 
      if ($lvl_2 ne "Unknown"){
            if ($links_ok{$lvl_2}{$lvl_1}==1){}
            else{
                  print "We link $lvl_2 to $lvl_1, this is not allowed\n";
                  print $out_c $lvl_2."\t".$lvl_1."\n";
                  $lvl_2="Unknown";
                  # <STDIN>;
            }
      }
      if ($biofilm eq ""){$biofilm="Unknown";}
      elsif($biofilm eq "yes"){$biofilm="Biofilm";}
      my $eco=$lvl_1.";".$lvl_2.";".$biofilm;
      print "Final result is $eco\n";
      print "for reasons $reason\n";
      # # <STDIN>;
      return ($eco,$reason);
}
