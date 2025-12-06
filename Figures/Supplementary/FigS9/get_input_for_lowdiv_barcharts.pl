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

### NOTE - we use "frozen" sometimes in there, but it's really "low diversity" (becomes frozen if it's low diversity across samples, that's a next step)
my $file_array_info="Data/Spacer_db/Array_info_filtered_for_db-Oct24-25.tsv";
my $file_sample_info="Data/Spacer_db/Runs_to_ecosystem_and_sequencing_and_study_for_db-Jul28-24.tsv";
my $in_alpha_base="Analyses/Spacer_database/spacer_sets_alphadiv.tsv";

my $out_file="input_lowdiv_barcharts.Jan3.tsv";
my $out_file_bis="all_counts.Jan3.tsv";
my $out_file_nosing="input_lowdiv_barcharts.nosingleton.Jan3.tsv";
my $out_file_nosing_bis="all_counts.nosingleton.Jan3.tsv";

my %info_array;
open my $tsv,"<",$file_array_info;
print "Reading $file_array_info\n";
while(<$tsv>){
      chomp($_);
      my @tab=split("\t",$_);
      if ($tab[0] eq "repeat_cluster"){next;}
      # $info_array{$tab[0]}{"type"}=$tab[1];
      my @t=split("-",$tab[1]);
      $info_array{$tab[0]}{"type"}=$t[0];
      $info_array{$tab[0]}{"taxo_source"}=$tab[2];
      $info_array{$tab[0]}{"taxo_class"}=$tab[4];
      $info_array{$tab[0]}{"taxo_genus"}=$tab[6];
      my @t=split((";",$tab[6]));
      $info_array{$tab[0]}{"taxo_phylum"}=$t[0].";".$t[1];
      if ($t[2] eq "c__unclassified"){
           # print $tab[4]." is actually not a class\n";
           $info_array{$tab[0]}{"taxo_class"}="NA";
      }
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
my %store_set_type;
my %store_set_type_nosing;
my %eco_to_set;
my %taxon_to_set;
open my $tsv,"<",$in_alpha_base;
print "Reading $in_alpha_base\n";
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "array"){
        # print join("\t",@tab[0..6])."\t".$tab[12]."\t".$tab[13]."\tP_min50\tP_max5\ttype\n";
        next;
    }
    my $set=$tab[0]."_".$tab[1];
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
    if ($type ne "low_coverage"){
	    ## Also compute type without singletons
	    my $t_ns=$tab[2]-$tab[6];
	    my $p_min50_nosing=$tab[12]/$t_ns;
	    my $type_nosing="NA";
	    if ($p_min50_nosing<=0.25){
		    $type_nosing="standard";
	    }
	    elsif($p_min50_nosing>=0.50){
		    $type_nosing="frozen";
	    }
	    else{
		    $type_nosing="other";
	    }
         ## Since we are mostly interested in high coverage anyway
         # print join("\t",@tab[0..6])."\t".$tab[12]."\t".$tab[13]."\t".$p_min50."\t".$p_max5."\t".$type."\n";
         $store_set_type{$set}=$type;
	    $store_set_type_nosing{$set}=$type_nosing;
         if (defined($info_sample{$tab[1]}{"eco"}) && $info_sample{$tab[1]}{"eco"} ne "Unknown"){
              $eco_to_set{$info_sample{$tab[1]}{"eco"}}{$set}=1;
         }
         if (defined($info_array{$tab[0]}{"taxo_class"}) && $info_array{$tab[0]}{"taxo_class"} ne "NA"){
              $taxon_to_set{$info_array{$tab[0]}{"taxo_class"}}{$set}=1;
         }
     }
     # <STDIN>;
     $i++;
     if ($i % 100000 == 0){print $i." ... \n";}
     # if ($i>=500000){last;}
}
close $tsv;
open my $s1,">",$out_file;
print $s1 "category\tgroup\tn_obs\taverage_lowdiv\taverage_sub\tstdev_sub\tstderr_sub\n";
open my $s2,">",$out_file_bis;
print $s2 "category\tgroup\tsampled\tn_obs\tn_lowdiv\tn_standard\n";
print "## General stats\n";
my $prefix="all\tall";
my ($n_obs,$global,$avg,$stdev,$stderr)=&process(\%store_set_type,50000,$prefix,$s2,\%store_set_type);
print $s1 "all\tall\t".$n_obs."\t".$global."\t".$avg."\t".$stdev."\t".$stderr."\n";
print "## By ecosystem\n";
foreach my $eco (sort keys %eco_to_set){
     print "\t$eco\n";
	my $prefix="ecosystem\t".$eco;
     my ($n_obs,$global,$avg,$stdev,$stderr)=&process($eco_to_set{$eco},5000,$prefix,$s2,\%store_set_type);
     print $s1 "ecosystem\t".$eco."\t".$n_obs."\t".$global."\t".$avg."\t".$stdev."\t".$stderr."\n";
}
print "## By taxon\n";
foreach my $taxon (sort keys %taxon_to_set){
     print "\t$taxon\n";
	my $prefix="taxon\t".$taxon;
     my ($n_obs,$global,$avg,$stdev,$stderr)=&process($taxon_to_set{$taxon},500,$prefix,$s2,\%store_set_type);
     if ($global ne "NA"){
          print $s1 "taxon\t".$taxon."\t".$n_obs."\t".$global."\t".$avg."\t".$stdev."\t".$stderr."\n";
     }
}
close $s1;

print "#### NOW CALCULATING THE SAME THING WHILE EXCLUDING ALL SINGLETONS ####\n";

open my $s1,">",$out_file_nosing;
print $s1 "category\tgroup\tn_obs\taverage_lowdiv\taverage_sub\tstdev_sub\tstderr_sub\n";
open my $s2,">",$out_file_nosing_bis;
print $s2 "category\tgroup\tsampled\tn_obs\tn_lowdiv\tn_standard\n";
print "## General stats\n";
my $prefix="all\tall";
my ($n_obs,$global,$avg,$stdev,$stderr)=&process(\%store_set_type,50000,$prefix,$s2,\%store_set_type_nosing);
print $s1 "all\tall\t".$n_obs."\t".$global."\t".$avg."\t".$stdev."\t".$stderr."\n";
print "## By ecosystem\n";
foreach my $eco (sort keys %eco_to_set){
     print "\t$eco\n";
	my $prefix="ecosystem\t".$eco;
     my ($n_obs,$global,$avg,$stdev,$stderr)=&process($eco_to_set{$eco},5000,$prefix,$s2,\%store_set_type_nosing);
     print $s1 "ecosystem\t".$eco."\t".$n_obs."\t".$global."\t".$avg."\t".$stdev."\t".$stderr."\n";
}
print "## By taxon\n";
foreach my $taxon (sort keys %taxon_to_set){
     print "\t$taxon\n";
	my $prefix="taxon\t".$taxon;
     my ($n_obs,$global,$avg,$stdev,$stderr)=&process($taxon_to_set{$taxon},500,$prefix,$s2,\%store_set_type_nosing);
     if ($global ne "NA"){
          print $s1 "taxon\t".$taxon."\t".$n_obs."\t".$global."\t".$avg."\t".$stdev."\t".$stderr."\n";
     }
}
close $s1;



sub process(){
     my $hash=$_[0];
     my $n_to_sample=$_[1];
	my $prefix=$_[2];
	my $sout=$_[3];
	my $hash_type=$_[4];
     my @list_sets=keys %{$hash};
     my $n_obs=scalar(@list_sets);
     if ($n_to_sample > ($n_obs/2)){
          print "We return NA because we only have ".$n_obs." observations and we asked for ".$n_to_sample." to subsample\n";
          return($n_obs,"NA","NA","NA");
     }
     else{
          print "We have $n_obs observations\n";
          print "Global stats:\n";
          my %tmp;
          foreach my $set (@list_sets){
               $tmp{"total"}++;
               # $tmp{$store_set_type{$set}}++;
			$tmp{$$hash_type{$set}}++;
			# print $set."\t".$store_set_type{$set}."\t".$$hash_type{$set}."\n";
          }
          my $global_frozen_ratio=$tmp{"frozen"}/($tmp{"frozen"}+$tmp{"standard"}+$tmp{"other"});
		print $sout $prefix."\tcomplete\t".($tmp{"frozen"}+$tmp{"standard"}+$tmp{"other"})."\t".$tmp{"frozen"}."\t".$tmp{"standard"}."\n";
          print "Average frozen in this category: ".$global_frozen_ratio."\n";
          print "We will do 10 subsamples of 10000 sets to get some standard error\n";
          my @stats;
          my $n=10;
          while ($n>0){
               ## subsample 10000 sets
               %tmp=();
               &fisher_yates_shuffle(\@list_sets);
               for (my $i=0;$i<$n_to_sample;$i++){
                    $tmp{"total"}++;
                    # $tmp{$store_set_type{$list_sets[$i]}}++;
				$tmp{$$hash_type{$list_sets[$i]}}++;
               }
               foreach my $type (keys %tmp){
                    print $type."\t".$tmp{$type}."\t".($tmp{$type}/$tmp{"total"}*100)."\n";
               }
               my $frozen_ratio=$tmp{"frozen"}/($tmp{"frozen"}+$tmp{"standard"}+$tmp{"other"});
			print $sout $prefix."\tsubsample\t".($tmp{"frozen"}+$tmp{"standard"}+$tmp{"other"})."\t".$tmp{"frozen"}."\t".$tmp{"standard"}."\n";
               push(@stats,$frozen_ratio);
               $n--;
          }
          my $average=&average(\@stats);
          my $stdev=&stdev(\@stats);
          my $stderr=$stdev/scalar(@stats);
          print "For full dataset: average is ".$average." and standard error is ".$stderr." for 10 samples of ".$n_to_sample." sets\n";
          print "#####################\n";
          return($n_obs,$global_frozen_ratio,$average,$stdev,$stderr);
     }
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

sub fisher_yates_shuffle {
        my $deck = shift;  # $deck is a reference to an array
        return unless @$deck; # must not be empty!
        my $i = @$deck;
        while (--$i) {
                my $j = int rand ($i+1);
                @$deck[$i,$j] = @$deck[$j,$i];
        }
}

sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty arrayn");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}

sub stdev{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}
