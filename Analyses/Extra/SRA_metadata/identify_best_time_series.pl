#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
GetOptions ('help' => \$h, 'h' => \$h); 
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to look for time series based on host subject id or lat/longq
# Arguments :
# none
";
    die "\n";
}

my $in_file="Data/Additional_data/Sample_table_all_subject_and_time.tsv.gz";
my $in_runs="Data/Spacer_db/Runs_to_ecosystem_and_sequencing_and_study_for_db-Jul28-24.tsv";
my $out_file="Overview_time_series.txt";
my $out_dir="Analyses/Spacer_database/Individual_time_series/";

my $min_date=5;

my %bs_to_lib;
print "Reading $in_runs\n";
open my $tsv,"<",$in_runs;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "library"){next;}
    $bs_to_lib{$tab[7]}{$tab[0]}=1;
}
close $tsv;

my %store;
print "Reading $in_file\n";
open my $tsv,"gunzip -c $in_file |";
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "BioSample"){next;}
    my $biosample=$tab[0];
    my $bioproject=$tab[1];
    my $eco=$tab[2];
    my $subject_id=$tab[6];
    my $subject=$tab[6];
    if ($subject eq "Blank" || $subject eq "Mock"){next;}
    my $date=$tab[7];
    if ($bioproject ne ""){
        my $tag_sub=&check_sub($subject);
        my $tag_date=&check_date($date);
        if ($tag_date==1 && $tag_sub==1){
            my $date_c=&clean_date($date);
            if ($date_c ne ""){
                if (!defined($bs_to_lib{$biosample})){
                    print "$biosample was maybe interesting but we did not process any library from there\n";
                }
                else{
                    foreach my $library (keys %{$bs_to_lib{$biosample}}){
                        $store{$bioproject}{$subject_id}{$eco}{$date_c}{$biosample}{$library}{"line"}=$_;
                    }
                }
            }
        }
    }
}
close $tsv;


my $n_series=0;
open my $s1,">",$out_file;
print $s1 "series_name\tbioproject\tsubject\tecosystem\tdate\tn_sample\tlist_dates\n";
foreach my $bp (sort keys %store){
    $n_series=0;
    foreach my $sub (sort keys %{$store{$bp}}){
        foreach my $eco (sort keys %{$store{$bp}{$sub}}){
            my @list_dates=sort keys %{$store{$bp}{$sub}{$eco}};
            my $n_date=scalar(@list_dates);
            my $n_sample=0;
            foreach my $date (@list_dates){
                $n_sample+=scalar(keys %{$store{$bp}{$sub}{$eco}{$date}});
            }
            if ($n_date>=$min_date){
                $n_series++;
                my $series_name=$bp."_".sprintf("%03d",$n_series);
                print $s1 $series_name."\t".$bp."\t".$sub."\t".$eco."\t".$n_date."\t".$n_sample."\t".join("|",@list_dates)."\n";
                # 
                my $series_dir=$out_dir.$series_name."/";
                if (!(-d $series_dir)){&run_cmd("mkdir $series_dir");}
                my $series_file=$series_dir.$series_name."_samples.tsv";
                open my $s2,">",$series_file;
                print $s2 "bioproject\tsubject_id\tecosystem\tdate\tlibrary\tbiosample\torganism\tloc_name\tmedium\tori_date\n";
                foreach my $date (@list_dates){
                    foreach my $sample (sort keys %{$store{$bp}{$sub}{$eco}{$date}}){
                        foreach my $library (sort keys %{$store{$bp}{$sub}{$eco}{$date}{$sample}}){
                            my @t=split("\t",$store{$bp}{$sub}{$eco}{$date}{$sample}{$library}{"line"});
                            print $s2 $bp."\t".$sub."\t".$eco."\t".$date."\t".$library."\t".$sample."\t".$t[3]."\t".$t[4]."\t".$t[5]."\t".$t[7]."\n";
                        }
                    }
                }
                close $s2;
            }
        }
    }
}
close $s1;


sub check_sub(){
    my $sub=$_[0];
    my $tag=0;
    if (&is_empty($sub)==0){}
    else{$tag=1;}
    return $tag;
}

sub check_date(){
    my $date=$_[0];
    my $tag=0;
    if (&is_empty($date)==0){}
    elsif ($date=~/^\d\d\d\d$/){$tag=0;}
    else{
        if ($date=~/\// || $date=~/-/ || $date=~/\./){
            $tag=1;
        }
        elsif($date=~/20\d{6}/){
            $tag=1;
        }
        else{
            print "not sure about $date ?\n";
            <STDIN>;
        }
    }
    return $tag;
}

sub clean_date(){
    my $date=$_[0];
    my $result="";
    if ($date=~/(\d{4})-(\d{2})-(\d{2})/){
        $result=$date;
    }
    elsif ($date=~/(\d{4})-(\d{2})/){
        $result=$date."-01";
    }
    elsif ($date=~/10-(\d{2})-(\d{2})/){
        $result="2010-".$1."-".$2;
    }
    elsif($date=~/(\d{2})\/(\d{2})\/(\d{4})/){
        $result=$3."-".$1."-".$2;
        if ($1>12){
            die("Pblm - month is $1 ? -- $date\n");
        }
    }
    elsif($date=~/(20\d{2})(\d{2})(\d{2})/){
        $result=$1."-".$2."-".$3;
        if ($2>12){
            die("Pblm - month is $2 ? -- $date\n");
        }
    }
    elsif($date=~/(\d{1,2})\.(\d{1,2})\.(\d{4})/){
        $result=$3."-".$1."-".$2;
        if ($1>12){
            die("Pblm - month is $1 ? -- $date\n");
        }
    }
    elsif($date=~/(\d{1,2})\/(\d{1,2})\/(\d{2})/){
        $result="20".$3."-".$1."-".$2;
        if ($1>12){
            die("Pblm - month is $1 ? -- $date\n");
        }
    }
    elsif($date=~/([A-z]{3})-(\d{4})/){
        $result=$2."-".&tr_month($1)."-01";
    }
    elsif($date=~/\d{4}\/\d{4}/){
        ## This happens, sometimes they just give a year range, not very useful
    }
    elsif($date eq "01-002-2007"){
        ## This happens, sometimes this is just weird, we ignore
    }
    else{
        print "Unexpected format $date\n";
        <STDIN>;
    }
    return $result;
}

sub tr_month(){
    my $code=$_[0];
    my $n=-1;
    if ($code eq "Jan"){$n="01";}
    elsif ($code eq "Feb"){$n="02";}
    elsif ($code eq "Mar"){$n="03";}
    elsif ($code eq "Apr"){$n="04";}
    elsif ($code eq "May"){$n="05";}
    elsif ($code eq "Jun"){$n="06";}
    elsif ($code eq "Jul"){$n="07";}
    elsif ($code eq "Aug"){$n="08";}
    elsif ($code eq "Sept"){$n="09";}
    elsif ($code eq "Oct"){$n="10";}
    elsif ($code eq "Nov"){$n="11";}
    else{
        die("Unknown month $code\n");
    }
}

sub is_empty(){
    my $text=$_[0];
    my $tag=0;
    if ($text eq "" || $text eq "not provided" || $text eq "NA" || $text eq "Not provided" || $text eq "missing" || $text eq "not collected" || $text eq "Not Applicable" || $text eq "not available" || $text eq "Not available" || $text eq "not applicable" || $text eq "Not applicable" || $text eq "unknown" || $text eq "none" || $text eq "not recorded"){}
    elsif($text=~/missing/i){}
    else{$tag=1;}
    return $tag;
}

sub run_cmd {
	my $cmd=$_[0];
	if ($_[1] ne "veryquiet"){print "$cmd\n";}
	if ($_[1] ne "dummy"){
		my $out=`$cmd`;
		if ($_[1] ne "quiet" && $_[1] ne "veryquiet"){
			if ($_[1] eq "stderr"){print STDERR "$out\n";}
			else{print "$out\n";}
		}
		return($out);
	}
	else{
		print " ### dummy run\n";
	}
}



