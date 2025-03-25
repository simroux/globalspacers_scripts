#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;

GetOptions ('help' => \$h, 'h' => \$h);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
    print "# Script to compare the results obtained on all simulations
# Arguments :
# none
";
    die "\n";
}

my $result="Compare_spacer_recovery_sub200.tsv";
my $out_file=$result;
$out_file=~s/\.tsv$//;
$out_file=$out_file."_metrics.tsv";

my %store;
my %header;
open my $tsv,"<",$result;
while(<$tsv>){
     chomp($_);
     my @tab=split("\t",$_);
     if ($tab[0] eq "code"){
          for (my $i=4;$i<=$#tab;$i++){
               $header{$i}=$tab[$i];
          }
          next;
     }
     my $code=$tab[0];
     my $exp=$tab[3];
     foreach my $i (keys %header){
          my @t=split(";",$tab[$i]);
          my $total_found=$t[0]+$t[1];
          my $recall=$t[0]/$exp; ## Found over expected
          if ($recall>1){
               die("Pblm recall $recall for $code / $exp / $header{$i}\n");
          }
          print $code."\t".$tab[1]."\t".$header{$i}."\t".$recall."\n";
          $store{$header{$i}}{$code}{"recall"}.=$recall.";";
          my $fdr=0;
          if ($total_found > 0){
               $fdr=$t[1]/($t[0]+$t[1]); ## Extra over Found + Extra
               ## Only compute FDR if there are some found, otherwise it's 0
          }
          $store{$header{$i}}{$code}{"fdr"}.=$fdr.";";
          my $f1=(2*$t[0])/(2*$t[0]+$t[1]+($exp-$t[0]));
          ## 2*Found over 2*Found + Wrong + (Total - Found)
          $store{$header{$i}}{$code}{"f1"}.=$f1.";";
     }
}
close $tsv;

my $den=sqrt(200);
print "###### LONG #######\n";
my %store_compute;
my %check_code;
open my $s1,">",$out_file;
print $s1 "tool\tcode\tcoverage\terror_profile\trecall\tfdr\tf1\trecall_avg\trecall_stderr\tfdr_avg\tfdr_stderr\tf1_avg\tf1_stderr\n";
foreach my $tool (sort keys %store){
     foreach my $code (sort keys %{$store{$tool}}){
          $check_code{$code}=1;
          my @t=split("_",$code);
          chop($store{$tool}{$code}{"recall"});
          my @tab_recall=split(";",$store{$tool}{$code}{"recall"});
          my $avg_recall=&average(\@tab_recall);
          my $median_recall=&median(\@tab_recall);
          my $stdev_recall=&stdev(\@tab_recall);
          my $stderr_recall=$stdev_recall/$den;
          my @tab_fdr=split(";",$store{$tool}{$code}{"fdr"});
          my $avg_fdr=&average(\@tab_fdr);
          my $median_fdr=&median(\@tab_fdr);
          my $stdev_fdr=&stdev(\@tab_fdr);
          my $stderr_fdr=$stdev_fdr/$den;
          my @tab_f1=split(";",$store{$tool}{$code}{"f1"});
          my $avg_f1=&average(\@tab_f1);
          my $median_f1=&median(\@tab_f1);
          my $stdev_f1=&stdev(\@tab_f1);
          my $stderr_f1=$stdev_f1/$den;
          print $tool."\t".$code."\t".$t[0]."\t".$t[1]."\t".$avg_recall.";".$median_recall.";".$stdev_recall.
          "\t".$avg_fdr.";".$median_fdr.";".$stdev_fdr."\t".$avg_f1.";".$median_f1.";".$stdev_f1."\n";
          print $s1 $tool."\t".$code."\t".$t[0]."\t".$t[1]."\t".$avg_recall.";".$median_recall.";".$stdev_recall.
          "\t".$avg_fdr.";".$median_fdr.";".$stdev_fdr."\t".$avg_f1.";".$median_f1.";".$stdev_f1."\t".
          $avg_recall."\t".$stderr_recall."\t".$avg_fdr."\t".$stderr_fdr."\t".$avg_f1."\t".$stderr_f1."\n";
          $store_compute{$tool}{$code}{"recall"}=$avg_recall;
          $store_compute{$tool}{$code}{"fdr"}=$avg_fdr;
          $store_compute{$tool}{$code}{"f1"}=$avg_f1;
     }
}
close $s1;

print "##### WIDE #######\n";
my @list_tools=sort keys %store;
my @list_codes=sort keys %check_code;
print "## Recall ##\n";
print "Recall\t".join("\t",@list_tools)."\n";
foreach my $code (@list_codes){
     my $line=$code;
     foreach my $tool (@list_tools){
          $line.="\t".sprintf("%.02f",$store_compute{$tool}{$code}{"recall"}*100)."%";
     }
     print $line."\n";
}
print "## FDR ##\n";
print "FDR\t".join("\t",@list_tools)."\n";
foreach my $code (@list_codes){
     my $line=$code;
     foreach my $tool (@list_tools){
          $line.="\t".sprintf("%.02f",$store_compute{$tool}{$code}{"fdr"}*100)."%";
     }
     print $line."\n";
}

print "## F1 ##\n";
print "F1\t".join("\t",@list_tools)."\n";
foreach my $code (@list_codes){
     my $line=$code;
     foreach my $tool (@list_tools){
          $line.="\t".sprintf("%.02f",$store_compute{$tool}{$code}{"f1"});
     }
     print $line."\n";
}


sub median{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $med="NA";
        @$data=sort {$a <=> $b} (@$data);
	if (scalar(@$data) % 2 ==0 ){
		$med=(@{$data}[scalar(@$data)/2]+@{$data}[scalar(@$data)/2-1])/2;
	}
	else{
		$med=@{$data}[int(scalar(@$data)/2)];
	}
        return $med;
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
