#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h=0;
my $tag_md10=0;
GetOptions ('help' => \$h, 'h' => \$h, 'm'=>\$tag_md10); # , 'i=s'=>\$in_file, 'r=s'=>\$out_dir_root);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
	print "# Script to calculate alpha and beta diversity per array
# Arguments :
# -m: to work with the spacer sets including less than 10 spacers
";
	die "\n";
}


my $in_file="spacer_clusters_and_metadata_for_alphadiv.tsv";
my $out_file="spacer_sets_alphadiv_up10.tsv";

if ($tag_md10==1){
    $in_file="spacer_clusters_and_metadata_for_alphadiv_md10.tsv";
    $out_file="spacer_sets_alphadiv_md10.tsv";
}

my $last;
my %seen;
my %tmp;
my $i=0;
open my $s1,">",$out_file;
print $s1 "array\tsample\tn_spacers\ttotal_cover\tmax_cover\tmedian_cover\tn_singleton\tshannon_diversity\tshannon_evenness\tsimpson\tsimpson_diversity\tn_spacer_min80\tn_spacer_min50\tn_spacer_max5\n";
open my $tsv,"<",$in_file;
while(<$tsv>){
    chomp($_);
    my @tab=split("\t",$_);
    if ($tab[0] eq "cluster_id"){next;}
    my $code=$tab[2].",".$tab[3];
    if ($last ne $code){
        if (defined($seen{$code})){
            print "Pblm, we already saw $tab[1] ??\n";
            die("\n");
        }
        if ($last ne ""){
            &process(\%tmp,$last,$s1);
            $i++;
            if ($i % 100 == 0){
                print ".. $i -- $last ..\n";
            }
            # if ($i>=1000){last;}
        }
        %tmp=();
        $last=$code;
        $seen{$code}=1;
    }
    if (defined($tmp{$tab[0]})){
        print "pblm, already saw $tab[0] for $code ?\n";
        die("\n");
    }
    $tmp{"by_spacer"}{$tab[0]}=$tab[1];
    $tmp{"total_cover"}+=$tab[1];
}
close $tsv;
if ($last ne ""){
    &process(\%tmp,$last,$s1);
    # last;
}
close $s1;

sub process(){
    my $hash=$_[0];
    my @t=split(",",$_[1]);
    # print "## Processing $t[0] -- $t[1]\n";
    my $s1=$_[2];
    ## Alpha diversity
    my %stats;
    my @ordered_list=sort {$$hash{"by_spacer"}{$b} <=> $$hash{"by_spacer"}{$a} or $a cmp $b} keys %{$$hash{"by_spacer"}};
    my $total=scalar(@ordered_list);
    my $max_cover=$$hash{"by_spacer"}{$ordered_list[0]};
    my $th_eighty=0.8*$max_cover;
    my $th_fifty=0.5*$max_cover;
    my $th_one=0.05*$max_cover;
    my @tab_cover=();
    my $n_singleton=0;
    my $n_eighty=0;
    my $n_fifty=0;
    my $n_one=0;
    foreach my $spacer (@ordered_list){
        my $p=$$hash{"by_spacer"}{$spacer}/$$hash{"total_cover"};
        # print $spacer."\t".$$hash{"by_spacer"}{$spacer}."\t".$p."\n";
        $stats{"shannon"}+=$p*log($p);
        $stats{"simpson"}+=$p*$p;
        $stats{"simpson_div"}+=$$hash{"by_spacer"}{$spacer}*($$hash{"by_spacer"}{$spacer}-1);
        $stats{"total"}++;
        if ($$hash{"by_spacer"}{$spacer}==1){$n_singleton++;}
        if ($$hash{"by_spacer"}{$spacer}>=$th_eighty){$n_eighty++;}
        if ($$hash{"by_spacer"}{$spacer}>=$th_fifty){$n_fifty++;}
        if ($$hash{"by_spacer"}{$spacer}<=$th_one){$n_one++;}
        push(@tab_cover,$$hash{"by_spacer"}{$spacer});
    }
    my $shannon_index=-1*$stats{"shannon"};
    my $shannon_even=$shannon_index/log($total);
    my $simpson_index=1/$stats{"simpson"};
    my $simpson_diversity=1-($stats{"simpson_div"}/($$hash{"total_cover"}*($$hash{"total_cover"}-1)));
    my $median=&median(\@tab_cover);
    # print $s1 "array\tsample\tn_spacers\ttotal_cover\tmax_cover\tmedian_cover\tn_singleton\tshannon\tsimpson\tn_spacer_min50\tn_spacer_max5\n";
    print $s1 $t[0]."\t".$t[1]."\t".$total."\t".$$hash{"total_cover"}."\t".$max_cover."\t".$median."\t".$n_singleton."\t".
    $shannon_index."\t".$shannon_even."\t".$simpson_index."\t".$simpson_diversity."\t".$n_eighty."\t".$n_fifty."\t".$n_one."\n";
    if ($median>$max_cover){
        die("pblm $median vs $max_cover\n");
    }
}


sub median(){
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
