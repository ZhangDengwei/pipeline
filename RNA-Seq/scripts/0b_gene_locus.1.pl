#!/usr/bin/perl -w

die ("Argument error!") if (@ARGV!=2);
open (R_1,"< $ARGV[0]") || die ("Cannot open $ARGV[0]: $!\n");
open (W_1,"> $ARGV[1]");

select W_1;

while ($_=<R_1>) {
 if ($_=~/^(.+?)\t.+?\texon\t(\d+)\t(\d+)\t.+?\t([\+\-])\t.+?gene_id\s+"(.+?)";/) {
  ($chr,$start,$end,$strand,$gene_id)=($1,$2,$3,$4,$5);
  $hash_chr{$gene_id}=$chr;
  $hash_1{$gene_id}.="$start,$end,";
  $hash_strand{$gene_id}=$strand;
 }
}

sub mysort {
 $a<=>$b;
}

while (($key,$value)=each(%hash_1)) {
 chop($value);
 @array_1_unsort=split(/,/,$value);
 @array_1=sort mysort @array_1_unsort;
 ($start,$end)=($array_1[0],$array_1[-1]);
 print "$key\t$hash_chr{$key}:$start\-$end\[$hash_strand{$key}\]\n";
}
close R_1;
close W_1;

