#!/usr/bin/perl -w

die ("Argument error!") if (@ARGV!=3);
open (R_1,"< $ARGV[0]") || die ("Cannot open $ARGV[0]: $!\n");
open (R_2,"< $ARGV[1]") || die ("Cannot open $ARGV[1]: $!\n");
open (W_1,"> $ARGV[2]");

select W_1;

while ($_=<R_1>) {
 if ($_=~/gene_id\s+"(.+?)";\s+transcript_id\s+"(.+?)";/) {
  $hash_1{$2}=$1;
 }
}

while ($_=<R_2>) {
 if ($_=~/^(.+?)\s+gene.+?\t(.+?)\n/) {
  print "$1\t$hash_1{$1}\t$2\n";
 }
}

close R_1;
close W_1;

