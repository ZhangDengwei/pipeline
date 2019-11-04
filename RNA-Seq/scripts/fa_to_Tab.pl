#!/usr/bin/perl -w

die ("Argument error!") if (@ARGV!=2);
open (R_1,"< $ARGV[0]") || die ("Cannot open $ARGV[0]: $!\n");
open (W_1,"> $ARGV[1]");

select W_1;

while ($_=<R_1>) {
 if ($_=~/^>(.+?)\n/) {
  print "$1\t";
 }
 else{
  print $_;
 }
}

close R_1;
close W_1;

