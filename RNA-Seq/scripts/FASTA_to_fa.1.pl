#!/usr/bin/perl -w

die ("Argument error!") if (@ARGV!=2);
open (R_1,"< $ARGV[0]") || die ("Cannot open $ARGV[0]: $!\n");
open (W_1,"> $ARGV[1]");

select W_1;
$n=0;

while ($_=<R_1>) {
 if ($n==0) {
  print $_;
 }
 else {
  if ($_=~/^>/) {
   print "\n$_";
  }
  else {
   chomp ($_);
   print $_;
  }
 }
 $n=1;
}

print "\n";

close R_1;
close W_1;

