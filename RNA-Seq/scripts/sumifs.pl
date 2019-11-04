#!/usr/bin/perl -w
die ("Argument error!") if (@ARGV!=2);
open (R_1,"< $ARGV[0]") || die ("Cannot open $ARGV[0]: $!\n");
open (W_1,"> $ARGV[1]");

select W_1;

while ($_=<R_1>) {
	unless ($_=~/^Gene ID\t/) {
		chomp($_);
		@array_1=split(/\t/,$_);
		$id=$array_1[0];
		$fpkm=$array_1[7];
		$hash_1{$id}+=$fpkm;
	}
}

open (R_1,"< $ARGV[0]") || die ("Cannot open $ARGV[0]: $!\n");

while ($_=<R_1>) {
	if ($_=~/^(.+?)\t/) {
		$id=$1;
		if ($id ne "Gene ID") {
			unless (exists($hash_2{$id})) {
				print $id."\t".$hash_1{$id}."\n";
				$hash_2{$id}=$id;
			}
		}
	}
}

close R_1;
close W_1;

