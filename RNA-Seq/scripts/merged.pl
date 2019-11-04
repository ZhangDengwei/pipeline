#!/usr/bin/perl -w

die ("Argument error!") if (@ARGV!=3);
open (R_1,"< $ARGV[0]") || die ("Cannot open $ARGV[0]: $!\n");
open (R_2,"< $ARGV[1]") || die ("Cannot open $ARGV[1]: $!\n");
open (W_1,"> $ARGV[2]");

select W_1;

print "TAIR_ID"."\t";

while ($_=<R_1>) {
	chomp($_);
	$head.=$_."\t";
}

chop($head);

print $head."\n";

while ($_=<R_2>) {
	chomp($_);
	@array_1=split(/\t/,$_);
	print $array_1[0]."\t";
	$num=@array_1;#print $num."\t".$array_1[1]."\n";
	for ($i=1;$i<=$num-1;$i+=2) {
		$tmp.=$array_1[$i]."\t";
	}
	chop($tmp);
	print $tmp."\n";
	$tmp="";
}


close R_1;
close R_2;
close W_1;

