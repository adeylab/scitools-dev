#!/usr/bin/perl

$die = "

ARGV0 = index 1 list (i.e. i7 index / I1)
ARGV1 = index 2 list (i.e. i5 index / I2)
ARGV2 = index 3 list (i.e. Tn5 in-line index)
ARGV3 = annot name (for 2nd column of output)

ARGV4 = optional; column number for index sequence
        default = column 1 / single column

";

if (!defined $ARGV[3]) {die $die};
if (!defined $ARGV[4]) {$ARGV[4] = 0} else {$ARGV[4]--};

open IN, "$ARGV[0]";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$I1{$P[$ARGV[4]]} = 1;
} close IN;

open IN, "$ARGV[1]";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$I2{$P[$ARGV[4]]} = 1;
} close IN;

open IN, "$ARGV[2]";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$I3{$P[$ARGV[4]]} = 1;
} close IN;

foreach $i1 (keys %I1) {
	foreach $i2 (keys %I2) {
		foreach $i3 (keys %I3) {
			print "$i1"."$i2"."$i3\t$ARGV[3]\n";
		}
	}
}

exit;