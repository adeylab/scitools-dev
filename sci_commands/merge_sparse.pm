package sci_commands::merge_sparse;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("merge_sparse");

sub merge_sparse {

@ARGV = @_;

$die2 = "
scitools merge-sparse [Output_prefix] [sample ID],[columns file],[rows file],[values file] [4 fields for 2nd matrix] (etc...)

Will merge sparse count matrixes that share the same peak set.
Requires an output prefix and the four comma separated items for each matrix.

";

if (!defined $ARGV[2]) {die $die2};
$out = $ARGV[0];

# assign the input matrixes

for ($matrixID = 1; $matrixID < @ARGV; $matrixID++) {
	($ID,$cols,$rows,$vals) = split(/,/, $ARGV[$matrixID]);
	$COLS{$ID} = $cols;
	$ROWS{$ID} = $rows;
	$VALS{$ID} = $vals;
}

# 1) get the unified rows file

%ALL_ROWS = ();
foreach $ID (keys %ROWS) {
	if ($ROWS{$ID} =~ /\.gz/) {
		open IN, "zcat $ROWS{$ID} |";
	} else {
		open IN, "$ROWS{$ID}";
	}
	$rowID = 1;
	while ($l = <IN>) {
		chomp $l;
		if (!defined $ALL_ROWS{$l}) {
			$ALL_ROWS{$l} = 1;
		}
		$ROW_ID_rowID{$ID}{$rowID} = $l;
		$rowID++;
	} close IN;
}

$rowID = 1;
open OUT, "| gzip > $out.rows.gz";
foreach $rowName (sort keys %ALL_ROWS) {
	print OUT "$rowName\n";
	$ROWname_ID{$rowName} = $rowID;
	$rowID++;
} close OUT;

# 2) make column file

open OUT, "gzip > $out.cols.gz";
$global_colID = 1;
foreach $ID (keys %COLS) {
	if ($COLS{$ID} =~ /\.gz/) {
		open IN, "zcat $COLS{$ID} |";
	} else {
		open IN, "$COLS{$ID}";
	}
	$colID = 1;
	while ($l = <IN>) {
		chomp $l;
		$COL_ID_colID_globalColID{$ID}{$colID} = $global_colID;
		$colID++;
		$global_colID++;
		print OUT "$ID\_$l\n";
	} close IN;
}
close OUT;

# 3) values file

open OUT, "| gzip > $out.values.gz";
foreach $ID (keys %VALS) {
	if ($VALS{$ID} =~ /\.gz/) {
		open IN, "zcat $VALS{$ID} |";
	} else {
		open IN, "$VALS{$ID}";
	}
	while ($l = <IN>) {
		chomp $l;
		($rowID,$colID,$value) = split(/\t/, $l);
		print OUT "$ROWname_ID{$ROW_ID_rowID{$ID}{$rowID}}\t$COL_ID_colID_globalColID{$ID}{$colID}\t$value\n";
	} close IN;
}
close OUT;


}
1;
