package commands::matrix_filter;

use lib "/home/adey/GitHub/scitools-dev"; #LIB#
use commands::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_filter");

sub matrix_filter {

@ARGV = @_;

# Defaults
$colMin = 1;
$rowMin = 1;

getopts("O:C:R:c:r:A:a:", \%opt);

$die2 = "
scitools matrix-filter [options] [counts matrix]
   or    filter-matrix
   
Note: The current version is not memory-efficient and is set to be updated.

Options:
   -O   [STR]   Output prefix (default is [input].filt_[colMin]_[rowMin].matrix)
   -C   [INT]   Number of nonZero sites per column to retain (def = $colMin)
   -R   [INT]   Number of nonZero sites per row to retain (def = $rowMin)
   -c   [STR]   List of CellIDs to retain
   -r   [STR]   List of RowIDs to retain
   -A   [STR]   Annotation file
   -a   [STR]   Comma separated list of annotations to include (requires -A)

Note: -C and -R filters are applied after all other filtering.

";

if (!defined $ARGV[0]) {die $die2};
if (defined $opt{'C'}) {$colMin = $opt{'C'}};
if (defined $opt{'R'}) {$rowMin = $opt{'R'}};
if (!defined $opt{'O'}) {
	$opt{'O'} = $ARGV[0];
	$opt{'O'} =~ s/\.matrix$//;
	$opt{'O'} .= ".filt_$colMin\_$rowMin";
}
if (defined $opt{'a'} && !defined $opt{'A'}) {die "\nMust provide an annotaiton file (-A) if specifying annotations to filter (-a)!\n$die2"};

if (defined $opt{'A'}) {read_annot($opt{'A'})};

if (defined $opt{'a'}) {
	@ANNOT_LIST = split(/,/, $opt{'a'});
	foreach $annot (@ANNOT_LIST) {
		$ANNOT_include{$annot} = 1;
	}
}

if (defined $opt{'c'}) {
	open IN, "$opt{'c'}";
	while ($cellID = <IN>) {
		chomp $cellID;
		$CELLID_list_include{$cellID} = 1;
	} close IN;
}

if (defined $opt{'r'}) {
	open IN, "$opt{'r'}";
	while ($rowID = <IN>) {
		chomp $rowID;
		$ROWID_list_include{$rowID} = 1;
	} close IN;
}

read_matrix($ARGV[0]);

foreach $cellID (keys %CELLID_FEATURE_value) {
	$includeCell = 1;
	if (defined $opt{'c'} && !defined $CELLID_list_include{$cellID}) {$includeCell = 0};
	if (defined $opt{'a'} && !defined $ANNOT_include{$CELLID_annot{$cellID}}) {$includeCell = 0};
	if ($includeCell>0) {
		foreach $rowID (keys %{$CELLID_FEATURE_value{$cellID}}) {
			$includeRow = 1;
			if (defined $opt{'r'} && !defined $ROWID_list_include{$rowID}) {$includeRow = 0};
			if ($includeRow>0 && abs($CELLID_FEATURE_value{$cellID}{$rowID})>0) {
				$CELLID_count{$cellID}++;
			}
		}
		if ($CELLID_count{$cellID}>=$colMin) {
			$CELLID_include{$cellID} = 1;
			$included_cells++;
		}
	}
}

foreach $rowID (keys %MATRIX_feature_nonZero) {
	if (!defined $opt{'r'} || defined $ROWID_list_include{$rowID}) {
		$ROWID_count{$rowID} = 0;
		foreach $cellID (keys %CELLID_include) {
			if (abs($CELLID_FEATURE_value{$cellID}{$rowID})>0) {
				$ROWID_count{$rowID}++;
			}
		}
	}
	if ($ROWID_count{$rowID}>=$rowMin) {
		$ROWID_include{$rowID} = 1;
		$included_rows++;
	}
}

open OUT, ">$opt{'O'}.matrix";

$out_header = "";
for ($cellNum = 0; $cellNum < @MATRIX_COLNAMES; $cellNum++) {
	if (defined $CELLID_include{$MATRIX_COLNAMES[$cellNum]}) {
		$out_header .= "$MATRIX_COLNAMES[$cellNum]\t";
	}
} $out_header =~ s/\t$//;
print OUT "$out_header\n";

for ($rowNum = 0; $rowNum < @MATRIX_ROWNAMES; $rowNum++) {
	$rowID = $MATRIX_ROWNAMES[$rowNum];
	if (defined $ROWID_include{$rowID}) {
		print OUT "$rowID";
		for ($cellNum = 0; $cellNum < @MATRIX_COLNAMES; $cellNum++) {
			if (defined $CELLID_include{$MATRIX_COLNAMES[$cellNum]}) {
				print OUT "\t$CELLID_FEATURE_value{$MATRIX_COLNAMES[$cellNum]}{$rowID}";
			}
		} print OUT "\n";
	}
} close OUT;

open LOG, ">$opt{'O'}.log";
$ts = localtime(time);
print LOG "$ts scitools atac-filter
Matrix file = $ARGV[0]
Options:
";
foreach $option (keys %opt) {
	print LOG "   $option   $opt{$option}\n";
}
print LOG "Total cells in input: $matrix_colNum
Total cells retained in output: $included_cells
Total rows in input: $matrix_rowNum
Total rows retained in output: $included_rows\n";
close LOG;

}
1;
