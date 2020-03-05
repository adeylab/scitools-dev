package sci_commands::matrix_subsample;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_subsample");

sub matrix_subsample {

@ARGV = @_;
getopts("O:r:L:W:", \%opt);

$sample_frac = 0.1;

$die2 = "
scitools matrix-subsample [options] [input sparse matrix]

Note: only sparse matrix files are currently accepted.

Will downsample the number of cells to be profiled.

Options:
   -O   [STR]  Output prefix (default is [input].cistopic.dims)
   -L   [STR]  Columns file for sparseMatrix (will try to auto-detect)
   -W   [STR]  Rows file for sparseMatrix (will try to auto-detect)
   -r   [FLT]  Approximate fraction of cells to retain (def = $sample_frac)

";

if (!defined $ARGV[0]) {die $die2};
if (defined $opt{'r'}) {$sample_frac = $opt{'r'}};

$prefix = $ARGV[0];
$prefix =~ s/\.gz$//;
$prefix =~ s/\.matrix$//;
$prefix =~ s/\.values$//;
$prefix =~ s/\.sparseMatrix$//;

if (!defined $opt{'O'}) {
   $opt{'O'} = $prefix;
}

if (defined $opt{'L'}) {
      $col_file = $opt{'L'};
} else {
	if (-e "$prefix.sparseMatrix.cols") {
		$col_file = "$prefix.sparseMatrix.cols";
	} elsif (-e "$prefix.sparseMatrix.cols.gz") {
		$col_file = "$prefix.sparseMatrix.cols.gz";
	} else {
		die "ERROR: Cannot detect cols file (e.g. $prefix.sparseMatrix.cols), please provide as -C\n";
	}
}

if (defined $opt{'W'}) {
   $row_file = $opt{'W'};
} else {
   if (-e "$prefix.sparseMatrix.rows") {
      $row_file = "$prefix.sparseMatrix.rows";
   } elsif (-e "$prefix.sparseMatrix.rows.gz") {
      $row_file = "$prefix.sparseMatrix.rows.gz";
   } else {
      die "ERROR: Cannot detect rows file (e.g. $prefix.sparseMatrix.rows), please provide as -C\n";
   }
}

$prefix .= ".subsampled";

open COL_OUT, "| gzip > $prefix.sparseMatrix.cols.gz";

$cellNum = 0; $newNum = 0;
if ($col_file =~ /\.gz$/) {
	open COLS, "zcat $col_file |";
} else {
	open COLS, "$col_file";
}
while ($cellID = <COLS>) {
	$cellNum++;
	$check = rand(1);
	if ($check <= $sample_frac) {
		$newNum++;
		$CELLNUM_newNum{$cellNum} = $newNum;
		print COL_OUT "$cellID";
	}
} close COLS; close COL_OUT;

open VALS_OUT, "| gzip > $prefix.sparseMatrix.values.gz";

if ($ARGV[0] =~ /\.gz$/) {
	open VALS, "zcat $ARGV[0] |";
} else {
	open VALS, "$ARGV[0]";
}

while ($l = <VALS>) {
	chomp $l;
	($siteNum,$cellNum,$val) = split(/\s/, $l);
	if (defined $CELLNUM_newNum{$cellNum}) {
		print VALS_OUT "$siteNum\t$CELLNUM_newNum{$cellNum}\t$val\n";
	}
} close VALS; close VALS_OUT;

if ($row_file =~ /\.gz$/) {
	system("cp $row_file $prefix.sparseMatrix.rows.gz");
} else {
	system("cp $row_file $prefix.sparseMatrix.rows");
}

}
1;
