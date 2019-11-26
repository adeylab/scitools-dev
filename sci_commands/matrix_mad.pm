package sci_commands::matrix_mad;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_mad");

sub matrix_mad {

@ARGV = @_;

getopts("O:", \%opt);

$die2 = "
scitools matrix-MAD [options] [window counts matrix]
   or    MAD

Calculates the MAD (median absolute deviation) score
for each cell in a matrix.

Currently only supports dense matrix format.

Options:
   -O   [STR]   Output prefix (default is [input].MAD.values)

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};

if ($ARGV[0] =~ /\.gz/) {
	open IN, "zcat $ARGV[0] |";
} else {
	open IN, "$ARGV[0]";
}

# load matrix and calculate the sum of the cell signal
$h = <IN>; chomp $h; @CELLIDS = split(/\t/, $h);
while ($l = <IN>) {
	chomp $l;
	@P = split (/\t/, $l);
	$rowID = shift(@P);
	for ($i = 0; $i < @P; $i++) {
		$CELLID_sum{$CELLIDS[$i]} += $P[$i];
		push @{$CELLID_values{$CELLIDS[$i]}}, $P[$i];
	} $rowCT++;
} close IN;

# calculate differences on normalized values (norm to 1)
foreach $cellID (keys %CELLID_sum) {
	for ($i = 1; $i < $rowCT; $i++) {
		push @{$CELLID_diffs{$cellID}}, (abs( (($CELLID_values{$cellID}[$i]/$CELLID_sum{$cellID})*$rowCT) - (($CELLID_values{$cellID}[$i-1]/$CELLID_sum{$cellID})*$rowCT) ));
	}
}

# print MAD
open OUT, ">$opt{'O'}.MAD.values";
foreach $cellID (keys %CELLID_sum) {
	$MAD = $CELLID_diffs{$cellID}[int($rowCT/2)];
	print OUT "$cellID\t$MAD\n";
} close OUT;

}
1;
