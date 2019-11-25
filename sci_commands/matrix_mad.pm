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
$h = <IN>; @CELLIDS = split(/\t/, $h);
while ($l = <IN>) {
	chomp $l;
	@P = split (/\t/, $l);
	$rowID = shift(@P);
	for ($i = 0; $i < @P; $i++) {
		$CELLID_sum{$CELLIDS[$i]} += $P[$i];
	} $rowCT++;
} close IN;

foreach $cellID (keys %CELLID_sum) {
	$CELLID_mean{$cellID} = $CELLID_sum{$cellID}/$rowCT;
}

if ($ARGV[0] =~ /\.gz/) {
	open IN, "zcat $ARGV[0] |";
} else {
	open IN, "$ARGV[0]";
}
$h = <IN>;
while ($l = <IN>) {
	chomp $l;
	@P = split (/\t/, $l);
	$rowID = shift(@P);
	for ($i = 0; $i < @P; $i++) {
		push @{$CELLID_DEVS{$CELLIDS[$i]}}, abs($CELLID_mean{$CELLIDS[$i]} - $P[$i])
	}
} close IN;

open OUT, ">$opt{'O'}.MAD.values";
foreach $cellID (keys %CELLID_DEVS) {
	$MAD = $CELLID_DEVS{$cellID}[($rowCT-1)];
	print OUT "$cellID\t$MAD\n";
}

}
1;
