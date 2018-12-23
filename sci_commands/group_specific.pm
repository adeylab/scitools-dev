package sci_commands::group_specific;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("group_specific");

sub group_specific {

@ARGV = @_;

getopts("b:", \%opt);

$die2 = "

scitools group-specific [options] [output prefix] [bed1] [bed2] ... (bedN)

Will take speaks called for each group and intersect to pull those specific to each group.

Options:
   -b   [STR]   bedtools call (def = $bedtools)
   -X           retain intermediate files (def = no)

";

if (!defined $ARGV[2]) {die $die2};
if (defined $opt{'b'}) {$bedtools = $opt{'b'}};

for ($i = 1; $i < @ARGV; $i++) {
	$files = "";
	for ($j = 1; $j < @ARGV; $j++) {
		if ($i != $j) {
			$files .= "$ARGV[$j] ";
		}
	}
	$name = $ARGV[$i]; $name =~ s/\.bed$//;
	system("cat $files | $bedtools sort -i - | $bedtools merge -i - > $ARGV[0].$i.tmp.bed");
	system("$bedtools intersect -v -b $ARGV[0].$i.tmp.bed -a $ARGV[$i] > $ARGV[0].$name.specific.bed");
}


}
1;
