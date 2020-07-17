package sci_commands::annot_compare;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("annot_compare");

sub annot_compare {

@ARGV = @_;
getopts("O:", \%opt);

$die2 = "
scitools annot-compare [options] [annot1] [annot2]
   or    compare-annot

Will produce a table of the number of cells present
in each annotaiton combination. Only cells in both
annot files will be included.

Rows are annot1, columns are annot2.

Options:
   -O   [STR]   Output prefix (def = annot1.annot2)
                Adds .compare.txt suffix

";

if (!defined $ARGV[1]) {die $die2};

if (!defined $opt{'O'}) {
	$pfx1 = $ARGV[0]; $pfx1 =~ s/\.annot$//;
	$pfx2 = $ARGV[1]; $pfx2 =~ s/\.annot$//;
	$opt{'O'} = "$pfx1.$pfx2";
}
$opt{'O'} =~ s/\.txt$//; $opt{'O'} =~ s/\.compare$//;

$out_header = "#"; @ANNOT1 = ();
open IN, "$ARGV[0]";
while ($l = <IN>) {
	chomp $l;
	($cellID,$annot) = split(/\t/, $l);
	$CELLID_annot1{$cellID} = $annot;
	if (!defined $ANNOT1_check{$annot}) {
		$out_header .= "\t$annot";
		push @ANNOT1, $annot;
		$ANNOT1_check{$annot} = 1;
	}
} close IN;

$cell_match_ct = 0; @ANNOT2 = ();
open IN, "$ARGV[1]";
while ($l = <IN>) {
	chomp $l;
	($cellID,$annot2) = split(/\t/, $l);
	if (!defined $ANNOT2_check{$annot2}) {
		push @ANNOT2, $annot2;
		$ANNOT2_check{$annot2} = 1;
	}
	if (defined $CELLID_annot1{$cellID}) {
		$annot1 = $CELLID_annot1{$cellID};
		$cell_match_ct++;
		$ANNOT1_ANNOT2_ct{$annot1}{$annot2}++;
	}
} close IN;

open OUT, ">$opt{'O'}.compare.txt";
print OUT "$out_header\n";
foreach $annot2 (@ANNOT2) {
	print OUT "$annot2";
	for ($i = 0; $i < @ANNOT1; $i++) {
		if (!defined $ANNOT1_ANNOT2_ct{$ANNOT1[$i]}{$annot2}) {$ANNOT1_ANNOT2_ct{$ANNOT1[$i]}{$annot2}=0};
		print OUT "\t$ANNOT1_ANNOT2_ct{$ANNOT1[$i]}{$annot2}";
	}
	print OUT "\n";
}
close OUT;

exit;

}
1;