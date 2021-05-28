package sci_commands::bam_coverage;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("bam_coverage");

sub bam_coverage {

    @ARGV = @_;

# Defaults:
    $bedtools = "bedtools";

    getopts("O:b:", \%opt);

$die2 = "
scitools bam-coverage [options] [duplicate removed and filtered bam file] [bed file]

Calculates per-cell coverage of the target bed.

Options:
   -O   [STR]   Output prefix (default is bam file prefix)
   -b   [STR]   Bedtools call (def = $bedtools)

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $ARGV[1]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.bam$//};
if (defined $opt{'b'}) {$bedtools = $opt{'b'}};

$total_bases = 0;
open IN, "$ARGV[1]";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$total_bases += ($P[2]-$P[1]);
} close IN;

print STDERR "\nTarget Size = $total_bases bp.\n";

open IN, "$bedtools intersect -bed -wo -abam $ARGV[0] -b $ARGV[1] |";
#chr10   8055246 8055298 GACTAGCAGCCAAGGCGGCATTCT:599818#0/2     60      +       8055246 8055298 0,0,0   1       52,     0,      chr10   8055286 8055896 12
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$cellID = $P[3]; $cellID =~ s/:.+$//;
	if ($P[1]<=$P[13]) {
		$start = $P[13];
	} else{
		$start = $P[1];
	}
	for ($i = 0; $i < $P[15]; $i++) {
		$pos = $i+$start;
		if (!defined $CELLID_CHROM_POS_cov{$cellID}{$P[0]}{$pos}) {
			$CELLID_CHROM_POS_cov{$cellID}{$P[0]}{$pos}++;
			$CELLID_bases_covered{$cellID}++;
		}
	}
} close IN;

open OUT, ">$opt{'O'}.frac_covered.txt";
%CELLID_CHROM_POS_cov = ();
foreach $cellID (keys %CELLID_bases_covered) {
	$frac = sprintf("%.3f", $CELLID_bases_covered{$cellID}/$total_bases);
	print OUT "$cellID\t$CELLID_bases_covered{$cellID}\t$frac\n";
} close OUT;
	
}

1;
