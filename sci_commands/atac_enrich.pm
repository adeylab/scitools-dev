package sci_commands::atac_enrich;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("atac_enrich");

sub atac_enrich {

@ARGV = @_;
getopts("O:b:", \%opt);

$die2 = "
scitools atac-enrich [options] [peaks.bed] [peak_sets.bed] [feature_sets.bed]
   or    atac-enrichment

This tools will perform a test for the enrichment of features (eg. motifs) within
a set, or sets of peaks. For example, peak_sets.bed may be cicero-linked CCANs.

Each bed must have 4 columns: chr, start, end, annotation; and the peak_sets
must be peaks that exactly match up with peaks in the peaks.bed file.

Output is columns = features and rows = peak sets, p-val for each.

Options:
   -O   [STR]   Output prefix (default is [peak_sets].[feature_sets])
   -b   [STR]   Bedtools call (def = $bedtools)
   
";

if (!defined $ARGV[2]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[1]; $opt{'O'} =~ s/\.bed$//; $opt{'O'} .= ".".$ARGV[2]; $opt{'O'} =~ s/\.bed$//};

# count peaks
open IN, "$ARGV[0]";
while ($l = <IN>) {$peakCT++};
close IN;

# annotate all peaks for features
open IN, "$bedtools intersect -a $ARGV[0] -b $ARGV[2] -wa -wb |";
while ($l = <IN>) {
	chomp $l;
	($peak_chr,$peak_start,$peak_end,$peakID,$feature_chr,$feature_start,$feature_end,$featureID) = split(/\t/, $l);
	$peakID = $peak_chr."_".$peak_start."_".$peak_end; # ensure that peakID format is this
	push @{$PEAKID_features{$peakID}}, $featureID;
	$FEATUREID_peak_count{$featureID}++;
} close IN;

# load in peak sets
open IN, "$ARGV[1]";
while ($l = <IN>) {
	chomp $l;
	($peak_chr,$peak_start,$peak_end,$setID) = split(/\t/, $l);
	$peakID = $peak_chr."_".$peak_start."_".$peak_end;
	$SETID_peakCT{$setID}++;
	push @{$SETID_peakList{$setID}}, $peakID;
} close IN;

# precompute logfact
build_logfact($peakCT);

# perform all by all
open OUT, ">$opt{'O'}.pvals.matrix";
$header = "";
foreach $featureID (sort keys %FEATUREID_peak_count) {
	$header .= "$featureID\t";
} $header =~ s/\t$//;
print OUT "$header\n";

foreach $setID (%SETID_peakCT) {
	print OUT "$setID";
	foreach $featureID (sort keys %FEATUREID_peak_count) {
		$peaks_with_feature = $FEATUREID_peak_count{$featureID};
		$peaks_without_feature = $peakCT - $FEATUREID_peak_count{$featureID};
		$set_peaks_with_feature = 0;
		foreach $peakID (@{$SETID_peakList{$setID}}) {
			foreach $check_feature (@{$PEAKID_features{$peakID}}) {
				if ($check_feature eq $featureID) {$set_peaks_with_feature++};
			}
		}
		$set_peaks_total = $SETID_peakCT{$setID};
		if ($set_peaks_with_feature > 0) {
			$hypergeom_cumulative_pval = hypergeom($peaks_with_feature,$peaks_without_feature,$set_peaks_total,$set_peaks_with_feature);
		} else {
			$hypergeom_cumulative_pval = 1;
		}
		print "\t$hypergeom_cumulative_pval";
		push @ALL_PVALS, $hypergeom_cumulative_pval;
	}
	print OUT "\n";
}

close OUT;

}

sub hypergeom {
    ($n, $m, $N, $i) = @_;
	
	$max_possible = ($n + $m + $N + $i);
	if (!defined $logfact[ $max_possible ]) {add_logfact($max_possible)};
	
    $loghyp1 = $logfact[ $m ] + $logfact[ $n ]
             + $logfact[ $N ] + $logfact[ $m + $n - $N ];
    $loghyp2 = $logfact[ $i ] + $logfact[ $n - $i ] 
             + $logfact[ $m + $i - $N ] + $logfact[ $N - $i ] 
             + $logfact[ $m + $n ];

    return exp($loghyp1 - $loghyp2);
}

sub build_logfact {
	$log_max = $_[0];
	$value = 0;
	for ($current = 1; $current <= $log_max; $current++) {
		$value += log($current);
		$logfact[$current] = $value;
	}
}

sub add_logfact {
	$new_max = $_[0];
	$value = $logfact[$log_max];
	for ($current = $log_max+1; $current <= $new_max; $current++) {
		$value += log($current);
		$logfact[$current] = $value;
	}
	$log_max = $new_max;
}

1;