package sci_commands::atac_classify;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("atac_classify");

sub atac_classify {

require "sci_commands/atac_deviation.pm";
import "sci_commands::atac_deviation", "atac_deviation";

@ARGV = @_;
$args = join("\t", @ARGV);

# Defaults:
$permCT = 100;
$binCT = 100;
$minTFCT = 10;
$TSS_flanking = 20000;
$corr_cutoff = 0.15;

getopts("O:Xb:P:B:F:G:C:c:T:t", \%opt);

$die2 = "
scitools atac-classify [options, -G required] [counts matrix] [specification file]
   or    classify-cells
         cell-classifier
   
Will perform atac-deviation on TF sets and gene sets and utilize up and down weights.
It is currently designed to utilize marker genes / TFs that should be on or off, but
future implementations will incorporate weights.

Specification file has '>' lines that have the cell type name (no whitespace),
subsequent lines have three tab-delimited columns:
   [gene/g or motif/m]   [gene or motif name]   [up/u or down/dn/d]
Motifs or gene names must match the names in the -G and -T files.

Required Options:
   -G   [STR]   Gene info (refGene.txt formats), Shortcut eg: hg38, hg19, mm10

Other Options:
   -O   [STR]   Output prefix directory (def = matrix_prefix.specification_prefix.classify)
   -T   [STR]   TF binding site bed
   -P   [INT]   Permutation count (def = $permCT)
   -B   [INT]   Bin count (def = $binCT)
   -F   [INT]   Min peaks for a feature to include it (def = $minTFCT)
   -C   [STR]   Cicero all-cons file. Will also use all linked sites for genes.
   -c   [FLT]   Cicero correlation cutoff (0-1, def = $corr_cutoff)
   -t           Also use cicero linked peaks for TF deviation (def = just motif sites)
   -S   [INT]   Flanking size (out from TSS in bp, def = $TSS_flanking)
   -b   [STR]   Bedtools call (def = $bedtools)
   -X           Retain intermediate files (def = remove)

";

if (!defined $ARGV[1] || !defined $opt{'G'}) {die $die2};
if (!defined $opt{'O'}) {
	$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix//;
	$opt{'O'} .= ".".$ARGV[1]; $opt{'O'} =~ s/\.txt//;
} else {$opt{'O'} =~ s/\.classify//};
if (defined $opt{'F'}) {$minTFCT = $opt{'F'}};
if (defined $opt{'P'}) {$permCT = $opt{'P'}};
if (defined $opt{'B'}) {$binCT = $opt{'B'}};
if (defined $opt{'c'}) {$corr_cutoff = $opt{'c'}};
if (defined $opt{'b'}) {$bedtools = $opt{'b'}};

system("mkdir $opt{'O'}.classify");

open LOG, ">$opt{'O'}.classify/cell_classifier.log";
$ts = localtime(time);
print LOG "$ts\tatac-classify called\n\t\t\t\t$args\n\t\t\t\tReading in specification file ... ";

open IN, "$ARGV[1]";
while ($l = <IN>) {
	chomp $l;
	if ($l =~ /^>/) {
		$classID = $l;
		$classID =~ s/^>//;
		$classID =~ s/\s.+$//;
		$CLASS_gene_up_list{$classID} = "";
		$CLASS_gene_dn_list{$classID} = "";
		$CLASS_gene_upCT{$classID} = 0;
		$CLASS_gene_dnCT{$classID} = 0;
		$classCT++;
	} else {
		($type,$name,$dir) = split(/\t/, $l);
		if ($type =~ /^g/i) { # gene
			$CLASS_geneCT{$classID}++;
			if ($dir =~ /^u/i) { # up
				$CLASS_GENE_UP{$classID}{$name} = 1;
				$CLASS_gene_up_list{$classID} .= "$name\n";
				$CLASS_gene_upCT{$classID}++;
			} elsif ($dir =~ /^d/i) { # down
				$CLASS_GENE_DN{$classID}{$name} = 1;
				$CLASS_gene_dn_list{$classID} .= "$name\n";
				$CLASS_gene_dnCT{$classID}++;
			}
			if (!defined $GENE_include{$name}) {
				$GENE_include{$name} = 1;
				$geneCT++;
			}
		} elsif ($type =~ /^m/i) { # motif
			if (!defined $opt{'T'}) {die "ERROR: Motifs are specified but no motif bed file provided (-T)\n";
			$CLASS_tfCT{$classID}++;
			if ($dir =~ /^u/i) { # up
				$CLASS_TF_UP{$classID}{$name} = 1;
			} elsif ($dir =~ /^d/i) { # down
				$CLASS_TF_DN{$classID}{$name} = 1;
			}
			if (!defined $TF_include{$name}) {
				$TF_include{$name} = 1;
				$tfCT++;
			}
		}
	}
} close IN;

$ts = localtime(time);
print LOG "$classCT total classes.\n";

if ($tfCT>0) {
	print LOG "$ts\tBuilding TF set bed file for $tfCT motifs.\n";

	open IN, "$opt{'T'}";
	open OUT, ">$opt{'O'}.classify/motif_set.bed";
	
	while ($l = <IN>) {
		chomp $l;
		($chr,$start,$end,$name) = split(/\t/, $l);
		if (defined $TF_include{$name}) {
			print OUT "$l\n";
			$TF_include{$name} = 2;
		}
	} close IN; close OUT;
	
	$missingTFct = 0; $foundTFct = 0; $missing_TFs = "";
	foreach $name (keys %TF_include) {
		if ($TF_include{$name} > 1) {
			$foundTFct++;
		} else {
			$missing_TFs .= "$name,";
			$missingTFct++;
		}
	}
	
	if ($missingTFct>0) {
		$missing_TFs =~ s/,$//;
		print LOG "\t\t\t\tWarning: $missingTFct missing motifs ($missing_TFs)\n";
	}
	
	if ($foundTFct>0) {
		print LOG "\t\t\t\tRunning deviation analysis for $foundTFct motifs (of $tfCT specified)\n";
		$dev_opts = "-O $opt{'O'}.classify/motif_set -b $bedtools -P $permCT -B $binCT -F $minTFCT";
		if (defined $opt{'t'} && defined $opt{'C'}) {
			$dev_opts .= " -C $opt{'C'} -c $corr_cutoff";
		}
		if (defined $opt{'X'}) {$dev_opts .= " -X"};
		print LOG "\t\t\t\tatac-deviation $dev_opts $ARGV[0] $opt{'O'}.classify/motif_set.bed\n";
		atac_deviation("$dev_opts $ARGV[0] $opt{'O'}.classify/motif_set.bed");
	} else {
		print LOG "\t\t\t\tNo specified motifs were found! - skipping motif portion of the analysis.\n";
	}
}

$ts = localtime(time);
print LOG "$ts\tBuilding gene sets.\n";

open OUT, ">$opt{'O'}.classify/gene_set.txt";
$gene_setCT = 0;
foreach $classID (keys %CLASS_gene_up_list) {
	if ($CLASS_gene_upCT{$classID}>0) {
		print OUT "\>$classID\_up\n$CLASS_gene_up_list{$classID}";
		$gene_setCT++;
	}
	if ($CLASS_gene_dnCT{$classID}>0) {
		print OUT "\>$classID\_dn\n$CLASS_gene_dn_list{$classID}";
		$gene_setCT++;
	}
} close OUT;

if ($gene_setCT<1) {
	print LOG "\t\t\t\tERROR: No gene sets were established. Exiting!\n";
	die "ERROR: No gene sets were established. Exiting!\n";
}

print LOG "\t\t\t\tRunning deviation analysis for $gene_setCT gene sets.\n";
$dev_opts = "-O $opt{'O'}.classify/motif_set -b $bedtools -P $permCT -B $binCT -F $minTFCT -G $opt{'G'} -S $TSS_flanking";
if (defined $opt{'C'}) {
	$dev_opts .= " -C $opt{'C'} -c $corr_cutoff";
}
if (defined $opt{'X'}) {$dev_opts .= " -X"};
print LOG "\t\t\t\tatac-deviation $dev_opts $ARGV[0] $opt{'O'}.classify/gene_set.txt\n";
atac_deviation("$dev_opts $ARGV[0] $opt{'O'}.classify/gene_set.txt");






}