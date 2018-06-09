#!/usr/bin/perl

# VERSION INFO
$version = "0.1.2";
%version_info = (
	"0.1.0" => "180215, alpha - initial development",
	"0.1.1d" => "180418, alpha - dev split",
	"0.1.2" => "180608, modulization"
);

# LOAD MODULES
use Getopt::Std; %opt = ();
use lib "/home/adey/scitools-dev"; #LIB#
use commands::general;

# LOAD DEFAULTS
$SCITOOLS_DEFAULTS = "/home/adey/scitools-dev/scitools.cfg"; #CONFIG#
if (-e "$ENV{'HOME'}/.scitools.cfg") {
	$SCITOOLS_DEFAULTS = "$ENV{'HOME'}/.scitools.cfg";
} elsif (-e "$ENV{'HOME'}/scitools.cfg") {
	$SCITOOLS_DEFAULTS = "$ENV{'HOME'}/scitools.cfg";
} elsif ($ARGV[0] =~ /\.cfg$/) {
	# special config specification as first argument - advanced use only
	$SCITOOLS_DEFAULTS = shift(@ARGV);
	print STDERR "INFO: Special config file detected as first argument. ($ARGV[0])\n";
}
load_defaults($SCITOOLS_DEFAULTS);

# LOAD GRADIENTS
load_gradient_defaults();

# COMMAND DIRECTORY
%COMMANDS = (
	"help" => 0,
	"merge" => 0, "filter" => 0, "split" => 0,
	"dependencies" => 0, "depend" => 0,
	"gradient" => 0, "gradients" => 0,
	
	"fastq-dump" => 0, "dump-fastq" => 0,
	"fastq-dump-new" => 0,
	"fastq-split" => 0, "split-fastq" => 0,
	"fastq-merge" => 0, "merge-fastq" => 0,
	"fastq-align" => 0, "align-fastq" => 0, "align" => 0,
	
	"bam-bulk2sci" => 0, "bulk2sci" => 0,
	"bam-addrg" => 0, "addrg" => 0,
	"bam-rmdup" => 0, "rmdup" => 0,
	"bam-split" => 0, "split-bam" => 0,
	"bam-filter" => 0, "filter-bam" => 0,
	"bam-merge" => 0, "merge-bam" => 0,
	"bam-project" => 0, "project-bam" => 0, "project" => 0,
	"bam-aggregate" => 0, "aggregate-bam" => 0,
	
	"signal-make" => 0, "make-signal" => 0,
	"plot-signal" => 0, "signal-plot" => 0,
	
	"annot-make" => 0, "make-annot" => 0,
	"annot-merge" => 0, "merge-annot" => 0,
	"rename-cells" => 0,
	
	"atac-callpeaks" => 0, "atac-callpeak" => 0, "callpeak" => 0, "callpeaks" => 0,
	"atac-mergepeaks" => 0, "atac-mergepeak" => 0, "mergepeak" => 0, "mergepeaks" => 0,
	"atac-counts" => 0, "atac-count" => 0, "count" => 0, "counts" => 0,
	"atac-chromvar" => 0, "chromvar" => 0,
	"atac-cicero" => 0, "cicero" => 0,
	
	"matrix-filter" => 0, "filter-matrix" => 0,
	"matrix-naomit" => 0, "naomit-matrix" => 0,
	"matrix-summarize" => 0, "summarize-matrix" => 0,
	"matrix-tf" => 0, "tf" => 0,
	"matrix-tfidf" => 0, "tfidf" => 0,
	"matrix-lsi" => 0, "lsi" => 0,
	"matrix-zscore" => 0, "zscore-matrix" => 0,
	"matrix-tsne" => 0, "tsne" => 0,
	"matrix-pca" => 0, "pca" => 0,
	"matrix-nmf" => -1, "nmf" => -1,
	"matrix-bicluster" => -1, "bicluster" => -1,
	"matrix-aggregate" => 0, "aggregate-matrix" => 0,
	"matrix-merge" => 0, "merge-matrix" => 0,
	"matrix-approx-factors" => 0, "matrix-factors" => 0, "factors" => 0,
	"matrix-swne" => 0, "swne" => 0, "piglet" => 0,
	"umap" => 0, "owl" => 0, "matrix-umap" => 0,
	
	"dims-kmeans" => 0, "kmeans" => 0,
	"dims-dbscan" => 0, "dbscan" => 0,
	"dims-pcurve" => 0, "pcurve" => 0,
	"pcurve-center" => 0, "center-pcurve" => 0, "lambda-center" => 0, "center-lambda" => 0,
	"prune-pcurve" => 0, "pcurve-prune" => 0,
	
	"aggregate-cells" => 0, "aggregate" => 0,

	"plot-complexity" => 0,
	"plot-dims" => 0,
	"plot-pcurve" => 0,
	"plot-reads" => 0,
	"plot-factors" => 0,
	
	"index-performance" => 0, "index-perform" => 0,
	"combine-data" => 0, "data-combine" => 0,
	"split-data" => 0, "data-split" => 0

);

# HELP TEXT
$die = "
scitools [command] [options] [arguments]

Version: $version ($version_info{$version})
    adeylab.org & github.com/adeylab

scitools is a set of commands for general processing of single-
cell combinatorial indexing data. It is predominantly a wrapper
for generating R or other scripts that will be stored and executed.

scitools commands are in the form of [class]-[operation], calling
a tool with just the [class] option will give a more detailed
description and a refined set of commands. Many commands can be
called using only the [operation] if it is unique to the class.

Command:               Description:

   help                Display additional scitools descriptions
   dependencies        Check dependencies
   gradient            Print out color gradient pre-sets and details

   fastq-dump          Go from illumina fastqs to SCI fastq format
   fastq-split         Split SCI fastq files using annotation file
   fastq-merge         Merge fastq files that have the same barcodes
   fastq-align         Align fastq files and sort resulting bam file

   bam-bulk2sci        Merge multiple bulk bam files to one SCI bam
   bam-addrg           Add RG lines to bam
   bam-rmdup           Barcode-based duplicate removal
   bam-filter          Filter bam based on a variety of options
   bam-split           Split bam by annotation file
   bam-merge           Merges one or more bam files
   bam-project         Use complexity to project additional sequence
   bam-aggregate       Aggregate cells in bam by annotation file
   
   annot-make          Make annotation file
   annot-merge         Merge annotation files
   rename-cells        Rename cells
   
   atac-callpeak       Call peaks on bam file using macs2
   atac-mergepeak      Merge ATAC-seq peak files
   atac-counts         Bam and peak file to a counts matrix
   atac-chromvar       Run chromVAR wrapper on sci-ATAC-seq data
   atac-cicero         Run cicero wrapper on sci-ATAC-seq data
   
   signal-make         Generate windowed signal over features from bam
   signal-plot         Plot windowed signal views
   
   matrix-summarize    Generate a summary and plots on matrix properties
   matrix-filter       Filter a sci-ATAC-seq counts matrix
   matrix-naomit       Filter out rows with NA values
   matrix-zscore       Z-scores matrix by rows, columns, or globally
   matrix-tf           Normalize only by term frequency
   matrix-tfidf        Perform tf-idf on counts matrix
   matrix-lsi          Perform Latent Semantic Indexing on matrix
   matrix-tsne         tSNE on matrix
   matrix-pca          PCA on matrix
   matrix-factors      Calculates reconstruction error for NMF and SWNE
   matrix-nmf          Non-negative Matrix Factorization of matrix
   matrix-swne 	       SWNE on matrix 
   matrix-umap         UMAP on matrix
   matrix-bicluster    Bicluster and plot a heatmap
   matrix-aggregate    Aggregate cells in counts matrix by annotation
   matrix-merge        Merge matrices, no overlap in cell names assumed
   
   dims-kmeans         Kmeans clustering on dims file
   dims-dbscan         Density-base (dbscan) clustering on dims file
   dims-pcurve         Project a principle curve through dims file
   pcurve-center       Centers and normalizes a pcurve lambda
   pcurve-prune        Prunes cells distant from the pcurve
   aggregate-cells     Aggregate cells in proximity with one another

   plot-complexity     Plot complexity data
   plot-dims           Plot tSNE or other dimensions file
   plot-pcurve         Genrate multiple princurve plots
   plot-factors        Generates factor vs reconstruction error plot
   plot-reads          Plot reads as points in genomic window
   
   index-perform       Index performance on fastq or bam
   combine-data        Combine matrixes, annotations, dims to a table
   split-data          Breaks a combined data file into component files

";

# PULL COMMAND OR DIE
$command = lc(shift(@ARGV));
if (!defined $COMMANDS{$command}) {die $die};

# IN DEVELOPMENT WARNING
if ($COMMANDS{$command}==0) {
	print STDERR "\nWARNING: The command you have specified ($command) is either currently being developed or unverified!\n";
} elsif ($COMMANDS{$command}<0) {
	die "\nERROR: The command you have specified ($command) is very much an aspirational command and has no actual code written for it. Feel free to volunteer.\n";
}

########## GENERAL FUNCTIONS ##########
if ($command eq "dependencies" || $command eq "depend") {
	use commands::dependencies; dependencies(@ARGV); exit;
}

if ($command eq "index-perform" || $command eq "index-performance") {
	use commands::index_performance; index_performance(@ARGV); exit;
}


########## FASTQ FUNCTIONS ##########
if ($command eq "fastq-dump" || $command eq "dump-fastq") {
	use commands::fastq_dump; fastq_dump(@ARGV); exit;	
}

if ($command eq "fastq-split" || $command eq "split-fastq") {
	use commands::fastq_split; fastq_split(@ARGV); exit;	
}

if ($command eq "fastq-merge" || $command eq "merge-fastq") {
	use commands::fastq_merge; fastq_merge(@ARGV); exit;	
}

if ($command eq "fastq-align" || $command eq "align-fastq" || $command eq "align") {
	use commands::fastq_align; fastq_align(@ARGV); exit;
}


########## ANNOTATION FUNCTIONS ##########
if ($command eq "annot-make" || $command eq "make-annot") {

getopts("O:I:P:ph", \%opt);

# DEFAULTS
@LETTERS = ("0", "A", "B", "C", "D", "E", "F", "G", "H");
%LETTER_NUM = ("A"=>"1", "B"=>"2", "C"=>"3", "D"=>"4", "E"=>"5", "F"=>"6", "G"=>"7", "H"=>"8");

$die2 = "
scitools annot-make [options] [annotation_description_1] [annotaiton_description_2] ...
   or    make-annot

Options:
   -O   [STR]   Output annotation file (default = STDOUT)
   -I   [STR]   Index file
         (default = $VAR{'SCI_index_file'})
         (Index names must be in form of: [Tier]_[set]_[i5/i7]_[A-H/1-12])
   -P   [STR]   Plate descriptor file (instead of written descriptors)
   -p           Print out sample plate file for modificaiton and exit
                (ExamplePlateDescriptor.csv)
   -h           More detailed description of plate / combo specification

";

$die3 = "
scitools annot-make [options] [annotation_description_1] [annotaiton_description_2] ...

Options:
   -O   [STR]   Output annotation file (default = STDOUT)
   -I   [STR]   Index file
         (default = $VAR{'SCI_index_file'})
         (Index names must be in form of: [Tier]_[set]_[i5/i7]_[A-H/1-12])
   -P   [STR]   Plate descriptor file (instead of written descriptors)
   -p           Print out sample plate file for modificaiton and exit
                (ExamplePlateDescriptor.csv)
   -h           More detailed description of plate / combo specification

Annotation descriptors are provided as:

[Annotation_Name]+[Transposase or PCR descriptor]+[Transposase or PCR descriptor]+[etc...]
  must provide at least 1 transposase descriptor and at least 1 PCR descriptor for
  each annotation.

Transposase descriptor:
[TN5],[TN5 index set combination]=[row],[row]:[columns],[row],etc...
 or [NEX]

PCR Descriptor:
[PCR],[PCR index set combination]=[row]:[columns],[row],[row]:[columns],etc...

Before the \"=\" are two comma separated fields. The first is TN5/NEX or PCR to
define what stage of indexing is specified. The second is the index set IDs in
the order of i5 and then i7, (e.g. AA).

After the \"=\" are a series of comma separated fields, where each field
corresponds to a column of a plate, where columns are numbered 1-12. The column
field can be further specified with a subset of rows after a colon, where rows
are A-H.

Note: each index stage can be specified multiple times. The result is an all by
all of the transposase and PCR index sets that are specified.

Example:
My_Sample_1+NEX,AA=ALL+PCR,AC=1-8+PCR,AD=1:A-D,2-5,6:ACDFH,7-12 My_Sample_2+NEX,AA=ALL+NEX,BB=ALL+PCR,GE=1-8

   row specifications can be listed as a range OR letters with no spacing (e.g. ABCGH)
   commas in the column specification are ONLY for separating out columns NOT rows
      (i.e. 1,2,3:A,B would NOT be OK because of the comma between row letters)

";

if (defined $opt{'h'}) {die $die3};

if (defined $opt{'p'}) {
open OUT, ">ExamplePlateDescriptor.csv";
	print OUT "#NEX,MySampleID1,AA,Partial
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
#NEX,MySampleID1,BB,All
#PCR,MySampleID1,CE,Partial
1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,1,1,1,1,1,1
0,0,0,0,0,0,0,0,0,0,0,0
#NEX,MySampleID2,AA,Partial
0,0,0,0,0,0,1,1,1,1,1,1
0,0,0,0,0,0,1,1,1,1,1,1
0,0,0,0,0,0,1,1,1,1,1,1
0,0,0,0,0,0,1,1,1,1,1,1
0,0,0,0,0,0,1,1,1,1,1,1
0,0,0,0,0,0,1,1,1,1,1,1
0,0,0,0,0,0,1,1,1,1,1,1
0,0,0,0,0,0,1,1,1,1,1,1
#PCR,MySampleID2,EE,All
#PCR,MySampleID2,DF,Partial
1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,1,1,1,1,1,1
0,0,0,0,0,0,0,0,0,0,0,0
0,0,0,0,0,0,0,0,0,0,0,0
0,0,0,0,0,0,0,0,0,0,0,0
0,0,0,0,0,0,0,0,0,0,0,0\n";
close OUT;
exit;
}

if (!defined $ARGV[0] && !defined $opt{'P'}) {die $die2};

# Read in index file
open IN, $VAR{'SCI_index_file'};
while ($l = <IN>) {
	chomp $l;
	($ID,$pos,$seq) = split(/\t/, $l);
	($tier,$set,$side,$wells) = split(/_/, $ID);
	if ($tier =~ /(Tn5|Nex)/i) {
		if ($side =~ /i5/) {
			$TN5SET_i5WELLS_seq{$set}{$wells} = $seq;
		} else {
			$TN5SET_i7WELLS_seq{$set}{$wells} = $seq;
		}
		$tier = "Tn5";
	} else {
		if ($side =~ /i5/) {
			$PCRSET_i5WELLS_seq{$set}{$wells} = $seq;
		} else {
			$PCRSET_i7WELLS_seq{$set}{$wells} = $seq;
		}
	}
} close IN;

if (defined $opt{'O'}) {open OUT, ">$opt{'O'}"};

if (defined $opt{'P'}) {
	%NEX_ID_i5_i7_pair = ();
	%PCR_ID_i5_i7_pair = ();
	open IN, "$opt{'P'}";
	while ($l = <IN>) {
		chomp $l;
		if ($l =~ /^#/) {
			($class,$annot,$combo,$subset) = split(/,/, $l);;
			$class =~ s/^#//; $annot =~ s/ /_/g;
			if (!defined $ANNOT_flag{$annot}) {
				$ANNOT_flag{$annot} = $class;
			} else {
				$ANNOT_flag{$annot} .= ",$class";
			}
			($i5_set,$i7_set) = split(//, $combo);
			if ($subset =~ /all/i) {
				for ($rowNum = 1; $rowNum <= 8; $rowNum++) {
					$rowLetter = $LETTERS[$rowNum];
					for ($colNum = 1; $colNum <= 12; $colNum++) {
						if ($class =~ /(Tn5|Nex)/i) {
							$pair = "$TN5SET_i5WELLS_seq{$i5_set}{$rowLetter},$TN5SET_i7WELLS_seq{$i7_set}{$colNum}";
							$NEX_ID_i5_i7_pair{$annot}{$pair} = 1;
						} else {
							$pair = "$PCRSET_i5WELLS_seq{$i5_set}{$rowLetter},$PCRSET_i7WELLS_seq{$i7_set}{$colNum}";
							$PCR_ID_i5_i7_pair{$annot}{$pair} = 1;
						}
					}
				}
			} else {
				for ($rowNum = 1; $rowNum <= 8; $rowNum++) {
					$row = <IN>; chomp $row; $rowLetter = $LETTERS[$rowNum];
					@ROW_COLS = split(/,/, $row); unshift @ROW_COLS, "0";
					for ($colNum = 1; $colNum <= 12; $colNum++) {
						if ($ROW_COLS[$colNum]>0) {
							if ($class =~ /(Tn5|Nex)/i) {
								$pair = "$TN5SET_i5WELLS_seq{$i5_set}{$rowLetter},$TN5SET_i7WELLS_seq{$i7_set}{$colNum}";
								$NEX_ID_i5_i7_pair{$annot}{$pair} = 1;
							} else {
								$pair = "$PCRSET_i5WELLS_seq{$i5_set}{$rowLetter},$PCRSET_i7WELLS_seq{$i7_set}{$colNum}";
								$PCR_ID_i5_i7_pair{$annot}{$pair} = 1;
							}
						}
					}
				}
			}
		}
	} close IN;
	
	foreach $annot (keys %ANNOT_flag) {
		if ($ANNOT_flag{$annot} =~ /(Tn5|Nex)/i && $ANNOT_flag{$annot} =~ /pcr/i) {
			print STDERR "Printing $annot index combinations.\n";
			foreach $NEX_pair (keys %{$NEX_ID_i5_i7_pair{$annot}}) {
				($ix3,$ix1) = split(/,/, $NEX_pair);
				foreach $PCR_pair (keys %{$PCR_ID_i5_i7_pair{$annot}}) {
					($ix4,$ix2) = split(/,/, $PCR_pair);
					if (defined $opt{'O'}) {
						print OUT "$ix1$ix2$ix3$ix4\t$annot\n";
					} else {
						print "$ix1$ix2$ix3$ix4\t$annot\n";
					}
				}
			}
		} else {
			print STDERR "\nWARNING: Transposase (Nex/Tn5) AND PCR specifications must both be included for each annotation!\nBoth were not found for $annot! - SKIPPING!\n";
		}
	}
	
} else {
	foreach $annot_descriptor (@ARGV) {

		@DESCRIPTORS = split(/\+/, $annot_descriptor);
		$annot = shift(@DESCRIPTORS);
		%NEX_i5_i7_pair = ();
		%PCR_i5_i7_pair = ();
		
		print STDERR "\n##### PARSING ANNOT $annot #####\n";
		
		for ($i = 0; $i < @DESCRIPTORS; $i++) {
		
			print STDERR "\n##### Input Parse Descriptor: $DESCRIPTORS[$i] #####\n";
			($class,$columns) = split(/=/, $DESCRIPTORS[$i]);
			print STDERR "\tClass=$class, columns=$columns\n";
			($tier,$combo) = split(/,/, $class);
			print STDERR "\tTier=$tier, combo=$combo\n";
			($i5_set,$i7_set) = split(//, $combo);
			print STDERR "\ti5=$i5_set, i7=$i7_set\n";
			@COLUMN_DESCRIPTORS = split(/,/, $columns);
			for ($j = 0; $j < @COLUMN_DESCRIPTORS; $j++) {
				
				if ($COLUMN_DESCRIPTORS[$j] =~ /ALL/i) {$COLUMN_DESCRIPTORS[$j] = "1-12"};
				
				# determine if a subset of rows are specified (after :)
				if ($COLUMN_DESCRIPTORS[$j] =~ /:/) {
					($col,$rows) = split(/:/, $COLUMN_DESCRIPTORS[$j]);
				} else {
					$col = $COLUMN_DESCRIPTORS[$j];
					$rows = "A-H";
				}
				
				print STDERR "\tRows=$rows, which is:";
				# Add rows to the row list for the column(s)
				@ROW_LIST = split(//, $rows);
				@ROW_INCLUDE = ();
				for ($k = 0; $k < @ROW_LIST; $k++) {
					if ($ROW_LIST[($k+1)] eq "-") {
						$rowStart = $LETTER_NUM{$ROW_LIST[$k]};
						$k++; $k++;
						$rowEnd = $LETTER_NUM{$ROW_LIST[$k]};
						for ($l = $rowStart; $l <= $rowEnd; $l++) {
							push @ROW_INCLUDE, $LETTERS[$l];
							print STDERR " $LETTERS[$l]";
						}
					} else {
						push @ROW_INCLUDE, $ROW_LIST[$k];
						print STDERR " $ROW_LIST[$k]";
					}
				}
				
				# Go through the column(s)
				print STDERR "\n\tCols=$col, which is:";
				@COL_INCLUDE = ();
				if ($col =~ /-/) {
					($startCol,$endCol) = split(/-/, $col);
					for ($l = $startCol; $l <= $endCol; $l++) {
						push @COL_INCLUDE, $l;
						print STDERR " $l";
					}
				} else {
					push @COL_INCLUDE, $col;
					print STDERR " $col";
				}
				
				print STDERR "\n\t  --> Adding all pairs\n";
				# Now add all relevant barcodes to their respective groupings
				foreach $add_col (@COL_INCLUDE) {
					foreach $add_row (@ROW_INCLUDE) {
						if ($class =~ /(Tn5|Nex)/i) {
							$pair = "$TN5SET_i5WELLS_seq{$i5_set}{$add_row},$TN5SET_i7WELLS_seq{$i7_set}{$add_col}";
							$NEX_i5_i7_pair{$pair} = 1;
						} else {
							$pair = "$PCRSET_i5WELLS_seq{$i5_set}{$add_row},$PCRSET_i7WELLS_seq{$i7_set}{$add_col}";
							$PCR_i5_i7_pair{$pair} = 1;
						}
					}
				}
			}
		}

		foreach $NEX_pair (keys %NEX_i5_i7_pair) {
			($ix3,$ix1) = split(/,/, $NEX_pair);
			foreach $PCR_pair (keys %PCR_i5_i7_pair) {
				($ix4,$ix2) = split(/,/, $PCR_pair);
				if (defined $opt{'O'}) {
					print OUT "$ix1$ix2$ix3$ix4\t$annot\n";
				} else {
					print "$ix1$ix2$ix3$ix4\t$annot\n";
				}
			}
		}
	}
}

if (defined $opt{'O'}) {close OUT};

exit;
}

if ($command eq "annot-merge" || $command eq "merge-annot") {

$joiner = "_";

getopts("O:xj:", \%opt);

$die2 = "
scitools annot-merge [options] [annot1] [annot2] (optional annot N etc...)
   or    merge-annot

Will merge annotaitons and include only the intersect

Options:
   -O   [STR]   Output file name / prefix (def = annot1 prefix w/ merge)
   -x           Include cellIDs not present in all files
   -j   [STR]   Character to join annotaitons (def = $joiner)

";

if (!defined $ARGV[1]) {die $die2};
if (defined $opt{'j'}) {$joiner = $opt{'j'}};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.annot$//; $opt{'O'} .= ".merged"};
$opt{'O'} =~ s/\.annot$//;

for ($annotID = 0; $annotID < @ARGV; $annotID++) {
	read_annot($ARGV[$annotID]);
	%{$ANNOTID_CELLID_annot{$annotID}} = %CELLID_annot;
	%CELLID_annot = (); %ANNOT_count = ();
	open IN, "$ARGV[$annotID]";
	while ($l = <IN>) {
		chomp $l; ($cellID,$annot) = split(/\t/, $l);
		$CELLID_annotCT{$cellID}++;
	}
}

open OUT, ">$opt{'O'}.annot";
foreach $cellID (keys %CELLID_annotCT) {
	if ($CELLID_annotCT{$cellID} == @ARGV || defined $opt{'x'}) {
		$new_annot = "";
		for ($annotID = 0; $annotID < @ARGV; $annotID++) {
			if (defined $ANNOTID_CELLID_annot{$annotID}{$cellID}) {
				$new_annot .= $ANNOTID_CELLID_annot{$annotID}{$cellID}.$joiner;
			}
		}
		$new_annot =~ s/$joiner$//;
		print OUT "$cellID\t$new_annot\n";
	}
}
close OUT;

exit;
}

if ($command eq "rename-cells") {

# Defaults
$naming_scheme = "Cell_[number]";

getopts("O:A:N:R:Dx", \%opt);

$die2 = "
scitools rename-cells [options] [input file] (additional input file) etc...
   
Will rename cells according to an annotation, or create an annotation.

Options:
   -O   [STR]   Output prefix (default is input file without file
                extension and 'renamed')
   -R   [STR]   Annotation file for renaming (uses the provided names)
                  Note: if a cell is not in this file, it will be excluded!

   -N   [STR]   New naming schema, includes specific variables:
                  text       = includes in name
                  [number]   = number of the barcode (required)
                  [annot]    = annotation from option -A
                  [orig]     = the original cell ID
                  Example:     MySample_[annot]_[number]
                  Default:     $naming_scheme
   -A   [STR]   Annotation file to include as [annot] in new names
   -x           Exclude cells not in the annot file (-A)
                  (def = annot: Cell)
   -D           Do not create new files, just an annotaiton file to later
                  be used as the -R option

";

if (!defined $ARGV[0]) {die $die2};

if (defined $opt{'O'}) {
	$out = $opt{'O'};
} else {
	$out = $ARGV[0];
	$out =~ s/(\.fq\.gz|\.fq|\.fastq|\.fastq\.gz|\.bam|\.sam|\.matrix|\.tfidf|\.tf|\.values|\.annot|\.annotation|\.dims|\.LSI)$//i;
}

if (defined $opt{'R'} && defined $opt{'A'}) {die "ERROR: Both -R and -A cannot be defined. -R should be used as the only option (other than -O if used).\n"};

# subroutine specific for this command:
sub rename_cell { # note - include an error message if -R is specified and a cellID is found that is not in the -R annot
	if (!defined $opt{'R'}) {
		$newID_number++;
		$rename_origID = $_[0];
		if (defined $opt{'x'}) {
			$rename_annot = "00EXCL00";
		} else {
			$rename_annot = "Cell";
		}
		if (defined $CELLID_annot{$rename_origID}) {
			$rename_annot = $CELLID_annot{$rename_origID};
		}
		$rename_newID = $naming_scheme;
		$rename_newID =~ s/\[annot\]/$rename_annot/;
		$rename_newID =~ s/\[orig\]/$rename_origID/;
		$rename_newID =~ s/\[number\]/$newID_number/;
		return $rename_newID;
	} else {
		return "00EXCL00";
	}
}

%ORIGINAL_newID = ();
if (defined $opt{'R'}) {
	print STDERR "\nUsing existing renaming annotaiton file: $opt{'R'}\n";
	read_annot($opt{'R'});
	foreach $cellID (keys %CELLID_annot) {
		$ORIGINAL_newID{$cellID} = $CELLID_annot{$cellID};
	}
} else {
	# parse scheme
	if (defined $opt{'N'}) {$naming_scheme = $opt{'N'}};
	if ($naming_scheme !~ /\[number\]/) {die "ERROR: Naming scheme must contain the [number] field.\n"};
	if ($naming_scheme =~ /\[annot\]/ && !defined $opt{'A'}) {die "\nERROR: When specifying [annot] in naming scheme, must provide an annot file as -A.\n"};
	if (defined $opt{'A'}) {read_annot($opt{'A'})};
}

$newID_number = 0;
for ($in_file = 0; $in_file < @ARGV; $in_file++) {

	# figure out file input
	if ($ARGV[$in_file] =~ /(\.fq\.gz|\.fq|\.fastq|\.fastq\.gz)$/i) { # fastq
		if ($ARGV[$in_file] =~ /\.gz$/) {
			open IN, "$zcat $ARGV[$in_file] |";
		} else {
			open IN, "$ARGV[$in_file]";
		}
		if (!defined $opt{'D'}) {
			open OUT, "| $gzip > $out.renamed.fq.gz";
		}
		while ($tag = <IN>) {
			chomp $tag; $seq = <IN>; chomp $seq; $null = <IN>; $qual = <IN>; chomp $qual;
			$tag =~ s/^\@//; ($origID,$tail_info) = split(/:/, $tag);
			if (!defined $ORIGINAL_newID{$origID}) {
				$ORIGINAL_newID{$origID} = rename_cell($origID);
			}
			if (!defined $opt{'D'} && $ORIGINAL_newID{$origID} !~ /00EXCL00/) {
				$newID = $ORIGINAL_newID{$origID};
				print OUT "\@$newID:$tail_info\n$seq\n\+\n$qual\n";
			}
		} close IN;
		if (!defined $opt{'D'}) {close OUT};
	} elsif ($ARGV[$in_file] =~ /(\.bam|\.sam)$/i) { # bam/sam
		if ($ARGV[$in_file] =~ /\.bam$/) {
			open IN, "$samtools view -h $ARGV[$in_file] |";
		} else {
			open IN, "$ARGV[$in_file]";
		}
		if (!defined $opt{'D'}) {
			open OUT, "| $samtools view -bS - > $out.renamed.bam 2>/dev/null";
		}
		while ($l = <IN>) {
			chomp $l;
			@P = split(/\t/, $l);
			if ($P[0] =~ /^\@/) {
				if ($P[0] =~ /\@RG/) {
					$RG_lines = "TRUE";
					$origID = $P[1]; $origID =~ s/^ID://;
					$ORIGINAL_newID{$origID} = rename_cell($origID);
					$newID = $ORIGINAL_newID{$origID};
					$out_line .= "\@RG\tID:$newID\tSM:$newID\tLB:$newID\tPL:SCI\n";
				} else {
					$out_line = $l;
				}
			} else {
				($origID,$tail_info) = split(/:/, $P[0]);
				if ($RG_lines eq "TRUE" && !defined $ORIGINAL_newID{$origID}) {
					$out_line = "00EXCL00";
				} else {
					if (!defined $ORIGINAL_newID{$origID}) {
						$ORIGINAL_newID{$origID} = rename_cell($origID);
					}
					$newID = $ORIGINAL_newID{$origID};
					$P[0] = "$newID:$tail_info";
					for ($i = 10; $i < @P; $i++) {
						if ($P[$i] =~ /^RG:Z:/) {
							$P[$i] = "RG:Z:$newID";
						}
					}
					$out_line = join("\t", @P);
				}
			}
			if (!defined $opt{'D'} && $out_line !~ /00EXCL00/) {
				print OUT "$out_line\n";
			}
		} close IN;
		if (!defined $opt{'D'}) {close OUT};
	} elsif ($ARGV[$in_file] =~ /(\.matrix|\.tfidf|\.tf|\.LSI)$/i) { # matrix
		@FP = split(/\./, $ARGV[$in_file]); $extension = pop(@FP);
		open IN, "$ARGV[$in_file]";
		if (!defined $opt{'D'}) {
			open OUT, ">$out.renamed.$extension";
		}
		$header = <IN>; chomp $header;
		@OH = split(/\t/, $header);
		$new_header = ""; @NH = ();
		for ($i = 0; $i < @OH; $i++) {
			$origID = $OH[$i];
			if (!defined $ORIGINAL_newID{$origID}) {
				$ORIGINAL_newID{$origID} = rename_cell($origID);
			}
			$newID = $ORIGINAL_newID{$origID};
			if ($newID !~ /00EXCL00/) {
				$new_header .= "$newID\t";
				push @NH, $newID;
			} else {
				$OH[$i] = "00EXCL00";
			}
		} $new_header =~ s/\t$//;
		if (!defined $opt{'D'}) {
			print OUT "$new_header\n";
			while ($l = <IN>) {
				chomp $l;
				@P = split(/\t/, $l);
				$siteID = shift(@P);
				$out_line = "$siteID";
				for ($i = 0; $i < @OH; $i++) {
					if ($OH[$i] !~ /00EXCL00/) {
						$out_line .= "\t$P[$i]";
					}
				}
				print OUT "$new_header\n";
			}
			close OUT;
		} close IN;
	} elsif ($ARGV[$in_file] =~ /(\.dims|\.values|\.annot|\.annotation)$/i) { # dims/values/annot
		@FP = split(/\./, $ARGV[$in_file]); $extension = pop(@FP);
		open IN, "$ARGV[$in_file]";
		if (!defined $opt{'D'}) {
			open OUT, ">$out.renamed.$extension";
		}
		while ($l = <IN>) {
			chomp $l;
			@P = split(/\t/, $l);
			$origID = $P[0];
			if (!defined $ORIGINAL_newID{$origID}) {
				$ORIGINAL_newID{$origID} = rename_cell($origID);
			}
			$newID = $ORIGINAL_newID{$origID};
			if (!defined $opt{'D'} && $newID !~ /00EXCL00/) {
				$P[0] = $newID;
				$out_line = join("\t", @P);
				print OUT "$out_line\n";
			}
		} close IN;
		if (!defined $opt{'D'}) {close OUT};
	} else {
		print STDERR "ERROR: Cannot determine input file type of $ARGV[$in_file] based ont he file extension.\n";
	}
}

if (!defined $opt{'R'}) { # print out cross-mapping annot files
	open DIR, ">$out.rename.annot";
	open REV, ">$out.rename_reverse.annot";
	foreach $origID (keys %ORIGINAL_newID) {
		if ($ORIGINAL_newID{$origID} !~ /00EXCL00/) {
			print DIR "$origID\t$ORIGINAL_newID{$origID}\n";
			print REV "$ORIGINAL_newID{$origID}\t$origID\n";
		}
	}
	close DIR; close REV;
}

exit;
}


########## BAM FUNCTIONS ##########
if ($command eq "bam-bulk2sci" || $command eq "bulk2sci") {

# Defaults
$memory = "2G";

getopts("s:O:Gm:", \%opt);

$die2 = "
scitools bam-bulk2sci [options] [output bam] [name1]=[bam1] [name2]=[bam2] ...
   or    bulk2sci

Will produce a SCI bam file where each input bam is treated
as a separate cell / barcode combo.

Will exclude /(M|Y|L|K|G|Un|Random|Alt|_)/i chroms

Options:
   -G           Do not add RG lines (def = add them)
   -m   [MEM]   Samtools sort mex memory per thread, K/M/G (def = $memory)
   -s   [STR]   Samtools call (def = $samtools)

";

if (defined $opt{'O'}) {unshift @ARGV, $opt{'O'}};
if (!defined $ARGV[1]) {die $die2};
if (defined $opt{'s'}) {$samtools = $opt{'s'}};
if (!defined $opt{'m'}) {$opt{'m'} = $memory};

$ARGV[0] =~ s/\.bam$//;

$RG_header = "";
for ($i = 1; $i < @ARGV; $i++) {
	($name,$bam) = split(/=/, $ARGV[$i]);
	$NAME_bam{$name} = $bam;
	$RG_header .= "\@RG\tID:$name\tSM:$name\tLB:$name\tPL:bulk2sci\n";
}

$header = "";
open H, "$samtools view -H $bam 2>/dev/null |";
while ($l = <H>) {
	chomp $l;
	$header .= "$l\n";
} close H;

open OUT, "| $samtools view -bSu - 2>/dev/null | $samtools sort -m $memory -T $ARGV[0].TMP - > $ARGV[0].bam 2>/dev/null";
print OUT "$header";
if (!defined $opt{'G'}) {
	print OUT "$RG_header";
}
foreach $name (keys %NAME_bam) {
	$num = 0;
	%READ_num = ();
	open IN, "$samtools view -q 10 $NAME_bam{$name} 2>/dev/null |";
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		if ($P[2] !~ /(M|Y|L|K|G|Un|Random|Alt|_)/i) {
			$tag = shift(@P); $tag =~ s/#.+$//;
			$line = join("\t", @P);
			if (defined $READ_num{$tag}) {
				$read_num = $READ_num{$tag};
			} else {
				$READ_num{$tag} = $num;
				$read_num = $num;
				$num++;
			}
			print OUT "$name:$read_num\t$line\n";
		}
	} close IN;
} close OUT;

exit;
}

if ($command eq "bam-rmdup" || $command eq "rmdup") {

getopts("s:O:", \%opt);

$die2 = "
scitools bam-rmdup [options] [sorted bam file]
   or    rmdup

Will produce a barcode-based duplicate removed bam
and associated complexity file.

Will exclude /(M|Y|L|K|G|Un|Random|Alt)/i chroms

Options:
   -O   [STR]   Output prefix (default is bam file prefix)
   -s   [STR]   Samtools call (def = $samtools)

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.bam$//};
if (defined $opt{'s'}) {$samtools = $opt{'s'}};

open OUT, "| $samtools view -bS - > $opt{'O'}.bbrd.q10.bam 2>/dev/null";

open H, "$samtools view -H $ARGV[0] |";
while ($l = <H>) {print OUT $l};
close H;

open IN, "$samtools view -q 10 $ARGV[0] |";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	($barc,$null) = split(/:/, $P[0]);
	if (defined $KEEP{$P[0]}) {
		print OUT "$l\n";
		$BARC_total{$barc}++;
		$BARC_kept{$barc}++;
	} elsif ($P[1] & 4) {} else {
		if ($P[2] !~ /(M|Y|L|K|G|Un|Random|Alt)/i) {
			$BARC_total{$barc}++;
			if (!defined $BARC_POS_ISIZE{$barc}{"$P[2]:$P[3]"} && !defined $OBSERVED{$P[0]}) {
				$BARC_POS_ISIZE{$barc}{"$P[2]:$P[3]"} = 1;
				$KEEP{$P[0]} = 1;
				print OUT "$l\n";
				$BARC_kept{$barc}++;
			}
			$OBSERVED{$P[0]} = 1;
		}
	}
} close IN; close OUT;


open OUT, ">$opt{'O'}.complexity.txt";
$rank = 1;
foreach $barc (sort {$BARC_kept{$b}<=>$BARC_kept{$a}} keys %BARC_kept) {
	$pct = sprintf("%.2f", ($BARC_kept{$barc}/$BARC_total{$barc})*100);
	print OUT "$rank\t$barc\t$BARC_total{$barc}\t$BARC_kept{$barc}\t$pct\n";
	$rank++;
} close OUT;

exit;
}

if ($command eq "bam-filter" || $command eq "filter-bam") {

getopts("s:O:A:a:N:C:c:L:v", \%opt);

$die2 = "
scitools bam-filter [options] [bam file]
   or    filter-bam

Filters a SCI bam file by a variety of parameters

Options:
   -O   [STR]   Output prefix (default is bam file prefix)
   -A   [STR]   Annotation file
   -a   [STR]   Comma separated list of annotations to include (requires -A)
   -N   [INT]   Minimum read count for each barcode to keep
                (performs much faster if -C is provided)
   -C   [STR]   Complexity file (can be comma separated)
   -c   [STR]   min,max percent complexity to include (requires -C)
                if just one number, will assume max
   -L   [STR]   File listing indexes to include
   -v           Keep the reciprocal of specified filters
   -s   [STR]   Samtools call (def = $samtools)
   
";

if (!defined $ARGV[0]) {die $die2};
$filter_params = "";
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.bam$//; $filter_params .= "O=$opt{'O'}_"};
if (defined $opt{'c'} && !defined $opt{'C'}) {die "\nMust provide a complexity file (-C) if specifying min and max complexity to filter (-c)!\n$die2"};
if (defined $opt{'a'} && !defined $opt{'A'}) {die "\nMust provide an annotaiton file (-A) if specifying annotations to filter (-a)!\n$die2"};
if (defined $opt{'c'}) {
	if ($opt{'c'} =~ /,/) {
		($min_compl,$max_compl) = split(/,/, $opt{'c'});
	} else {
		$min_compl = 0;
		$max_compl = $opt{'c'};
	}
	$filter_params .= "c=$opt{'c'}_";
}
if (defined $opt{'s'}) {$samtools = $opt{'s'}};
if (defined $opt{'A'}) {read_annot($opt{'A'}); $filter_params .= "A=$opt{'A'}_"};
if (defined $opt{'a'}) {
	@ANNOT_LIST = split(/,/, $opt{'a'});
	foreach $annot (@ANNOT_LIST) {
		$ANNOT_include{$annot} = 1;
	}
	$filter_params .= "a=$opt{'a'}_";
}

if (defined $opt{'C'}) {read_complexity($opt{'C'}); $filter_params .= "C=$opt{'C'}_"};
if (defined $opt{'L'}) {
	open LIST, "$opt{'L'}";
	while ($list_line = <LIST>) {chomp $list_line; $CELLID_inList{$list_line} = 1};
	close LIST;
	$filter_params .= "L=$opt{'L'}_";
}

$filter_params =~ s/_$//;

open LOG, ">$opt{'O'}.filt.log";
$ts = localtime(time);
print LOG "$ts\tscitools bam-filter called:\n";
foreach $option (keys %opt) {print LOG "   $option   $opt{$option}\n"};

$included_reads = 0; $total_reads = 0;
if (defined $opt{'N'} && !defined $opt{'C'}) {
	%CELLID_uniq_reads = ();
	open IN, "$samtools view $ARGV[0] 2>/dev/null |";
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		if ($P[2] !~ /(M|Y|L|K|G|Un|Random|Alt|_)/i) {
			$cellID = $P[0];
			$cellID =~ s/:.*$//;
			if (!defined $CELLID_uniq_reads{$cellID}) {
				$CELLID_uniq_reads{$cellID} = 1;
			} else {
				$CELLID_uniq_reads{$cellID}++;
			}
		}
	} close IN;
}

$out_header = "";
open IN, "$samtools view -h $ARGV[0] 2>/dev/null |";
open OUT, "| $samtools view -bS - > $opt{'O'}.filt.bam 2>/dev/null";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$lineType = "READ";
	if ($P[0] =~ /^\@/) {
		if ($P[0] =~ /\@RG/) {
			$lineType = "RG";
		} else {
			$out_header .= "$l\n";
			$lineType = "HEADER";
		}
	}
	
	if ($lineType eq "READ" || $lineType eq "RG") {
	
		if ($lineType eq "READ" && $out_header ne "done") {
			print OUT "$out_header\@PG\tID:scitools_bam-filter_$filter_params\tVN:$version\n";
			$out_header = "done";
		}
	
		if ($lineType eq "READ") {
			$cellID = $P[0];
			$cellID =~ s/:.*$//;
		} else {
			($null,$cellID) = split(/:/, $P[1]);
		}
		
		@INCLUDE_FLAGS = ();
		
		if (defined $opt{'a'}) {
			if (defined $CELLID_annot{$cellID}) {
				$annot = $CELLID_annot{$cellID};
				if (defined $ANNOT_include{$annot}) {
					push @INCLUDE_FLAGS, 1;
				} else {push @INCLUDE_FLAGS, 0};
			} else {push @INCLUDE_FLAGS, 0};
		}
		
		if (defined $opt{'N'}) {
			if ($CELLID_uniq_reads{$cellID}>=$opt{'N'}) {
				push @INCLUDE_FLAGS, 1;
			} else {push @INCLUDE_FLAGS, 0};
		}
		
		if (defined $opt{'L'}) {
			if (defined $CELLID_inList{$cellID}) {
				push @INCLUDE_FLAGS, 1;
			} else {push @INCLUDE_FLAGS, 0};
		}
		
		if (defined $opt{'c'}) {
			if ($CELLID_complexity{$cellID}<=$max_compl&&$CELLID_complexity{$cellID}>=$min_compl) {
				push @INCLUDE_FLAGS, 1;
			} else {push @INCLUDE_FLAGS, 0};
		}
		
		if (defined $opt{'v'}) {
			$print = 1;
			foreach $flag (@INCLUDE_FLAGS) {if ($flag>0) {$print = 0}};
			if ($print>0) {
				if ($lineType eq "READ") {
					print OUT "$l\n";
					$included_reads++;
				} else {
					$out_header .= "$l\n";
				}
			}
		} else {
			$print = 1;
			foreach $flag (@INCLUDE_FLAGS) {if ($flag<1) {$print = 0}};
			if ($print>0) {
				if ($lineType eq "READ") {
					print OUT "$l\n";
					$included_reads++;
				} else {
					$out_header .= "$l\n";
				}
			}
		}
		if ($lineType eq "READ") {$total_reads++};
	}
} close IN; close OUT;

$ts = localtime(time);
print LOG "$ts\tTotal reads: $total_reads, Included reads: $included_reads\n";
close LOG;

exit;
}

if ($command eq "bam-split" || $command eq "split-bam") {

getopts("s:O:A:a:", \%opt);

$die2 = "
scitools bam-split [options] [input bam]
   or    split-bam

Options:
   -A   [STR]   Annotation file (required)
   -O   [STR]   Output prefix (default is bam file prefix)
   -a   [STR]   Comma separated list of annotations to include (requires -A)
                (Default = a seaprate bam for each annotation)
   -s   [STR]   Samtools call (def = $samtools)

";

if (!defined $ARGV[0] || !defined $opt{'A'}) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.bam$//};
if (defined $opt{'s'}) {$samtools = $opt{'s'}};

read_annot($opt{'A'});

$header = "";
open IN, "$samtools view -H $ARGV[0] 2>/dev/null |";
while ($l = <IN>) {
	chomp $l;
	$header .= "$l\n";
} close IN;

if (defined $opt{'a'}) {
	@ANNOT_LIST = split(/,/, $opt{'a'});
	foreach $annot (@ANNOT_LIST) {$ANNOT_include{$annot} = 1};
} else {
	foreach $annot (keys %ANNOT_count) {$ANNOT_include{$annot} = 1};
}

foreach $annot (keys %ANNOT_include) { 
	$file = "$annot\_out";
	open $file, "| $samtools view -bS - 2>/dev/null > $opt{'O'}.$annot.bam";
	print $file "$header";
}

open IN, "$samtools view $ARGV[0] 2>/dev/null |";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	($cellID,$null) = split(/:/, $P[0]);
	$annot = $CELLID_annot{$cellID};
	if (defined $ANNOT_include{$annot}) {
		$file = "$annot\_out";
		print $file "$l\n";
	}
} close IN;

foreach $annot (keys %ANNOT_include) {
	$file = $ANNOT_file{$annot};
	close $file;
}

exit;
}

if ($command eq "bam-merge" || $command eq "merge-bam") { # TODO: Add in runID awareness

getopts("s:", \%opt);

$die2 = "
scitools bam-merge [options] [output bam] [input bam1] [input bam2] ...
   or    merge-bam

Simple wrapper for samtools merge
(use samtools merge for more complicated operations)

Options:
   -s   [STR]   Samtools call (def = $samtools)

";

if (defined $opt{'O'}) {unshift @ARGV, $opt{'O'}};
if (!defined $ARGV[1]) {die $die2};
if (defined $opt{'s'}) {$samtools = $opt{'s'}};

$ARGV[0] =~ s/\.bam$//;
$args = "";
for ($i = 1; $i < @ARGV; $i++) {$args .= "$ARGV[$i] "};

system("$samtools merge $ARGV[0].bam $args");

exit;
}

if ($command eq "bam-addrg" || $command eq "addrg") {

getopts("s:A:L:O:", \%opt);

$die2 = "
scitools bam-addrg [options] [input bam]
   or    addrg

Produces a sorted bam file with read name and RG barcodes.

Options:
   -O   [STR]   Output (default is bam file prefix, adds .RG.bam)
   -A   [STR]   Annotation file (optional, for faster RG pre-processing)
                (Note: will exclude barcodes not in the annot file)
   -L   [STR]   List of CellIDs (optional, for faster RG pre-processing)
                (Note: will exclude barodes not in the list file)
   -s   [STR]   Samtools call (def = $samtools)

";

if (!defined $ARGV[0]) {die $die2};
if (defined $opt{'s'}) {$samtools = $opt{'s'}};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.bam$//; $opt{'O'} =~ s/\.RG//;

if (defined $opt{'L'}) {
	%CELLID_annot = ();
	open IN, "$opt{'L'}";
	while ($cellID = <IN>) {
		chomp $cellID;
		$CELLID_annot{$cellID} = "Cell";
	} close IN;
}
if (defined $opt{'A'}) {read_annot($opt{'A'})};

$RG_header = "";
if (defined $opt{'A'} || defined $opt{'L'}) {
	foreach $cellID (keys %CELLID_annot) {
		$RG_header .= "\@RG\tID:$cellID\tSM:$cellID\tLB:$cellID\tPL:SCI\n";
		$CELLID_include{$cellID} = 1;
	}
} else {
	open IN, "$samtools view $ARGV[0] |";
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		$cellID = $P[0]; $cellID =~ s/:.*$//;
		$RG_header .= "\@RG\tID:$cellID\tSM:$cellID\tLB:$cellID\tPL:SCI\n";
		$CELLID_include{$cellID} = 1;
	} close IN;
}

$out_header = "";
open IN, "$samtools view -h $ARGV[0] 2>/dev/null |";
open OUT, "| $samtools view -bS - 2>/dev/null > $opt{'O'}.RG.bam";
while ($l = <IN>) {
	chomp $l;
	if ($l =~ /^\@/) {
		$out_header .= "$l\n";
	} else {
		if ($out_header ne "done") {
			print OUT $out_header.$RG_header."\@PG\tID:scitools_bam-addrg\tVN:$version\n";
			$out_header = "done";
		} else {
			@P = split(/\t/, $l);
			$cellID = $P[0]; $cellID =~ s/:.*$//;
			print OUT "$l\tRG:Z:$cellID\n";
		}
	}
} close IN; close OUT;

exit;
}

if ($command eq "bam-aggregate" || $command eq "aggregate-bam") {

getopts("O:s:r", \%opt);

$die2 = "
scitools bam-aggregate [options] [bam file] [annotation file]
   or    aggregate-bam

Note: if a non counts matrix is provided, it will still sum the values.
Support for comma-separated annotations (i.e. cell aggregation with
oversampling) has been added.

Options:
   -O   [STR]   Output prefix (default is [input annot].aggregate.bam)
   -s   [STR]   Samtools call (def = $samtools)
   -r           Add RG fields (def = no)

";

if (!defined $ARGV[1]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[1]; $opt{'O'} =~ s/\.annot$//};
if (defined $opt{'s'}) {$samtools = $opt{'s'}};
read_annot($ARGV[1]);

open OUT, "| $samtools view -bS - 2>/dev/null >$opt{'O'}.aggregate.bam";
open HEAD, "$samtools view -H $ARGV[0] 2>/dev/null |";
while ($l = <HEAD>) {
	chomp $l;
	@P = split(/\t/, $l);
	if ($P[0] !~ /RG/) {
		print OUT "$l\n";
	}
} close HEAD;

if (defined $opt{'r'}) {
	foreach $annot_field (keys %ANNOT_count) {
		@ANNOTS = split(/,/, $annot_field);
		foreach $annot (@ANNOTS) {
			if (!defined $ANNOT_inRG{$annot}) {
				print OUT "\@RG\tID:$annot\tSM:$annot\tLB:$annot\tPL:SCI_aggregate\n";
				$ANNOT_inRG{$annot}++;
			}
		}
	}
}

open IN, "$samtools view $ARGV[0] 2>/dev/null |";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	($cellID,$other) = split(/:/, $P[0]);
	if (defined $CELLID_annot{$cellID}) {
		@ANNOTS = split(/,/, $CELLID_annot{$cellID});
		foreach $annot (@ANNOTS) {
			print OUT "$annot:$cellID.$other";
			for ($i = 1; $i < @P; $i++) {
				if ($P[$i] !~ /^RG:Z:/) {
					print OUT "\t$P[$i]";
				}
			}
			if (defined $opt{'r'}) {
				print OUT "\tRG:Z:$annot";
			}
			print OUT "\n";
		}
	}
} close IN; close OUT;

exit;
}

if ($command eq "bam-project" || $command eq "project-bam" || $command eq "project") { # TODO: Add in memory efficiency modification from standalone version

# Defaults
$training_increment = 0.01;
$complexity_increment = 0.05;
$min_unique = 1000;
$read_increments = 10000000;
$gradient_def = "BuYlRd";

getopts("s:O:n:t:c:R:r:XPfG:", \%opt);

$die2 = "
scitools bam-project [options] [non rmdup bam OR rand_reads.txt]
   or    project-bam
   or    project

Will use the current sequencing to build a projection curve
and project out the expected unique reads at additional levels
of read depth. (use only a NON RMDUP bam). Projects to 0.05
complexity (median of cells).

If run with a bam, will create a [output].rand_reads.txt file
which can be used for subsequent runs for a faster runtime.

This model typically predicts within 1-2% of the actual unique
read counts; however, its intended use is just to get a rough
estimate of the target reads desired to continue sequencing a
library to the desired complexity point.

Note: It is also important to consider the bam input file has
already been barcode matched - therefore to predict raw read
counts for additional sequencing, you need to account for the
unmatching reads. (ie divide by the fraction matching)

Options:
   -O   [STR]   Output dir prefix (default is bam file prefix)
   -n   [INT]   Minimum unique read count to include cellID
                (def = $min_unique)
   -t   [FLT]   Fraction of reads to use as each input data point
                (def = $training_increment)
   -c   [FLT]   Complexity increments (median) to print at
                (def = $complexity_increment)
   -r   [INT]   Read counts to increment projections
                (def = $read_increments)
   -G   [GRD]   Color gradient in plot (def = $gradient_def)
                For all available gradients, run 'scitools gradient'
   -s   [STR]   Samtools call (def = $samtools)
   -R   [STR]   Rscript call (def = $Rscript)
   -X           Keep intermediate files (def = delete)
   -f           Filter out models
                (agressive; may fail, but if it works, it may
                 be a more accurate projection)

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.bam$//; $opt{'O'} =~ s/\.rand_reads.txt$//;
if (defined $opt{'s'}) {$samtools = $opt{'s'}};
if (defined $opt{'t'}) {$training_increment = $opt{'t'}};
if (defined $opt{'n'}) {$min_unique = $opt{'n'}};
if (defined $opt{'c'}) {$complexity_increment = $opt{'c'}};
if (!defined $opt{'G'}) {$opt{'G'} = $gradient_def};
$gradient_function = get_gradient($opt{'G'});

if (-e "$opt{'O'}.read_projections") {
	print STDERR "\n\nWARNING: $opt{'O'}.read_projections Directory Exists! Will Rewrite Contents!\n";
} else {
	system("mkdir $opt{'O'}.read_projections");
}

$readCT_for_proj = 0;
if ($ARGV[0] =~ /\.bam$/) {
	open IN, "$samtools view $ARGV[0] 2>/dev/null |";
	while ($l = <IN>) {
		chomp $l;
		$readCT_for_proj++;
		@P = split(/\t/, $l);
		($cellID,$null) = split(/:/, $P[0]);
		$seed = 1;
		while (defined $SEED_info{$seed}) {$seed = rand(1e20)};
		if ($P[1] & 4) {$P[4] = 0};
		if ($P[2] =~ /(M|Y|L|K|G|Un|Random|Alt)/i) {$P[4] = 0};
		$P[0] =~ s/#.+$//;
		$SEED_info{$seed} = "$P[0]\t$P[4]\t$P[2]\t$P[3]";
	} close IN;
	open OUT, ">$opt{'O'}.rand_reads.txt";
	foreach $seed (sort keys %SEED_info) {
		print OUT "$SEED_info{$seed}\n";
	} close OUT;
	%SEED_info = ();
	open IN, "$opt{'O'}.rand_reads.txt";
} elsif ($ARGV[0] =~ /\.rand_reads\.txt$/) {
	open IN, "$ARGV[0]";
	while ($l = <IN>) {$readCT_for_proj++};
	close IN;
	open IN, "$ARGV[0]";
} else {
	die "\nCannot determine file type. Must be .bam or .rand_reads.txt\n$die2";
}

$totalCT = 0; $allKept = 0;
$increment = int(($readCT_for_proj*$training_increment)+1);
$report = $increment; $report_num = 0;
%CELLID_total = (); %CELLID_kept = ();
open OUT, ">$opt{'O'}.read_projections/model_all.txt";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	($cellID,$null) = split(/:/, $P[0]);
	$totalCT++;
	if (defined $KEEP{$P[0]}) {
		$CELLID_total{$cellID}++;
		$CELLID_kept{$cellID}++;
		$allKept++;
	} elsif ($P[1] < 10) {} else {
		$CELLID_total{$cellID}++;
		if (!defined $CELLID_POS_ISIZE{$cellID}{"$P[2]:$P[3]"} && !defined $OBSERVED{$P[0]}) {
			$CELLID_POS_ISIZE{$cellID}{"$P[2]:$P[3]"} = 1;
			$KEEP{$P[0]} = 1;
			$CELLID_kept{$cellID}++;
			$allKept++;
		}
		$OBSERVED{$P[0]} = 1;
	}
	if ($totalCT>=$report) {
		$report_num++;
		foreach $cellID (keys %CELLID_total) {
			$frac = sprintf("%.3f", $CELLID_kept{$cellID}/$CELLID_total{$cellID});
			print OUT "$cellID\t$report_num\t$report\t$CELLID_total{$cellID}\t$CELLID_kept{$cellID}\t$frac\t".(1/$CELLID_total{$cellID})."\t".(1/$CELLID_kept{$cellID})."\n";
		}
		$report += $increment;
	}
} close IN; close OUT;

open IN, "$opt{'O'}.read_projections/model_all.txt";
open OUT, ">$opt{'O'}.read_projections/model_increments.txt";
while ($l = <IN>) {
	chomp $l; @P = split(/\t/, $l);
	$cellID = $P[0];
	if ($CELLID_kept{$cellID} >= $min_unique) {
		print OUT "$l\n";
	}
} close IN; close OUT;
if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.read_projections/model_all.txt");
}

open R, ">$opt{'O'}.read_projections/model.r";
print R "
INCREMENTS<-read.table(\"$opt{'O'}.read_projections/model_increments.txt\")
colnames(INCREMENTS) <- c(\"cellID\",\"ReportID\",\"ReportCT\",\"Total\",\"Uniq\",\"Complexity\",\"X\",\"Y\")\n";

$cellIDs_included = 0;
foreach $cellID (keys %CELLID_kept) {
	if ($CELLID_kept{$cellID} >= $min_unique) {
		$CELLID_fraction{$cellID} = $CELLID_total{$cellID}/$totalCT;
		if ($cellIDs_included<1) {
			print R "MODEL<-rbind(c(\"$cellID\",t(as.matrix(lm(subset(INCREMENTS,cellID==\"$cellID\")\$Y ~ subset(INCREMENTS,cellID==\"$cellID\")\$X)\$coefficients))))\n";
		} else {
			print R "MODEL<-rbind(MODEL,c(\"$cellID\",t(as.matrix(lm(subset(INCREMENTS,cellID==\"$cellID\")\$Y ~ subset(INCREMENTS,cellID==\"$cellID\")\$X)\$coefficients))))\n";
		}
		$cellIDs_included++;
	}
}
print R "write.table(MODEL,file=\"$opt{'O'}.read_projections/model.txt\",row.names=FALSE,col.names=FALSE,quote=FALSE,sep=\"\\t\")\n";
close R;

system("$Rscript $opt{'O'}.read_projections/model.r");

$test_count = 10000; $failed_models = 0; $passed_models = 0;
open MODEL, "$opt{'O'}.read_projections/model.txt";
while ($l = <MODEL>) {
	chomp $l;
	($cellID,$int,$slope) = split(/\t/, $l);
	$Vmax = (1/$int); $Km = ($slope/$int);
	$test_uniq = (($Vmax*$test_count)/($Km+$test_count));
	if ($test_count<$test_uniq || !defined $opt{'f'}) {
		$CELLID_Vmax{$cellID} = $Vmax;
		$CELLID_Km{$cellID} = $Km;
		$passed_models++;
	} else {
		$failed_models++;
	}
} close MODEL;

if ($passed_models>0) {
	$fail_frac = sprintf("%.3f", $failed_models/$passed_models);
} else {$fail_frac=1};
if ($fail_frac > 0.1) {
	print STDERR "\nWARNING: $failed_models cell models failed, $passed_models passed ($fail_frac)\n";
	if ($fail_frac > 0.5) {
		die "ERROR: failed model fraction less than 0.5, not enough data to proceed. Exiting.\n";
	}
}

open CELL_PROJ, ">$opt{'O'}.read_projections/cell_projections.txt";
open MED_PROJ, ">$opt{'O'}.read_projections/summary_projections.txt";

$projected_complexity = 1;
$previous_complexity = 1.1;
$projected_read_total = $read_increments;
while ($projected_complexity > 0.05) {
	@COMPLEXITY = (); @UNIQ_RDS = (); $uniq_read_sum = 0;
	foreach $cellID (keys %CELLID_Vmax) {
		$cell_total_read_count = int($projected_read_total*$CELLID_fraction{$cellID});
		$cell_unique_read_count = int(($CELLID_Vmax{$cellID}*$cell_total_read_count)/($CELLID_Km{$cellID}+$cell_total_read_count));
		$cell_unique_read_frac = sprintf("%.4f", $cell_unique_read_count/$cell_total_read_count);
		print CELL_PROJ "$cellID\t$projected_read_total\t$cell_total_read_count\t$cell_unique_read_count\t$cell_unique_read_frac\n";
		push @COMPLEXITY, $cell_unique_read_frac;
		push @UNIQ_RDS, $cell_unique_read_count;
		$uniq_read_sum += $cell_unique_read_count;
	}
	@COMPLEXITYs = sort {$a<=>$b} @COMPLEXITY; $len = @COMPLEXITYs;
	$projected_complexity = sprintf("%.2f", $COMPLEXITYs[int($len/2)]);
	@UNIQ_RDSs = sort {$a<=>$b} @UNIQ_RDS;
	$projected_unique = $UNIQ_RDSs[int($len/2)];
	$projected_mean = int($uniq_read_sum/$len);
	if ($projected_complexity != $previous_complexity) {
		print MED_PROJ "$projected_complexity\t$projected_read_total\t$projected_unique\t$projected_mean\n";
		$previous_complexity = $projected_complexity;
	}
	$projected_read_total += $read_increments;
}
close CELL_PROJ; close MED_PROJ;

open R, ">$opt{'O'}.read_projections/plot_projections.r";
print R "
library(ggplot2)
$gradient_function
CELLS<-read.table(\"$opt{'O'}.read_projections/cell_projections.txt\")
SUMMARY<-read.table(\"$opt{'O'}.read_projections/summary_projections.txt\")
PLT<-ggplot() + theme_bw() +
	geom_line(aes((CELLS\$V5*100),log10(CELLS\$V4),group=CELLS\$V1),alpha=0.15,size=0.15,color=\"lightsteelblue4\") +
	geom_line(aes((SUMMARY\$V1*100),log10(SUMMARY\$V4)),color=\"black\",size=1,linetype=\"dashed\") +
	geom_line(aes((SUMMARY\$V1*100),log10(SUMMARY\$V3)),color=\"black\",size=1) +
	geom_point(aes((SUMMARY\$V1*100),log10(SUMMARY\$V3)),color=\"black\",size=3) +
	geom_point(aes((SUMMARY\$V1*100),log10(SUMMARY\$V3),color=log10(SUMMARY\$V2)),size=2) +
	scale_color_gradientn(colours=gradient_funct(99)) +
	scale_x_continuous(limits=c(0,100)) +
	scale_y_continuous(limits=c(2,6)) +
	xlab(\"Complexity\") +
	ylab(\"log10 Unique Reads\") +
	labs(color=\"Log10 Total\\nReads\")
ggsave(plot=PLT,filename=\"$opt{'O'}.read_projections/projected_complexity.png\",width=6,height=5)
ggsave(plot=PLT,filename=\"$opt{'O'}.read_projections/projected_complexity.pdf\",width=6,height=5)
"; close R;

system("$Rscript $opt{'O'}.read_projections/plot_projections.r 2>/dev/null");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.read_projections/plot_projections.r $opt{'O'}.read_projections/model.r $opt{'O'}.read_projections/model.txt $opt{'O'}.read_projections/model_increments.txt");
}

exit;
}


########## ATAC FUNCTIONS ##########
if ($command eq "atac-callpeak" || $command eq "atac-callpeaks" || $command eq "callpeak" || $command eq "callpeaks") {

# Defaults:
$min_feature_size = 500;

getopts("s:m:b:O:l:f:X", \%opt);

$die2 = "
scitools atac-callpeak [options] [duplicate removed and filtered bam file]
   or    callpeak(s)

Options:
   -O   [STR]   Output prefix (default is bam file prefix)
   -l   [INT]   Minimum feature size (def = $min_feature_size);
   -b   [STR]   Bedtools call (def = $bedtools)
   -m   [STR]   Macs2 call (def = $macs2)
   -s   [STR]   Samtools call (def = $samtools)
   -f   [STR]   Fai file for chr lengths (shorcut examples: hg19, hg38, and mm10 if in .cfg)
                If toggled will ensure n peaks extend beyond
   -X           Retain intermediate files (def = remove)

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.bam$//};
if (defined $opt{'l'}) {$min_feature_size = $opt{'l'}};
if (defined $opt{'s'}) {$samtools = $opt{'s'}};

if (defined $opt{'f'}) {
	if (defined $REF{$opt{'f'}}) {
		open FAI, "$REF{$opt{'f'}}.fai";
	} else {
		open FAI, "$opt{'f'}";
	}
	while ($l = <FAI>) {
		chomp $l;
		@P = split(/\t/, $l);
		$CHR_length{$P[0]} = $P[1];
	} close FAI;
}

system("$macs2 callpeak -t $ARGV[0] -n $opt{'O'} >> $opt{'O'}.macs2.log 2>> $opt{'O'}.macs2.log");

open IN, "$bedtools merge -i $opt{'O'}_peaks.narrowPeak 2>/dev/null |";
open OUT, "| $bedtools sort -i - 2>/dev/null | $bedtools merge -i - > $opt{'O'}.$min_feature_size.tmp 2>/dev/null";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	if ($P[0] !~ /(M|Y|L|K|G|Un|Random|Alt|_)/i) {
		if (($P[2]-$P[1])<$min_feature_size) {
			$mid = ($P[2]+$P[1])/2;
			$start = int($mid-($min_feature_size/2));
			$end = int($mid+($min_feature_size/2));
			print OUT "$P[0]\t$start\t$end\n";
		} else {
			print OUT "$P[0]\t$P[1]\t$P[2]\n";
		}
	}
} close IN; close OUT;

open IN, "$opt{'O'}.$min_feature_size.tmp";
open OUT, ">$opt{'O'}.$min_feature_size.bed";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	if (defined $opt{'f'}) {
		if ($P[2] < $CHR_length{$P[0]}) {
			print OUT "$P[0]\t$P[1]\t$P[2]\t$P[0]_$P[1]_$P[2]\n";
		}
	} else {
		print OUT "$P[0]\t$P[1]\t$P[2]\t$P[0]_$P[1]_$P[2]\n";
	}
} close IN; close OUT;

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.$min_feature_size.tmp $opt{'O'}_peaks.narrowPeak $opt{'O'}_model.r $opt{'O'}_peaks.xls $opt{'O'}_summits.bed");
}

exit;
}

if ($command eq "atac-mergepeak" || $command eq "atac-mergepeaks" || $command eq "mergepeak" || $command eq "mergepeaks") {

getopts("b:O:", \%opt);

$die2 = "
scitools atac-mergepeak [options] [peak bed file] [peak bed file 2] ...
   or    mergepeak(s)

Options:
   -O   [STR]   Output prefix (default is stdout, adds .bed)
   -b   [STR]   Bedtools call (def = $bedtools)

";

if (!defined $ARGV[1]) {die $die2};
if (defined $opt{'O'}) {$opt{'O'} =~ s/\.bed//; open OUT, ">$opt{'O'}.bed"};

$bed_list = ""; for ($i = 0; $i < @ARGV; $i++) {$bed_list .= "$ARGV[$i] "};
open IN, "cat $bed_list | $bedtools sort -i - 2>/dev/null | $bedtools merge -i - 2>/dev/null |";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	if (defined $opt{'O'}) {
		print OUT "$P[0]\t$P[1]\t$P[2]\t$P[0]_$P[1]_$P[2]\n";
	} else {
		print "$P[0]\t$P[1]\t$P[2]\t$P[0]_$P[1]_$P[2]\n";
	}
} close IN;

if (defined $opt{'O'}) {close OUT};

exit;
}

if ($command eq "atac-count" || $command eq "atac-counts" || $command eq "counts" || $command eq "count") {

getopts("s:b:O:BC:X:", \%opt);

$die2 = "
scitools atac-count [options] [bam file] [peaks bed file or matrix]
   or    count(s)

Creates a counts or binary matrix file and a values file
with the fraction of reads on target for each cell. In addition if matrix defined instead of bed file
can add cells to counts matrix if cells are new (at moment merge bam if cell names do not overlap).

Options:
   -O   [STR]   Output prefix (default is peaks prefix)
                (adds .counts.matrix)
   -b   [STR]   Bedtools call (def = $bedtools)
   -s   [STR]   Samtools call (def = $samtools)
   -B           Make it a binary matrix instead of counts matrix
                (file name will end in .binary.matrix)
   -C   [STR]   Complexity file (speeds up on-target calc)
   -X 		Remove temp files

";

if (!defined $ARGV[1]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[1]; $opt{'O'} =~ s/\.bed//};

#if ends in .matrix then uses peaks from matrix to do counts or binary
if ($ARGV[1] =~ /\.matrix$/)
{
print "USING MATRIX FILE, ADDING CELLS: \n";
#read in matrix
read_matrix($ARGV[1]);
$opt{'O'} = $ARGV[1]; $opt{'O'} =~ s/\.binary.matrix//; $opt{'O'} =~ s/\.counts.matrix//; $opt{'O'} =~ s/\.matrix//;
# create temporary bed file from matrix pre defined peaks
$tembedfilename=$ARGV[0];
$tembedfilename=~ s/\.bam//;
$tembedfilename=~ s/.*\///;
open OUT, ">$tembedfilename.temp.bed";
foreach $peak (@MATRIX_ROWNAMES) 
{
	print OUT join("\t",split("_",$peak))."\n"}; 
close(OUT);

#create options string
	$common_opts = "";
	if (defined $opt{'b'}) {$common_opts .= "-b $opt{'b'} "};
	if (defined $opt{'s'}) {$common_opts .= "-s $opt{'s'} "};
	if (defined $opt{'B'}) {$common_opts .= "-B $opt{'B'} "};
	if (defined $opt{'C'}) {$common_opts .= "-C $opt{'C'} "};
	$common_opts =~ s/\s$//;
	#call same script but on temp bed
	system("scitools atac-count $common_opts $ARGV[0] $tembedfilename.temp.bed");
	#join matrixes and remove temp matrix and bed file if X is not defined
	#DEV: merging matrix files if overlap in cell names see notes in merge-matrix
	if (defined $opt{'B'}) {
	system("scitools matrix-merge -O $opt{'O'}_$tembedfilename $ARGV[1] $tembedfilename.temp.binary.matrix");
	if (!defined $opt{'X'}) {system("rm -f $tembedfilename.temp.bed $tembedfilename.binary.temp.matrix")};
	}
	else
	{
	system("scitools matrix-merge -O $opt{'O'}_$tembedfilename $ARGV[1] $tembedfilename.temp.counts.matrix");
	if (!defined $opt{'X'}) {system("rm -f $tembedfilename.temp.bed $tembedfilename.counts.temp.matrix")};
	}
}
else
{

if (defined $opt{'C'}) {
	read_complexity($opt{'C'});
} else {
	open IN, "$samtools view $ARGV[0] |";
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		if ($P[2] !~ /(M|Y|L|K|G|Un|Random|Alt)/i) {
			($cellID,$null) = split(/:/, $P[0]);
			$CELLID_uniq_reads{$cellID}++;
		}
	} close IN;
}

open IN, "$bedtools intersect -abam $ARGV[0] -b $ARGV[1] -bed -wa -wb |";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	($cellID,$null) = split(/:/, $P[3]);
	$siteID = $P[12]."_".$P[13]."_".$P[14];
	if (!defined $SITE_CELL_STATUS{$siteID}{$cellID}) {
		$SITE_CELL_STATUS{$siteID}{$cellID} = 1;
	} else {
		if (defined $opt{'B'}) {
			$SITE_CELL_STATUS{$siteID}{$cellID} = 1;
		} else {
			$SITE_CELL_STATUS{$siteID}{$cellID}++;
		}
	}
	$CELLID_onTarget{$cellID}++;
	if (!defined $CELL_IDS{$cellID}) {
		push @CELL_ID_LIST, $cellID;
		$CELL_IDS{$cellID} = 1;
	}
} close IN;

if (defined $opt{'B'}) {
	open OUT, ">$opt{'O'}.binary.matrix";
} else {
	open OUT, ">$opt{'O'}.counts.matrix";
}

print OUT "$CELL_ID_LIST[0]";
for ($i = 1; $i < @CELL_ID_LIST; $i++) {
	print OUT "\t$CELL_ID_LIST[$i]";
} print OUT "\n";

foreach $siteID (sort keys %SITE_CELL_STATUS) {
	print OUT "$siteID";
	for ($i = 0; $i < @CELL_ID_LIST; $i++) {
		if ($SITE_CELL_STATUS{$siteID}{$CELL_ID_LIST[$i]}>0.5) {
			print OUT "\t$SITE_CELL_STATUS{$siteID}{$CELL_ID_LIST[$i]}";
		} else {
			print OUT "\t0";
		}
	} print OUT "\n";
} close OUT;

open OUT, ">$opt{'O'}.fracOnTarget.values";
for ($i = 0; $i < @CELL_ID_LIST; $i++) {
	if ($CELLID_uniq_reads{$cellID}>0) {
		$cellID = $CELL_ID_LIST[$i];
		$frac = sprintf("%.3f", $CELLID_onTarget{$cellID}/$CELLID_uniq_reads{$cellID});
		print OUT "$cellID\t$frac\n";
	}
} close OUT;
}
exit;
}

if ($command eq "atac-chromvar" || $command eq "chromvar") {

# Defaults
$count_cutoff = 1000;
$frac_in_peaks = 0.15;
$motif_set = "human_pwms_v2";
%MOTIF_SETS = (
   "human_pwms_v1" => 1,
   "mouse_pwms_v1" => 1,
   "human_pwms_v2" => 1,
   "mouse_pwms_v2" => 1,
   "homer_pwms" => 1,
   "encode_pwms" => 1
);
%GENOMES = (
   "hg19" => "BSgenome.Hsapiens.UCSC.hg19",
   "hg38" => "BSgenome.Hsapiens.UCSC.hg38",
   "mm10" => "BSgenome.Mmusculus.UCSC.mm10"
);

getopts("O:R:Xg:p:M:c:", \%opt);

$die2 = "
scitools atac-chromvar [options] [bam file] [peaks bed file]
   or    chromvar

Bam file MUST have RG lines present. If they are not, add them using
'scitools bam-addrg'. This wrapper wille xecute chromVAR and output
deviation and deviation z-score matrix files, along with a variaiton
file for each motif or peak set assessed.

Options:
   -O   [STR]   Output prefix (default is bam prefix)
                (creates prefix.chromVAR directory)
   -g   [STR]   Genome (hg38, hg19, or mm10; def = hg38)
   -M   [STR]   Motif set (def = $motif_set)
   -c   [INT,FRAC]   count_cutoff,frac_in_peaks (def = $count_cutoff,$frac_in_peaks)
   -R   [STR]   Rscript call (def = $Rscript)
   -X           Rewrite existing directory

Note: Requires the chromVAR R package to be installed and loadable.
      It can be found here: github.com/GreenleafLab/chromVAR
      It is also recommended to install chromVARmotifs vaialble here:
	  github.com/GreenleafLab/chromVARmotifs

Motif Sets: human_pwms_v1, mouse_pwms_v1, human_pwms_v2, mouse_pwms_v2,
            homer_pwms, encode_pwms

";

if (!defined $ARGV[1]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[1]; $opt{'O'} =~ s/\.bam//};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (!defined $opt{'g'}) {$opt{'g'} = "hg38"};
if (!defined $GENOMES{$opt{'g'}}) {die "\n\nERROR: Genome $opt{'g'} is not a prper genome reference! Exiting!\n"} else {$genome = $GENOMES{$opt{'g'}}};
if (!defined $opt{'M'}) {$opt{'M'} = "human_pwms_v2"};
if (!defined $MOTIF_SETS{$opt{'M'}}) {die "\n\nERROR: The motif set provided ($opt{'M'}) does not exist in chromVARmotifs\n"};
if (defined $opt{'c'}) {($count_cutoff,$frac_in_peaks) = split(/,/, $opt{'c'})};

if (-e "$opt{'O'}.chromVAR") {
	if (defined $opt{'X'}) {
		print STDERR "\n\nWARNING: $opt{'O'}.chromVAR exists, will rewrite contents!\n";
	} else {
		die "\n\nERROR: $opt{'O'}.chromVAR exists! Exiting!\n";
	}
}

system("mkdir $opt{'O'}.chromVAR");

$peakfile = $ARGV[1];
@BAMS = split(/,/, $ARGV[0]);
$bam_list = "c(";
foreach $bam (@BAMS) {$bam_list .= "$bam,"};
$bam_list =~ s/,$/\)/;

open R, ">$opt{'O'}.chromVAR/chromVAR.r";
$ts = localtime(time);
print R "# chromVAR script generated: $ts
# bam: $ARGV[0]
# peak bed: $ARGV[1]
# genome: $opt{'g'}
# motif set: $opt{'M'}
# count_cutoff,frac_cutoff = $count_cutoff,$frac_in_peaks

# load libraries
library(chromVAR)
library(chromVARmotifs)
library(motifmatchr)
library(SummarizedExperiment)
library(Matrix)
library(BiocParallel)
register(MulticoreParam(8))
library(JASPAR2016)
library(ggplot2)

# load motifs and genome
library($genome)
motifs <- data(\"$opt{'M'}\")

# read in peaks and filter them
peaks <- getPeaks(\"$ARGV[1]\",sort=TRUE)
peaks <- resize(peaks, width = 500, fix = \"center\")

# make counst matrix and normalize / filter
counts<-getCounts(\"$ARGV[0]\", peaks, paired = TRUE, by_rg = TRUE, format = \"bam\", colData = DataFrame(celltype = \"cells\"))
counts <- addGCBias(counts, genome = $genome)
counts <- filterSamples(counts, min_depth = $count_cutoff, min_in_peaks = $frac_in_peaks, shiny = FALSE)
counts <- filterPeaks(counts, non_overlapping = TRUE)

# make motif index
motifs <- getJasparMotifs()
motif_ix <- matchMotifs(motifs, counts, genome = $genome)

# calculate & print deviations
dev <- computeDeviations(object = counts, annotations = motif_ix)
write.table(as.matrix(deviations(dev)),file = \"$opt{'O'}.chromVAR/deviations.matrix\", col.names = TRUE, row.names = TRUE, sep = \"\\t\", quote = FALSE)
write.table(as.matrix(deviationScores(dev)),file = \"$opt{'O'}.chromVAR/deviation_scores.matrix\", col.names = TRUE, row.names = TRUE, sep = \"\\t\", quote = FALSE)

# calculate & print variabilities
var <- computeVariability(dev)
write.table(as.matrix(var),file = \"$opt{'O'}.chromVAR/variability.txt\", col.names = TRUE, row.names = TRUE, sep = \"\\t\", quote = FALSE)
plot<-plotVariability(var,use_plotly=FALSE)
ggsave(plot,file = \"$opt{'O'}.chromVAR/variability.png\")

# generate tSNE on the deviations
#############tsne causes segmentation faults############
###################needs to be corrected###############
#tsne <- deviationsTsne(dev, threshold = 1.5, perplexity = 10, shiny = FALSE)
#write.table(as.matrix(tsne),file = \"$opt{'O'}.chromVAR/tsne.dims\", col.names = TRUE, row.names = TRUE, sep = \"\\t\", quote = FALSE)

"; close R;

system("$Rscript $opt{'O'}.chromVAR/chromVAR.r >> $opt{'O'}.chromVAR/chromVAR.log 2>> $opt{'O'}.chromVAR/chromVAR.log");

exit;
}


########## SIGNAL FUNCTIONS ##########
if ($command eq "signal-make" || $command eq "make-signal") {

# Defaults
$span_size = 5000;
$sub_size = 500;
$slide_size = 250;
$matching_annot = "TRUE";

getopts("O:A:a:W:w:l:b:s:Xn:r", \%opt);

$die2 = "
scitools signal-make [options] [bam file] [bed file of features]
   or    make-signal

This tool will generate a specialized matrix file with multiple entries
over each window provided. Columns are the annotaiton and the sub-window
and rows are the bed file features. It also produces a z-scored singal
matrix where the z-score is produced for all sub windows within the
annotation.

Options:
   -O   [STR]   Output prefix (default is [bam prefix].[bed prefix].signal)
   -A   [STR]   Annotation file (will aggregate signal over the annotations)
                  Note: if not provided, will assume signal for every cellID
                  wich is not recommended - for plotting individual cells, use
                  'scitools plot-reads'
   -a   [STR]   Comma separated list of annotations to include (requires -A)
                  Will position annotations in the specified order
   -r           Order rows by the order of the input bed file (overrides -n)
   -n   [STR]   Sort by the signal for this annotation.
                  Default is to sort by all (if no -A), by subsets of matching
                  bed and cell annotations (ie bed names are the same as annot
                  names), or by the first annotation listed.
   -W   [INT]   Window size - the span of the signal view (def = $span_size)
   -w   [INT]   Sub-window size for signal (def = $sub_size)
   -l   [INT]   Sub-window slide size (def = $slide_size)
   -b   [STR]   Bedtools call (def = $bedtools)
   -s   [STR]   Samtools call (def = $samtools)
   -X           Do not delete intermediate files (def = delete)

";

if (!defined $ARGV[1]) {die $die2};
if (defined $opt{'a'} && !defined $opt{'A'}) {die "\nMust provide an annotaiton file (-A) if specifying annotations to filter (-a)!\n$die2"};
if (!defined $opt{'O'}) {
	$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/bam$//;
	$opt{'O'} .= $ARGV[1]; $opt{'O'} =~ s/\.bed$//;
}
if (defined $opt{'W'}) {$span_size = $opt{'W'}};
if (defined $opt{'w'}) {$sub_size = $opt{'w'}};
if (defined $opt{'l'}) {$slide_size = $opt{'l'}};

if (!defined $opt{'A'}) {
	print STDERR "
WARNING: No annotation specified - will treat as a single annotation.\n";
	open IN, "$samtools view $ARGV[0] |";
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		($cellID,$other) = split(/:/, $P[0]);
		$CELLID_annot{$cellID} = $cellID;
		$ANNOT_read_count{'all'}++;
	} close IN;
} else {
	read_annot($opt{'A'});
	open IN, "$samtools view $ARGV[0] |";
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		($cellID,$other) = split(/:/, $P[0]);
		if (defined $CELLID_annot{$cellID}) {
			$ANNOT_read_count{$CELLID_annot{$cellID}}++;
		}
	} close IN;
}

if (defined $opt{'a'}) {
	@ANNOT_LIST = split(/,/, $opt{'a'});
	foreach $annot (@ANNOT_LIST) {
		$ANNOT_include{$annot} = 1;
	}
} else {
	foreach $annot (sort keys %ANNOT_read_count) {
		$ANNOT_include{$annot} = 1;
		push @ANNOT_LIST, $annot;
	}
}

$flank_size = int($span_size/2);

# make temporary bed file with subwindows
open IN, "$ARGV[1]";
open OUT, "| $bedtools sort -i - > $opt{'O'}.signal.subwin.bed";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$chr = $P[0];
	$start = $P[1];
	$end = $P[2];
	$winID = "$chr\_$start\_$end";
	if (defined $P[3]) {
		$fgroup = $P[3];
		if (!defined $ANNOT_include{$P[3]}) {
			$matching_annot = "FALSE";
		}
	} else {
		$fgroup = "Feature";
	}
	$FGROUP_count{$fgroup}++;
	$WINID_fgroup{$winID} = $fgroup;
	$FGROUP_WINID_list{$fgroup}{$winID} = 1;
	$WINID_list{$winID} = 1;
	$midPt = int(($start+$end)/2);
	if (($midPt-$flank_size)>0) {
		$subNum = 0;
		for ($i = ($midPt-$flank_size); $i <= ($midPt+$flank_size)-$sub_size; $i += $slide_size) {
			print OUT "$chr\t$i\t".($i+$sub_size)."\t$winID\t$subNum\n";
			$subNum++;
		}
		$totwin = $subNum-1;
	}
	$featureCT++;
} close IN; close OUT;

# intitalize counts
foreach $winID (keys %WINID_list) {
	for ($i = 0; $i <= $totwin; $i++) {
		foreach $annot (keys %ANNOT_include) {
			$WINID_SUB_ANNOT_count{$winID}{"$annot\_window_$i"} = 0;
		}
	}
}

# intersect to get counts
open IN, "$bedtools intersect -abam $ARGV[0] -b $opt{'O'}.signal.subwin.bed -bed -wa -wb 2>/dev/null |";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	($cellID,$other) = split(/:/, $P[3]);
	if (defined $CELLID_annot{$cellID}) {
		$annot = $CELLID_annot{$cellID};
		$winID = $P[15];
		$subWin = $P[16];
		$subID = "$annot\_window_$subWin";
		$SUBID_annot{$subID} = $annot;
		$WINID_SUB_ANNOT_count{$winID}{$subID}++;
	}
} close IN;

# print the signal matrix
open OUT, ">$opt{'O'}.raw.signal";
open NRM, ">$opt{'O'}.norm.signal";

$header = ""; @COLNAMES = ();
for ($annot_pos = 0; $annot_pos < @ANNOT_LIST; $annot_pos++) {
	$annot = $ANNOT_LIST[$annot_pos];
	for ($subWin = 0; $subWin <= $totwin; $subWin++) {
		$header .= "$annot\_window_$subWin\t";
		push @COLNAMES, "$annot\_window_$subWin";
	}
} $header =~ s/\t$//;
print OUT "$header\n";
print NRM "$header\n";

# order rows
$midWin = int($totwin/2); @ROW_ORDER = ();
if (defined $opt{'r'}) {
	open IN, "$ARGV[1]";
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		$winID = "$P[0]\_$P[1]\_$P[2]";
		if (defined $WINID_list{$winID}) {
			push @ROW_ORDER, $winID;
		}
	} close IN;
} elsif (!defined $opt{'A'}) {
	$midWinID = "all_window_$midWin";
	foreach $fgroup (sort keys %FGROUP_WINID_list) {
		foreach $winID (sort {$WINID_SUB_ANNOT_count{$b}{$midWinID}<=>$WINID_SUB_ANNOT_count{$a}{$midWinID}} keys %{$FGROUP_WINID_list{$fgroup}}) {
			push @ROW_ORDER, $winID;
		}
	}
} else {
	if (defined $opt{'n'}) {
		$midWinID = "$opt{'n'}\_window_$midWin";
		foreach $winID (sort {$WINID_SUB_ANNOT_count{$b}{$midWinID}<=>$WINID_SUB_ANNOT_count{$a}{$midWinID}} keys %WINID_SUB_ANNOT_count) {
			push @ROW_ORDER, $winID;
		}
	} elsif ($matching_annot = "TRUE") {
		for ($annot_pos = 0; $annot_pos < @ANNOT_LIST; $annot_pos++) {
			$annot = $ANNOT_LIST[$annot_pos];
			$midWinID = "$annot\_window_$midWin";
			foreach $winID (sort {$WINID_SUB_ANNOT_count{$b}{$midWinID}<=>$WINID_SUB_ANNOT_count{$a}{$midWinID}} keys %{$FGROUP_WINID_list{$annot}}) {
				push @ROW_ORDER, $winID;
			}
		}
	} else {
		($annot,$null) = split(/_window_/, $COLNAMES[0]);
		print STDERR "INFO: Annotaiton for cells provided, but annotaitons of peaks do not all match cell annotations.
      If this is incorrect - double check that annotaitons and bed feature names match.
      Proceeding with sorting on the first listed annotation: $annot\n";
		$midWinID = "$annot\_window_$midWin";
		foreach $fgroup (sort keys %FGROUP_WINID_list) {
			foreach $winID (sort {$WINID_SUB_ANNOT_count{$b}{$midWinID}<=>$WINID_SUB_ANNOT_count{$a}{$midWinID}} keys %{$FGROUP_WINID_list{$fgroup}}) {
				push @ROW_ORDER, $winID;
			}
		}
	}
}


# print raw and norm
for ($row_pos = 0; $row_pos < @ROW_ORDER; $row_pos++) {
	$winID = $ROW_ORDER[$row_pos];
	print OUT "$winID"; print NRM "$winID";
	for ($colNum = 0; $colNum < @COLNAMES; $colNum++) {
		$subID = $COLNAMES[$colNum];
		$annot = $SUBID_annot{$subID};
		if (!defined $WINID_SUB_ANNOT_count{$winID}{$subID}) {
			$WINID_SUB_ANNOT_count{$winID}{$subID} = 0;
		}
		print OUT "\t$WINID_SUB_ANNOT_count{$winID}{$subID}";
		$nrm = sprintf("%.4f", ($WINID_SUB_ANNOT_count{$winID}{$subID}/$ANNOT_read_count{$annot})*($featureCT*$totwin));
		print NRM "\t$nrm";
		$ANNOT_posCT{$annot}++;
		$ANNOT_sum{$annot}+=$WINID_SUB_ANNOT_count{$winID}{$subID};
		push @{$ANNOT_values{$annot}}, $WINID_SUB_ANNOT_count{$winID}{$subID};
	}
	print OUT "\n"; print NRM "\n";
}
close OUT; close NRM;

# calculate each annot mean
foreach $annot (keys %ANNOT_posCT) {
	$ANNOT_mean{$annot} = $ANNOT_sum{$annot}/$ANNOT_posCT{$annot};
	$stdev_sum = 0;
	foreach $value (@{$ANNOT_values{$annot}}) {
		$stdev_sum += ($ANNOT_mean{$annot} - $value)**2;
	}
	$ANNOT_stdev{$annot} = sqrt($stdev_sum/$ANNOT_posCT{$annot});
}

open OUT, ">$opt{'O'}.zscore.signal";
print OUT "$header\n";

for ($row_pos = 0; $row_pos < @ROW_ORDER; $row_pos++) {
	$winID = $ROW_ORDER[$row_pos];
	print OUT "$winID";
	for ($colNum = 0; $colNum < @COLNAMES; $colNum++) {
		$subID = $COLNAMES[$colNum];
		$annot = $SUBID_annot{$subID};
		$zscore = sprintf("%.4f", ($WINID_SUB_ANNOT_count{$winID}{$subID}-$ANNOT_mean{$annot})/$ANNOT_stdev{$annot});
		print OUT "\t$zscore";
	}
	print OUT "\n";
}
close OUT;


if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.signal.subwin.bed");
}

exit;
}


########## MATRIX FUNCTIONS ##########
if ($command eq "matrix-summarize" || $command eq "summarize-matrix") { # TODO: EXPAND ON THIS!

getopts("O:R:XA:a:", \%opt);

$die2 = "
scitools matrix-summarize [options] [matrix of any kind]
   or    summarize-matrix
   
Generates summaries on the matrix file and plots some of the properties.

Options:
   -O   [STR]   Output prefix (default is [input_matrix].suffix)
   -R   [STR]   Rscript call (def = $Rscript)
   -X           Retain intermediate files (def = delete)

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};

# initialize stats
$ncol = 0; $nrow = 0;
$maxVal = "na"; $minVal = "na";
$sum = 0;
$nonzero = 0; $nrow_nonzero = 0; $ncol_nonzero = 0;
$mean = 0; $stdev = 0;

open STATS, ">$opt{'O'}.stats";
open ROW_DATA, ">$opt{'O'}.row.stats";
print ROW_DATA "mean\tstdev\tmin\tmax\tsum\tnonzero\n";
open COL_DATA, ">$opt{'O'}.col.stats";
print COL_DATA "mean\tstdev\tmin\tmax\tsum\tnonzero\n";
open COL_DIST, ">$opt{'O'}.col.vals";

open IN, "$ARGV[0]";
$h = <IN>; chomp $h; @H = split(/\t/, $h);
$ncol = @H;

for ($i = 0; $i < @H; $i++) {
	@{$COL_values[$i]} = ();
}

while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$rowID = shift(@P);
	$row_sum = 0; $row_mean = 0; $row_stdev_sum = 0; $row_stdev = 0; $row_nonzero = 0; $row_min = $P[0]; $row_max = $P[0];
	for ($i = 0; $i < @H; $i++) {
		$sum+=$P[0];
		if ($P[$i] != 0) {$nonzero++; $row_nonzero++};
		if ($minVal eq "na" || $P[$i] < $minVal) {$minVal = $P[$i]};
		if ($maxVal eq "na" || $P[$i] > $maxVal) {$maxVal = $P[$i]};
		$row_sum += $P[$i];
		if ($P[$i] < $row_min) {$row_min = $P[$i]};
		if ($P[$i] > $row_max) {$row_max = $P[$i]};
		push @{$COL_values[$i]}, $P[$i];
	}
	$row_mean = $row_sum/$ncol;
	for ($i = 0; $i < @H; $i++) {
		$row_stdev_sum += ($P[$i] - $row_mean)**2;
	}
	$row_stdev = sqrt($row_stdev_sum/$ncol);
	print ROW_DATA "$rowID\t$row_mean\t$row_stdev\t$row_min\t$row_max\t$row_sum\t$row_nonzero\n";
	if ($row_nonzero>0) {$nrow_nonzero++};
	$nrow++;
} close IN; close ROW_DATA;

$mean = $sum/($ncol*$nrow);
$stdevSum = 0;

for ($i = 0; $i < @H; $i++) {
	$col_sum = 0; $col_mean = 0; $col_stdev_sum = 0; $col_stdev = 0; $col_nonzero = 0; $col_min = $COL_values[$i][0]; $col_max = $COL_values[$i][0];
	for ($j = 0; $j < @{$COL_values[$i]}; $j++) {
		$val = $COL_values[$i][$j];
		$col_sum += $val;
		if ($val > $col_max) {$col_max = $val};
		if ($val < $col_min) {$col_min = $val};
		if ($val != 0) {$col_nonzero++};
	}
	$col_mean = $col_sum/@{$COL_values[$i]};
	for ($j = 0; $j < @{$COL_values[$i]}; $j++) {
		$val = $COL_values[$i][$j];
		$col_stdev_sum += ($val - $col_mean)**2;
		$stdev_sum += ($val - $mean)**2;
		print COL_DIST "$H[$i]\t$val\n";
	}
	$col_stdev = sqrt($col_stdev_sum/@{$COL_values[$i]});
	print COL_DATA "$H[$i]\t$col_mean\t$col_stdev\t$col_min\t$col_max\t$col_sum\t$col_nonzero\n";
	if ($col_nonzero>0) {$ncol_nonzero++};
} close COL_DATA; close COL_DIST;

$stdev = sqrt($stdev_sum/($ncol*$nrow));

print STATS "Matrix        $ARGV[0]
Columns       $ncol
Rows          $nrow
Mean          $mean
Stdev         $stdev
Min           $minVal
Max           $maxVal
Nonzero_rows  $nrow_nonzero
Nonzero_cols  $ncol_nonzero
Nonzero_vals  $nonzero
Matrix_sum    $sum
"; close STATS;

exit;
}

if ($command eq "matrix-filter" || $command eq "filter-matrix") { # TODO: Mem efficiency (at some point)

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

exit;
}

if ($command eq "matrix-naomit" || $command eq "naomit-matrix") { 


getopts("O:", \%opt);

$die2 = "
scitools matrix-naomit [options] [matrix]
   or    naomit-matrix
   
Note: Print out number of NA's in matrix then remove NA's from matrix

Options:
   -O   [STR]   Output prefix (default is [input].naomit.matrix)


";

if (!defined $ARGV[0]) {die $die2};

if (!defined $opt{'O'}) {
	$opt{'O'} = $ARGV[0];
	$opt{'O'} =~ s/\.matrix$//;
}
#count how many
$matrix_NA=0;
open IN, $ARGV[0];
open OUT, ">$opt{'O'}.naomit.matrix";
$head = <IN>; print OUT "$head";
while ($line = <IN>) {
	if ($line=~ m/NA/) {
		$matrix_NA++;
	} else {
		print OUT "$line";
	}
} close IN; close OUT;

open LOG, ">$opt{'O'}.log";
$ts = localtime(time);
print LOG "$ts scitools naomit
Matrix file = $ARGV[0]
Options:
";
foreach $option (keys %opt) {
	print LOG "   $option   $opt{$option}\n";
}
print LOG "
Total rows with NA: $matrix_NA. \n
";
close LOG;

exit;
}

if ($command eq "matrix-tf" || $command eq "tf") {

getopts("O:N:", \%opt);

$die2 = "
scitools matrix-tf [options] [counts matrix]
   or    tf

Term frequency normalization of matrix

Options:
   -O   [STR]   Output prefix (default is [input].tf)
   -N   [VAL]   Normalization constant (def = row number)

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};

read_matrix($ARGV[0]);
if (!defined $opt{'N'}) {$opt{'N'} = $matrix_rowNum};

open IN, "$ARGV[0]";
open OUT, ">$opt{'O'}.tf";
$h = <IN>; print OUT "$h";
while ($l = <IN>) {
	chomp $l;
	@P = split (/\t/, $l);
	$rowID = shift(@P);
	print OUT "$rowID";
	for ($i = 0; $i < @P; $i++) {
		$tf = ($P[$i]/$MATRIX_CellID_signal{$MATRIX_COLNAMES[$i]});
		$score = $tf*$opt{'N'};
		print OUT "\t$score";
	}
	print OUT "\n";
} close IN; close OUT;

exit;
}

if ($command eq "matrix-tfidf" || $command eq "tfidf") {

getopts("O:", \%opt);

$die2 = "
scitools matrix-tfidf [options] [counts matrix]
   or    tfidf

Options:
   -O   [STR]   Output prefix (default is [input].tfidf)

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};

read_matrix($ARGV[0]);

open IN, "$ARGV[0]";
open OUT, ">$opt{'O'}.tfidf";
$h = <IN>; print OUT "$h";
while ($l = <IN>) {
	chomp $l;
	@P = split (/\t/, $l);
	$rowID = shift(@P);
	print OUT "$rowID";
	for ($i = 0; $i < @P; $i++) {
		$tf = ($P[$i]/$MATRIX_CellID_signal{$MATRIX_COLNAMES[$i]});
		$idf = (log(1+($matrix_colNum/($MATRIX_feature_signal{$rowID}+1))));
		$score = ($tf*$idf);
		print OUT "\t$score";
	}
	print OUT "\n";
} close IN; close OUT;

exit;
}

if ($command eq "matrix-lsi" || $command eq "lsi") {

# Defaults
$range_default = "1-15";

getopts("O:D:d:XR:", \%opt);

$die2 = "
scitools matrix-lsi [options] [tfidf matrix]
   or    lsi

Options:
   -O   [STR]   Output prefix (def = [input].tfidf.LSI.matrix)
   -D   [RNG]   Range of dimensions to include (range format)
                (e.g. 1-5,6,8; def = $range_default)
   -X           Retain intermediate files (def = delete)
   -R   [STR]   Rscript call (def = $Rscript)

Note: Requires svd R package

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
if (!defined $opt{'D'}) {$opt{'D'} = $range_default};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
read_ranges($opt{'D'});

open PARAMS, ">$opt{'O'}.LSI.params";
print PARAMS "scitools matrix-lsi ON $ARGV[0]
Dimensions: $range_R_set\n";
close PARAMS;

open R, ">$opt{'O'}.LSI.SVD.r";
print R "
library(svd)
IN<-read.table(\"$ARGV[0]\")
SVD<-svd(as.matrix(IN))
ds<-diag(SVD\$d[$range_R_set])
us<-as.matrix(SVD\$u[,$range_R_set])
vs<-as.matrix(SVD\$v[,$range_R_set])
write.table(ds,file=\"$opt{'O'}.LSI.SVD.d\",sep=\"\\t\",row.names=FALSE,col.names=FALSE)
write.table(us,file=\"$opt{'O'}.LSI.SVD.u\",sep=\"\\t\",row.names=FALSE,col.names=FALSE)
write.table(vs,file=\"$opt{'O'}.LSI.SVD.v\",sep=\"\\t\",row.names=FALSE,col.names=FALSE)
";
close R;

system("$Rscript $opt{'O'}.LSI.SVD.r");

open DIAG, "$opt{'O'}.LSI.SVD.d";
$rowNum = 0;
while ($l = <DIAG>) {
	chomp $l;
	@P = split(/\t/, $l);
	$diagSum+=$P[$rowNum];
	$rowNum++;
} close DIAG;

open DIAG, "$opt{'O'}.LSI.SVD.d";
open OUT, ">$opt{'O'}.LSI.SVD.d.norm";

$rowNum = 0;
while ($l = <DIAG>) {
	chomp $l;
	@P = split(/\t/, $l);
	$P[$rowNum] = ($P[$rowNum]/$diagSum);
	$out = join("\t", @P);
	print OUT "$out\n";
	$rowNum++;
} close DIAG; close OUT;

open R, ">$opt{'O'}.LSI.SVD.reconstruct.r";
print R "
D<-as.matrix(read.table(\"$opt{'O'}.LSI.SVD.d.norm\"))
U<-as.matrix(read.table(\"$opt{'O'}.LSI.SVD.u\"))
V<-as.matrix(read.table(\"$opt{'O'}.LSI.SVD.v\"))
LSI<-t(scale(V %*% D %*% t(U)))
write.table(LSI,file=\"$opt{'O'}.LSI.SVD.reconstructed\",sep=\"\\t\",row.names=FALSE,col.names=FALSE)
";
close R;

system("$Rscript $opt{'O'}.LSI.SVD.reconstruct.r");

open IN, "$ARGV[0]";
open LSI, "$opt{'O'}.LSI.SVD.reconstructed";
open OUT, ">$opt{'O'}.LSI.matrix";
$h = <IN>; chomp $h;
print OUT "$h\n";
while ($l = <IN>) {
	$v = <LSI>;
	chomp $l; chomp $v;
	@P = split(/\t/, $l);
	$siteID = shift(@P);
	print OUT "$siteID\t$v\n";
} close IN; close OUT; close LSI;

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.LSI.SVD.r $opt{'O'}.LSI.SVD.d $opt{'O'}.LSI.SVD.d.norm $opt{'O'}.LSI.SVD.u $opt{'O'}.LSI.SVD.v $opt{'O'}.LSI.SVD.reconstructed $opt{'O'}.LSI.SVD.reconstruct.r");
}

exit;
}

if ($command eq "matrix-zscore" || $command eq "zscore-matrix") {

getopts("O:CRG", \%opt);

$die2 = "
scitools matrix-zscore [options] [matrix]
   or    zscore-matrix
   
Performs z-scoring on the matrix

Options:
   -O   [STR]   Output prefix (default is [input].zscore.matrix)
   -R           Z-score rows
   -C           Z-score columns
   -G           Z-score globally

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};
if ((defined $opt{'R'} && defined $opt{'C'}) ||
    (defined $opt{'R'} && defined $opt{'G'}) ||
	(defined $opt{'C'} && defined $opt{'G'}) ||
	(!defined $opt{'R'} && !defined $opt{'C'} && !defined $opt{'G'})) {die "ERROR: Please perform row, column, OR global z-scoring.\n"};

if (defined $opt{'R'}) {
	open OUT, ">$opt{'O'}.R_zscore.matrix";
} elsif (defined $opt{'C'}) {
	open OUT, ">$opt{'O'}.C_zscore.matrix";
	@SUMS = (); @COUNTS = (); @MEANS = (); @STDEV_SUM = (); @STDEV = ();
} elsif (defined $opt{'G'}) {
	open OUT, ">$opt{'O'}.G_zscore.matrix";
	$sum = 0; $count = 0; $stdev_sum = 0; $stdev = 0;
}

open IN, "$ARGV[0]";
$h = <IN>; chomp $h; @H = split(/\t/, $h);
print OUT "$h\n";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$siteID = shift(@P);
	if (defined $opt{'R'}) {
		$sum = 0;
		for ($i = 0; $i < @P; $i++) {
			$sum += $P[$i];
		}
		$mean = $sum/(@P);
		$stdev_sum = 0;
		for ($i = 0; $i < @P; $i++) {
			$stdev_sum += ($mean-$P[$i])**2;
		}
		$stdev = sqrt($stdev_sum/(@P));
		if ($stdev != 0) {
			print OUT "$siteID";
			for ($i = 0; $i < @P; $i++) {
				$zscore = ($P[$i]-$mean)/$stdev;
				print OUT "\t$zscore";
			}
			print OUT "\n";
		} else {
			print STDERR "WARNING: row $siteID has a stdev of 0.\n";
		}
	} elsif (defined $opt{'C'}) {
		for ($i = 0; $i < @P; $i++) {
			$SUMS[$i] += $P[$i];
			$COUNTS[$i]++;
		}
	} elsif (defined $opt{'G'}) {
		for ($i = 0; $i < @P; $i++) {
			$sum += $P[$i]; $count++;
		}
	}
} close IN;

if (defined $opt{'R'}) {close OUT; exit};

if (defined $opt{'C'}) {
	for ($i = 0; $i < @SUMS; $i++) {
		$MEANS[$i] = $SUMS[$i]/$COUNTS[$i];
	}
} else {
	$mean = $sum/$count;
}

open IN, "$ARGV[0]";
$h = <IN>;
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$siteID = shift(@P);
	if (defined $opt{'C'}) {
		for ($i = 0; $i < @P; $i++) {
			$STDEV_SUM[$i] += ($MEANS[$i]-$P[$i])**2;
		}
	} elsif (defined $opt{'G'}) {
		for ($i = 0; $i < @P; $i++) {
			$stdev_sum += ($mean-$P[$i])**2;
		}
	}
} close IN;

if (defined $opt{'C'}) {
	for ($i = 0; $i < @SUMS; $i++) {
		$STDEV[$i] = sqrt($STDEV_SUM[$i]/$COUNTS[$i]);
	}
} else {
	$stdev = sqrt($stdev_sum/$count);
}

open IN, "$ARGV[0]";
$h = <IN>;
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$siteID = shift(@P);
	print OUT "$siteID";
	if (defined $opt{'C'}) {
		for ($i = 0; $i < @P; $i++) {
			$zscore = ($P[$i]-$MEANS[$i])/$STDEV[$i];
			print OUT "\t$zscore";
		}
	} elsif (defined $opt{'G'}) {
		for ($i = 0; $i < @P; $i++) {
			$zscore = ($P[$i]-$mean)/$stdev;
			print OUT "\t$zscore";
		}
	}
	print OUT "\n";
} close IN;

exit;
}

if ($command eq "matrix-tsne" || $command eq "tsne") {

# Deafults
$dims = 2;
$perp = 30;

getopts("O:XD:R:P:", \%opt);

$die2 = "
scitools matrix-tsne [options] [input matrix]
   or    tsne

Produces a dims file with tSNE coordinates

Options:
   -O   [STR]   Output prefix (default is [input].tSNE)
   -D   [INT]   Dimensions to embed tSNE in (def = $dims)
   -P   [INT]   Perplexity (def = $perp)
   -X           Retain intermediate files (def = delete)
   -R   [STR]   Rscript call (def = $Rscript)

Note: Requires Rtsne R package

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};
if (defined $opt{'D'}) {$dims = $opt{'D'}};
if (defined $opt{'P'}) {$perp = $opt{'P'}};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};

open R, ">$opt{'O'}.p$perp.tSNE.r";
print R "
library(Rtsne)
tIN<-t(read.table(\"$ARGV[0]\"))
TSNE<-Rtsne(tIN,dims=$dims,perplexity=$perp,check_duplicates=FALSE)
write.table(TSNE\$Y,file=\"$opt{'O'}.p$perp.tSNE.raw\",col.names=FALSE,row.names=FALSE)
"; close R;

system("$Rscript $opt{'O'}.p$perp.tSNE.r 2>/dev/null");

open IN, "$ARGV[0]";
$h = <IN>; close IN;
chomp $h; $h =~ s/\r//;
@CELLID_list = split(/\t/, $h);

open IN, "$opt{'O'}.p$perp.tSNE.raw";
open OUT, ">$opt{'O'}.p$perp.tSNE.dims";
while ($l = <IN>) {
	chomp $l;
	$cellID = shift(@CELLID_list);
	print OUT "$cellID\t$l\n";
} close IN; close OUT;

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.p$perp.tSNE.r $opt{'O'}.p$perp.tSNE.raw");
}

exit;
}

if ($command eq "matrix-pca" || $command eq "pca") {

getopts("O:XR:", \%opt);

$die2 = "
scitools matrix-pca [options] [input matrix]
   or    pca

Produces a dims file with PCA coordinates

Options:
   -O   [STR]   Output prefix (default is [input].PCA.dims)
   -X           Retain intermediate files (def = delete)
   -R   [STR]   Rscript call (def = $Rscript)

Note: Requires methods R package

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};

open R, ">$opt{'O'}.PCA.r";
print R "
library(methods)
tM<-t(read.table(\"$ARGV[0]\"))
D<-dist(tM,method=\"euclidean\",diag=FALSE,upper=FALSE)
PCA<-prcomp(D,center=TRUE)
write.table(PCA\$x,file=\"$opt{'O'}.PCA.dims\",sep=\"\\t\",row.names=TRUE,col.names=FALSE,quote=FALSE)
";
close R;

system("$Rscript $opt{'O'}.PCA.r");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.PCA.r");
}

exit;
}

if ($command eq "matrix-approx-factors" || $command eq "matrix-factors" || $command eq "factors") {

getopts("O:XR:", \%opt);

$die2 = "
scitools matrix-approx_factors [options] [input tf-idf matrix]
   or    factors

Produces a k-range vs error txt of NMF

Options:
   -O   [STR]   Output prefix (default is [input].factors.txt)
   -X           Retain intermediate files (def = delete)
   -R   [STR]   Rscript call (def = $Rscript)

Note: Requires Seurat and swne R packages so dependencies need to look for those

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};

open R, ">$opt{'O'}.factors.r";
print R "
library(Seurat)
library(swne)
norm<-as.matrix(read.table(file=\"$ARGV[0]\",row.names=1), \"dgCMatrix\")

#read in tf-idf matrix into dgc format, might want to do original count matrix
norm.counts<-as.matrix(norm,\"dgCMatrix\")

#or do original matrix and maybe frequency transform
#norm.counts <- ScaleCounts(counts, batch = NULL, method = \"ft\", adj.var = T)

## Unguided NMF
loss <- \"mse\" ## Loss function
#loss <- \"mkl\" ## Loss function
k.range <- seq(1,100,1) ## Range of factors to iterate over
n.cores <- 25 ## Number of cores to use
seed <- 32566 ## Set seed for 

## Identify optimal number of factors

n.comp.res <- FindNumFactors(norm.counts, k.range = k.range, n.cores = n.cores, do.plot = F, loss = loss,max.iter = 10000)
output=data.frame(\"k\"=k.range,\"recon.err\"=n.comp.res\$err[loss,])
output2=data.frame(\"k\"=k.range,n.comp.res\$err)

#output k vs error txt to plot with plot-k
write.table(output,file=\"$opt{'O'}.factors.txt\",sep=\"\\t\",row.names=FALSE,col.names=TRUE,quote=FALSE)

#output k vs error txt to plot with plot-k
write.table(output2,file=\"$opt{'O'}.allerr.factors.txt\",sep=\"\\t\",row.names=FALSE,col.names=TRUE,quote=FALSE)


";
close R;

system("$Rscript $opt{'O'}.factors.r");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.factors.r");
}

exit;
}

if ($command eq "matrix-swne" || $command eq "swne" || $command eq "piglet") {
# Defaults
$range_default = "1-15";

getopts("O:XRAC:D:", \%opt);

$die2 = "
scitools matrix-swne [options] [input tf-idf matrix] [k value based on the reconstruction error]
   or    swne
   or    piglet
 Applies SWNE based dim reduction and clustering on given TF-IDF matrix 

Options:
   -O   [STR]   Output prefix (default is [input].k.)
   -X           Retain intermediate files (def = delete)
   -D           Dims to use for pca def: 1-15
   -R   [STR]   Rscript call (def = $Rscript)
   -A           Annotation file, when provided [input].cells.factors.annot is output where factors are added 
   -C			Color Annotation file, when provided [input].colors.annot is output where factor with color back is added 

Note: Requires Seurat and swne R packages so dependencies need to look for those

";


if (!defined $ARGV[0]) {die $die2};
if (!defined $ARGV[1]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (!defined $opt{'D'}) {$opt{'D'} = $range_default};
read_ranges($opt{'D'});
$k=$ARGV[1];


open R, ">$opt{'O'}.SWNE.r";
print R "


#packages
library(Seurat)
library(swne)

norm.matrix<-read.table(file=\"$ARGV[0]\",row.names=1)

#use seurat package for PCA
S_matrix  <- CreateSeuratObject(raw.data = norm.matrix, project = \"SWNE\")

#Create scaled data, which is TF-IDF matrix in our case
S_matrix\@scale.data<-S_matrix\@raw.data
#run SVD much like in normal analysis, this is redundant in future
max_PC=max($range_R_set)+5


#first calculate PCA, calculate first 20 PC usuallym, this is done by the irlba package which computes fast truncated svd (very similar to what we do normally
S_matrix <- RunPCA(object = S_matrix,pcs.compute=max_PC,pc.genes = rownames(S_matrix\@data),reduction.name=\"svd\", do.print = F)

#png(\"$opt{'O'}.PCA.bow.plot.png\")
#PCElbowPlot(S_matrix, num.pc = max_PC)
#dev.off()

#calculate SNN based on PCA/svd and find internal clusters
S_matrix <- FindClusters(object = S_matrix, reduction.type = \"svd\", dims.use = $range_R_set, resolution = 0.6, print.output = 0, save.SNN = TRUE)
#or use this if you already know cell type, add case when cell types read in
#S_matrix <- BuildSNN(S_matrix, dims.use = $range_R_set, k.param = 30, k.scale = 10, prune.SNN = 1/20, force.recalc = T)

#if find clusters run to build SNN
clusters <- S_matrix\@ident; names(clusters) <- S_matrix\@cell.names;

#write out SNN clusters
write.table(as.matrix(clusters),file=\"$opt{'O'}.SNN.Clusters.$k.SWNE.annot\",quote=FALSE,sep=\"\\t\",row.names=TRUE,col.names=FALSE)

norm.counts<-as.matrix(S_matrix\@scale.data,\"dgCMatrix\")

## Unguided NMF
loss <- \"mse\" ## Loss function
n.cores <- 30 ## Number of cores to use
seed <- 32566 ## Set seed for 

# do NMF with k provided
nmf.res <- RunNMF(norm.counts, k = $ARGV[1], alpha = 0, init = \"nnsvd\", n.cores = n.cores, loss = loss)

#calc weight matrix
nmf.res\$W <- ProjectFeatures(norm.counts, nmf.res\$H, loss = \"mse\", n.cores = n.cores)
nmf.scores <- nmf.res\$H

## Run SWNE embedding
snn.matrix <- S_matrix\@snn[colnames(nmf.scores), colnames(nmf.scores)]
swne.embedding <- EmbedSWNE(nmf.scores, snn.matrix, alpha.exp = 1.0, snn.exp = 1, n_pull = 4, dist.use = \"IC\")



#write out dims and factors then combine and write combined dims
write.table(as.matrix(swne.embedding\$sample.coords),file=\"$opt{'O'}.Embedded.cells.$k.SWNE.dims\",quote=FALSE,sep=\"\\t\",row.names=TRUE,col.names=FALSE)

reform<-as.data.frame(cbind(swne.embedding\$H.coords\$x,swne.embedding\$H.coords\$y))
annot_form<-as.data.frame(cbind(swne.embedding\$H.coords\$name,rep(\"factor\",times=length(swne.embedding\$H.coords\$name))))

names(reform)<-c(\"x\",\"y\")
row.names(reform)<-swne.embedding\$H.coords\$name
out_comb<-rbind(as.matrix(swne.embedding\$sample.coords),as.matrix(reform))
write.table(as.matrix(reform),file=\"$opt{'O'}.Embedded.factors.$k.SWNE.dims\",quote=FALSE,sep=\"\\t\",row.names=TRUE,col.names=FALSE)
write.table(out_comb,file=\"$opt{'O'}.Embedded.cell.factors.$k.SWNE.dims\",quote=FALSE,sep=\"\\t\",row.names=TRUE,col.names=FALSE)
write.table(as.matrix(annot_form),file=\"$opt{'O'}.Embedded.factors.$k.SWNE.annot\",quote=FALSE,sep=\"\\t\",row.names=FALSE,col.names=FALSE)

## Associate factors with genes using the gene loadings (W) matrix
gene.loadings <- nmf.res\$W
gene.loadings <- t(apply(nmf.res\$W, 1, function(x) (x - min(x))/(max(x) - min(x))))



write.table(as.matrix(gene.loadings),file=\"$opt{'O'}.Embedded.factors.$k.loadings.SWNE.txt\",quote=FALSE,sep=\"\\t\",row.names=TRUE,col.names=TRUE)
top.factor.genes.df <- SummarizeAssocFeatures(gene.loadings, features.return = 1000)
write.table(as.matrix(top.factor.genes.df),file=\"$opt{'O'}.Embedded.factors.loadings.$k.top1000.loadings.SWNE.txt\",quote=FALSE,sep=\"\\t\",row.names=TRUE,col.names=TRUE)

## Make gene factor association heatmaps
gene.factor.heat <- gene.loadings[unique(top.factor.genes.df\$feature),]

pdf(\"$opt{'O'}.gene_factor_heatmap.$k.SWNE.pdf\", width = 7.5, height = 7.0)
ggHeat(gene.factor.heat, clustering = \"both\", x.lab.size = 14, y.lab.size = 5)
dev.off()


## Associate factors with cell clusters
clusters.list <- UnflattenGroups(clusters)
clusters.matrix <- t(swne:::.genesets_indicator(clusters.list, inv = F, return.numeric = T))
cluster.nmf.assoc <- FactorAssociation(clusters.matrix, nmf.scores, n.cores = n.cores, metric = \"IC\")

pdf(\"$opt{'O'}cluster_factor_heatmap.$k.SWNE.pdf\", width = 7.5, height = 4.5)
ggHeat(cluster.nmf.assoc, clustering = \"both\", x.lab.size = 14, y.lab.size = 14)
dev.off()

";
close R;

system("$Rscript $opt{'O'}.SWNE.r");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.SWNE.r");
}

#create new annot files with factors
if (defined $opt{'A'}) 
{

system("cat $opt{'A'} $opt{'O'}.Embedded.factors.$k.SWNE.annot > $opt{'O'}.cells.factors.$k.SWNE.annot");

};

#create new color annot files with factors
if (defined $opt{'C'}) 
{

system("cat $opt{'C'} \"factor\t\"black\"\" > $opt{'O'}.Embedded.factors.$k.colors.SWNE.annot");

};


exit;
}

if ($command eq "matrix-umap" || $command eq "umap" || $command eq "owl") {

# Defaults
$dims = 2;
$neigh = 30;
$mdist=0.1;
$metric="euclidean";

getopts("O:n:d:m:D:XP:", \%opt);

$die2 = "
scitools umap [options] [matrix]

Options:
   -O   [STR]   Output prefix (default is matrix file prefix)
   -D   [INT]   Dimensions to embed UMAP in (def = $dims)
   -n   [INT]   number of neighbors to use for UMAP high dim embedding (def = $neigh)
   -d   [INT]   min distance for mapping (def = $mdist)
   -m   [STR]   metric to use (def = $metric)
   -X           Retain intermediate files (def = delete)
   -P   [STR]   python script call (def = $Pscript)
   
Note: Requires python, numpy, and umap to be installed and callable
      This command is a wrapper for executing the python code.

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.matrix$//;
if (defined $opt{'D'}) {$dims = $opt{'D'}};
if (defined $opt{'n'}) {$neigh = $opt{'n'}};
if (defined $opt{'d'}) {$mdist = $opt{'d'}};
if (defined $opt{'m'}) {$metric = $opt{'m'}};
if (defined $opt{'P'}) {$Pscript = $opt{'P'}};
read_matrix($ARGV[0]);

open OUT, ">$opt{'O'}.UMAP.py";
print OUT"
import numpy
import umap
data_matrix=numpy.loadtxt(\"$ARGV[0]\",skiprows=1,usecols=range(1,$matrix_colNum))

fit = umap.UMAP(
       n_neighbors=$neigh,
       min_dist=$mdist,
       n_components=$dims,
       metric=\'$metric\'
   )
data=data_matrix.T
u = fit.fit_transform(data)
numpy.savetxt(\"$opt{'O'}.temp.UMAP.dims\",u,delimiter=\"\\t\")
";
close OUT;
system("$Pscript $opt{'O'}.UMAP.py");


open OUT, ">$opt{'O'}.UMAP.dims";
$counter=0;
open IN, "$opt{'O'}.temp.UMAP.dims"; 
while($l=<IN>)
{
$l =~ s/"//g;   
print OUT $MATRIX_COLNAMES[$counter]."\t".$l;
$counter++;
}
close(IN);
close OUT;



if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.temp.UMAP.dims $opt{'O'}.UMAP.py");
}


exit;
}

if ($command eq "matrix-aggregate" || $command eq "aggregate-matrix") {

use Getopt::Std; %opt = ();
getopts("O:r:c:a:b", \%opt);

$die2 = "
scitools matrix-aggregate [options] [counts matrix] [annotation file]
   or    aggregate-matrix

Note: if a non counts matrix is provided, it will still sum the values.
Support for comma-separated annotations (i.e. cell aggregation with
oversampling) has been added.

Options:
   -O   [STR]   Output prefix (default is [input annot].matrix)
   -c   [STR]   File with list of CellIDs to retain
   -r   [STR]   File with list of RowIDs to retain
   -a   [STR]   Comma separated list of annotations to include
   -b           Treat each cell as a binary value (resulting matrix is
                  the number of cells at the peak. def = total counts)
   
";

#name output 
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[1]; $opt{'O'} =~ s/\.matrix$//};


#read in annotation, scitools approach
if (!defined $ARGV[1]) {die $die2}
else {read_annot($ARGV[1])};

#define on your own which annot to include, scitools approach
if (defined $opt{'a'}) {
	@ANNOT_LIST = split(/,/, $opt{'a'});
	foreach $annot (@ANNOT_LIST) {
		$ANNOT_include{$annot} = 1;
	}
}

#define cells to keep
if (defined $opt{'c'}) {
	open IN, "$opt{'c'}";
	while ($cellID = <IN>) {
		chomp $cellID;
		$CELLID_list_include{$cellID} = 1;
	} close IN;
}

#define feature to keep
if (defined $opt{'r'}) {
	open IN, "$opt{'r'}";
	while ($rowID = <IN>) {
		chomp $rowID;
		$ROWID_list_include{$rowID} = 1;
	} close IN;
}

#read in matrix scitools standard
read_matrix($ARGV[0]);

#define your matrix based on which cells and features are allowed to be kept
#same style as CELL_FEATURE_VALUE
%ANNOT_FEATURE_aggregate_value;
#annot that we aggreagate by
@MATRIX_aggregate_COLNAMES;

#then the cells and aggregate signal
foreach $cellID (keys %CELLID_FEATURE_value) {
	#define cells to include based -c and -a options
	#start out with all Cells being included
	$includeCell = 1;
	#remove cells that are not in include list provided
	if (defined $opt{'c'} && !defined $CELLID_list_include{$cellID}) {$includeCell = 0};
	
	# handle multi-annot
	@ANNOTS = split(/,/, $CELLID_annot{$cellID});
	foreach $annot (@ANNOTS) {
		#remove cells that are not in annot list provided
		if (defined $opt{'a'} && !defined $ANNOT_include{$annot}) {$includeCell = 0};
		if ($includeCell>0) {
			foreach $rowID (keys %{$CELLID_FEATURE_value{$cellID}}) {
				$includeRow = 1;
				#define which rows to include based on -r option
				if (defined $opt{'r'} && !defined $ROWID_list_include{$rowID}) {$includeRow = 0};
				if ($includeRow>0 ) {
					if(!defined $opt{'b'}) {
						if (defined $ANNOT_FEATURE_aggregate_value{$annot}{$rowID}) {
							$ANNOT_FEATURE_aggregate_value{$annot}{$rowID}+=$CELLID_FEATURE_value{$cellID}{$rowID};
						} else {
							$ANNOT_FEATURE_aggregate_value{$annot}{$rowID}=$CELLID_FEATURE_value{$cellID}{$rowID};
						}
					} else {
						if ($CELLID_FEATURE_value{$cellID}{$rowID}>0) {
							if (defined $ANNOT_FEATURE_aggregate_value{$annot}{$rowID}) {
								$ANNOT_FEATURE_aggregate_value{$annot}{$rowID}++;
							}
						} else {
							$ANNOT_FEATURE_aggregate_value{$annot}{$rowID}=1;
						}
					}
				}  
			}
			$included_cells++;
			$CELLID_include{$cellID} = 1;
			#define annotations to include
			$ANNOT_include{$annot} = 1;  
		}
	}
}

#create a matrix of included annot
foreach $annot (keys %ANNOT_include) {
	push(@MATRIX_aggregate_COLNAMES,$annot);
}

#then look which rows are allowed to be kept
foreach $rowID (keys %MATRIX_feature_nonZero) {
	if (!defined $opt{'r'} || defined $ROWID_list_include{$rowID}) {
		$ROWID_count{$rowID} = 0;
		foreach $cellID (keys %CELLID_include) {
			if ($CELLID_FEATURE_value{$cellID}{$rowID}>0) {
				$ROWID_count{$rowID}++;
			}
		}
		$ROWID_include{$rowID} = 1;
		$included_rows++;
	}       
}

#do write out depending on b option
if(!defined $opt{'b'}) {
	open OUT, ">$opt{'O'}.aggregate.matrix";
	open LOG, ">$opt{'O'}.aggregate.matrix.log";
} else {
	open OUT, ">$opt{'O'}.binary.aggregate.matrix";
	open LOG, ">$opt{'O'}.binary.aggregate.matrix.log";
}

#first define header
$out_header = join("\t",@MATRIX_aggregate_COLNAMES);

print OUT "$out_header\n";

for ($rowNum = 0; $rowNum < @MATRIX_ROWNAMES; $rowNum++) {
	$rowID = $MATRIX_ROWNAMES[$rowNum];
	if (defined $ROWID_include{$rowID}) {
		print OUT "$rowID";
		for ($aggrNum = 0; $aggrNum < @MATRIX_aggregate_COLNAMES; $aggrNum++) {
				print OUT "\t$ANNOT_FEATURE_aggregate_value{$MATRIX_aggregate_COLNAMES[$aggrNum]}{$rowID}";
		} print OUT "\n";
	}
} close OUT;

$ts = localtime(time);
print LOG "$ts scitools annot-aggr
Matrix file = $ARGV[1]
Annotation file = $ARGV[0]

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

exit;
}

if ($command eq "matrix-merge" || $command eq "merge-matrix") {

$joiner = "_";

getopts("O:U", \%opt);

$die2 = "
scitools matrix-merge [options] [matrix1] [matrix2] (optional matrix N etc...)
   or    merge-matrix

Will merge matrices if names of cells do not match (others can be further developed).
Features will be the intersect across the matrices.

Options:
   -O   [STR]   Output file name / prefix (def = matrix1 prefix w/ merge)
   -U           Print the union of all features (those not present in a
                 matrix will be set to 0)

";

if (!defined $ARGV[1]) {die $die2};

if (!defined $opt{'O'}) {
$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.counts.matrix$//; $opt{'O'} =~ s/\.binary.matrix$//; $opt{'O'} =~ s/\.matrix$//;
};
#
if ($ARGV[0] =~ /\.matrix$/) 
{
if ($ARGV[0] =~ /\.counts.matrix$/) {$opt{'O'} .= ".merged.counts";} elsif ($ARGV[0] =~ /\.binary.matrix$/) {$opt{'O'} .= ".merged.binary";} else {print "Warning: Type of matrix not in name\n"; $opt{'O'} .= ".merged";};
} else {print "ERROR: This is not a matrix, rename\n"; die};

#create a merged matrix with CELLID_FEATURE_value structure
print "reading matrix $ARGV[0]\n";
read_matrix($ARGV[0]);
%CELLID_FEATURE_merged_value=%CELLID_FEATURE_value;
#create an included rowid hash 
foreach $cellID (keys %CELLID_FEATURE_value) 
{
	foreach $rowID (keys %{$CELLID_FEATURE_value{$cellID}}) 
	{	
	$All_included_rowid{$rowID}=1;
	}
}
#create an all cells Array to preserve order. 
@MATRIX_ALL_COLNAMES=@MATRIX_COLNAMES;

for ($matrixID = 1; $matrixID < @ARGV; $matrixID++) 
{
	print "reading matrix $ARGV[$matrixID]\n";
	read_matrix($ARGV[$matrixID]);
	
#quick check for cellID overlaps (faster and can change to warning later)
	foreach $cellID (keys %CELLID_FEATURE_value) 
	{
		if (defined $CELLID_FEATURE_merged_value{$cellID})
		{
		print "ERROR: Overlapping cellIDs, renaming not yet developed\nThis is the cell name $cellID\n"; die
		}
	}
	#if no overlap continue
	foreach $cellID (keys %CELLID_FEATURE_value) 
	{
		if (!defined $CELLID_FEATURE_merged_value{$cellID})
		{
			foreach $rowID (keys %{$CELLID_FEATURE_value{$cellID}}) 
			{
				if (defined $All_included_rowid{$rowID})
				{# add values to the merged matrix
				$CELLID_FEATURE_merged_value{$cellID}{$rowID}=$CELLID_FEATURE_value{$cellID}{$rowID};
				}
				else
				{# if not defined, define and add values
				$CELLID_FEATURE_merged_value{$cellID}{$rowID}=$CELLID_FEATURE_value{$cellID}{$rowID};
				$All_included_rowid{$rowID}=1;
				}
			}	
		}
		else {print "ERROR: Overlapping cellIDs, renaming not yet developed\nThis is the cell name $cellID\n"; die};
		#can add part that handles overlapping IDs
		#need to consider two cases: 1. cellIDs are same because different capacity run on same cells, 2. CellIDs are same but different cells, need to rename
	}
	#add new cellnames
	push(@MATRIX_ALL_COLNAMES,@MATRIX_COLNAMES);
	
}

open OUT, ">$opt{'O'}.matrix";



print OUT join("\t",@MATRIX_ALL_COLNAMES)."\n";
foreach $rowID (sort keys %All_included_rowid) 
	{
	#UNION vs INTERSECT
	$in_all = 1;
	#OUTPUT hash
	@OUTPUT_line=();
	push(@OUTPUT_line,$rowID);
		foreach $cellID (@MATRIX_ALL_COLNAMES) 
		{
			if (defined $CELLID_FEATURE_merged_value{$cellID}{$rowID})
			{#if defined input value
			push(@OUTPUT_line,$CELLID_FEATURE_merged_value{$cellID}{$rowID});
			}
			else
			{#if not defined input  0
			push(@OUTPUT_line,0);
			$in_all = -1;
			}
		}
		if (defined $opt{'U'} || $in_all>0)
		{
			print OUT join("\t",@OUTPUT_line)."\n";
		}
	}
close OUT;

exit;
}


########## DIMENSION FILE FUNCTIONS ##########
if ($command eq "dims-kmeans" || $command eq "kmeans") {

# Defaults
$range_default = "1-15";
$xdim = 1;
$ydim = 2;

getopts("O:XR:K:D:s:P:x:y:p", \%opt);

$die2 = "
scitools dims-kmeans [options] [input dims]
   or    kmeans

Performs kmeans clustering on a dims file

Options:
   -K   [INT]   K-value (number of clusters, Required)
   -O   [STR]   Output prefix (default is [input].Kmeans_K[K].annot)
   -D   [RNG]   Range of dimensions to include (range format)
                (e.g. 1-5,6,8; def = $range_default)
   -X           Retain intermediate files (def = delete)
   -R   [STR]   Rscript call (def = $Rscript)
   -p           Plot the output using specified dims file
   -P   [DIMS]  Plot the resulting Kmeans clustering using the
                specified dims file (for more plot options use
                scitools plot-dims)
   -x   [INT]   X-dimension to plot (def = $xdim)
   -y   [INT]   Y-dimension to plot (def = $ydim)
   -s   [STR]   scitools call (def = $scitools)

";

if (!defined $ARGV[0] | !defined $opt{'K'}) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (defined $opt{'x'}) {$xdim = $opt{'x'}};
if (defined $opt{'y'}) {$ydim = $opt{'y'}};
if (!defined $opt{'D'}) {$opt{'D'} = $range_default};
if (defined $opt{'p'} && !defined $opt{'P'}) {$opt{'P'} = $ARGV[0]};
read_ranges($opt{'D'});

open IN, "$ARGV[0]";
$dim_line = <IN>; close IN;
chomp $dim_line; @DL = split(/\t/, $dim_line);
if (@DL < $RANGE_VALUES[@RANGE_VALUES-1]) {
	die "ERROR: The dimension ranges (max = $RANGE_VALUES[@RANGE_VALUES-1]) specified are beyond the number of dimensions in $ARGV[0] (".@DL.")!\n";
}

open R, ">$opt{'O'}.Kmeans_K$opt{'K'}.r";
print R "
D<-read.table(\"$ARGV[0]\",row.names=1)[,$range_R_set]
K<-kmeans(D,$opt{'K'})
ANN<-as.matrix(K\$cluster)
ANN[,1]<-sub(\"\^\", \"K\", ANN[,1])
write.table(ANN,file=\"$opt{'O'}.Kmeans_K$opt{'K'}.annot\",col.names=FALSE,row.names=TRUE,quote=FALSE,sep=\"\\t\")
"; close R;

system("$Rscript $opt{'O'}.Kmeans_K$opt{'K'}.r");

if (defined $opt{'P'}) {
	print STDERR "PLOTTING: $scitools plot-dims -O $opt{'O'}.Kmeans_K$opt{'K'} -A $opt{'O'}.Kmeans_K$opt{'K'}.annot -x $xdim -y $ydim $opt{'P'}\n";
	system("$scitools plot-dims -O $opt{'O'}.Kmeans_K$opt{'K'} -A $opt{'O'}.Kmeans_K$opt{'K'}.annot -x $xdim -y $ydim $opt{'P'}");
}

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.Kmeans_K$opt{'K'}.r");
}

exit;
}

if ($command eq "dims-dbscan" || $command eq "dbscan") {

# Defaults
$range_default = "1-2";
$xdim = 1;
$ydim = 2;
$epsilon = 2;
$minPts = 10;

getopts("O:XR:D:s:P:x:y:e:m:p", \%opt);

$die2 = "
scitools dims-dbscan [options] [input dims]
   or    dbscan

Performs dbscan clustering on a dims file

Options:
   -O   [STR]   Output prefix (default is [input].dbscan.annot)
   -e   [FLT]   Epsilon value, akin to distance (def = $epsilon)
   -m   [INT]   Minimum points to seed a cluster (def = $minPts)
   -D   [RNG]   Range of dimensions to include (range format)
                (e.g. 1-5,6,8; def = $range_default)
   -X           Retain intermediate files (def = delete)
   -R   [STR]   Rscript call (def = $Rscript)
   -p           Plot on the dims file (-P will toggle)
   -P   [DIMS]  Plot the resulting dbscan clustering using the
                specified dims file (for more plot options use
                scitools plot-dims)
   -x   [INT]   X-dimension to plot (def = $xdim)
   -y   [INT]   Y-dimension to plot (def = $ydim)
   -s   [STR]   scitools call (def = $scitools)

Requires the dbscan R package.

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (defined $opt{'x'}) {$xdim = $opt{'x'}};
if (defined $opt{'y'}) {$ydim = $opt{'y'}};
if (defined $opt{'e'}) {$epsilon = $opt{'e'}};
if (defined $opt{'m'}) {$minPts = $opt{'m'}};
if (!defined $opt{'D'}) {$opt{'D'} = $range_default};
if (defined $opt{'p'} && !defined $opt{'P'}) {$opt{'P'} = $ARGV[0]};
read_ranges($opt{'D'});

open IN, "$ARGV[0]";
$dim_line = <IN>; close IN;
chomp $dim_line; @DL = split(/\t/, $dim_line);
if (@DL < $RANGE_VALUES[@RANGE_VALUES-1]) {
	die "ERROR: The dimension ranges (max = $RANGE_VALUES[@RANGE_VALUES-1]) specified are beyond the number of dimensions in $ARGV[0] (".@DL.")!\n";
}

open R, ">$opt{'O'}.dbscan.r";
print R "
library(dbscan)
DIMS<-as.matrix(read.table(\"$ARGV[0]\",row.names=1)[,$range_R_set])
DIST<-dist(DIMS)
DBS<-dbscan(DIST,$epsilon,minPts = $minPts)
ANN<-as.matrix(DBS\$cluster)
rownames(ANN)<-rownames(DIMS)
ANN[,1]<-sub(\"\^\", \"D\", ANN[,1])
write.table(ANN,file=\"$opt{'O'}.dbscan.annot\",col.names=FALSE,row.names=TRUE,sep=\"\\t\",quote=FALSE)
"; close R;

system("$Rscript $opt{'O'}.dbscan.r");

if (defined $opt{'P'}) {
	print STDERR "PLOTTING: $scitools plot-dims -O $opt{'O'}.dbscan -A $opt{'O'}.dbscan.annot -x $xdim -y $ydim $opt{'P'}\n";
	system("$scitools plot-dims -O $opt{'O'}.dbscan -A $opt{'O'}.dbscan.annot -x $xdim -y $ydim $opt{'P'}");
}

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.dbscan.r");
}

exit;
}

if ($command eq "aggregate-cells" || $command eq "aggregate") {

# Defaults
$range_default = "1-2";
$oversample = 1;
$xdim = 1;
$ydim = 2;
$aggN = 15;
$minN = 5;
$theme = "Clean";

getopts("O:D:N:A:a:C:c:x:y:R:XK:T:n:v:", \%opt);

$die2 = "
scitools dims-aggregate [options] [input dims/lambda/values]
   or    aggregate

Aggregates cells via k-means or similar. Produces a merged dims for the
new aggregate cells, and a merge annotation file (can merge matrix
using this new annotation file)

Options:
   -O   [STR]   Output prefix (default is [input].aggregate.annot)
   -D   [RNG]   Range of dimensions to aggregate from (range format)
                (e.g. 1-5,6,8; def = $range_default; 1 if lambda input)
   -N   [INT]   Target number of cells to aggregate (def = $aggN)
   -n   [INT]   Min number of cells in cluster (def = $minN)
   -v   [FLT]   Oversampling value (ie mean assignments for each cell)
                Note: this will create an annot file with comma separated
                annotation memberships if >1. (def = $oversample)
   -K   [INT]   Number of clusters (overrides -N to be the min/aggregate)
                This will utilize standard k-means to aggregate cells
                For lambda clustering it will evenly space them
   -A   [STR]   Annotation file (will only merge within annot)
   -a   [STR]   Comma separated list of annotations to include
                  (requires -A to be specified)

Plotting Options (def is input dims file coordinates - only for dims files):
   -C   [STR]   Color coding file (annot (tab) #hexColor)
   -c   [STR]   Color coding string
                  Annot=#hexColor,Annot2=#hexColor
   -x   [INT]   X-dimension to plot (def = $xdim)
   -y   [INT]   Y-dimension to plot (def = $ydim)
   -T   [STR]   Theme: (def = $theme)
                  Clean = no axis lines (for tSNE)
                  Reg = regular, include axes and grid (for PCA)
   -R   [STR]   Rscript call (def = $Rscript)
   -X           Retain intermediate files (def = delete)

";

#die "ERROR: This function is under development and cannot be used in its current form!\n";

if (!defined $ARGV[0]) {die $die2};
if (defined $opt{'a'} && !defined $opt{'A'}) {die "\nMust provide an annotaiton file (-A) if specifying annotations to plot (-a)!\n$die2"};
if (defined $opt{'C'} && defined $opt{'c'}) {die "\nSpecify either a color string (-c) or a color coding file (-C), not both!\n$die2"};
if (!defined $opt{'D'}) {$opt{'D'} = $range_default};
if (defined $opt{'N'}) {$aggN = $opt{'N'}};
if (defined $opt{'n'}) {$minN = $opt{'n'}};
if (defined $opt{'T'}) {$theme = $opt{'T'}};
if (defined $opt{'v'}) {$oversample = $opt{'v'}};

if ($ARGV[0] =~ /\.(lambda|values)$/) {
	print STDERR "INFO: A lambda or values file was detected - will aggregate in one dimension\n";
	if (!defined $opt{'D'}) {$opt{'D'} = 1};
	read_values($ARGV[0]);
	if (defined $opt{'A'}) {
		read_annot($opt{'A'});
	} else {
		foreach $cellID (keys %CELLID_value) {
			$CELLID_annot{$cellID} = "Agg";
		}
		$ANNOT_include{"Agg"} = 1;
	}
} else {
	read_dims($ARGV[0]);
	if (defined $opt{'A'}) {
		read_annot($opt{'A'});
	} else {
		foreach $cellID (keys %CELLID_DIMS) {
			$CELLID_annot{$cellID} = "Agg";
		}
		$ANNOT_include{"Agg"} = 1;
	}
}

if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.(dims|lambda|values)$//;

if (defined $opt{'R'}) {$Rscript = $opt{'R'}};

if (defined $opt{'a'}) {
	@ANNOT_LIST = split(/,/, $opt{'a'});
	foreach $annot (@ANNOT_LIST) {
		$ANNOT_include{$annot} = 1;
	}
} else {
	foreach $annot (keys %ANNOT_count) {
		$ANNOT_include{$annot} = 1;
	}
}

if (defined $opt{'C'}) {read_color_file($opt{'C'})};
if (defined $opt{'c'}) {read_color_string($opt{'c'})};

if (defined $opt{'x'}) {$xdim = $opt{'x'}};
if (defined $opt{'y'}) {$ydim = $opt{'y'}};
$xpos = $xdim-1; $ypos = $ydim-1;

read_ranges($opt{'D'});

if ($ARGV[0] =~ /\.(lambda|values)$/) {
	# figure out even spaced vs. even cell N
	# if -K, then will be evenly spaced on lambda - uniform for all annotations
	# elsif -N, then will be even number of cells - varable for annotations
	if (defined $opt{'K'}) { # find centroids that will be constant across all annotations
		%CLUST_global_centroid = ();
		$clust_span = $value_range/$opt{'K'};
		$center = $clust_span/2;
		$clustID = 0;
		$CLUST_global_centroid{$clustID} = $center;
		for ($clustID = 1; $clustID < $opt{'K'}; $clustID++) {
			$center += $clust_span;
			$CLUST_global_centroid{$clustID} = $center;
		}
	}
}

open CNT, ">$opt{'O'}.centroids.dims";
open OUT, ">$opt{'O'}.aggregate.annot";

# loop through annotations individually
foreach $annot (keys %ANNOT_include) {
if ($ARGV[0] =~ /\.(lambda|values)$/) {
	# lambda file
	# setup new CLUST center hash for the annot - all can be annot-only, no need for global to be kept
	%CLUST_center = (); $assignment_count = 0; %CLUST_assignments = (); %CELLID_initial_cluster = (); %CLUST_cellIDs = ();
	# setup cluster centroids based on K or N
	if (defined $opt{'K'}) { # copy over the clusters and add in the annot
		foreach $clustID (keys %CLUST_global_centroid) {
			$CLUST_center{$annot."_".$clustID} = $CLUST_global_centroid{$clustID};
		}
		# now do initial assignments based on closest center
		foreach $cellID (keys %CELLID_value) {
			if ($CELLID_annot{$cellID} eq $annot) {
				$minDist = 1e9; $winAsn = "NA";
				foreach $clustID (keys %CLUST_center) {
					$dist = abs($CELLID_value{$cellID} - $CLUST_center{$clustID});
					if ($dist<$minDist) {
						$minDist = $dist;
						$winAsn = $clustID;
					}
				}
				$CELLID_initial_cluster{$cellID} = $winAsn;
				push @{$CLUST_cellIDs{$winAsn}}, $cellID;
				$CLUST_assignments{$winAsn}++;
				$assignment_count++;
			}
		}
		# now check for reaching min n, remove other clusters and do orphan assignment
		%CELLID_orphan = ();
		foreach $clustID (keys %CLUST_center) {
			if ($CLUST_assignments{$clustID} < $minN) {
				print STDERR "INFO: Cluster $clustID has $CLUST_assignments{$clustID} cells assigned which is < $minN, and will be removed and have cells re-assigned.\n";
				delete $CLUST_center{$clustID};
				foreach $cellID (@{$CLUST_cellIDs{$clustID}}) {
					$CELLID_orphan{$cellID} = 1;
				}
				delete $CLUST_cellIDs{$clustID};
			}
		}
		foreach $cellID (keys %CELLID_orphan) {
			$minDist = 1e9; $winAsn = "NA";
			foreach $clustID (keys %CLUST_center) {
				$dist = abs($CELLID_value{$cellID} - $CLUST_center{$clustID});
				if ($dist<$minDist) {
					$minDist = $dist;
					$winAsn = $clustID;
				}
			}
			$CELLID_initial_cluster{$cellID} = $winAsn;
			push @{$CLUST_cellIDs{$winAsn}}, $cellID;
			$CLUST_assignments{$winAsn}++;
		}
	} else {
		# go through cells and calculate the center for each grouping - also do initial assignment
		$clustNum = 0; $clust_sum = 0; $clust_memberCT = 0; $clustID = $annot."_".$clustNum; $annot_cellCT_in_lambda = 0;
		foreach $cellID (keys %CELLID_value) {
			if ($CELLID_annot{$cellID} eq $annot) {
				$annot_cellCT_in_lambda++;
			}
		}
		$cluster_ct = int($annot_cellCT_in_lambda/$aggN);
		if ($cluster_ct < 2) {
			$cluster_ct = 2;
			print STDERR "WARNING: The N value specified, $aggN, results in cluster numbers less than 2 (cells in $annot = $annot_cellCT_in_lambda), setting to 2\n";
		}
		$target_assignments = int($annot_cellCT_in_lambda/$cluster_ct);
		foreach $cellID (sort {$CELLID_value{$a}<=>$CELLID_value{$b}} keys %CELLID_value) {
			if ($CELLID_annot{$cellID} eq $annot) {
				$clust_memberCT++;
				$clust_sum += $CELLID_value{$cellID};
				$CELLID_initial_cluster{$cellID} = $clustID;
				push @{$CLUST_cellIDs{$clustID}}, $cellID;
				$assignment_count++;
				if ($clust_memberCT >= $target_assignments) {
					$center = $clust_sum/$clust_memberCT;
					$CLUST_center{$clustID} = $center;
					$clustNum++; $clust_sum = 0; $clust_memberCT = 0; $clustID = $annot."_".$clustNum; 
				}
			}
		}
		# last one
		$center = $clust_sum/$clust_memberCT;
		$CLUST_center{$clustID} = $center;
	}
	
	# now check if there is oversampling
	if ($oversample>1) {
		%CELLID_cluster = ();
		if (defined $opt{'K'}) { # oversample for each cluster
			foreach $clustID (keys %CLUST_center) {
				$target_assignments = int(($CLUST_assignments{$clustID}*$oversample)+1);
				%CELLID_dist = ();
				foreach $cellID (keys %CELLID_value) {
					if ($CELLID_annot{$cellID} eq $annot && $CELLID_initial_cluster{$cellID} ne $clustID) {
						$CELLID_dist{$cellID} = abs($CELLID_value{$cellID} - $CLUST_center{$clustID});
					} elsif ($CELLID_annot{$cellID} eq $annot && $CELLID_initial_cluster{$cellID} eq $clustID) {
						$CELLID_cluster{$cellID} .= "$clustID,";
					}
				}
				foreach $cellID (sort {$CELLID_dist{$a}<=>$CELLID_dist{$b}} keys %CELLID_dist) {
					if ($CLUST_assignments{$clustID}<$target_assignments) {
						$CELLID_cluster{$cellID} .= "$clustID,";
						$CLUST_assignments{$clustID}++;
						push @{$CLUST_cellIDs{$clustID}}, $cellID;
					}
				}
			}
			foreach $clustID (keys %CLUST_center) { # no center changes, so print centroids
				print CNT "$clustID\t$CLUST_center{$clustID}\n";
			}
		} else {
			$added_assignments = ($annot_cellCT_in_lambda*($oversample-1));
			$target_assignments += int(($added_assignments/$clustNum)+1);
			foreach $clustID (keys %CLUST_center) {
				%CELLID_dist = ();
				foreach $cellID (keys %CELLID_value) {
					if ($CELLID_annot{$cellID} eq $annot && $CELLID_initial_cluster{$cellID} ne $clustID) {
						$CELLID_dist{$cellID} = abs($CELLID_value{$cellID} - $CLUST_center{$clustID});
					} elsif ($CELLID_annot{$cellID} eq $annot && $CELLID_initial_cluster{$cellID} eq $clustID) {
						$CELLID_cluster{$cellID} .= "$clustID,";
					}
				}
				foreach $cellID (sort {$CELLID_dist{$a}<=>$CELLID_dist{$b}} keys %CELLID_dist) {
					if ($CLUST_assignments{$clustID}<$target_assignments) {
						$CELLID_cluster{$cellID} .= "$clustID,";
						$CLUST_assignments{$clustID}++;
						push @{$CLUST_cellIDs{$clustID}}, $cellID;
					}
				}
			}
			# remake centroid & print
			foreach $clustID (keys %CLUST_center) {
				$clust_sum = 0; $clust_memberCT = 0;
				foreach $cellID (@{$CLUST_cellIDs{$clustID}}) {
					$clust_sum+=$CELLID_value{$cellID};
					$clust_memberCT++;
				}
				$CLUST_center{$clustID} = $clust_sum/$clust_memberCT;
				print CNT "$clustID\t$CLUST_center{$clustID}\n";
			}
		}
		# print out annot
		foreach $cellID (keys %CELLID_cluster) {
			$cluster_set = $CELLID_cluster{$cellID};
			$cluster_set =~ s/,$//;
			print OUT "$cellID\t$cluster_set\n";
		}
	} else { # no oversampling
		foreach $clustID (keys %CLUST_center) {
			print CNT "$clustID\t$CLUST_center{$clustID}\n";
			foreach $cellID (@{$CLUST_cellIDs{$clustID}}) {
				print OUT "$cellID\t$clustID\n";
			}
		}
	}
	
} else {
	# Dimensions file
	
	%CLUST_include = ();
	%CELLID_initial_cluster = ();
	%CLUST_cellCT = ();
	$ANNOT_count{$annot} = 0;
	
	foreach $cellID (keys %CELLID_DIMS) {
		if ($CELLID_annot{$cellID} eq $annot) {$ANNOT_count{$annot}++};
	}
	
	if (defined $opt{'K'}) {
		$K = $opt{'K'};
		$aggN = int($ANNOT_count{$annot}/$K);
	} else {
		$K = int(($ANNOT_count{$annot}/$aggN)+0.5);
		if ($K<1) {$K=1};
	}

	open AO, ">$opt{'O'}.$annot.agg_dims";
	foreach $cellID (keys %CELLID_DIMS) {
		if ($CELLID_annot{$cellID} eq $annot) {
			print AO "$cellID";
			for ($dimID = 0; $dimID < @RANGE_VALUES; $dimID++) {
				print AO "\t$CELLID_DIMS{$cellID}[$RANGE_VALUES[$dimID]]";
			}
			print AO "\n";
		}
	} close AO;
	open R, ">$opt{'O'}.$annot.agg_dims.r";
	if ($ARGV[0] =~ /\.(lambda|values)$/) {
	print R "
D<-read.table(\"$opt{'O'}.$annot.agg_dims\",row.names=1)";
	} else {
	print R "
D<-read.table(\"$opt{'O'}.$annot.agg_dims\",row.names=1)[,$range_R_set]";
	}
	print R"
K<-kmeans(D,$K)
ANN<-as.matrix(K\$cluster)
ANN[,1]<-sub(\"\^\", \"K\", ANN[,1])
write.table(ANN,file=\"$opt{'O'}.$annot.agg_dims.annot\",col.names=FALSE,row.names=TRUE,quote=FALSE,sep=\"\\t\")
"; close R;
	system("$Rscript $opt{'O'}.$annot.agg_dims.r 2>/dev/null");

	# find passing clusters (needed for separation)
	open IN, "$opt{'O'}.$annot.agg_dims.annot";
	while ($l = <IN>) {
		chomp $l;
		($cellID,$clust_raw) = split(/\t/, $l);
		$cluster = "$annot\_$clust_raw";
		$CLUST_cellCT{$cluster}++;
	} close IN;
	foreach $cluster (keys %CLUST_cellCT) {
		if ($CLUST_cellCT{$cluster}>=$minN) {
			$CLUST_include{$cluster} = 1;
		}
	}

	# build assignments
	open IN, "$opt{'O'}.$annot.agg_dims.annot";
	while ($l = <IN>) {
		chomp $l;
		($cellID,$clust_raw) = split(/\t/, $l);
		$cluster = "$annot\_$clust_raw";
		if (defined $CLUST_include{$cluster}) {
			$CELLID_initial_cluster{$cellID} = "$cluster";
		}
		push @{$CLUSTER_CELLS{$cluster}}, $cellID;
	} close IN;
	
	# now compute centroids
	foreach $cluster (keys %CLUST_include) {
		@DIM_SUMS = ();
		@{$CLUSTER_centroidDims{$cluster}} = ();
		foreach $cellID (@{$CLUSTER_CELLS{$cluster}}) {
			for ($dimID = 0; $dimID < @RANGE_VALUES; $dimID++) {
				$DIM_SUMS[$dimID] += $CELLID_DIMS{$cellID}[$RANGE_VALUES[$dimID]];
			}
		}
		for ($dimID = 0; $dimID < @DIM_SUMS; $dimID++) {
			$mean_dim = $DIM_SUMS[$dimID]/$CLUST_cellCT{$cluster};
			$CLUSTER_centroidDims{$cluster}[$dimID] = $mean_dim;
		}
	}
	
	# now assign cells based on the centroids & factor in oversampling
	if ($oversample>1) {
		foreach $cluster (keys %CLUST_include) {
			# calculate cell to centroid distances
			%CELLID_distance = ();
			foreach $cellID (keys %CELLID_initial_cluster) {
				$dist_sum = 0;
				for ($dimID = 0; $dimID < @{$CLUSTER_centroidDims{$cluster}}; $dimID++) {
					$dist_sum += ($CLUSTER_centroidDims{$cluster}[$dimID] - $CELLID_DIMS{$cellID}[$RANGE_VALUES[$dimID]])**2;
				}
				$CELLID_distance{$cellID} = sqrt($dist_sum);
			}
			
			# now calculate target number of cells in cluster
			if ($CLUST_cellCT{$cluster}>=($oversample*$aggN)) {
				$CLUST_memberCT{$cluster} = $CLUST_cellCT{$cluster};
			} else {
				$CLUST_memberCT{$cluster} = int($oversample*$aggN);
			}
			
			# now assign & recompute centroids
			$added = 0; @DIM_SUMS = ();
			foreach $cellID (sort {$CELLID_distance{$a}<=>$CELLID_distance{$b}} keys %CELLID_distance) {
				if ($added < $CLUST_memberCT{$cluster}) {
					for ($dimID = 0; $dimID < @RANGE_VALUES; $dimID++) {
						$DIM_SUMS[$dimID] += $CELLID_DIMS{$cellID}[$RANGE_VALUES[$dimID]];
					}
					$CELLID_cluster{$cellID} .= "$cluster,";
					$added++;
				}
			}
			for ($dimID = 0; $dimID < @DIM_SUMS; $dimID++) {
				$mean_dim = $DIM_SUMS[$dimID]/$CLUST_memberCT{$cluster};
				$CLUSTER_centroidDims{$cluster}[$dimID] = $mean_dim;
			}
		}
	} else {
		foreach $cellID (keys %CELLID_initial_cluster) {
			$CELLID_cluster{$cellID} = $CELLID_initial_cluster{$cellID};
		}
	}
	
	foreach $cluster (keys %CLUST_include) {
		print CNT "$cluster";
		for ($dimID = 0; $dimID < @{$CLUSTER_centroidDims{$cluster}}; $dimID++) {
			print CNT "\t$CLUSTER_centroidDims{$cluster}[$dimID]";
		} print CNT "\n";
	}
	foreach $cellID (keys %CELLID_initial_cluster) {
		$CELLID_cluster{$cellID} =~ s/,$//;
		if ($CELLID_cluster{$cellID} ne "") {
			print OUT "$cellID\t$CELLID_cluster{$cellID}\n";
		}
	}
	
	if (!defined $opt{'X'}) {
		system("rm -f $opt{'O'}.$annot.agg_dims.annot $opt{'O'}.$annot.agg_dims.r $opt{'O'}.$annot.agg_dims");
	}
}
}
close CNT; close OUT;

# Plot the projections of cells to their centroids
if ($ARGV[0] =~ /\.(lambda|values)$/) {
	# labmda data plots?
	# NOT CURRENTLY SUPPORTED
} else {
# Dimensions data & plot
open PD, ">$opt{'O'}.aggregate_data";
foreach $cellID (keys %CELLID_cluster) {
	@CLUSTER_SET = split(/,/, $CELLID_cluster{$cellID});
	foreach $cluster (@CLUSTER_SET) {
		print PD "$cellID:$cluster\t$CELLID_annot{$cellID}\t$CELLID_DIMS{$cellID}[$xdim]\t$CLUSTER_centroidDims{$cluster}[$xpos]\t$CELLID_DIMS{$cellID}[$ydim]\t$CLUSTER_centroidDims{$cluster}[$ypos]\n";
	}
} close PD;

open R, ">$opt{'O'}.aggregate_data.r";
print R "
library(ggplot2)
CELLS<-read.table(\"$opt{'O'}.aggregate_data\",row.names=1)
colnames(CELLS)<-c(\"annot\",\"xdim\",\"centx\",\"ydim\",\"centy\")
PLT<-ggplot() +";

if (!defined $opt{'c'} && !defined $opt{'C'} && !defined $opt{'A'}) {
	print R "
	geom_segment(aes(x=CELLS\$xdim,xend=CELLS\$centx,y=CELLS\$ydim,yend=CELLS\$centy),size=0.5,color=\"lightsteelblue4\",alpha=0.25) +
	geom_point(aes(CELLS\$xdim,CELLS\$ydim),color=\"lightsteelblue4\",size=0.5) +
	geom_point(aes(CELLS\$centx,CELLS\$centy),color=\"gray20\",size=1.75) +
	geom_point(aes(CELLS\$centx,CELLS\$centy),color=\"lightsteelblue4\",size=1) +";
} elsif (!defined $opt{'V'}) {
	print R "
	geom_segment(aes(x=CELLS\$xdim,xend=CELLS\$centx,y=CELLS\$ydim,yend=CELLS\$centy,color=CELLS\$annot),alpha=0.25,size=0.5) +
	geom_point(aes(CELLS\$xdim,CELLS\$ydim,color=CELLS\$annot),size=0.5) +
	geom_point(aes(CELLS\$centx,CELLS\$centy),color=\"gray20\",size=1.75) +
	geom_point(aes(CELLS\$centx,CELLS\$centy,color=CELLS\$annot),size=1) +
	guides(colour = guide_legend(override.aes = list(size=4))) +";
}

if ($color_mapping !~ /none/i) {
	print R "
	scale_colour_manual(values = c($color_mapping)) +";
}

if ($theme =~ /Clean/i) {
	print R "
	theme_bw() +
	theme(panel.border=element_blank(),
		  panel.grid=element_blank(),
		  axis.line=element_blank(),
		  axis.ticks=element_blank(),
		  legend.background=element_blank(),
		  legend.title=element_blank(),
		  panel.background=element_blank(),
		  axis.text=element_blank(),
		  axis.title.x=element_blank(),
		  axis.title.y=element_blank(),
		  plot.background=element_blank(),
		  plot.margin=unit(c(0,0,0,0),\"pt\"))\n";
} else {
	print R "
	theme_bw()\n";
}

print R "
ggsave(plot=PLT,filename=\"$opt{'O'}.aggregate.png\",width=5,height=4,dpi=900)
ggsave(plot=PLT,filename=\"$opt{'O'}.aggregate.pdf\",width=6,height=5)
";
close R;

system("$Rscript $opt{'O'}.aggregate_data.r");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.aggregate_data $opt{'O'}.aggregate_data.r");
}
}

exit;
}

if ($command eq "dims-pcurve" || $command eq "pcurve") {

# Defaults
$range_default = "1-2";

getopts("O:D:XR:", \%opt);

$die2 = "
scitools dims-pcurve [options] [input dims]
   or    pcurve

Projects a principle curve through the dimensions specified

Options:
   -O   [STR]   Output prefix (default is [input].pcurve)
   -D   [RNG]   Range of dimensions to include (range format)
                (e.g. 1-5,6,8; def = $range_default)
   -R   [STR]   Rscript call (def = $Rscript)
   -X           Retain intermediate files (def = delete)

Requires the princurve R package.

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.dims$//; $opt{'O'} =~ s/\.txt$//};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (!defined $opt{'D'}) {$opt{'D'} = $range_default};
read_ranges($opt{'D'});

open R, ">$opt{'O'}.pcurve.r";
print R "
library(princurve)
IN <- as.matrix(read.table(\"$ARGV[0]\",row.names=1))
PCURVE <- principal.curve(IN[,$range_R_set])
write.table(PCURVE\$s,file=\"$opt{'O'}.pcurve.proj.dims\",sep=\"\\t\",col.names=FALSE,row.names=TRUE,quote=FALSE)
write.table(IN[,$range_R_set],file=\"$opt{'O'}.pcurve.orig.dims\",sep=\"\\t\",col.names=FALSE,row.names=TRUE,quote=FALSE)
LAMBDA <- cbind(PCURVE\$lambda)
rownames(LAMBDA) <- rownames(PCURVE\$s)
write.table(LAMBDA,file=\"$opt{'O'}.pcurve.lambda\",sep=\"\\t\",col.names=FALSE,row.names=TRUE,quote=FALSE)
"; close R;

system("$Rscript $opt{'O'}.pcurve.r");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.pcurve.r");
}

exit;
}

if ($command eq "pcurve-center" || $command eq "center-pcurve" || $command eq "lambda-center" || $command eq "center-lambda") { # TODO: Verify

# Defaults
$max = 1;
$type = "centered";

getopts("O:A:a:M:P:r", \%opt);

$die2 = "
scitools pcurve-center [options] [input pcurve.lambda]
   or    center-pcurve
         lambda-center
         center-lambda

Centers a pcurve lambda and adjusts the scale.
Default is to center on all cells and extend lambda out.

Options:
   -O   [STR]   Output prefix (default is [prefix].centered.lambda)
   -A   [STR]   Annotation file (optional)
   -a   [STR]   Annotation that should be used to center the pcurve
                  (requires -A to be specified)
   -M   [FLT]   Maximum absolute value to adjust lambda to (either side of 0)
                  (def = $max)
   -P   [STR]   Progression type (centered/cen, progressive/pro)
                  (def = $type, if -a it will only center on that annotation)
   -r           Reverse order (defaults to current order)

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.lambda$//;
if (defined $opt{'M'}) {$max = $opt{'M'}};
if (defined $opt{'P'}) {$type = $opt{'P'}};
if ($type !~ /^(cen|pro)/i) {
	die "\nERROR: Cannot determine preogresstion type ($type) - specify as centered/cen or progressive/pro\n";
}
if (defined $opt{'A'}) {read_annot($opt{'A'})};
if (defined $opt{'r'}) {$r = -1} else {$r = 1};

read_values($ARGV[0]);

if (defined $opt{'a'}) {
	$type = "centered";
	if (defined $ANNOT_count{$opt{'a'}}) {
		$annot_value_sum = 0; $annot_cell_count = 0; $center_value = 0;
		foreach $cellID (keys %CELLID_value) {
			if ($CELLID_annot{$cellID} eq $opt{'a'}) {
				$annot_cell_count++;
				$annot_value_sum += $CELLID_value{$cellID};
			}
		}
		if ($annot_cell_count>0) {
			$center_value = $annot_value_sum/$annot_cell_count;
			if (abs($value_max - $center_value)>abs($value_min - $center_value)) {
				$progression_span = abs($value_max-$center_value);
			} else {
				$progression_span = abs($value_min - $center_value);
			}
		} else {
			die "\nERROR: No cells in labda file $ARGV[0] have the annotation $opt{'a'} that was provided.\n";
		}
	} else {
		die "\nERROR: Provided annotation $opt{'a'} cannot be found in the annotaiton file provided: $opt{'A'}, or the annotation file was not specified (-A)\n";
	}
} elsif ($type =~ /^cen/i) {
	$center_value = $value_mean;
	if (abs($value_max - $center_value)>abs($value_min - $center_value)) {
		$progression_span = abs($value_max-$center_value);
	} else {
		$progression_span = abs($value_min - $center_value);
	}
} else {
	$center_value = $value_min;
	$progression_span = $value_range;
}

open OUT, ">$opt{'O'}.centered.lambda";
foreach $cellID (sort {$CELLID_value{$a}<=>$CELLID_value{$b}} keys %CELLID_value) {
	$lambda = ((($CELLID_value{$cellID}-$center_value)/$progression_span)*$max*$r);
	print OUT "$cellID\t$lambda\n";
} close OUT;

exit;
}

if ($command eq "pcurve-prune" || $command eq "prune-pcurve") {

# Defaults
$prune = 0.1;

getopts("O:o:p:l:e:", \%opt);

$die2 = "
scitools pcurve-prune [options] [pcurve prefix]
   or    prune-pcurve

Prunes cells distant from the pcurve. It will print a new file for
the orig and proj dims, and lambda.

Will search for [prefix].orig.dims
                [prefix].proj.dims
                [prefix].lambda (or [prefix].centered.lambda)

Options:
   -O   [STR]   Output prefix (default is [prefix].pruned.lambda)
   -e   [FLT]   Fraction of cells to exclude by distance (def = $prune)
   -o   [STR]   Original dims file (def = auto find [prefix].orig.dims)
   -p   [STR]   Pcurve dims file (def = auto find [prefix].proj.dims)
   -l   [STR]   Lambda file (def = auto find [prefix].lambda)
                (if [prefix].centered.lambda is found, it will use it)

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.lambda$//;

if (!defined $opt{'o'}) {
	if (-e "$ARGV[0].orig.dims") {
		$opt{'o'} = "$ARGV[0].orig.dims";
	} else {
		die "\n\nERROR: Cannot find $ARGV[0].orig.dims, check your prefix ($ARGV[0])\n";
	}
}
if (!defined $opt{'p'}) {
	if (-e "$ARGV[0].proj.dims") {
		$opt{'p'} = "$ARGV[0].proj.dims";
	} else {
		die "\n\nERROR: Cannot find $ARGV[0].proj.dims, check your prefix ($ARGV[0])\n";
	}
}
if (!defined $opt{'l'}) {
	if (-e "$ARGV[0].centered.lambda") {
		$opt{'l'} = "$ARGV[0].centered.lambda";
	} elsif (-e "$ARGV[0].lambda") {
		$opt{'l'} = "$ARGV[0].lambda";
	} else {
		die "\n\nERROR: Cannot find $ARGV[0].centered.lambda OR $ARGV[0].lambda, check your prefix ($ARGV[0])\n";
	}
}

read_dims($opt{'o'});
read_pcurve_dims($opt{'p'});
read_values($opt{'l'});

# calculate the distances
%CELLID_pcurveDist = ();
$cells_in_set = 0;
foreach $cellID (keys %CELLID_DIMS) {
	if (defined $CELLID_PCURVE_DIMS{$cellID}[0]) {
		$dist_sum = 0;
		for ($dim = 1; $dim < @{$CELLID_DIMS{$cellID}}; $dim++) {
			if (defined $CELLID_PCURVE_DIMS{$cellID}[$dim]) {
				$dist_sum += ($CELLID_PCURVE_DIMS{$cellID}[$dim] - $CELLID_DIMS{$cellID}[$dim])**2;
			}
		}
		$CELLID_pcurveDist{$cellID} = sqrt($dist_sum);
		$cells_in_set++;
	}
}

# figure out which ones to keep & print to new lambda file
system("cp $opt{'o'} $opt{'O'}.pruned.orig.dims"); # this stays the same
open LAMBDA, ">$opt{'O'}.pruned.lambda";
open PROJ, ">$opt{'O'}.pruned.proj.dims";
$excluded_ct = 0; $checked_ct = 0;
foreach $cellID (sort {$CELLID_pcurveDist{$b}<=>$CELLID_pcurveDist{$a}} keys %CELLID_pcurveDist) {
	$checked_ct++;
	if (($checked_ct/$cells_in_set)<=$prune) {
		$excluded_ct++;
		$threshold = $CELLID_pcurveDist{$cellID};
	} else {
		print LAMBDA "$cellID\t$CELLID_value{$cellID}\n";
		print PROJ "$cellID";
		for ($dim = 1; $dim < @{$CELLID_PCURVE_DIMS{$cellID}}; $dim++) {
			print PROJ "\t$CELLID_PCURVE_DIMS{$cellID}[$dim]";
		} print PROJ "\n";
	}
} close LAMBDA; close PROJ;

print STDERR "INFO: $excluded_ct cells were excluded out of $cells_in_set, with a distance of $threshold or greater from the pcurve.\n";

exit;
}


########## PLOTTING FUNCTIONS ##########
if ($command eq "plot-complexity") {

getopts("O:A:a:C:c:R:M", \%opt);

$die2 = "
scitools plot-complexity [options] [complexity file(s) can be comma separated]

Options:
   -O   [STR]   Output prefix (default is complexity file 1 prefix)
   -M           Run mixed model to determine read count cutoff for cells (def = no)
   -A   [STR]   Annotation file (to color code points)
   -a   [STR]   Comma separated list of annoations to include in plot
                (requires -A to be specified)
   -C   [STR]   Color coding file (annot (tab) #hexColor)
   -c   [STR]   Color coding string
                Annot=#hexColor,Annot2=#hexColor
   -R   [STR]   Rscript call (def = $Rscript)

Note: Requires ggplot2 R package

";

if (!defined $ARGV[0]) {die $die2};
if (defined $opt{'a'} && !defined $opt{'A'}) {die "\nMust provide an annotaiton file (-A) if specifying annotations to plot (-a)!\n$die2"};
if (defined $opt{'C'} && defined $opt{'c'}) {die "\nSpecify either a color string (-c) or a color coding file (-C), not both!\n$die2"};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (defined $opt{'A'}) {read_annot($opt{'A'})};
if (defined $opt{'a'}) {
	@ANNOT_LIST = split(/,/, $opt{'a'});
	foreach $annot (@ANNOT_LIST) {
		$ANNOT_include{$annot} = 1;
	}
}
if (defined $opt{'C'}) {read_color_file($opt{'C'})};
if (defined $opt{'c'}) {read_color_string($opt{'c'})};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.txt$//;

read_complexity($ARGV[0]);

open OUT, ">$opt{'O'}.plot.txt";
foreach $cellID (keys %CELLID_complexity) {
	if (defined $opt{'a'}) {
		$annot = $CELLID_annot{$cellID};
		if (defined $ANNOT_include{$annot} && defined $CELLID_annot{$cellID}) {
			print OUT "$cellID\t$CELLID_annot{$cellID}\t$CELLID_uniq_reads{$cellID}\t$CELLID_complexity{$cellID}\n";
		}
	} elsif (defined $opt{'A'} && defined $CELLID_annot{$cellID}) {
		print OUT "$cellID\t$CELLID_annot{$cellID}\t$CELLID_uniq_reads{$cellID}\t$CELLID_complexity{$cellID}\n";
	} else {
		print OUT "$cellID\tCell\t$CELLID_uniq_reads{$cellID}\t$CELLID_complexity{$cellID}\n";
	}
} close OUT;

open R, ">$opt{'O'}.plot.r";
print R "
library(ggplot2)
IN<-read.table(\"$opt{'O'}.plot.txt\")
PLT<-ggplot(data=subset(IN,V4<100&V4>0)) + theme_bw() +
   geom_point(aes(V4,log10(V3),color=V2),size=1,alpha=0.3) +
   geom_density2d(aes(V4,log10(V3),color=V2),size=0.3) +
";
if (defined $opt{'C'} || defined $opt{'c'}) {
	print R "   scale_colour_manual(values = c($color_mapping)) +
	guides(colour = guide_legend(override.aes = list(size=4))) +";
}
print R "
   scale_x_continuous(limits=c(0,100)) +
   scale_y_continuous(limits=c(0,6)) +
   xlab(\"Complexity\") +
   ylab(\"log10 Unique Reads\") +";
if (defined $opt{'A'}) {
print R "
	theme(legend.background=element_blank(),legend.title=element_blank())";
} else {
print R "
	theme(legend.position=\"none\")";
}
print R "
ggsave(plot=PLT,filename=\"$opt{'O'}.png\",width=7,height=6)
ggsave(plot=PLT,filename=\"$opt{'O'}.pdf\",width=7,height=6)
";

if (defined $opt{'M'}) {
	print R "
IN_sub=subset(IN,V4<100&V4>0)

#take unique aligned 
IN_sub\$unique_aligned<-IN_sub\$V3 
x1 <- IN_sub\$unique_aligned[IN_sub\$unique_aligned != 0]

# trimodal fit
km <- kmeans(log10(x1),centers=3)
clustr <- as.factor(km\$cluster)
data_1<-data.frame(\"val\"=log10(x1),\"km\"=clustr)
highest_cluster<-which.max(km\$centers)
second<-km\$center[-highest_cluster,1]
second_highest_cluster<-as.numeric(names(which.max(second)))
cluster_top<-data_1\$val[which(data_1\$km==highest_cluster)]
border1=mean(max(data_1\$val[which(data_1\$km==second_highest_cluster)]),min(data_1\$val[which(data_1\$km==highest_cluster)]))

#in case the groups seperate we do normal fit
sink(\"$opt{'O'}.log\")
fit <- fitdistr(cluster_top, \"normal\")
para <- fit\$estimate
threshold_final_1=para[1]-1.96*para[2]

#in case 3 groups do not sepearte well do a mixed model instead of normal
MMS<-normalmixEM(log10(cluster_top),maxrestarts=1000,maxit = 10000)
threshold_final_2 =10**(MMS\$mu[2]-1.96*MMS\$sigma[2])

p1<-ggplot(data_1, aes(x=data_1\$val)) + geom_histogram(aes(fill=data_1\$km,y=..count../sum(..count..)),color=\"grey50\")+ylab(\"Density\")+stat_density(geom=\"line\",color=\"red\")+geom_vline(xintercept = threshold_final_1)+xlab(\"Log10 Unique Reads\")+theme_bw()+theme(legend.background=element_blank(),text = element_text(size=10),legend.title=element_blank())
ggsave(p1,filename = \"$opt{'O'}.dist.threshold_normal.pdf\")
ggsave(p1,filename = \"$opt{'O'}.dist.threshold_normal.png\")

p2<-ggplot(data_1, aes(x=data_1\$val)) + geom_histogram(aes(fill=data_1\$km,y=..count../sum(..count..)),color=\"grey50\")+ylab(\"Density\")+stat_density(geom=\"line\",color=\"red\")+geom_vline(xintercept = threshold_final_2)+xlab(\"Log10 Unique Reads\")+theme_bw()+theme(legend.background=element_blank(),text = element_text(size=10),legend.title=element_blank())
ggsave(p2,filename = \"$opt{'O'}.dist.threshold_mixed.pdf\")
ggsave(p2,filename = \"$opt{'O'}.dist.threshold_mixed.png\")

print(paste(\"the number of total cells: \",length(IN_sub\$unique_aligned)))
print(paste(\"threshold 1: \",10**threshold_final_1))
print(paste(\"the number of cells above this threshold 1: \",length(which(IN_sub\$unique_aligned>round(10**threshold_final_1)))))

print(paste(\"threshold 2: \",10**threshold_final_2))
print(paste(\"the number of cells above this threshold 2: \",length(which(IN_sub\$unique_aligned>round(10**threshold_final_2)))))
sink()
";
}

print R "
PLT<-ggplot(data=subset(IN,V4<100&V4>0)) + theme_bw() +
	geom_histogram(aes(log10(V3),fill=V2)) +";
if (defined $opt{'C'} || defined $opt{'c'}) {
	print R "   scale_colour_manual(values = c($color_mapping)) +";
}
print R "
	xlab(\"log10 Unique Reads\") +
	ylab(\"Counts\") +
	scale_x_continuous(limits=c(0,6)) +";
if (defined $opt{'A'}) {
print R "
	theme(legend.background=element_blank(),legend.title=element_blank())";
} else {
print R "
	theme(legend.position=\"none\")";
}
print R "
ggsave(plot=PLT,filename=\"$opt{'O'}.hist.png\",width=7,height=6)
ggsave(plot=PLT,filename=\"$opt{'O'}.hist.pdf\",width=7,height=6)
";

close R;

system("$Rscript $opt{'O'}.plot.r");

exit;

}

if ($command eq "plot-dims") {

# Defaults
$xdim = 1;
$ydim = 2;
$minV = -10;
$maxV = 10;
$theme = "Clean";
$gradient_def = "BuG90Rd";
$ptSize = 1;

getopts("O:A:a:C:c:R:x:y:T:V:M:XS:s:G:p:", \%opt);

$die2 = "
scitools plot-dims [options] [dimensions file(s), comma sep]

Options:
   -O   [STR]   Output prefix (default is dims file 1 prefix)
   -A   [STR]   Annotation file (to color code points)
   -a   [STR]   Comma separated list of annotations to include in plot
                  (requires -A to be specified)
   -x   [INT]   X-dimension to plot (def = $xdim)
   -y   [INT]   Y-dimension to plot (def = $ydim)
   -T   [STR]   Theme: (def = $theme)
                  Clean = no axis lines (for tSNE)
                  Reg = regular, include axes and grid (for PCA)
   -p   [FLT]   Point size (def = $ptSize)
   -C   [STR]   Color coding file (annot (tab) #hexColor)
   -c   [STR]   Color coding string
                  Annot=#hexColor,Annot2=#hexColor
   -V   [STR]   Values file. tab-delimited, cellID (tab) value
                  Will plot as teh color of points (overrides
                  annotation colors)
   -S   [MIN,MAX]   Min and max values for scaling (if -V specified)
                  (default = -10,10)
   -G   [GRD]   Color gradient (def = $gradient_def)
                  For all available gradients, run 'scitools gradient'
   -M   [STR]   Matrix file (e.g. deviation z-scores from chromVAR)
                  Will create a folder as -O option and produce
                  a plot for each row of the matrix. Overrides -V.
                  Only includes rows with at least one non-zero value.
   -R   [STR]   Rscript call (def = $Rscript)
   -s   [STR]   scitools call (def = $scitools)
   -X           Do not delete intermediate files (def = delete)

Note: Requires ggplot2 R package

";

if (!defined $ARGV[0]) {die $die2};
if (defined $opt{'a'} && !defined $opt{'A'}) {die "\nMust provide an annotaiton file (-A) if specifying annotations to plot (-a)!\n$die2"};
if (defined $opt{'C'} && defined $opt{'c'}) {die "\nSpecify either a color string (-c) or a color coding file (-C), not both!\n$die2"};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.txt$//; $opt{'O'} =~ s/\.dims$//;
if (!defined $opt{'G'}) {$opt{'G'} = $gradient_def};
$gradient_function = get_gradient($opt{'G'});
if (defined $opt{'p'}) {$ptSize = $opt{'p'}};

read_dims($ARGV[0]);

if (defined $opt{'M'}) {
	print STDERR "SCITOOLS: Matrix file plotting detected! Will plot a separate value-based plot for each row entry of the matrix.\n";
	
	$matrix_out = "matrix_plots";
	
	if (-e "$opt{'O'}.$matrix_out") {
		die "\nFATAL: $opt{'O'}.$matrix_out directory already exists! Exiting!\n$die2";
	}
	system("mkdir $opt{'O'}.$matrix_out");
	
	read_matrix($opt{'M'});
	
	$common_opts = "";
	if (defined $opt{'A'}) {$common_opts .= "-A $opt{'A'} "};
	if (defined $opt{'a'}) {$common_opts .= "-a $opt{'a'} "};
	if (defined $opt{'x'}) {$common_opts .= "-x $opt{'x'} "};
	if (defined $opt{'y'}) {$common_opts .= "-y $opt{'y'} "};
	if (defined $opt{'T'}) {$common_opts .= "-T $opt{'T'} "};
	if (defined $opt{'S'}) {$common_opts .= "-S $opt{'S'} "};
	if (defined $opt{'R'}) {$common_opts .= "-R $opt{'R'} "};
	if (defined $opt{'X'}) {$common_opts .= "-X "};
	$common_opts =~ s/\s$//;
	
	foreach $feature (keys %MATRIX_feature_nonZero) {
		$feature_polished = $feature;
		$feature_polished =~ s/(:|;|'|\.|\(|\)|\{|\}|\[|\]|\s|\|)/_/;
		open VALS, ">$opt{'O'}.$matrix_out/$feature_polished.values";
		foreach $cellID (keys %CELLID_FEATURE_value) {
			if (defined $CELLID_FEATURE_value{$cellID}{$feature}) {
				print VALS "$cellID\t$CELLID_FEATURE_value{$cellID}{$feature}\n";
			}
		} close VALS;
		system("$scitools plot-dims $common_opts -O $opt{'O'}.$matrix_out/$feature_polished -V $opt{'O'}.$matrix_out/$feature_polished.values -G $opt{'G'} $ARGV[0]");
		if (!defined $opt{'X'}) {system("rm -f $opt{'O'}.$matrix_out/$feature_polished.values")};
	}
	
	exit;
}

if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (defined $opt{'s'}) {$scitools = $opt{'s'}};
if (defined $opt{'A'}) {read_annot($opt{'A'})};
if (defined $opt{'a'}) {
	@ANNOT_LIST = split(/,/, $opt{'a'});
	foreach $annot (@ANNOT_LIST) {
		$ANNOT_include{$annot} = 1;
	}
}
if (defined $opt{'C'}) {read_color_file($opt{'C'})};
if (defined $opt{'c'}) {read_color_string($opt{'c'})};
if (defined $opt{'V'}) {read_values($opt{'V'})};
if (defined $opt{'x'}) {$xdim = $opt{'x'}};
if (defined $opt{'y'}) {$ydim = $opt{'y'}};
if (defined $opt{'S'}) {($minV,$maxV) = split(/,/, $opt{'S'})};

open DATA, ">$opt{'O'}.plot.txt";
foreach $cellID (keys %CELLID_DIMS) {

	if (defined $opt{'A'}) {
		if (defined $opt{'a'}) {
			if (defined $CELLID_annot{$cellID}) {
				if (defined $ANNOT_include{$CELLID_annot{$cellID}}) {
					$annot = $CELLID_annot{$cellID}
				} else {$annot = "Exclude"};
			} else {$annot = "Exclude"};
		} else {
			if (defined $CELLID_annot{$cellID}) {
				$annot = $CELLID_annot{$cellID};
			} else {$annot = "Exclude"};
		}
	} else {
		$annot = "Cell";
	}
	
	if ($annot !~ /Exclude/i) {
		if (!defined $opt{'V'}) { # qualitative annotations
			print DATA "$cellID\t$annot\t$CELLID_DIMS{$cellID}[$xdim]\t$CELLID_DIMS{$cellID}[$ydim]\n";
		} else { # value annotations
			if (defined $CELLID_value{$cellID}) {
				print DATA "$cellID\t$CELLID_value{$cellID}\t$CELLID_DIMS{$cellID}[$xdim]\t$CELLID_DIMS{$cellID}[$ydim]\n";
			} else {
				print STDERR "WARNING: $cellID does not have values specified in $opt{'V'}, skipping.\n";
			}
		}
	}
	
} close DATA;

open R, ">$opt{'O'}.plot.r";

print R "
library(ggplot2)
IN<-read.table(\"$opt{'O'}.plot.txt\")
$gradient_function
PLT<-ggplot() +";

if (!defined $opt{'c'} && !defined $opt{'C'} && !defined $opt{'A'} && !defined $opt{'V'}) {
	print R "
	geom_point(aes(IN\$V3,IN\$V4),color=\"lightsteelblue4\",size=$ptSize) +";
} elsif (!defined $opt{'V'}) {
	print R "
	geom_point(aes(IN\$V3,IN\$V4,color=IN\$V2),size=$ptSize) +
	guides(colour = guide_legend(override.aes = list(size=4))) +";
} else {
	print R "
	geom_point(aes(IN\$V3,IN\$V4,color=IN\$V2),size=$ptSize) +";
}

if ($color_mapping !~ /none/i && !defined $opt{'V'}) {
	print R "
	scale_colour_manual(values = c($color_mapping)) +";
} elsif (defined $opt{'V'}) {
	print R "
	scale_colour_gradientn(colours=gradient_funct(21),limits=c($minV,$maxV)) +";
}

if ($theme =~ /Clean/i) {
	if (defined $opt{'A'}) {
	print R "
	theme_bw() +
	theme(panel.border=element_blank(),
		  panel.grid=element_blank(),
		  axis.line=element_blank(),
		  axis.ticks=element_blank(),
		  legend.background=element_blank(),
		  legend.title=element_blank(),
		  panel.background=element_blank(),
		  axis.text=element_blank(),
		  axis.title.x=element_blank(),
		  axis.title.y=element_blank(),
		  plot.background=element_blank(),
		  plot.margin=unit(c(0,0,0,0),\"pt\"))\n";
	} else {
	print R "
	theme_bw() +
	theme(panel.border=element_blank(),
		  panel.grid=element_blank(),
		  axis.line=element_blank(),
		  axis.ticks=element_blank(),
		  legend.position=\"none\",
		  panel.background=element_blank(),
		  axis.text=element_blank(),
		  axis.title.x=element_blank(),
		  axis.title.y=element_blank(),
		  plot.background=element_blank(),
		  plot.margin=unit(c(0,0,0,0),\"pt\"))\n";
	}
} else {
	if (defined $opt{'A'}) {
	print R "
	theme_bw()\n";
	} else {
	print R "
	theme_bw(legend.position=\"none\")\n";
	}
}

print R "
ggsave(plot=PLT,filename=\"$opt{'O'}.plot.png\",width=5,height=4,dpi=900)
ggsave(plot=PLT,filename=\"$opt{'O'}.plot.pdf\",width=6,height=5)
";

close R;

system("$Rscript $opt{'O'}.plot.r");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.plot.r $opt{'O'}.plot.txt");
}

exit;
}

if ($command eq "plot-pcurve") {

# Defaults
$xdim = 1;
$ydim = 2;
$alpha = 1;
$ptSize = 1;
$theme = "Clean";
$gradient_def = "YlOrRd";

getopts("O:A:a:C:c:R:x:y:T:Xo:p:l:G:M:f:s:", \%opt);

$die2 = "
scitools plot-dims [options] [pcurve prefix]

Will search for [prefix].orig.dims
                [prefix].proj.dims
                [prefix].lambda (or [prefix].centered.lambda)
                                (or [prefix].pruned.lambda)

Options:
   -O   [STR]   Output prefix (default is input prefix)
   -o   [STR]   Original dims file (def = auto find [prefix].orig.dims)
   -p   [STR]   Pcurve dims file (def = auto find [prefix].proj.dims)
   -l   [STR]   Lambda file (def = auto find [prefix].lambda)
                (if [prefix].centered.lambda is found, it will it)
   -M   [STR]   Matrix file - will plot each value set in a matrix
                along the lambda values of the pcurve
   -f   [FLT]   Alpha of points (for -M plotting; def = $alpha)
   -s   [FLT]   Point size (def = $ptSize)
   -A   [STR]   Annotation file (to color code points)
   -a   [STR]   Comma separated list of annoations to include in plot
                  (requires -A to be specified)
   -x   [INT]   X-dimension to plot (def = $xdim)
   -y   [INT]   Y-dimension to plot (def = $ydim)
   -T   [STR]   Theme: (def = $theme)
                  Clean = no axis lines (for tSNE)
                  Reg = regular, include axes and grid (for PCA)
   -C   [STR]   Color coding file (annot (tab) #hexColor)
   -c   [STR]   Color coding string
                  Annot=#hexColor,Annot2=#hexColor
   -G   [GRD]   Color gradient for labmda (def = $gradient_def)
                  For all available gradients, run 'scitools gradient'
   -R   [STR]   Rscript call (def = $Rscript)
   -X           Do not delete intermediate files (def = delete)

Note: Requires ggplot2 R package

If -O -o -p and -l are all defined [pcurve prefix] does not need to be
specified to execute this command.

";

if (!defined $ARGV[0]) { 
    if(!defined $opt{'O'} &&
	   !defined $opt{'o'} &&
	   !defined $opt{'p'} &&
	   !defined $opt{'l'}) {die $die2};
}
if (defined $opt{'a'} && !defined $opt{'A'}) {die "\nMust provide an annotaiton file (-A) if specifying annotations to plot (-a)!\n$die2"};
if (defined $opt{'C'} && defined $opt{'c'}) {die "\nSpecify either a color string (-c) or a color coding file (-C), not both!\n$die2"};

if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (defined $opt{'s'}) {$ptSize = $opt{'s'}};
if (defined $opt{'f'}) {$alpha = $opt{'f'}};
if (defined $opt{'A'}) {read_annot($opt{'A'})};
if (defined $opt{'a'}) {
	@ANNOT_LIST = split(/,/, $opt{'a'});
	foreach $annot (@ANNOT_LIST) {
		$ANNOT_include{$annot} = 1;
	}
}
if (defined $opt{'C'}) {read_color_file($opt{'C'})};
if (defined $opt{'c'}) {read_color_string($opt{'c'})};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.txt$//; $opt{'O'} =~ s/\.dims$//;
if (defined $opt{'x'}) {$xdim = $opt{'x'}};
if (defined $opt{'y'}) {$ydim = $opt{'y'}};
if (!defined $opt{'G'}) {$opt{'G'} = $gradient_def};
$gradient_function = get_gradient($opt{'G'});

if (!defined $opt{'o'}) {
	if (-e "$ARGV[0].orig.dims") {
		$opt{'o'} = "$ARGV[0].orig.dims";
	} else {
		die "\n\nERROR: Cannot find $ARGV[0].orig.dims, check your prefix ($ARGV[0])\n";
	}
}
if (!defined $opt{'p'}) {
	if (-e "$ARGV[0].proj.dims") {
		$opt{'p'} = "$ARGV[0].proj.dims";
	} else {
		die "\n\nERROR: Cannot find $ARGV[0].proj.dims, check your prefix ($ARGV[0])\n";
	}
}
if (!defined $opt{'l'}) {
	if (-e "$ARGV[0].centered.lambda") {
		$opt{'l'} = "$ARGV[0].centered.lambda";
	} elsif (-e "$ARGV[0].lambda") {
		$opt{'l'} = "$ARGV[0].lambda";
	} else {
		die "\n\nERROR: Cannot find $ARGV[0].centered.lambda OR $ARGV[0].lambda, check your prefix ($ARGV[0])\n";
	}
}

read_dims($opt{'o'});
read_pcurve_dims($opt{'p'});
read_values($opt{'l'});


if (defined $opt{'M'}) {
	print STDERR "SCITOOLS: Matrix file plotting detected! Will plot a separate value-based plot for each row entry of the matrix.\n";
	
	$matrix_out = "matrix_plots";
	
	if (-e "$opt{'O'}.$matrix_out\pcurve") {
		die "\nFATAL: $opt{'O'}.$matrix_out\pcurve directory already exists! Exiting!\n$die2";
	}
	system("mkdir $opt{'O'}.$matrix_out\pcurve");
	
	read_matrix($opt{'M'});
	
	foreach $feature (keys %MATRIX_feature_nonZero) {
		$feature_polished = $feature;
		$feature_polished =~ s/(:|;|'|\.|\(|\)|\{|\}|\[|\]|\s|\|)/_/;
		open VALS, ">$opt{'O'}.$matrix_out\pcurve/$feature_polished.values";

		foreach $cellID (keys %CELLID_FEATURE_value) {
			if (defined $opt{'A'}) {
				if (defined $opt{'a'}) {
					if (defined $CELLID_annot{$cellID}) {
						if (defined $ANNOT_include{$CELLID_annot{$cellID}}) {
							$annot = $CELLID_annot{$cellID}
						} else {$annot = "Exclude"};
					} else {$annot = "Exclude"};
				} else {
					if (defined $CELLID_annot{$cellID}) {
						$annot = $CELLID_annot{$cellID};
					} else {$annot = "Cell"};
				}
			} else {
				$annot = "Cell";
			}
			if (defined $CELLID_FEATURE_value{$cellID}{$feature} && $CELLID_value{$cellID} && ($annot !~ /Exclude/i)) {
				print VALS "$cellID\t$CELLID_value{$cellID}\t$CELLID_FEATURE_value{$cellID}{$feature}\t$annot\n";
			}
		} close VALS;
		
		open R, ">$opt{'O'}.plot.r";

# lambda vals plot
	print R "
library(ggplot2)
IN<-read.table(\"$opt{'O'}.$matrix_out\pcurve/$feature_polished.values\",header=FALSE)
# Make the lambda and matrix val plot
PointPlot<-ggplot() + ";

# plot the original dims w/ respective color
if (!defined $opt{'c'} && !defined $opt{'C'} && !defined $opt{'A'}) {
	print R "
	geom_point(aes(x=IN\$V2,y=IN\$V3),color=\"lightsteelblue4\",size=$ptSize,alpha=$alpha) +";
} else {
	print R "
	geom_point(aes(x=IN\$V2,y=IN\$V3,color=IN\$V4),size=$ptSize,alpha=$alpha) +
	guides(colour = guide_legend(override.aes = list(size=4))) +";
}

if ($color_mapping !~ /none/i) {
	print R "
	scale_colour_manual(values = c($color_mapping)) +";
}

print R "
geom_smooth(aes(x=IN\$V2,y=IN\$V3)) +
theme_bw() + xlab(\"Lambda\") + ylab(\"Feature value\") +";
if (defined $opt{'A'}) {
print R "
	theme(legend.background=element_blank(),
		legend.title=element_blank(),
		panel.background=element_blank(),
		plot.background=element_blank())";
} else {
print R "
	theme(legend.position=\"none\",
		panel.background=element_blank(),
		plot.background=element_blank())";
}
print R "

ggsave(plot=PointPlot,filename=\"$opt{'O'}.$matrix_out\pcurve/$feature_polished.lambda.png\",width=5,height=4,dpi=900)
ggsave(plot=PointPlot,filename=\"$opt{'O'}.$matrix_out\pcurve/$feature_polished.lambda.plot.pdf\",width=5,height=4)

a=cor(IN\$V2,IN\$V3,method =\"pearson\")
b=cor(IN\$V2,IN\$V3,method =\"spearman\")
output<-c(\"$feature_polished\",a)
output2<-c(\"$feature_polished\",b)
write(output,file=\"$opt{'O'}.$matrix_out\pcurve/corr.pearson.txt\",append=TRUE,ncolumns=2,sep = \"\\t\")
write(output2,file=\"$opt{'O'}.$matrix_out\pcurve/corr.spearman.txt\",append=TRUE,ncolumns=2,sep = \"\\t\")

";
	close R;
	
		system("$Rscript $opt{'O'}.plot.r");
		if (!defined $opt{'X'}) {system("rm -f $opt{'O'}.$matrix_out\pcurve/$feature_polished.values $opt{'O'}.plot.r")};
	}
	
	
	open R, ">$opt{'O'}.plot.r";
	print R "
	library(ggplot2)
	IN<-read.table(\"$opt{'O'}.$matrix_out\pcurve/corr.pearson.txt\",header=FALSE)
	# Make the corr distribution plot
	
	HISTPlot<-ggplot(IN,aes(x=V2)) + geom_histogram()+ 
	theme_bw() + xlab(\"Correlation\") + ylab(\"Counts\") +";
if (defined $opt{'A'}) {
print R "
	theme(legend.background=element_blank(),
		legend.title=element_blank(),
		panel.background=element_blank(),
		plot.background=element_blank())";
} else {
print R "
	theme(legend.position=\"none\",
		panel.background=element_blank(),
		plot.background=element_blank())";
}
print R "

	ggsave(plot=HISTPlot,filename=\"$opt{'O'}.$matrix_out\pcurve/Hist.corr.Pearson.lambda.png\",width=5,height=4,dpi=900)
	ggsave(plot=HISTPlot,filename=\"$opt{'O'}.$matrix_out\pcurve/Hist.corr.Pearson.lambda.plot.pdf\",width=5,height=4)
	
		IN<-read.table(\"$opt{'O'}.$matrix_out\pcurve/corr.spearman.txt\",header=FALSE)
	# Make the corr distribution plot
	
	HISTPlot<-ggplot(IN,aes(x=V2)) + geom_histogram()+ 
	theme_bw() + xlab(\"Correlation\") + ylab(\"Counts\") +";
if (defined $opt{'A'}) {
print R "
	theme(legend.background=element_blank(),
		legend.title=element_blank(),
		panel.background=element_blank(),
		plot.background=element_blank())";
} else {
print R "
	theme(legend.position=\"none\",
		panel.background=element_blank(),
		plot.background=element_blank())";
}
print R "

	ggsave(plot=HISTPlot,filename=\"$opt{'O'}.$matrix_out\pcurve/Hist.corr.Spearman.lambda.png\",width=5,height=4,dpi=900)
	ggsave(plot=HISTPlot,filename=\"$opt{'O'}.$matrix_out\pcurve/Hist.corr.Spearman.lambda.plot.pdf\",width=5,height=4)
	";
	close(R);
	system("$Rscript $opt{'O'}.plot.r");
	if (!defined $opt{'X'}) {system("rm -f $opt{'O'}.plot.r")};
	
	exit;
}


open DATA, ">$opt{'O'}.plot.txt";
print DATA "cellID\tannot\tlambda\todim1\todim2\tpdim1\tpdim2\n";
foreach $cellID (keys %CELLID_DIMS) {

	if (defined $opt{'A'}) {
		if (defined $opt{'a'}) {
			if (defined $CELLID_annot{$cellID}) {
				if (defined $ANNOT_include{$CELLID_annot{$cellID}}) {
					$annot = $CELLID_annot{$cellID}
				} else {$annot = "Exclude"};
			} else {$annot = "Exclude"};
		} else {
			if (defined $CELLID_annot{$cellID}) {
				$annot = $CELLID_annot{$cellID};
			} else {
				$annot = "Cell";
			}
		}
	} else {
		$annot = "Cell";
	}
	
	if ($annot !~ /Exclude/i && defined $CELLID_PCURVE_DIMS{$cellID}[$xdim]) {
		print DATA "$cellID\t$annot\t$CELLID_value{$cellID}\t$CELLID_DIMS{$cellID}[$xdim]\t$CELLID_DIMS{$cellID}[$ydim]\t$CELLID_PCURVE_DIMS{$cellID}[$xdim]\t$CELLID_PCURVE_DIMS{$cellID}[$ydim]\n";
	}
	
} close DATA;

open R, ">$opt{'O'}.plot.r";

# DIMS PLOT
print R "
library(ggplot2)
IN<-read.table(\"$opt{'O'}.plot.txt\",header=TRUE)

# Make the odim and pdim plot
PointPlot<-ggplot() +";

# plot the original dims w/ respective color
if (!defined $opt{'c'} && !defined $opt{'C'} && !defined $opt{'A'}) {
	print R "
	geom_point(aes(IN\$odim1,IN\$odim2),color=\"lightsteelblue4\",size=1) +";
} else {
	print R "
	geom_point(aes(IN\$odim1,IN\$odim2,color=IN\$annot),size=1) +
	guides(colour = guide_legend(override.aes = list(size=4))) +";
}

if ($color_mapping !~ /none/i) {
	print R "
	scale_colour_manual(values = c($color_mapping)) +";
}

# plot the projected points and color by lambda
print R "
	geom_point(aes(IN\$pdim1,IN\$pdim2),color=\"gray30\",size=1) +";

# set theme
if ($theme =~ /Clean/i) {
	if (defined $opt{'A'}) {
	print R "
	theme_bw() +
	theme(panel.border=element_blank(),
		  panel.grid=element_blank(),
		  axis.line=element_blank(),
		  axis.ticks=element_blank(),
		  legend.background=element_blank(),
		  legend.title=element_blank(),
		  panel.background=element_blank(),
		  axis.text=element_blank(),
		  axis.title.x=element_blank(),
		  axis.title.y=element_blank(),
		  plot.background=element_blank(),
		  plot.margin=unit(c(0,0,0,0),\"pt\"))\n";
	} else {
	print R "
	theme_bw() +
	theme(panel.border=element_blank(),
		  panel.grid=element_blank(),
		  axis.line=element_blank(),
		  axis.ticks=element_blank(),
		  legend.position=\"none\",
		  panel.background=element_blank(),
		  axis.text=element_blank(),
		  axis.title.x=element_blank(),
		  axis.title.y=element_blank(),
		  plot.background=element_blank(),
		  plot.margin=unit(c(0,0,0,0),\"pt\"))\n";
	}
} else {
	if (defined $opt{'A'}) {
	print R "
	theme_bw()\n";
	} else {
	print R "
	theme_bw(legend.position=\"none\")\n";
	}
}

# plot them
print R "
ggsave(plot=PointPlot,filename=\"$opt{'O'}.plot.png\",width=5,height=4,dpi=900)
ggsave(plot=PointPlot,filename=\"$opt{'O'}.plot.pdf\",width=5,height=4)";

# DIMS PLOT - "CENTIPEDE"
print R "
# Make the odim and pdim \"centipede\" style plot
$gradient_function
CentPlot<-ggplot() +";

# plot the original dims w/ respective color
print R "geom_segment(aes(x=IN\$odim1,xend=IN\$pdim1,y=IN\$odim2,yend=IN\$pdim2),size=0.2,color=\"lightsteelblue4\") +
	geom_point(aes(IN\$odim1,IN\$odim2),color=\"lightsteelblue4\",size=0.5) +
	geom_point(aes(IN\$pdim1,IN\$pdim2,color=IN\$lambda),size=0.75) +
	scale_color_gradientn(colours=gradient_funct(15)) +
	labs(color=\"Lambda\") +";

# set theme
if ($theme =~ /Clean/i) {
	if (defined $opt{'A'}) {
	print R "
	theme_bw() +
	theme(panel.border=element_blank(),
		  panel.grid=element_blank(),
		  axis.line=element_blank(),
		  axis.ticks=element_blank(),
		  legend.background=element_blank(),
		  legend.title=element_blank(),
		  panel.background=element_blank(),
		  axis.text=element_blank(),
		  axis.title.x=element_blank(),
		  axis.title.y=element_blank(),
		  plot.background=element_blank(),
		  plot.margin=unit(c(0,0,0,0),\"pt\"))\n";
	} else {
	print R "
	theme_bw() +
	theme(panel.border=element_blank(),
		  panel.grid=element_blank(),
		  axis.line=element_blank(),
		  axis.ticks=element_blank(),
		  legend.position=\"none\",
		  panel.background=element_blank(),
		  axis.text=element_blank(),
		  axis.title.x=element_blank(),
		  axis.title.y=element_blank(),
		  plot.background=element_blank(),
		  plot.margin=unit(c(0,0,0,0),\"pt\"))\n";
	}
} else {
	if (defined $opt{'A'}) {
	print R "
	theme_bw()\n";
	} else {
	print R "
	theme_bw(legend.position=\"none\")\n";
	}
}

# plot them
print R "
ggsave(plot=CentPlot,filename=\"$opt{'O'}.vector_plot.png\",width=5,height=4,dpi=900)
ggsave(plot=CentPlot,filename=\"$opt{'O'}.vector_plot.pdf\",width=5,height=4)";

# DENSITY LAMBDA PLOT
print R "
# Plot the density over lambda

DensityPlot<-ggplot() +";
if (!defined $opt{'c'} && !defined $opt{'C'} && !defined $opt{'A'}) {
	print R "
	geom_density(aes(IN\$lambda),color=\"lightsteelblue4\",size=1) +";
} else {
	print R "
	geom_density(aes(IN\$lambda,color=IN\$annot),size=1) +
	guides(colour = guide_legend(override.aes = list(size=4))) +";
}
if ($color_mapping !~ /none/i) {
	print R "
	scale_colour_manual(values = c($color_mapping)) +";
}
print R "
	theme_bw() + xlab(\"Lambda\") +	ylab(\"Density\") +";
if (defined $opt{'A'}) {
print R "
	theme(legend.background=element_blank(),
		legend.title=element_blank(),
		panel.background=element_blank(),
		plot.background=element_blank())";
} else {
print R "
	theme(legend.position=\"none\",
		panel.background=element_blank(),
		plot.background=element_blank())";
}
print R "

ggsave(plot=DensityPlot,filename=\"$opt{'O'}.lambda.png\",width=5,height=3,dpi=900)
ggsave(plot=DensityPlot,filename=\"$opt{'O'}.lambda.pdf\",width=5,height=3)";

# VIOLIN LAMBDA PLOT
if (defined $opt{'A'}) {
print R "
# Horizontal violin plot over lambda

Violin<-ggplot() + theme_bw() +
	geom_violin(aes(IN\$annot,IN\$lambda,fill=IN\$annot),alpha=0.5,color=\"gray30\",size=0.5) +
	geom_jitter(aes(IN\$annot,IN\$lambda),color=\"gray30\",size=0.15) +";
if ($color_mapping !~ /none/i) {
	print R "
	scale_fill_manual(values = c($color_mapping)) +";
}
print R "
	scale_x_discrete(limits = rev(levels(IN\$V2))) +
	guides(fill=FALSE,colour=FALSE) +
	xlab(\"Annot\") +
	ylab(\"Lambda\") +
	coord_flip()
ggsave(plot=Violin,filename=\"$opt{'O'}.violin.png\",width=7,height=3,dpi=900)
ggsave(plot=Violin,filename=\"$opt{'O'}.violin.pdf\",width=7,height=3)";
}

close R;

system("$Rscript $opt{'O'}.plot.r");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.plot.r $opt{'O'}.plot.txt");
}

exit;
}

if ($command eq "plot-factors") {

getopts("O:XR:", \%opt);

$die2 = "
scitools plot-factors [options] [factors.txt file]

plot-factors  plots  k-range vs error txt

Options:
   -O   [STR]   Output prefix (default is [input].plot.png and [input].plot.pdf)
   -X           Retain intermediate files (def = delete)
   -R   [STR]   Rscript call (def = $Rscript)

Note: Requires ggplot package 
";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};

open R, ">$opt{'O'}.plot.factors.r";
print R "
#read packages
library(ggplot2)
#input
k_input<-read.table(file=\"$ARGV[0]\",header=TRUE)

PLT<-ggplot(data=k_input,aes(x=k,y=recon.err)) + theme_bw() +
	geom_point() + geom_line() +
	xlab(\"Number of factors\") + ylab(\"Reconstruction error\") +
	theme(strip.background=element_rect(fill=\"transparent\")) +
	theme(strip.background=element_blank(),
		panel.grid=element_blank(),
		panel.background=element_blank(),
		plot.background=element_blank())
ggsave(plot=PLT,filename=\"$opt{'O'}.plot.png\",width=5,height=4,dpi=900)
ggsave(plot=PLT,filename=\"$opt{'O'}.plot.pdf\",width=5,height=4);
";

close R;

system("$Rscript $opt{'O'}.plot.factors.r");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.plot.factors.r");
}

exit;
}

if ($command eq "plot-reads") {

# Defaults
$height = 6;
$width = 6;
$pt_size = 0.5;
$flanking_size = 100000;
$gene_scale_factor = 1;
$gene_text_size = 1.5;

getopts("O:A:a:C:c:R:Xs:Dh:w:rp:B:G:S:f:t:V:v:", \%opt);

$die2 = "
scitools plot-reads [options] [rmdup sci bam file] [chrN:start-end] [region 2] ...

Options:
   -O   [STR]   Output prefix (default is input prefix)
   -D           Output to a directory (will add to it if it exists)
   -B   [BED]   Bed file of peaks (optional)
   -G   {STR]   Gene info (refGene.txt formats)
                Shortcut eg: hg38, hg19, mm10 (see scitools -h for more details)
   -S   [INT]   Flanking size if gene names are specified
                (reuires -G, def = $flanking_size)
   -A   [STR]   Annotation file (to color code reads)
   -a   [STR]   Comma separated list of annotations to include in plot
                  (requires -A, will plot annots in the specified order)
   -C   [STR]   Color coding file (annot (tab) #hexColor)
   -c   [STR]   Color coding string
                  Annot=#hexColor,Annot2=#hexColor
   -h   [IN]    Height (inches, def = $height)
   -w   [IN]    Width (inches, def = $width)
   -p   [FLT]   Point size (def = $pt_size)
   -f   [FLT]   Gene plot spacing factor (def = $gene_scale_factor)
                  (for -G, larger values = more vertical spread)
   -t   [FLT]   Gene name text size (def = $gene_text_size)
   -r           Order cells by start of first read (def = randomize)
   -V   [STR]   Values or lambda file for cell ordering
                  (overrides -r, only cells in file will plot)
   -v           Reverse order of the provided values/lambda file
   -R   [STR]   Rscript call (def = $Rscript)
   -s   [STR]   Samtools call (def = $samtools)
   -X           Do not delete intermediate files (def = delete)

Note: Requires ggplot2 R package

";

if (!defined $ARGV[1]) {die $die2};
if (defined $opt{'a'} && !defined $opt{'A'}) {die "\nMust provide an annotaiton file (-A) if specifying annotations to plot (-a)!\n$die2"};
if (defined $opt{'C'} && defined $opt{'c'}) {die "\nSpecify either a color string (-c) or a color coding file (-C), not both!\n$die2"};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (defined $opt{'s'}) {$samtools = $opt{'s'}};
if (defined $opt{'h'}) {$height = $opt{'h'}};
if (defined $opt{'w'}) {$width = $opt{'w'}};
if (defined $opt{'S'}) {$flanking_size = $opt{'S'}};
if (defined $opt{'f'}) {$gene_scale_factor = $opt{'f'}};
if (defined $opt{'t'}) {$gene_text_size = $opt{'t'}};
if (defined $opt{'A'}) {
	read_annot($opt{'A'});
	foreach $annot (keys %ANNOT_count) {
		$ANNOT_ids{$annot} = 0;
	}
} else {
	$ANNOT_ids{"Cell"} = 0;
}
if (defined $opt{'a'}) {
	@ANNOT_LIST = split(/,/, $opt{'a'});
	foreach $annot (@ANNOT_LIST) {
		$ANNOT_include{$annot} = 1;
	}
}
if (defined $opt{'C'}) {read_color_file($opt{'C'})};
if (defined $opt{'c'}) {read_color_string($opt{'c'})};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.bam$//;

if (defined $opt{'G'}) {
	if (defined $REF{$opt{'G'}}) {
		$ref_file = $REF{$opt{'G'}};
		$opt{'G'} = $ref_file;
		$opt{'G'} =~ s/\.fa$/\.refGene.txt/;
	}
	read_refgene($opt{'G'});
}

if (defined $opt{'D'}) {
	if (-e "$opt{'O'}.read_plots") {
		print STDERR "INFO: Output directory exists - adding new regions to $opt{'O'}.read_plots.\n";
	} else {
		system("mkdir $opt{'O'}.read_plots");
	}
	$opt{'O'} .= ".read_plots/";
} else {$opt{'O'} .= "read_plot\."};

$bam = shift(@ARGV);
if (-e "$bam.bai") {
#	print STDERR "\n\nINFO: bai file found.\n";
} else {
	print STDERR "\n\nINFO: No bai file found - generating.\n";
	system("$samtools index $bam");
}

# LOOP THROUGH SPECIFIED REGIONS
foreach $region (@ARGV) {
print STDERR "INFO: Processing $region\n";
if ($region =~ /.+:.+-.+/) { # region
	($chr,$start,$end) = split(/[:-]/, $region);
	$region_name = "$chr\_$start\_$end";
} else { # gene?
	if (!defined $opt{'G'}) {
		print STDERR "WARNING: Cannot interpret coordinates $region - is it a gene name? If so, specify -G\n";
		$region = "skip";
	} else {
		if (defined $GENENAME_coords{$region}) {
			$region_name = $region; $region_name =~ s/[-:,\.\+\|\(\)]/_/g;
			$region = $GENENAME_coords{$region};
			($chr,$start,$end) = split(/[:-]/, $region);
			$region = "$chr:".($start-$flanking_size)."-".($end+$flanking_size);
		} elsif (defined $GENEID_coords{$region}) {
			$region_name = $region; $region_name =~ s/[-:,\.\+\|\(\)]/_/g;
			$region = $GENEID_coords{$region};
			($chr,$start,$end) = split(/[:-]/, $region);
			$region = "$chr:".($start-$flanking_size)."-".($end+$flanking_size);
		} else {
			print STDERR "WARNING: Cannot interpret coordinates $region - if it is a gene name or accession, it is not in $opt{'G'}\n";
			$region = "skip";
		}
	}
}

if ($region ne "skip") {
($chr,$start,$end) = split(/[:-]/, $region);
open OUT, ">$opt{'O'}$region_name.data";
open IN, "$samtools view $bam $region |";
$reads_in_region = 0; $cell_ct = 0; %CELLS_present = ();
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$cellID = $P[0]; $cellID =~ s/:.+$//; 
	if (!defined $CELLS_present{$cellID}) {
		$CELLS_present{$cellID} = 1;
		$cell_ct++;
	}
	$annot = "null";
	if (!defined $opt{'A'}) {
		$annot = "Cell"; $CELLID_annot{$cellID} = "Cell";
	} elsif (defined $CELLID_annot{$cellID}) {
		if ((defined $opt{'a'} && defined $ANNOT_include{$CELLID_annot{$cellID}}) || !defined $opt{'a'}) {
			$annot = $CELLID_annot{$cellID};
		}
	}
	if ($annot ne "null") {
		$posn = int($P[3]+(length($P[9])/2));
		if (defined $CELLID_id{$cellID}) {
			$id = $CELLID_id{$cellID};
		} else {
			$id = $ANNOT_ids{$annot}; $ANNOT_ids{$annot}++;
			$CELLID_id{$cellID} = $id;
			if (!defined $opt{'r'}) {$ID_randVal{$cellID} = rand(1e20)};
		}
		print OUT "$annot\t$id\t$posn\n";
		$reads_in_region++;
	}
} close IN; close OUT;

if ($reads_in_region>0) {

# order annot if -a specified, otherwise order however
@ANNOT_ORDER = ();
if (!defined $opt{'a'}) {
	foreach $annot (sort keys %ANNOT_ids) {
		push @ANNOT_ORDER, $annot;
	}
} else {
	for ($i = 0; $i < @ANNOT_LIST; $i++) {
		if (defined $ANNOT_ids{$ANNOT_LIST[$i]}) {
			push @ANNOT_ORDER, $annot;
		}
	}
}

# DO reordering stuff here
$newID = 0;
if (!defined $opt{'V'}) {
	for ($i = 0; $i < @ANNOT_ORDER; $i++) {
		$annot = $ANNOT_ORDER[$i];
		if (!defined $opt{'r'}) {
			foreach $cellID (sort {$ID_randVal{$a}<=>$ID_randVal{$b}} keys %ID_randVal) {
				if ($CELLID_annot{$cellID} eq $annot) {
					$ANNOT_CELLID_newID{$annot}{$CELLID_id{$cellID}} = $newID;
					$newID++;
				}
			}
		} else {
			foreach $cellID (sort {$CELLID_id{$a}<=>$CELLID_id{$b}} keys %CELLID_id) {
				if ($CELLID_annot{$cellID} eq $annot) {
					$ANNOT_CELLID_newID{$annot}{$CELLID_id{$cellID}} = $newID;
					$newID++;
				}
			}
		}
	}
} else {
	read_values($opt{'V'});
	if (!defined $opt{'v'}) {
		foreach $cellID (sort {$CELLID_value{$a}<=>$CELLID_value{$b}} keys %CELLID_value) {
			if (defined $CELLID_id{$cellID} && defined $CELLID_annot{$cellID}) {
				$ANNOT_CELLID_newID{$CELLID_annot{$cellID}}{$CELLID_id{$cellID}} = $newID;
				$newID++;
			}
		}
	} else {
		foreach $cellID (sort {$CELLID_value{$b}<=>$CELLID_value{$a}} keys %CELLID_value) {
			if (defined $CELLID_id{$cellID} && defined $CELLID_annot{$cellID}) {
				$ANNOT_CELLID_newID{$CELLID_annot{$cellID}}{$CELLID_id{$cellID}} = $newID;
				$newID++;
			}
		}
	}
}

open IN, "$opt{'O'}$region_name.data";
open OUT, ">$opt{'O'}$region_name.data.full";
while ($l = <IN>) {
	chomp $l;
	($annot,$id,$posn) = split(/\t/, $l);
	if (defined $ANNOT_CELLID_newID{$annot}{$id}) {
		$newID = $ANNOT_CELLID_newID{$annot}{$id};
		print OUT "READ\t$annot\t$newID\t$posn\n";
	}
} close IN;

if (defined $opt{'B'}) {
	$peaks_in_region = 0;
	open IN, "$opt{'B'}";
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		if ($P[0] eq $chr && $P[1] > $start && $P[2] < $end) {
			$posn = int(($P[1]+$P[2])/2);
			print OUT "PEAK\tpeak\t0\t$posn\n";
			$peaks_in_region++;
		}
	}
	if ($peaks_in_region==0) {
		print STDERR "WARNING: There are no peaks in $region_name, proceeding without peak plotting.\n";
	}
}

close OUT;
system("rm -f $opt{'O'}$region_name.data && mv $opt{'O'}$region_name.data.full $opt{'O'}$region_name.data");

if (defined $opt{'G'}) {
	$scale_increment = int(($cell_ct/80)*$gene_scale_factor);
	$gene_num = -1*$scale_increment;
	$genes_in_region = 0;
	open IN, "$opt{'G'}";
	open OUT, ">$opt{'O'}$region_name.gene_data";
	while ($l = <IN>) {
		chomp $l;
		if ($l !~ /^#/) {
			@P = split(/\t/, $l);
			if ($P[2] eq $chr && !defined $GENES{$P[12]}) {
				if ($P[4]>=$start && $P[5]<=$end) { # internal to window
					$gene_start = $P[4]; $gene_end = $P[5];
				} elsif ($P[4]<$start && $P[5]>=$end) { # starts before and ends after
					$gene_start = $start; $gene_end = $end;
				} elsif ($P[4]<$start && $P[5]>$start) { # starts before, ends internal
					$gene_start = $start; $gene_end = $P[5];
				} elsif ($P[4]>$start && $P[4]<$end) { # starts internal, ends after
					$gene_start = $P[4]; $gene_end = $end;
				} else { #external
					$gene_start = "null";
				}
				if ($gene_start ne "null") { # use it
					$genes_in_region++;
					$gene_num-=$scale_increment;
					$GENES{$P[12]} = 1;
					if ($P[3] eq "+") {
						print OUT "T\t$gene_start\t$gene_end\t".($gene_start-5000)."\t$P[12]\t$gene_num\n";
					} else {
						print OUT "T\t$gene_end\t$gene_start\t".($gene_start-5000)."\t$P[12]\t$gene_num\n";
					}
					@E_STARTS = split(/,/, $P[9]);
					@E_ENDS = split(/,/, $P[10]);
					for ($exon = 0; $exon < @E_STARTS; $exon++) {
						$e_start = $E_STARTS[$exon]; $e_end = $E_ENDS[$exon];
						if ($e_start>=$start && $e_end<=$end) { # internal to window
							# keep values
						} elsif ($e_start<$start && $e_end>=$end) { # starts before and ends after
							$e_start = $start; $e_end = $end;
						} elsif ($e_start<$start && $e_end>$start) { # starts before, ends internal
							$e_start = $start;
						} elsif ($e_start>$start && $e_start<$end) { # starts internal, ends after
							$e_end = $end;
						} else { #external
							$e_start = "null";
						}
						if ($e_start ne "null") {
							print OUT "E\t$e_start\t$e_end\t$exon\t$P[12]\t$gene_num\n";
						}
					}
				}
			}
		}
	} close IN; close OUT;
}


open R, ">$opt{'O'}$region_name.r";
print R "# plotting region $region
library(ggplot2)
IN <- subset(read.table(\"$opt{'O'}$region_name.data\"),V1==\"READ\")
colnames(IN) <- c(\"type\",\"annot\",\"y\",\"x\")";

if (defined $opt{'B'} && $peaks_in_region>0) {
	print R "
PEAKS <- subset(read.table(\"$opt{'O'}$region_name.data\"),V1==\"PEAK\")
colnames(PEAKS) <- c(\"type\",\"annot\",\"y\",\"x\")";
}

if (defined $opt{'G'} && $genes_in_region>0) {
	print R "
TRANS <- subset(read.table(\"$opt{'O'}$region_name.gene_data\"),V1==\"T\")
colnames(TRANS) <- c(\"type\",\"xstart\",\"xend\",\"namepos\",\"name\",\"y\")
EXONS <- subset(read.table(\"$opt{'O'}$region_name.gene_data\"),V1==\"E\")
colnames(EXONS) <- c(\"type\",\"start\",\"end\",\"exon\",\"name\",\"y\")";
}

print R "
PLT <- ggplot() + theme_bw() +";

if (defined $opt{'B'} && $peaks_in_region>0) {
	print R "
	geom_vline(aes(xintercept=PEAKS\$x),size=0.75,color=\"gray95\") +";
}

if (defined $opt{'G'} && $genes_in_region>0) {
	print R "
	geom_segment(aes(x=TRANS\$xstart,xend=TRANS\$xend,y=TRANS\$y,yend=TRANS\$y),size=0.3,color=\"darkblue\",arrow=arrow(length=unit(3,\"points\"))) +
	geom_segment(aes(x=EXONS\$start,xend=EXONS\$end,y=EXONS\$y,yend=EXONS\$y),size=1,color=\"darkblue\") +
	geom_text(aes(x=TRANS\$namepos,y=TRANS\$y,label=TRANS\$name),size=$gene_text_size,color=\"black\",hjust=1) +";
}

if (!defined $opt{'c'} && !defined $opt{'C'} && !defined $opt{'A'}) {
	print R "
	geom_point(aes(IN\$x,IN\$y),color=\"lightsteelblue4\",size=$pt_size,shape=15) +";
} else {
	print R "
	geom_point(aes(IN\$x,IN\$y,color=IN\$annot),size=$pt_size,shape=15) +"
}

if ($color_mapping !~ /none/i) {
	print R "
	scale_colour_manual(values = c($color_mapping)) +";
}

print R "
	scale_y_reverse(expand=c(0.02,0.02)) +
	scale_x_continuous(expand=c(0,0)) +
	xlab(\"$region_name\") +
	ylab(\"Cells\") +
	guides(colour = guide_legend(override.aes = list(size=4))) +
	theme(strip.background=element_blank(),";
if (defined $opt{'B'} && $peaks_in_region>0) {
print R "
		panel.grid.major.x=element_blank(),
		panel.grid.minor.x=element_blank(),";
}
if (defined $opt{'A'}) {
	print R "
		axis.ticks.y=element_blank(),
		axis.text.y=element_blank(),
		legend.title=element_blank(),
		panel.background=element_blank(),
		plot.background=element_blank())";
} else {
	print R "
		axis.ticks.y=element_blank(),
		axis.text.y=element_blank(),
		legend.position=\"none\",
		panel.background=element_blank(),
		plot.background=element_blank())";
}
print R "
ggsave(PLT,filename=\"$opt{'O'}$region_name.png\",width=$width,height=$height,dpi=900)
ggsave(PLT,filename=\"$opt{'O'}$region_name.pdf\",width=$width,height=$height)
";

system("$Rscript $opt{'O'}$region_name.r");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}$region_name.r $opt{'O'}$region_name.*data*");
}

} else {
	print STDERR "WARNING: $region_name has no reads in the provided bam file. Skipping!\n";
}
}
}

exit;
}

if ($command eq "plot-signal" || $command eq "signal-plot") {

# Defaults
$gradient_def = "WtPu";
$min = 0;
$max = 5;
$width_spec = "each=1";
$height=6;

getopts("O:r:G:R:Xh:w:", \%opt);

$die2 = "
scitools plot-signal [options] [(zscored) signal file]
   or    signal-plot

Options:
   -O   [STR]   Output prefix (default is signal file prefix)
   -r   [m,M]   Range: Min,Max values to impose for plotting (def = $min,$max)
   -G   [GRD]   Color gradient (def = $gradient_def)
                  For all available gradients, run 'scitools gradient'
   -h   [INS]   Height of plot (inches, def = $height)
   -w   [STR]   Width of plot (in, each=[for each annot] or all=[total width]
                  if a single number will assume total width def = $width_spec)
   -R   [STR]   Rscript call (def = $Rscript)
   -X           Do not delete intermediate files (def = delete)

Note: Requires ggplot2 R package

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
if (!defined $opt{'G'}) {$opt{'G'} = $gradient_def};
$gradient_function = get_gradient($opt{'G'});
if (defined $opt{'r'}) {($min,$max) = split(/,/, $opt{'r'})};
if (defined $opt{'h'}) {$height = $opt{'h'}};
if (defined $opt{'w'}) {$width_spec = $opt{'w'}};
if ($width_spec =~ /=/) {
	($width_type,$width_value) = split(/=/, $width_spec);
} else {
	$width_type = "all"; $width_value = $width_spec;
}

open IN, "$ARGV[0]";
open OUT, ">$opt{'O'}.plot";
$h = <IN>; chomp $h; @H = split(/\t/, $h);
$win_pos = 0;
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$winID = shift(@P);
	for ($i = 0; $i < @P; $i++) {
		if ($P[$i] < $min) {$P[$i] = $min};
		if ($P[$i] > $max) {$P[$i] = $max};
		print OUT "$win_pos\t$i\t$P[$i]\n";
	}
	$win_pos++;
} close OUT;

# determine plot width
%ANNOT_list = (); $annot_count = 0;
for ($i = 0; $i < @H; $i++) {
	($annot,$subWin) = split(/_window_/, $H[$i]);
	if (!defined $ANNOT_list{$annot}) {
		$ANNOT_list{$annot} = 1;
		$annot_count++;
	}
}
if ($width_type =~ /all/i) {
	$width = $width_value;
} elsif ($annot_count <= 1) {
	$width = $width_value+0.5;
} else {
	$width = ($width_value*$annot_count)+0.5;
}

open R, ">$opt{'O'}.plot.r";
print R "
library(ggplot2)
$gradient_function
IN<-read.table(\"$opt{'O'}.plot\")
colnames(IN)<-c(\"row\",\"col\",\"val\")
PLT<-ggplot() + theme_bw() +
	geom_tile(aes(IN\$col,IN\$row,fill=IN\$val)) +
	theme(axis.title.x=element_blank(),
		axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
		axis.title.y=element_blank(),
		axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(),
		axis.line=element_blank(),
		panel.background=element_blank()) +
	labs(fill=\"Signal\") +
	scale_y_reverse(expand=c(0,0)) +
	scale_x_continuous(expand=c(0,0)) +
	scale_fill_gradientn(colours=gradient_funct(99))
ggsave(plot=PLT,filename=\"$opt{'O'}.plot.png\",width=$width,height=$height,dpi=900)
ggsave(plot=PLT,filename=\"$opt{'O'}.plot.pdf\",width=$width,height=$height)
";
close R;

system("$Rscript $opt{'O'}.plot.r");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.plot.r $opt{'O'}.plot");
}

exit;
}


########## DATA ARCHIVE FUNCTIONS ##########
if ($command eq "combine-data" || $command eq "data-combine") { # TODO (low priority): Add in hdf5 version

# Defaults
$maxDim = 15;

getopts("O:D:A:R:z", \%opt);

$die2 = "
scitools combine-data [output_file] [data1_name],[data1_file] [data2_name],[data2_file] etc...
   or    data-combine

Will report all data from files into the specified output file. If the data name is not
specified, the file name will be used instead.

Will auto detect the following filetypes from names:
   annot, dims, matrix, lambda, values, complexity
   
For cells without information in a file, NA will be reported.

Options:
   -D   [INT]   Maximum dimension to include for dims files (def = $maxDim)
   -A   [STR]   Only include cells within specified annot file.
   -a   [STR]   Include only specified annotations (in -A), comma sep.
   -R   [STR]   Rename cells using rename.annot file
                (if absent will create as the provided file name)
   -z           Gzip the output

";

if (defined $opt{'O'}) {unshift @ARGV, $opt{'O'}};
if (!defined $ARGV[1]) {die $die2};
if (defined $opt{'D'}) {$maxDim = $opt{'D'}};
if (defined $opt{'a'} && !defined $opt{'A'}) {die "\nMust provide an annotaiton file (-A) if specifying annotations to filter (-a)!\n$die2"};

if (defined $opt{'A'}) {read_annot($opt{'A'})};

if (defined $opt{'a'}) {
	@ANNOT_LIST = split(/,/, $opt{'a'});
	foreach $annot (@ANNOT_LIST) {
		$ANNOT_include{$annot} = 1;
	}
} else {
	foreach $annot (keys %ANNOT_count) {
		$ANNOT_include{$annot} = 1;
	}
}

# get full file type info and cells present
$included_cell_ct = 0; $cellID_out = ""; @INCLUDED_CELLIDS = ();
for ($i = 1; $i < @ARGV; $i++) {
	if ($ARGV[$i] =~ /[,=]/) {
		($name,$file) = split(/[,=]/, $ARGV[$i]);
	} else {
		$name = $ARGV[$i]; $file = $ARGV[$i];
	}
	if (-e "$file") {
		if ($file =~ /\.annot$/) {$CLASS{$i} = "annot"}
		elsif ($file =~ /\.dims$/) {$CLASS{$i} = "dims"}
		elsif ($file =~ /\.complexity\.txt$/) {$CLASS{$i} = "complexity"}
		elsif ($file =~ /\.lambda$/) {$CLASS{$i} = "lambda"}
		elsif ($file =~ /\.values$/) {$CLASS{$i} = "values"}
		elsif ($file =~ /\.matrix$/) {$CLASS{$i} = "matrix"}
		else {
			print STDERR "\nWARNING: Cannot determine file type for file: $file, skipping!\n";
			$SKIP{$i} = 1;
		}
		if (!defined $SKIP{$i}) {
			$FILES{$i} = $file;
			$NAMES{$i} = $name;
			if ($CLASS{$i} =~ /(annot|dims|lambda|values|compelxity)/) {
				open IN, "$file";
				while ($l = <IN>) {
					chomp $l;
					@P = split(/\t/, $l);
					$cellID = $P[0];
					if (!defined $opt{'A'} || defined $ANNOT_include{$CELLID_annot{$cellID}}) {
						if (!defined $CELLID_include{$cellID}) {
							$CELLID_include{$cellID} = 1;
							$included_cell_ct++; $cellID_out .= "$cellID\t";
							push @INCLUDED_CELLIDS, $cellID;
						} else {
							$CELLID_include{$cellID}++;
						}
					}
				}
				close IN;
			} else {
				open IN, "$file"; $h = <IN>; close IN; chomp $h;
				@H = split(/\t/, $h);
				foreach $cellID (@H) {
					if (!defined $opt{'A'} || defined $ANNOT_include{$CELLID_annot{$cellID}}) {
						if (!defined $CELLID_include{$cellID}) {
							$CELLID_include{$cellID} = 1;
							$included_cell_ct++; $cellID_out .= "$cellID\t";
						push @INCLUDED_CELLIDS, $cellID;
						} else {
							$CELLID_include{$cellID}++;
						}
					}
				}
			}
		}
	} else {
		print STDERR "\nWARNING: Cannot find file: $file, skipping!\n";
	}
}
if ($included_cell_ct==0) {
	die "\nERROR: Included cell count is 0! Check your file cellID compatability!\n";
} elsif ($included_cell_ct<100) {
	print STDERR "\nWARNING: Included cell count is < 100, if this is less than expected, re-check your files and cellIDs!\n";
}

# rename cellIDs
%CELLID_cellName = ();
if (defined $opt{'R'}) {
	if (-e "$opt{'R'}") {
		open RENAME, "$opt{'R'}";
		while ($l = <RENAME>) {
			chomp $l;
			($cellID,$cellName) = split(/\t/, $l);
			$CELLID_cellName{$cellID} = $cellName;
		} close RENAME;
	} else {
		$opt{'R'} =~ s/\.annot$//; $opt{'R'} =~ s/\.rename$//;
		$cellNum = 1;
		open RENAME, ">$opt{'R'}.rename.annot";
		for ($i = 0; $i < @INCLUDED_CELLIDS; $i++) {
			$cellName = "CellID_".sprintf("%09s", $cellNum);
			$cellNum++;
			$cellID = $INCLUDED_CELLIDS[$i];
			$CELLID_cellName{$cellID} = $cellName;
			print RENAME "$cellID\t$cellName\n";
		} close RENAME;
	}
	$cellID_out = "";
	for ($i = 0; $i < @INCLUDED_CELLIDS; $i++) {
		$cellID = $INCLUDED_CELLIDS[$i];
		if (defined $CELLID_cellName{$cellID}) {
			$cellID_out .= "$CELLID_cellName{$cellID}\t";
		} else {
			die "\nERROR: Renaming file $opt{'R'} does not contain all cellIDs! Cannot proceed!\n";
		}
	}
}

$cellID_out =~ s/\t$//;

# open OUT and print cellIDs and some basic info
if (!defined $opt{'z'}) {
	open OUT, ">$ARGV[0]";
} else {
	$ARGV[0] =~ s/\.gz$//;
	open OUT, "| $gzip > $ARGV[0].gz";
}
if (defined $opt{'A'}) {
	print OUT "#CELL_IDS\tCELL_INCLUSION_ANNOTATION_FILE=$opt{'A'}";
	if (defined $opt{'a'}) {
		print OUT "\tINCLUDED_ANNOTATIONS=$opt{'a'}";
	}
} else {
	print OUT "#CELL_IDS";
}

for ($i = 1; $i < @ARGV; $i++) {
	print OUT "\tDATA=$ARGV[$i]";
} print OUT "\n";

print OUT "CELL_IDS\tUniqueIdentifier\t$cellID_out\n";
for ($i = 1; $i < @ARGV; $i++) {
	if (!defined $SKIP{$i}) {
		if ($CLASS{$i} eq "annot") {
			print OUT "#ANNOTATION_DATA\tNAME=$NAMES{$i}\tFILE=$FILES{$i}\n";
			read_annot($FILES{$i});
			print OUT "$NAMES{$i}\tANNOTATION";
			for ($j = 0; $j < @INCLUDED_CELLIDS; $j++) {
				$cellID = $INCLUDED_CELLIDS[$j];
				if (defined $CELLID_annot{$cellID}) {
					print OUT "\t$CELLID_annot{$cellID}";
				} else {
					print OUT "\tNA";
				}
			} print OUT "\n";
		} elsif ($CLASS{$i} eq "values") {
			print OUT "#VALUES_DATA\tNAME=$NAMES{$i}\tFILE=$FILES{$i}\n";
			read_values($FILES{$i});
			print OUT "$NAMES{$i}\tVALUES";
			for ($j = 0; $j < @INCLUDED_CELLIDS; $j++) {
				$cellID = $INCLUDED_CELLIDS[$j];
				if (defined $CELLID_value{$cellID}) {
					print OUT "\t$CELLID_value{$cellID}";
				} else {
					print OUT "\tNA";
				}
			} print OUT "\n";
		} elsif ($CLASS{$i} eq "lambda") {
			print OUT "#LAMBDA_DATA\tNAME=$NAMES{$i}\tFILE=$FILES{$i}\n";
			read_values($FILES{$i});
			print OUT "$NAMES{$i}\tLAMBDA";
			for ($j = 0; $j < @INCLUDED_CELLIDS; $j++) {
				$cellID = $INCLUDED_CELLIDS[$j];
				if (defined $CELLID_value{$cellID}) {
					print OUT "\t$CELLID_value{$cellID}";
				} else {
					print OUT "\tNA";
				}
			} print OUT "\n";
		} elsif ($CLASS{$i} eq "complexity") {
			print OUT "#COMPLEXITY_DATA\tNAME=$NAMES{$i}\tFILE=$FILES{$i}\n";
			read_complexity($FILES{$i});
			print OUT "$NAMES{$i}\tRAW_READS";
			for ($j = 0; $j < @INCLUDED_CELLIDS; $j++) {
				$cellID = $INCLUDED_CELLIDS[$j];
				if (defined $CELLID_raw_reads{$cellID}) {
					print OUT "\t$CELLID_raw_reads{$cellID}";
				} else {
					print OUT "\tNA";
				}
			}
			print OUT "\n$NAMES{$i}\tUNIQUE_READS";
			for ($j = 0; $j < @INCLUDED_CELLIDS; $j++) {
				$cellID = $INCLUDED_CELLIDS[$j];
				if (defined $CELLID_uniq_reads{$cellID}) {
					print OUT "\t$CELLID_uniq_reads{$cellID}";
				} else {
					print OUT "\tNA";
				}
			}
			print OUT "\n$NAMES{$i}\tCOMPLEXITY";
			for ($j = 0; $j < @INCLUDED_CELLIDS; $j++) {
				$cellID = $INCLUDED_CELLIDS[$j];
				if (defined $CELLID_complexity{$cellID}) {
					print OUT "\t$CELLID_complexity{$cellID}";
				} else {
					print OUT "\tNA";
				}
			} print OUT "\n$NAMES{$i}\nRANK";
			for ($j = 0; $j < @INCLUDED_CELLIDS; $j++) {
				$cellID = $INCLUDED_CELLIDS[$j];
				if (defined $CELLID_complexity_rank{$cellID}) {
					print OUT "\t$CELLID_complexity_rank{$cellID}";
				} else {
					print OUT "\tNA";
				}
			} print OUT "\n";
		} elsif ($CLASS{$i} eq "dims") {
			print OUT "#DIMENSIONS_DATA\tNAME=$NAMES{$i}\tFILE=$FILES{$i}\n";
			read_dims($FILES{$i});
			if ($Ndims<$maxDim) {$lastDim = $Ndims} else {$lastDim = $maxDim};
			for ($dim = 1; $dim < $lastDim; $dim++) {
				print OUT "$NAMES{$i}\tDIMENSION_$dim";
				for ($j = 0; $j < @INCLUDED_CELLIDS; $j++) {
					$cellID = $INCLUDED_CELLIDS[$j];
					if (defined $CELLID_DIMS{$cellID}[$dim]) {
						print OUT "\t$CELLID_DIMS{$cellID}[$dim]";
					} else {
						print OUT "\tNA";
					}
				} print OUT "\n";
			}
		} elsif ($CLASS{$i} eq "matrix") {
			print OUT "#MATRIX_DATA\tNAME=$NAMES{$i}\tFILE=$FILES{$i}\n";
			read_matrix($FILES{$i});
			foreach $feature (sort {$a cmp $b} @MATRIX_ROWNAMES) {
				print OUT "$NAMES{$i}\t$feature";
				for ($j = 0; $j < @INCLUDED_CELLIDS; $j++) {
					$cellID = $INCLUDED_CELLIDS[$j];
					if (defined $CELLID_FEATURE_value{$cellID}{$feature}) {
						print OUT "\t$CELLID_FEATURE_value{$cellID}{$feature}";
					} else {
						print OUT "\tNA";
					}
				} print OUT "\n";
			}
		}
	}
}
close OUT;

exit;
}

if ($command eq "split-data" || $command eq "data-split") {

getopts("O:", \%opt);

$die2 = "
scitools split-data [options] [input_data_file]

Will take a combined data file and split into scitools format files.
Can be gzipped or not.
   
For cells without information in a file, NA will be reported.

Options:
   -O   [STR]   Output prefix (def = [input - .data].[suffix])

";

if (!defined $ARGV[0]) {die $die2};

if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.gz//; $opt{'O'} =~ s/\.data//;

if ($ARGV[0] =~ /\.gz$/) {
	open IN, "$zcat $ARGV[0] |";
} else {
	open IN, "$ARGV[0]";
}

$info = <IN>;
$cellID_data = <IN>; chomp $cellID_data;
@CELLIDS = split(/\t/, $cellID_data);

$out_type = "null";
while ($l = <IN>) {
	chomp $l; @P = split(/\t/, $l);
	if ($l =~ /^#/) {
		if ($out_type eq "dims") {
			for ($i = 2; $i < @CELLIDS; $i++) {
				print OUT "$CELLIDS[$i]";
				foreach $dim (sort {$a<=>$b} @DIMS) {
					print OUT "\t$CELLID_DIMS{$CELLIDS[$i]}{$dim}";
				} print OUT "\n";
			} close OUT;
		} elsif ($out_type eq "matrix") {close OUT};
		$name = $P[1]; $name =~ s/NAME=//;
		if ($P[0] =~ /(ANNOTATION|LAMBDA|VALUES)/) {
			if ($P[0] eq "#ANNOTATION_DATA") {
				open OUT, ">$opt{'O'}.$name.annot";
			} elsif ($P[0] eq "#VALUES_DATA") {
				open OUT, ">$opt{'O'}.$name.values";
			} elsif ($P[0] eq "#LAMBDA_DATA") {
				open OUT, ">$opt{'O'}.$name.lambda";
			}
			$l = <IN>; chomp $l; @P = split(/\t/, $l);
			for ($i = 2; $i < @P; $i++) {
				print OUT "$CELLIDS[$i]\t$P[$i]\n";
			} close OUT;
		} elsif ($P[0] eq "#COMPLEXITY_DATA") {
			%COMP_DATA = ();
			for ($comp_line = 1; $comp_line <= 4; $comp_line++) {
				$l = <IN>; chomp $l; @P = split(/\t/, $l);
				for ($i = 2; $i < @P; $i++) {
					$COMP_DATA{$P[1]}{$CELLIDS[$i]} = $P[$i];
				}
			}
			open OUT, ">$opt{'O'}.$name.complexity.txt";
			for ($i = 2; $i < @P; $i++) {
				print OUT "$COMP_DATA{'RANK'}{$CELLIDS[$i]}\t$CELLIDS[$i]\t$COMP_DATA{'RAW_READS'}{$CELLIDS[$i]}\t$COMP_DATA{'UNIQUE_READS'}{$CELLIDS[$i]}\t$COMP_DATA{'COMPLEXITY'}{$CELLIDS[$i]}\n";
			} close OUT;
		} elsif ($P[0] eq "#DIMENSIONS_DATA") {
			open OUT, ">$opt{'O'}.$name.dims";
			$out_type = "dims";
			%CELLID_DIMS = (); @DIMS = ();
		} elsif ($P[0] eq "#MATRIX_DATA") {
			open OUT, ">$opt{'O'}.$name.matrix";
			print OUT "$CELLIDS[2]";
			for ($i = 3; $i < @CELLIDS; $i++) {
				print OUT "\t$CELLIDS[$i]";
			} print OUT "\n";
			$out_type = "matrix";
		}
	} else {
		if ($out_type eq "dims") {
			($null,$dim) = split(/_/, $P[1]);
			push @DIMS, $dim;
			for ($i = 2; $i < @P; $i++) {
				$CELLID_DIMS{$CELLIDS[$i]}{$dim} = $P[$i];
			}
		} elsif ($out_type eq "matrix") {
			print OUT "$P[1]";
			for ($i = 2; $i < @P; $i++) {
				print OUT "\t$P[$i]";
			} print OUT "\n";
		}
	}
} close IN;

exit;
}


########################################################################
#################### HELP TEXT & OTHER SUCH OPTIONS ####################
########################################################################

if ($command eq "help") {

$help = "

scitools Version: $version ($version_info{$version})
Adey Lab (www.adeylab.org)

scitools is a set of scripts designed for working with single-cell
combinatorial indexing data. It includes tools to go from fastq files
off of the sequencer (after bcl2fastq) to a processed dataset. It
includes a number of external tools and R packages that are called
by scitools. If scitools is used, be sure to cite those tools!

scitools commands are in the form of [class]-[operation]. Most can
be specified in reverse order, or if the operation is unique to the
class of files or the class type can be determined by the files in
the arguments then just the operation name can be specified.

Dependencies: (command-line callable, can be specified in options)
Run scitools dependencies to check software and R requirements.

Config file for defaults: scitools has a number of set defaults
that can be conveniently specified in a scitools.cfg file. There
is a hardcoded config file location, but the best way for personalized
defaults, include a '.scitools.cfg' file in your home directory.

Executables:
   gzip         For gzipped fastq files. Default: $gzip & $zcat
                (these are hardoded at the beginning of the scitools
                 code and not available as options)
   samtools     Bam-related commands. Default: $samtools
   bedtools     Bed-related commands. Default: $bedtools
   Rscript      Numerous R operations. Default: $Rscript
   python       For calling UMAP. Default: $Pscript
   bwa          For alignment only. Default: $bwa
   macs2        For atac-callpeak only. Default: $macs2
   scitools     Can call itself. Default: $scitools

R packages:
   ggplot2         For plotting commands
   svd             For Latent Semantic Indexing (LSI)
   Rtsne           For tSNE visualization
   methods         For PCA
   dbscan          For density-based clustering
   princurve       Principle curve projections
   chromVAR        For chromVAR only
   chromVARmotifs  For chromVAR only
   cicero          For cicero only

Some default locations / shortcuts:
   Fastq directory (where bcl2fastq outputs fastq files)
      $VAR{'fastq_input_directory'}
   Output fastq directory (for processed SCI fastq files)
      $VAR{'SCI_fastq_directory'}
   SCI index file (should comtain all barcodes in the proper format)
      $VAR{'SCI_index_file'}

   Default Refs / shortcuts ([reference.fa], [reference.fa].fai, [reference.fa].[bwa_index], and [reference].refGene.txt)
$ref_shortcuts

File Types used by scitools:

   fastq    Standard fastq files, input can be gzipped or not, output fastq
               files will always be gzipped.
   
   bam      Bam files. For scitools, it is expected the barcode sequence
               or cellID is in the read name ([barcode]:[read_number])
               Sam files are not supported to encourage (force) space saving
			   
   annot    Annotation file: tab delimited, 2 columns
               Col 1 = cellID (barcode), could also be feature name
               Col 2 = annotation information (e.g. experimental condition)
            values file A special case of an annotaiton file where the
               annotation is a continuous variable and not discrete.
			   
   matrix   Matrix file, tab delimited
               Row 1 = CellIDs, has 1 less column than all other rows
               Rows = field 1 is name, then 1 field for each cellID
			   
   dims     Dimensions file, tab delimited
               Col 1 = cellID
               Col 2-N = dimensions for each cell
            range    Range format is not a file type but a way to input a
               range of values (e.g. which dimensions to use from a dims file)
               This is a comma separated list and can include dashes for a
               range of values (e.g. 1,3-6,8-10,13)

Typical scitools analysis for sci-ATAC-seq:

   1) Make an annotation file for demultiplexing a sequencing run:
      scitools annot-make -p   # this will output an example annot csv file, edit in excel
      scitools annot-make -P [myExperimentDescriptorFile.csv] > [mySamples.annot]

   2) Demultiplex and reformat raw fastq files:
      scitools fastq-dump -A [mySamples.annot] -R [RUN_NAME, must be in fastq directory]

   3) Align fastq files (if reference defaults are set up use the reference shortcut)
      scitools align -t [n threads] [reference_prefix or shortcut] [sampleA_prefix] [reads.sampleA.1.fq.gz] [reads.sampleA.2.fq.gz]

   4) Remove duplicates and create complexity file:
      scitools rmdup [sampleA].bam

   5) Plot complexity file to assess library performance and determine bam filters:
      scitools plot-complexity [sampleA].complexity.txt
         Note: if multiple conditions are present in the sample, create a condition annot file as in step 1
               you can then add '-A [myConditions.annot]' to the above command to plot by conditions.
               It is also possible to specify colors for plotting conditions using -C or -c

   6) Determine performance and filter bam to remove noise reads:
      scitools bam-filter -N [read threshold, ~1000] -C [sampleA].complexity.txt -c [min compl, 0],[max compl, 60]

   7) Examine index-based perfoemance on filtered bam:
      scitools index-performance [sampleA].bbrd.q10.filt.bam

   8) Call ATAC peaks (macs2):
      scitools callpeak -f [reference.fai file] [sampleA].bbrd.q10.filt.bam

   9) Make a counts matrix:
      scitools atac-counts [sampleA].bbrd.q10.filt.bam [sampleA].bbrd.q10.filt.bed

   10) Filter the counts matrix:
      scitools matrix-filter [sampleA counts matrix]

   11) Perform term-frequency inverse-document-frequency tarnsformation:
      scitools tfidf [sampleA filtered matrix]

   12) Latent semantic indexing:
      scitools lsi [sampleA tfidf matrix]

   13) Visualize via tSNE:
      scitools tsne [sampleA LSI matrix]

   14) Plot tSNE:
      scitools plot-dims -A [myConditions.annot] [sampleA LSI].tsne.dims
";

die "$help";
}

if ($command eq "gradient" || $command eq "gradients") {

print "
scitools available gradients and gradient specification. For a number 
of plotting options, a gradient can be specified. Scitools has several
default gradients shown below, or a custom gradient can be specified.
These commands will add in a color ramp function in R.

Default Gradients: (name, R color function)
";

foreach $gradient (sort keys %COLOR_GRADIENT) {
	print "   $gradient\t($COLOR_GRADIENT{$gradient})\n";
}

print "
To specify a custom gradient, specify a comma-separated list of R
colors, or hex codes:
   eg.   red,yellow,blue   would be a gradient of those three colors.
   
";

exit;
}

################## NO COMMAND FOUND ##################

die "\nERROR: It appears you have specified a command that is registered, but does not actually exist. If you see this message, please inform us of the bug.\n";

#################### END SCITOOLS ####################
