#!/usr/bin/perl

# START GLOBAL DEFAULTS
#   To tailor the defaults to your own specified locations
#   to avoid having to specify them as options for each run, change the
#   variables below. A setup script may be included in future versions to
#   perform this for you via a series of prompts.

# Location Defaults
$fastq_input_directory = "/home/groups/oroaklab/fastq"; #DEFAULT=fastq_input_directory
$SCI_fastq_directory = "/home/groups/oroaklab/demultiplex"; #DEFAULT=SCI_fastq_directory
$SCI_index_file = "/home/groups/oroaklab/adey_lab/reference_files/index_files/SCI_ALL_current.txt"; #DEFAULT=SCI_index_file

# Reference genome shortcuts
# Note: scitools looks for the following files:
#  [ref.fa].[bwa_index]
#  [ref.fa].fai
#  [ref].refGene.txt (for plot-reads)
# To enable "hg38", "hg19", and "mm10" shorcut usage, ensure all files are present
$hg19_ref = "/home/groups/oroaklab/refs/hg19/hg19.fa"; #DEFAULT=hg19_ref
$hg38_ref = "/home/groups/oroaklab/refs/hg38/hg38.fa"; #DEFAULT=hg38_ref
$mm10_ref = "/home/groups/oroaklab/refs/mm10/mm10.fa"; #DEFAULT=mm10_ref

# Software Defaults
$gzip = "gzip"; #DEFAULT=gzip
$zcat = "zcat"; #DEFAULT=zcat
$bwa = "bwa"; #DEFAULT=bwa
$samtools = "samtools"; #DEFAULT=samtools
$scitools = "scitools"; #DEFAULT=scitools
$macs2 = "macs2"; #DEFAULT=macs2
$bedtools = "bedtools"; #DEFAULT=bedtools
$Rscript = "Rscript"; #DEFAULT=Rscript

# END GLOBAL DEFAULTS

# MODULES
# Note: to decrease dependencies, other modules are not used in favor of
#   including relevant code within this script tailored to the specific
#   functions
use Getopt::Std; %opt = ();
load_gradient_defaults();

# VERSION INFO
$version = "0.1.1";
%version_info = (
	"0.1.0" => "180215, alpha - initial development",
	"0.1.1" => "180418, alpha - dev split"
);

# COMMAND DIRECTORY
%COMMANDS = (
	"help" => 1,
	"merge" => 1, "filter" => 1, "split" => 1,
	"dependencies" => 1, "depend" => 1,
	"gradient" => 1, "gradients" => 1,
	
	"fastq-dump" => 1, "dump-fastq" => 1,
	"fastq-split" => 1, "split-fastq" => 1,
	"fastq-merge" => 1, "merge-fastq" => 1,
	"fastq-align" => 1, "align-fastq" => 1, "align" => 1,
	
	"bam-bulk2sci" => 1, "bulk2sci" => 1,
	"bam-addrg" => 1, "addrg" => 1,
	"bam-rmdup" => 1, "rmdup" => 1,
	"bam-split" => 1, "split-bam" => 1,
	"bam-filter" => 1, "filter-bam" => 1,
	"bam-merge" => 1, "merge-bam" => 1,
	"bam-project" => 1, "project-bam" => 1, "project" => 1,
	"bam-aggregate" => 0, "aggregate-bam" => 0,
	
	"signal-make" => 1, "make-signal" => 1,
	"plot-signal" => 1, "signal-plot" => 1,
	
	"annot-make" => 1, "make-annot" => 1,
	"annot-merge" => 1, "merge-annot" => 1,
	
	"atac-callpeaks" => 1, "atac-callpeak" => 1, "callpeak" => 1, "callpeaks" => 1,
	"atac-mergepeaks" => 1, "atac-mergepeak" => 1, "mergepeak" => 1, "mergepeaks" => 1,
	"atac-counts" => 1, "atac-count" => 1, "count" => 1, "counts" => 1,
	"atac-chromvar" => 0, "chromvar" => 0,
	"atac-cicero" => -1, "cicero" => -1,
	
	"matrix-filter" => 1, "filter-matrix" => 1,
	"matrix-summarize" => 0, "summarize-matrix" => 0,
	"matrix-tf" => 1, "tf" => 1,
	"matrix-tfidf" => 1, "tfidf" => 1,
	"matrix-lsi" => 1, "lsi" => 1,
	"matrix-tsne" => 1, "tsne" => 1,
	"matrix-pca" => 1, "pca" => 1,
	"matrix-nmf" => -1, "nmf" => -1,
	"matrix-bicluster" => -1, "bicluster" => -1,
	"matrix-aggregate" => 1, "aggregate-matrix" => 1,
	
	"dims-kmeans" => 1, "kmeans" => 1,
	"dims-dbscan" => 1, "dbscan" => 1,
	"dims-pcurve" => 1, "pcurve" => 1,
	"pcurve-center" => 1, "center-pcurve" => 1, "lambda-center" => 1, "center-lambda" => 1,
	
	"aggregate-cells" => 0, "aggregate" => 0,

	"plot-complexity" => 1,
	"plot-dims" => 1,
	"plot-pcurve" => 1,
	"plot-reads" => 1,
	
	"index-performance" => 1, "index-perform" => 1,
	"combine-data" => 1, "data-combine" => 1,
	"split-data" => 1, "data-split" => 1

);

# GLOBAL EMPTY VARIABLES
$color_mapping = "none";

# HELP TEXT
$die = "
scitools [command] [options] [arguments]

Version: $version ($version_info{$version}), www.adeylab.org

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
   
   atac-callpeak       Call peaks on bam file using macs2
   atac-mergepeak      Merge ATAC-seq peak files
   atac-counts         Bam and peak file to a counts matrix
   atac-chromvar       Run chromVAR wrapper on sci-ATAC-seq data
   atac-cicero         Run cicero wrapper on sci-ATAC-seq data
   
   signal-make         Generate windowed signal over features from bam
   signal-plot         Plot windowed signal views
   
   matrix-summarize    Generate a summary and plots on matrix properties
   matrix-filter       Filter a sci-ATAC-seq counts matrix
   matrix-tf           Normalize only by term frequency
   matrix-tfidf        Perform tf-idf on counts matrix
   matrix-lsi          Perform Latent Semantic Indexing on matrix
   matrix-tsne         tSNE on matrix
   matrix-pca          PCA on matrix
   matrix-nmf          Non-negative Matrix Factorization of matrix
   matrix-bicluster    Bicluster and plot a heatmap
   matrix-aggregate    Aggregate cells in counts matrix by annotation
   
   dims-kmeans         Kmeans clustering on dims file
   dims-dbscan         Density-base (dbscan) clustering on dims file
   dims-pcurve         Project a principle curve through dims file
   pcurve-center       Centers and normalizes a pcurve lambda
   aggregate-cells     Aggregate cells in proximity with one another

   plot-complexity     Plot complexity data
   plot-dims           Plot tSNE or other dimensions file
   plot-pcurve         Genrate multiple princurve plots
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

###############################################################
#################### SCITOOLS COMMAND CODE ####################
###############################################################

########## AUTO-DETECT COMMON OPERATIONS ##########
if ($command eq "merge") {
	if ($ARGV[(@ARGV-1)] =~ /\.bam$/) {$command eq "bam-merge"}
	elsif ($ARGV[(@ARGV-1)] =~ /\.fq.gz$/ || $ARGV[(@ARGV-1)] =~ /\.fq$/ || $ARGV[(@ARGV-1)] =~ /\.fastq.gz$/ || $ARGV[(@ARGV-1)] =~ /\.fastq$/) {$command eq "fastq-merge"}
	elsif ($ARGV[(@ARGV-1)] =~ /\.bed$/) {$command eq "atac-mergepeak"}
	elsif ($ARGV[(@ARGV-1)] =~ /\.annot$/) {$command eq "annot-merge"}
	else {die "\n\nscitools merge was specified, but the type of files to merge cannot be detected. Please specify: bam-merge, atac-mergepeak, OR annot-merge\n\n"}
}

if ($command eq "filter") {
	if ($ARGV[(@ARGV-1)] =~ /\.bam$/) {$command eq "bam-filter"}
	elsif ($ARGV[(@ARGV-1)] =~ /\.matrix$/) {$command eq "bam-matrix"}
	else {die "\n\nscitools filter was specified, but the type of file to filter cannot be detected. Please specify: bam-filter OR matrix-filter\n\n"}
}

if ($command eq "split") {
	if ($ARGV[(@ARGV-1)] =~ /\.bam$/) {$command eq "bam-split"}
	elsif ($ARGV[(@ARGV-1)] =~ /\.fq.gz$/ || $ARGV[(@ARGV-1)] =~ /\.fq$/ || $ARGV[(@ARGV-1)] =~ /\.fastq.gz$/ || $ARGV[(@ARGV-1)] =~ /\.fastq$/) {$command eq "fastq-split"}
	else {die "\n\nscitools split was specified, but the type of file to split cannot be detected. Please specify: bam-split OR fastq-split (or data-split)\n\n"}
}


########## GENERAL FUNCTIONS ##########
if ($command eq "dependencies" || $command eq "depend") { # TODO: Add in checking for default file locations (e.g. hg38 etc...)

getopts("R:F:O:o:I:1:2:A:B:", \%opt);

$die2 = "
scitools dependencies [options] [report output file]
  or     depend

This will run a check to make sure the required dependencies
are present and command-line callable. If defaults are not
used (shown below), then specify the executables using
options. It will then verify R packages are installed and
output all software & R package status in a report file.

Options:
   -s   [STR]   Samtools call (def = $samtools)
   -b   [STR]   Bedtools call (def = $bedtools)
   -B   [STR]   Bwa call (def = $bwa)
   -m   [STR]   Macs2 call (def = $macs2)
   -R   [STR]   Rscript call (def = $Rscript)
   -S   [STR]   Scitools call (def = $scitools)
                (for self-calling)

Hardcoded:      Gzip (def = $gzip)
                Zcat (def = $zcat)

";

if (!defined $ARGV[0]) {die $die2};
if (defined $opt{'s'}) {$samtools = $opt{'s'}};
if (defined $opt{'b'}) {$bedtools = $opt{'b'}};
if (defined $opt{'B'}) {$bwa = $opt{'B'}};
if (defined $opt{'m'}) {$macs2 = $opt{'m'}};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (defined $opt{'S'}) {$scitools = $opt{'S'}};

$executable_warnings = 0;
$R_package_warnings = 0;

open OUT, ">$ARGV[0]";
$ts = localtime(time);
print OUT "$ts scitools dependencies check.

Checking command-line executables:\n";

@DEPENDENCIES = ($samtools,$bedtools,$bwa,$macs2,$Rscript,$scitools,$gzip,$zcat);
foreach $dependency (@DEPENDENCIES) {
	$path = "";
	open WHICH, "which $dependency 2>/dev/null |";
	$path = <WHICH>; close WHICH;
	if ($path =~ /\//) {
		print OUT "  \"$dependency\" found. Path = $path";
	} else {
		print OUT "  WARNING: \"$dependency\" is not command-line callable!\n    Ensure that the executable is correct. If the tool is not installed, scitools functions that use this tool will not be useable.\n";
		$executable_warnings++;
	}
}

print OUT "\nChecking R packages:\n"; close OUT;

@R_PACKAGES = ("ggplot2","svd","Rtsne","methods","dbscan","chromVAR","chromVARmotifs","cicero","princurve");
open R, ">$ARGV[0].r";
foreach $package (@R_PACKAGES) {
	print R "if (!require($package)) {print(\"  R: WARNING: $package cannot be loaded! Some scitools functions will not be useable unless it is installed.\",quote=FALSE)} else {print(\"  $package found and loadable!\",quote=FALSE)}\n";
} close R;
system ("$Rscript $ARGV[0].r >> $ARGV[0] 2>/dev/null && rm -f $ARGV[0].r");

open IN, "$ARGV[0]";
while ($l = <IN>) {if ($l =~ /R: WARNING:/) {$R_package_warnings++}};
close IN;

open OUT, ">>$ARGV[0]";
print OUT "\nDependency check completed:
  $executable_warnings failed command-line executables
  $R_package_warnings failed R packages
"; close OUT;

print "\nDependency check completed with: $executable_warnings failed command-line executables, and $R_package_warnings failed R packages.\n\n";

exit;
}

if ($command eq "index-perform" || $command eq "index-performance") {

# Defaults
$gradient_def = "BuYlRd";

getopts("O:I:R:s:t:A:xG:", \%opt);

# DEFAULTS
@LETTERS = ("0", "A", "B", "C", "D", "E", "F", "G", "H");
%LETTER_NUM = ("A"=>"1", "B"=>"2", "C"=>"3", "D"=>"4", "E"=>"5", "F"=>"6", "G"=>"7", "H"=>"8");
%WELL_xy = ();
foreach $letter (keys %LETTER_NUM) {
	for ($number = 1; $number <= 12; $number++) {
		$well = $letter.$number;
		$WELL_xy{$well} = "$LETTER_NUM{$letter}\t$number";
	}
}
$threshold = 1000000;

$die2 = "
scitools index-perform [options] [fastq or bam]

Will generate the read counts from each well of each tier of indexing.
Should be on a pre-filetered since it will plot all barcode combos present.

Options:
   -O   [STR]   Output prefix, will create a folder (def = input prefix)
   -I   [STR]   Index file
         (default = $SCI_index_file)
         (Index names must be in form of: [Tier]_[set]_[i5/i7]_[A-H/1-12])
   -A   [STR]   Annotation file (only include cell IDs in the annot file)
   -t   [INT]   Threshold of reads for a plate to include (def = $threshold)
   -x           Do not report plate-plate-well stats (just plate-well)
   -G   [GRD]   Color gradient for plots (def = $gradient_def, biased 0.65)
                For all available gradients, run 'scitools gradient'
   -R   [STR]   Rscript call (def = $Rscript)
   -s   [STR]   samtools call (if bam; def = $samtools)

Note: Requires ggplot2 R package for plotting

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.gz$//; $opt{'O'} =~ s/\.fq$//; $opt{'O'} =~ s/\.bam$//;
if (-e "$opt{'O'}.index_performance") {die "$opt{'O'}.index_performance Directory already exists! Exiting!\n\n"};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (defined $opt{'s'}) {$samtools = $opt{'s'}};
if (defined $opt{'I'}) {$SCI_index_file = $opt{'I'}};
if (defined $opt{'t'}) {$threshold = $opt{'t'}};
if (defined $opt{'A'}) {read_annot($opt{'A'})};
if (!defined $opt{'G'}) {$opt{'G'} = $gradient_def};
$gradient_function = get_gradient($opt{'G'});
$gradient_function =~ s/\)$//; $gradient_function .= ",bias=0.65)";

read_indexes($SCI_index_file);

%CELLID_count = ();
if ($ARGV[0] =~ /\.bam$/) {
	open IN, "$samtools view $ARGV[0] 2>/dev/null |";
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		$cellID = $P[0]; $cellID =~ s/:.+$//;
		$CELLID_count{$cellID}++;
	} close IN;
} else {
	if ($ARGV[0] =~ /\.gz$/) {
		open IN, "$zcat $ARGV[0] |";
	} elsif ($ARGV[0] =~ /\.fq$/) {
		open IN, "$ARGV[0]";
	} else {die "\n\nCannot determine file input type! Provide either a fastq (can be gzipped) or a bam file!\n\n$die2"};
	while ($cellID = <IN>) {
		chomp $cellID; $null = <IN>; $null = <IN>; $null = <IN>;
		$cellID =~ s/:.+$//; $cellID =~ s/^\@//;
		$CELLID_count{$cellID}++;
	} close IN;
}

$NEX_set_count = 0; $PCR_set_count = 0;
foreach $cellID (keys %CELLID_count) {
	if (!defined $opt{'A'} || defined $CELLID_annot{$cellID}) {

		$ix1 = substr($cellID,0,8);
		$ix2 = substr($cellID,8,10);
		$ix3 = substr($cellID,18,8);
		$ix4 = substr($cellID,26,10);
		
		if (!defined $INDEX_POS_SEQ_id{'1'}{$ix1} ||
			!defined $INDEX_POS_SEQ_id{'2'}{$ix2} ||
			!defined $INDEX_POS_SEQ_id{'3'}{$ix3} ||
			!defined $INDEX_POS_SEQ_id{'4'}{$ix4}) {
			print STDERR "WARNING: Barcode combo: $cellID ($ix1)($ix2)($ix3)($ix4) has a barcode not present in the index file!\n";
		}
		
		$NEX_set = $INDEX_POS_SEQ_id{'3'}{$ix3}.$INDEX_POS_SEQ_id{'1'}{$ix1};
		$PCR_set = $INDEX_POS_SEQ_id{'4'}{$ix4}.$INDEX_POS_SEQ_id{'2'}{$ix2};
		
		$NEX_well = $INDEX_POS_SEQ_well{'3'}{$ix3}.$INDEX_POS_SEQ_well{'1'}{$ix1};
		$PCR_well = $INDEX_POS_SEQ_well{'4'}{$ix4}.$INDEX_POS_SEQ_well{'2'}{$ix2};
		
		if (!defined $NEX_SET_total{$NEX_set}) {$NEX_set_count++};
		if (!defined $PCR_SET_total{$PCR_set}) {$PCR_set_count++};
		
		$NEX_SET_total{$NEX_set}+=$CELLID_count{$cellID};
		$PCR_SET_total{$PCR_set}+=$CELLID_count{$cellID};
		
		$NEX_SET_WELL_total{$NEX_set}{$NEX_well}+=$CELLID_count{$cellID};
		$PCR_SET_WELL_total{$PCR_set}{$PCR_well}+=$CELLID_count{$cellID};
		
		$NEX_SET_PCR_SET_total{$NEX_set}{$PCR_set}+=$CELLID_count{$cellID};
		
		$NEX_SET_PCR_SET_NEX_WELL_total{$NEX_set}{$PCR_set}{$NEX_well}+=$CELLID_count{$cellID};
		$NEX_SET_PCR_SET_PCR_WELL_total{$NEX_set}{$PCR_set}{$PCR_well}+=$CELLID_count{$cellID};
		
	}
}

system("mkdir $opt{'O'}.index_performance");

open SUMMARY, ">$opt{'O'}.index_performance/summary.txt";
$ts = localtime(time);
print SUMMARY "$ts scitools index-performance on $ARGV[0]
Transposase-based index totals:\n";
foreach $NEX_set (keys %NEX_SET_total) {
	print SUMMARY "  $NEX_set\t$NEX_SET_total{$NEX_set}\n";
}
print SUMMARY "PCR-based index totals:\n";
foreach $PCR_set (keys %PCR_SET_total) {
	print SUMMARY "  $PCR_set\t$PCR_SET_total{$PCR_set}\n";
}

$plates_to_plot = 0;
open OUT, ">$opt{'O'}.index_performance/plate_performance.txt";
foreach $NEX_set (keys %NEX_SET_WELL_total) {
	if ($NEX_SET_total{$NEX_set}>=$threshold) {
		foreach $NEX_well (sort keys %WELL_xy) {
			$coord = $WELL_xy{$NEX_well};
			if (defined $NEX_SET_WELL_total{$NEX_set}{$NEX_well}) {
				print OUT "NEX\t$NEX_set\tALL\t$NEX_well\t$NEX_SET_WELL_total{$NEX_set}{$NEX_well}\t$coord\n";
			} else {
				print OUT "NEX\t$NEX_set\tALL\t$NEX_well\t0\t$coord\n";
			}
		}
		$plates_to_plot++;
	}
}
foreach $PCR_set (keys %PCR_SET_WELL_total) {
	if ($PCR_SET_total{$PCR_set}>=$threshold) {
		foreach $PCR_well (sort keys %WELL_xy) {
			$coord = $WELL_xy{$PCR_well};
			if (defined $PCR_SET_WELL_total{$PCR_set}{$PCR_well}) {
				print OUT "PCR\tALL\t$PCR_set\t$PCR_well\t$PCR_SET_WELL_total{$PCR_set}{$PCR_well}\t$coord\n";
			} else {
				print OUT "PCR\tALL\t$PCR_set\t$PCR_well\t0\t$coord\n";
			}
		}
		$plates_to_plot++;
	}
}
if ($PCR_set_count > 1 && $NEX_set_count > 1 && !defined $opt{'x'}) {
	foreach $NEX_set (keys %NEX_SET_WELL_total) {
		foreach $PCR_set (keys %PCR_SET_WELL_total) {
			if ($NEX_SET_PCR_SET_total{$NEX_set}{$PCR_set}>=$threshold) {
				foreach $NEX_well (sort keys %WELL_xy) {
					$coord = $WELL_xy{$NEX_well};
					if (defined $NEX_SET_PCR_SET_NEX_WELL_total{$NEX_set}{$PCR_set}{$NEX_well}) {
						print OUT "NEX\t$NEX_set\t$PCR_set\t$NEX_well\t$NEX_SET_PCR_SET_NEX_WELL_total{$NEX_set}{$PCR_set}{$NEX_well}\t$coord\n";
					} else {
						print OUT "NEX\t$NEX_set\t$PCR_set\t$NEX_well\t0\t$coord\n";
					}
				}
				$plates_to_plot++;
			}
		}
	}
	foreach $NEX_set (keys %NEX_SET_WELL_total) {
		foreach $PCR_set (keys %PCR_SET_WELL_total) {
			if ($NEX_SET_PCR_SET_total{$NEX_set}{$PCR_set}>=$threshold) {
				foreach $PCR_well (sort keys %WELL_xy) {
					$coord = $WELL_xy{$PCR_well};
					if (defined $NEX_SET_PCR_SET_PCR_WELL_total{$NEX_set}{$PCR_set}{$PCR_well}) {
						print OUT "PCR\t$NEX_set\t$PCR_set\t$PCR_well\t$NEX_SET_PCR_SET_PCR_WELL_total{$NEX_set}{$PCR_set}{$PCR_well}\t$coord\n";
					} else {
						print OUT "PCR\t$NEX_set\t$PCR_set\t$PCR_well\t0\t$coord\n";
					}
				}
				$plates_to_plot++;
			}
		}
	}
}
close OUT;

if ($plates_to_plot<2) {
	die "\nThe number of plates passing filters to plot is zero. Try reducing the threshold to plot. (currently $threshold)\n";
}

$plot_height = 0.5+(int(($plates_to_plot/3)+1)*2.5);

open R, ">$opt{'O'}.index_performance/plot_plates.r";
print R "
library(ggplot2)
$gradient_function
IN<-read.table(\"$opt{'O'}.index_performance/plate_performance.txt\")
colnames(IN)<-c(\"Plate_type\",\"NEX_set\",\"PCR_set\",\"Plate_well\",\"Well_count\",\"Row\",\"Column\")
PLT<-ggplot(data=IN) + theme_bw() +
	geom_tile(aes(Column,Row,fill=log10(Well_count+1))) +
	scale_y_reverse(breaks=c(1,2,3,4,5,6,7,8)) +
	xlab(\"Column\") + ylab(\"Row\") +
	facet_wrap(c(\"Plate_type\",\"NEX_set\",\"PCR_set\"),ncol=3,labeller=label_wrap_gen(multi_line=FALSE)) +
	theme(strip.background=element_rect(fill=\"transparent\")) +
	scale_fill_gradientn(colours=gradient_funct(99)) +
	scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12)) +
	labs(fill=\"Log10\nReads\") +
	theme(strip.background=element_blank(),
		panel.grid=element_blank(),
		axis.line=element_blank(),
		axis.ticks=element_blank(),
		panel.background=element_blank(),
		plot.background=element_blank())
ggsave(plot=PLT,filename=\"$opt{'O'}.index_performance/plate_performance.png\",height=$plot_height,width=12)
ggsave(plot=PLT,filename=\"$opt{'O'}.index_performance/plate_performance.pdf\",height=$plot_height,width=12)
"; close R;

system("$Rscript $opt{'O'}.index_performance/plot_plates.r");

exit;
}


########## FASTQ FUNCTIONS ##########
if ($command eq "fastq-dump" || $command eq "dump-fastq") {

getopts("R:F:O:o:I:1:2:A:i:j:r:", \%opt);

$die2 = "
scitools fastq-dump [options]
   or    dump-fastq

Takes sequencer fastq files (from bcl2fastq) and will format
them into fastq files with matched barcodes.

Options:
   -R   [STR]   Run name (preferred mode)
   -r   [STR]   Run ID (if additional sequencing for some
                libraries. Will add .[ID] to end of read
                names; optional)
   -A   [STR]   Annotation file (will split fastqs)

Defaults:
   -F   [STR]   Fastq directory
         ($fastq_input_directory)
   -O   [STR]   Output fastq directory
         ($SCI_fastq_directory)
   -o   [STR]   Output prefix
         (def = run name)
   -I   [STR]   SCI index master file
         ($SCI_index_file)

To specify specific fastq files instead of defaults:
   (can be comma sep for multiple)
   -1   [STR]   Read 1 fastq
   -2   [STR]   Read 2 fastq
   -i   [STR]   Index 1 fastq (opt)
   -j   [STR]   Index 2 fastq (opt)

Will split with a hamming distance of 2

";

if (!defined $opt{'R'}) {die $die2};
if (!defined $opt{'F'}) {$opt{'F'} = $fastq_input_directory};
if (!defined $opt{'O'}) {$opt{'O'} = $SCI_fastq_directory};
if (!defined $opt{'o'}) {$opt{'o'} = $opt{'R'}};
if (!defined $opt{'I'}) {$opt{'I'} = $SCI_index_file};

open IN, "$opt{'I'}";
while ($l = <IN>) {
	chomp $l;
	($id,$pos,$seq) = split(/\t/, $l);
	$POS_SEQ_seq{$pos}{$seq} = $seq;
	$POS_length{$pos} = length($seq);
} close IN;

# make all-1-away hash
foreach $pos (keys %POS_SEQ_seq) {
	foreach $seq (keys %{$POS_SEQ_seq{$pos}}) {
		@TRUE = split(//, $seq);
		for ($i = 0; $i < @TRUE; $i++) {
			if ($TRUE[$i] =~ /A/i) {
				@NEW = @TRUE; $NEW[$i] = "C"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
				@NEW = @TRUE; $NEW[$i] = "G"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
				@NEW = @TRUE; $NEW[$i] = "T"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
				@NEW = @TRUE; $NEW[$i] = "N"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
			} elsif ($TRUE[$i] =~ /C/i) {
				@NEW = @TRUE; $NEW[$i] = "A"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
				@NEW = @TRUE; $NEW[$i] = "G"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
				@NEW = @TRUE; $NEW[$i] = "T"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
				@NEW = @TRUE; $NEW[$i] = "N"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
			} elsif ($TRUE[$i] =~ /G/i) {
				@NEW = @TRUE; $NEW[$i] = "C"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
				@NEW = @TRUE; $NEW[$i] = "A"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
				@NEW = @TRUE; $NEW[$i] = "T"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
				@NEW = @TRUE; $NEW[$i] = "N"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
			} elsif ($TRUE[$i] =~ /T/i) {
				@NEW = @TRUE; $NEW[$i] = "C"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
				@NEW = @TRUE; $NEW[$i] = "G"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
				@NEW = @TRUE; $NEW[$i] = "A"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
				@NEW = @TRUE; $NEW[$i] = "N"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
			}
		}
	}
}

# make all-2-away hash
foreach $pos (keys %POS_SEQ_seq) {
	foreach $id_seq (keys %{$POS_SEQ_seq{$pos}}) {
		$seq = $POS_SEQ_seq{$pos}{$id_seq};
		@TRUE = split(//, $seq);
		for ($i = 0; $i < @TRUE; $i++) {
			if ($TRUE[$i] =~ /A/i) {
				@NEW = @TRUE; $NEW[$i] = "C"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
				@NEW = @TRUE; $NEW[$i] = "G"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
				@NEW = @TRUE; $NEW[$i] = "T"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
				@NEW = @TRUE; $NEW[$i] = "N"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
			} elsif ($TRUE[$i] =~ /C/i) {
				@NEW = @TRUE; $NEW[$i] = "A"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
				@NEW = @TRUE; $NEW[$i] = "G"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
				@NEW = @TRUE; $NEW[$i] = "T"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
				@NEW = @TRUE; $NEW[$i] = "N"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
			} elsif ($TRUE[$i] =~ /G/i) {
				@NEW = @TRUE; $NEW[$i] = "C"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
				@NEW = @TRUE; $NEW[$i] = "A"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
				@NEW = @TRUE; $NEW[$i] = "T"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
				@NEW = @TRUE; $NEW[$i] = "N"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
			} elsif ($TRUE[$i] =~ /T/i) {
				@NEW = @TRUE; $NEW[$i] = "C"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
				@NEW = @TRUE; $NEW[$i] = "G"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
				@NEW = @TRUE; $NEW[$i] = "A"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
				@NEW = @TRUE; $NEW[$i] = "N"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
			}
		}
	}
}


if (!defined $opt{'1'}) {
	if (-e "$opt{'F'}/$opt{'R'}/Undetermined_S0_L001_R1_001.fastq.gz") {
		$r1 = "$opt{'F'}/$opt{'R'}/Undetermined_S0_L001_R1_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L002_R1_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L003_R1_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L004_R1_001.fastq.gz";
		$r2 = "$opt{'F'}/$opt{'R'}/Undetermined_S0_L001_R2_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L002_R2_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L003_R2_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L004_R2_001.fastq.gz";
		$i1 = "$opt{'F'}/$opt{'R'}/Undetermined_S0_L001_I1_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L002_I1_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L003_I1_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L004_I1_001.fastq.gz";
		$i2 = "$opt{'F'}/$opt{'R'}/Undetermined_S0_L001_I2_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L002_I2_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L003_I2_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L004_I2_001.fastq.gz";
	} elsif (-e "$opt{'F'}/$opt{'R'}/Undetermined_S0_R1_001.fastq.gz") {
		$r1 = "$opt{'F'}/$opt{'R'}/Undetermined_S0_R1_001.fastq.gz";
		$r2 = "$opt{'F'}/$opt{'R'}/Undetermined_S0_R2_001.fastq.gz";
		$i1 = "$opt{'F'}/$opt{'R'}/Undetermined_S0_I1_001.fastq.gz";
		$i2 = "$opt{'F'}/$opt{'R'}/Undetermined_S0_I2_001.fastq.gz";
	}
} else {
	$r1 = "$opt{'1'}";
	$r2 = "$opt{'2'}";
	$i1 = "$opt{'i'}";
	$i2 = "$opt{'j'}";
}

system("mkdir $opt{'O'}/$opt{'o'}");

if (defined $opt{'A'}) {
	read_annot($opt{'A'});
	foreach $annot (keys %ANNOT_count) {
		$out1_handle = "$annot.1";
		open $out1_handle, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.$annot.1.fq.gz";
		$HANDLE1{$annot} = $out1_handle;
		$out2_handle = "$annot.2";
		open $out2_handle, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.$annot.2.fq.gz";
		$HANDLE2{$annot} = $out2_handle;
	}
	open O1, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.unassigned.1.fq.gz";
	open O2, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.unassigned.2.fq.gz";
	system("cp $opt{'A'} $opt{'O'}/$opt{'o'}/$opt{'o'}.split.annot");
} else {
	open R1OUT, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.1.fq.gz";
	open R2OUT, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.2.fq.gz";
}

open R1FAIL, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.fail.1.fq.gz";
open R2FAIL, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.fail.2.fq.gz";

$totalCT = 0; $failCT = 0;

open R1, "$zcat $r1 |";
open R2, "$zcat $r2 |";
open I1, "$zcat $i1 |";
open I2, "$zcat $i2 |";

	
while ($r1tag = <R1>) {

	$r2tag = <R2>; chomp $r1tag; chomp $r2tag;
	$r1seq = <R1>; $r2seq = <R2>; chomp $r1seq; chomp $r2seq;
	$null = <R1>; $null = <R2>;
	$r1qual = <R1>; $r2qual = <R2>; chomp $r1qual; chomp $r2qual;
	
	$i1tag = <I1>; chomp $i1tag; $i2tag = <I2>; chomp $i2tag;
	$i1seq = <I1>; chomp $i1seq; $i2seq = <I2>; chomp $i2seq;
	$null = <I1>; $null = <I1>;
	$null = <I2>; $null = <I2>;
	
	$ix1 = substr($i1seq,0,$POS_length{'1'});
	$ix2 = substr($i1seq,$POS_length{'1'},$POS_length{'2'});
	$ix3 = substr($i2seq,0,$POS_length{'3'});
	$ix4 = substr($i2seq,$POS_length{'3'},$POS_length{'4'});
	
	if (defined $POS_SEQ_seq{'1'}{$ix1} &&
		defined $POS_SEQ_seq{'2'}{$ix2} &&
		defined $POS_SEQ_seq{'3'}{$ix3} &&
		defined $POS_SEQ_seq{'4'}{$ix4}) {
		
		$barc = $POS_SEQ_seq{'1'}{$ix1}.$POS_SEQ_seq{'2'}{$ix2}.$POS_SEQ_seq{'3'}{$ix3}.$POS_SEQ_seq{'4'}{$ix4};
		
		$totalCT++;
		
		if (defined $opt{'r'}) {
			$r1out = "\@$barc:$totalCT.$opt{'r'}#0/1\n$r1seq\n\+\n$r1qual";
			$r2out = "\@$barc:$totalCT.$opt{'r'}#0/2\n$r2seq\n\+\n$r2qual";
		} else {
			$r1out = "\@$barc:$totalCT#0/1\n$r1seq\n\+\n$r1qual";
			$r2out = "\@$barc:$totalCT#0/2\n$r2seq\n\+\n$r2qual";
		}
		
		if (!defined $opt{'A'}) {
			print R1OUT "$r1out\n";
			print R2OUT "$r2out\n";
		} else {
			if (defined $CELLID_annot{$barc}) {
				$annot = $CELLID_annot{$barc};
				$ANNOT_count{$annot}++;
				$out1_handle = $HANDLE1{$annot}; $out2_handle = $HANDLE2{$annot};
				print $out1_handle "$r1out\n";
				print $out2_handle "$r2out\n";
			} else {
				$non_annot_count++;
				print O1 "$r1out\n";
				print O2 "$r2out\n";
			}
		}
		
	} else {
	
		$barc = $ix1.$ix2.$ix3.$ix4;
		
		if (defined $opt{'r'}) {
			print R1FAIL "\@$barc:F_$failCT.$opt{'r'}#0/1\n$r1seq\n\+\n$r1qual\n";
			print R2FAIL "\@$barc:F_$failCT.$opt{'r'}#0/1\n$r2seq\n\+\n$r2qual\n";
		} else {
			print R1FAIL "\@$barc:F_$failCT#0/1\n$r1seq\n\+\n$r1qual\n";
			print R2FAIL "\@$barc:F_$failCT#0/1\n$r2seq\n\+\n$r2qual\n";
		}
		
		$failCT++;
	}
}

close R1; close R2; close I1; close I2;

if (defined $opt{'A'}) {
	foreach $annot (keys %ANNOT_count) {
		$out1_handle = $HANDLE1{$annot}; $out2_handle = $HANDLE2{$annot};
		close $out1_handle; close $out2_handle;
		print STDERR "Annot: $annot, count = $ANNOT_count{$annot}\n";
	}
	close O1; close O2;
}

close R1FAIL; close R2FAIL;

exit;	
}

if ($command eq "fastq-split" || $command eq "split-fastq") {

getopts("A:O:X", \%opt);

$die2 = "
scitools fastq-split [options] read1.fq(.gz) read2.fq(.gz)
   or    split-fastq

Will split your fastq files by annotation. Reads that do not
match the annotastion will go to [out_prefix].unassigned.
Matching will go to [out_prefix].[annot]...

Options:
   -A   [STR]   Annotation file (comma separated for more than one)
                (If multiple, must be non-conflicting)
				(Required)
   -O   [STR]   Output prefix (def = annot file prefix)
   -X           Do not output unassigned reads

";

if (!defined $ARGV[0] || !defined $ARGV[1] || !defined $opt{'A'}) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $opt{'A'}; $opt{'O'} =~ s/\.annot$//};

$r1 = $ARGV[0]; $r2 = $ARGV[1];

read_annot($opt{'A'});
print STDERR "$annot_count total annotations found.\n";

# setup output files
foreach $annot (keys %ANNOT_count) {
	$out1_handle = "$annot.1";
	open $out1_handle, "| $gzip > $opt{'O'}.$annot.1.fq.gz";
	$HANDLE1{$annot} = $out1_handle;
	$out2_handle = "$annot.2";
	open $out2_handle, "| $gzip > $opt{'O'}.$annot.2.fq.gz";
	$HANDLE2{$annot} = $out2_handle;
}

if (!defined $opt{'X'}) {
	open O1, "| $gzip > $opt{'O'}.unassigned.1.fq.gz";
	open O2, "| $gzip > $opt{'O'}.unassigned.2.fq.gz";
}

if ($r1 =~ /\.gz$/) {open R1, "$zcat $r1 |"} else {open R1, "$r1"};
if ($r2 =~ /\.gz$/) {open R2, "$zcat $r2 |"} else {open R2, "$r2"};

while ($r1tag = <R1>) {
	$r2tag = <R2>; chomp $r1tag; chomp $r2tag;
	$r1seq = <R1>; chomp $r1seq; $null = <R1>; $r1qual = <R1>; chomp $r1qual;
	$r2seq = <R2>; chomp $r2seq; $null = <R2>; $r2qual = <R2>; chomp $r2qual;
	$cellID = $r1tag; $cellID =~ s/^\@//; $cellID =~ s/:.+$//;
	if (defined $CELLID_annot{$cellID}) {
		$annot = $CELLID_annot{$cellID};
		$ANNOT_count{$annot}++;
		$out1_handle = $HANDLE1{$annot}; $out2_handle = $HANDLE2{$annot};
		print $out1_handle "$r1tag\n$r1seq\n\+\n$r1qual\n";
		print $out2_handle "$r2tag\n$r2seq\n\+\n$r2qual\n";
	} else {
		$non_annot_count++;
		if (!defined $opt{'X'}) {
			print O1 "$r1tag\n$r1seq\n\+\n$r1qual\n";
			print O2 "$r2tag\n$r2seq\n\+\n$r2qual\n";
		}
	}
}

foreach $annot (keys %ANNOT_count) {
	$out1_handle = $HANDLE1{$annot}; $out2_handle = $HANDLE2{$annot};
	close $out1_handle; close $out2_handle;
	print STDERR "Annot: $annot, count = $ANNOT_count{$annot}\n";
}
if (!defined $opt{'X'}) {close O1; close O2};
print STDERR "Not in annotaiton file: $non_annot_count\n";

exit;

}

if ($command eq "fastq-merge" || $command eq "merge-fastq") {
getopts("O:r:", \%opt);

$die2 = "
scitools fastq-merge [options] [output prefix] [read1 fq comma sep] [read2 fq comma sep]
   or    merge-fastq

Will merge fastq files that have the same set of cellIDs present
Adds in a run identifier (after read number) based on the order of fastqs.

Options:
   -r   [STR]   Comma separated run IDs (def = 1,2,3... etc)

";

if (defined $opt{'O'}) {unshift @ARGV, $opt{'O'}};
if (!defined $ARGV[2]) {die $die2};
$ARGV[0] =~ s/\.gz$//; $ARGV[0] =~ s/\.fq$//; $ARGV[0] =~ s/\.$//; 
@R1_FQ = split(/,/, $ARGV[1]);
@R2_FQ = split(/,/, $ARGV[2]);
if (length(@R1_FQ) != length(@R2_FQ)) {die "\n\nERROR: must specify the same number of read 1 and read 2 fastq files!\n"};
if (defined $opt{'r'}) {
	@RUNIDS = split(/,/, $opt{'r'});
	if (length(@RUNIDS) != length(@R1_FQ)) {
		die "\n\nERROR: When specifying run IDs you must specify the same number as input fastq files\n";
	}
}
open O1, "| $gzip > $ARGV[0].1.fq.gz";
open O2, "| $gzip > $ARGV[0].2.fq.gz";
for ($id = 0; $id <= @R1_FQ; $id++) {
	if ($R1_FQ[$id] =~ /\.gz$/) {open R1, "$zcat $R1_FQ[$id] |"} else {open R1, "$R1_FQ[$id]"};
	if ($R2_FQ[$id] =~ /\.gz$/) {open R2, "$zcat $R2_FQ[$id] |"} else {open R2, "$R2_FQ[$id]"};
	if (defined $opt{'r'}) {$runID = $RUNIDS[$id]} else {$runID = $id};
	while ($r1tag = <R1>) {
		chomp $r1tag; $r1seq = <R1>; $null = <R1>; $r1qual = <R1>;
		$r2tag = <R2>; chomp $r2tag; $r2seq = <R2>; $null = <R2>; $r2qual = <R2>;
		$r1tag =~ s/\#.+$//; $r2tag =~ s/\#.+$//;
		print O1 "$r1tag.$runID#0/1\n$r1seq\+\n$r1qual";
		print O2 "$r2tag.$runID#0/2\n$r2seq\+\n$r2qual";
	} close R1; close R2;
} close O1; close O2;

exit;
}

if ($command eq "fastq-align" || $command eq "align-fastq" || $command eq "align") {

# Defaults
$threads = 1;
$memory = "2G";

getopts("t:b:s:A:L:rO:m:", \%opt);

$die2 = "
scitools fastq-align [options] [bwa reference] [read1.fq] [read2.fq] [output_prefix]
   or    align-fastq
   or    algn

Produces a sorted bam file.

Options:
   -t   [INT]   Threads for alignment (def = $threads)
   -O   [STR]   Output prefix (if not set in the other inputs)
   -b   [STR]   Bwa call (def = $bwa)
   -s   [STR]   Samtools call (def = $samtools)
   -m   [MEM]   Samtools sort mex memory per thread, K/M/G (def = $memory)

Bwa reference shortcuts:
   hg19   $hg19_ref
   hg38   $hg38_ref
   mm10   $mm10_ref

";

if (defined $opt{'O'}) {$ARGV[3] = $opt{'O'}};
$ARGV[3] =~ s/\.bam$//; $ARGV[3] =~ s/\.$//;
if (!defined $ARGV[3]) {die $die2};
if (defined $opt{'b'}) {$bwa = $opt{'b'}};
if (defined $opt{'s'}) {$samtools = $opt{'s'}};
if (!defined $opt{'t'}) {$opt{'t'} = $threads};
if (!defined $opt{'m'}) {$opt{'m'} = $memory};

if ($ARGV[0] eq "hg19") {$ARGV[0] = $hg19_ref};
if ($ARGV[0] eq "hg38") {$ARGV[0] = $hg38_ref};
if ($ARGV[0] eq "mm10") {$ARGV[0] = $mm10_ref};

system("$bwa mem -t $opt{'t'} $ARGV[0] $ARGV[1] $ARGV[2] 2>> $ARGV[3].align.log | $samtools view -bSu - 2>> $ARGV[3].align.log | $samtools sort -m $opt{'m'} -T $ARGV[3].TMP - > $ARGV[3].bam 2>> $ARGV[3].align.log");

exit;
}


########## ANNOTATION FUNCTIONS ##########
if ($command eq "annot-make" || $command eq "make-annot") {

getopts("O:I:P:ph", \%opt);

# DEFAULTS
@LETTERS = ("0", "A", "B", "C", "D", "E", "F", "G", "H");
%LETTER_NUM = ("A"=>"1", "B"=>"2", "C"=>"3", "D"=>"4", "E"=>"5", "F"=>"6", "G"=>"7", "H"=>"8");
$ALL = "1-8";

$die2 = "
scitools annot-make [options] [annotation_description_1] [annotaiton_description_2] ...
   or    make-annot

Options:
   -O   [STR]   Output annotation file (default = STDOUT)
   -I   [STR]   Index file
         (default = $SCI_index_file)
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
         (default = $SCI_index_file)
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
open IN, $SCI_index_file;
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
			$class =~ s/^#//;
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
   -f   [STR]   Fai file for chr lengths (hg19, hg38, and mm10 are shortcuts)
                If toggled will ensure n peaks extend beyond
   -X           Retain intermediate files (def = remove)

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.bam$//};
if (defined $opt{'l'}) {$min_feature_size = $opt{'l'}};
if (defined $opt{'s'}) {$samtools = $opt{'s'}};

if (defined $opt{'f'}) {
	if ($opt{'f'} eq "hg19") {
		open FAI, "$hg19_ref.fai";
	} elsif ($opt{'f'} eq "mm10") {
		open FAI, "$mm10_ref.fai";
	} elsif ($opt{'f'} eq "hg38") {
		open FAI, "$hg38_ref.fai";
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

getopts("s:b:O:BC:", \%opt);

$die2 = "
scitools atac-count [options] [bam file] [peaks bed file]
   or    count(s)

Creates a counts or binary matrix file and a values file
with the fraction of reads on target for each cell.

Options:
   -O   [STR]   Output prefix (default is peaks prefix)
                (adds .counts.matrix)
   -b   [STR]   Bedtools call (def = $bedtools)
   -s   [STR]   Samtools call (def = $samtools)
   -B           Make it a binary matrix instead of counts matrix
                (file name will end in .binary.matrix)
   -C   [STR]   Complexity file (speeds up on-target calc)

";

if (!defined $ARGV[1]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[1]; $opt{'O'} =~ s/\.bed//};

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
if (defined $opt{'M'}) {$motif_set = $opt{'M'}};
if (!defined $MOTIF_SETS{$motif_set}) {die "\n\nERROR: The motif set provided ($motif_set) does not exist in chromVARmotifs\n"};
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
motif_ix <- matchMotifs(motifs, counts, genome = $genome)

# calculate & print deviations
dev <- computeDeviations(object = counts, annotations = motif_ix)
write.table(as.matrix(deviations(dev)),file = \"$opt{'O'}.chromVAR/deviations.matrix\", col.names = TRUE, row.names = TRUE, sep = \"\\t\", quote = FALSE)
write.table(as.matrix(deviationScores(dev)),file = \"$opt{'O'}.chromVAR/deviation_scores.matrix\", col.names = TRUE, row.names = TRUE, sep = \"\\t\", quote = FALSE)

# calculate & print variabilities
var <- computeVariability(dev)
write.table(as.matrix(var),file = \"$opt{'O'}.chromVAR/variabiltiy.txt\", col.names = TRUE, row.names = TRUE, sep = \"\\t\", quote = FALSE)

# generate tSNE on the deviations
tsne <- deviationsTsne(dev, threshold = 1.5, perplexity = 30, shiny = FALSE)
write.table(as.matrix(tsne),file = \"$opt{'O'}.chromVAR/tsne.txt\", col.names = TRUE, row.names = TRUE, sep = \"\\t\", quote = FALSE)

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

if ($command eq "aggregate-cells" || $command eq "aggregate") { # TODO: LABDA BASED! Also - reassign orphans.

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
}

if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.(dims|lambda|values)$//;

if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
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
if (defined $opt{'C'}) {read_color_file($opt{'C'})};
if (defined $opt{'c'}) {read_color_string($opt{'c'})};

if (defined $opt{'x'}) {$xdim = $opt{'x'}};
if (defined $opt{'y'}) {$ydim = $opt{'y'}};
$xpos = $xdim-1; $ypos = $ydim-1;

read_ranges($opt{'D'});

if (!defined $opt{'A'}) {
	foreach $cellID (keys %CELLID_DIMS) {
		$CELLID_annot{$cellID} = "Agg";
	}
	$ANNOT_include{"Agg"} = 1;
}

open CNT, ">$opt{'O'}.centroids.dims";
open OUT, ">$opt{'O'}.aggregate.annot";

# loop through annotaitons individually
foreach $annot (keys %ANNOT_include) {
if ($ARGV[0] =~ /\.(lambda|values)$/) {
	# lambda file
	read_values($ARGV[0]);
	
} else {
	# Dimensions file
	read_dims($ARGV[0]);
	
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
	# labmda data plot
	
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
   ylab(\"log10 Unique Reads\") +
   theme(legend.background=element_blank(),legend.title=element_blank())
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
	scale_x_continuous(limits=c(0,6)) +
	theme(legend.background=element_blank(),legend.title=element_blank())
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
	scale_colour_gradientn(colours=gradient_funct(21),limits=c($minV,$maxV)) +"
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

Options:
   -O   [STR]   Output prefix (default is input prefix)
   -o   [STR]   Original dims file (def = auto find [prefix].orig.dims)
   -p   [STR]   Pcurve dims file (def = auto find [prefix].proj.dims)
   -l   [STR]   Lambda file (def = auto find [prefix].lambda)
                (if [prefix].centered.lambda is found, it will use that)
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
			if (defined $CELLID_FEATURE_value{$cellID}{$feature} && ($annot !~ /Exclude/i)) {
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
theme_bw() + xlab(\"Lambda\") + ylab(\"Feature value\") +
theme(legend.background=element_blank(),
legend.title=element_blank(),
panel.background=element_blank(),
plot.background=element_blank())

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
	theme_bw() + xlab(\"Correlation\") + ylab(\"Counts\") +
	theme(legend.background=element_blank(),
	legend.title=element_blank(),
	panel.background=element_blank(),
	plot.background=element_blank())

	ggsave(plot=HISTPlot,filename=\"$opt{'O'}.$matrix_out\pcurve/Hist.corr.Pearson.lambda.png\",width=5,height=4,dpi=900)
	ggsave(plot=HISTPlot,filename=\"$opt{'O'}.$matrix_out\pcurve/Hist.corr.Pearson.lambda.plot.pdf\",width=5,height=4)
	
		IN<-read.table(\"$opt{'O'}.$matrix_out\pcurve/corr.spearman.txt\",header=FALSE)
	# Make the corr distribution plot
	
	HISTPlot<-ggplot(IN,aes(x=V2)) + geom_histogram()+ 
	theme_bw() + xlab(\"Correlation\") + ylab(\"Counts\") +
	theme(legend.background=element_blank(),
	legend.title=element_blank(),
	panel.background=element_blank(),
	plot.background=element_blank())

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
	
	if ($annot !~ /Exclude/i) {
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
		  plot.background=element_blank())\n";
} else {
	print R "
	theme_bw()\n";
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
	print R "
	theme_bw() +
	theme(panel.border=element_blank(),
		  panel.grid=element_blank(),
		  axis.line=element_blank(),
		  axis.ticks=element_blank(),
		  legend.background=element_blank(),
		  panel.background=element_blank(),
		  axis.text=element_blank(),
		  axis.title.x=element_blank(),
		  axis.title.y=element_blank(),
		  plot.background=element_blank())\n";
} else {
	print R "
	theme_bw()\n";
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
	theme_bw() + xlab(\"Lambda\") +	ylab(\"Density\") + theme(legend.background=element_blank(),legend.title=element_blank())";
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

if ($command eq "plot-reads") {

# Defaults
$height = 6;
$width = 6;
$pt_size = 0.5;
$flanking_size = 100000;
$gene_scale_factor = 1;
$gene_text_size = 1.5;

getopts("O:A:a:C:c:R:Xs:Dh:w:rp:B:G:S:f:t:", \%opt);

$die2 = "
scitools plot-reads [options] [rmdup sci bam file] [chrN:start-end] [region 2] ...

Options:
   -O   [STR]   Output prefix (default is input prefix)
   -D           Output to a directory (will add to it if it exists)
   -B   [BED]   Bed file of peaks (optional)
   -G   [BED]   Gene info (refGene.txt formats)
                Defaults: hg38, hg19, mm10 (see scitools -h for more details)
   -S   [INT]   Flanking size if gene names are specified
                (reuires -G, def = $flanking_size)
   -A   [STR]   Annotation file (to color code reads)
   -a   [STR]   Comma separated list of annoations to include in plot
                  (requires -A to be specified)
   -C   [STR]   Color coding file (annot (tab) #hexColor)
   -c   [STR]   Color coding string
                  Annot=#hexColor,Annot2=#hexColor
   -h   [IN]    Height (inches, def = $height)
   -w   [IN]    Width (inches, def = $width)
   -p   [FLT]   Point size (def = $pt_size)
   -f   [FLT]   Gene plot spacing factor (def = $gene_scale_factor)
                  (for -G, larger values = more vertical spread)
   -t   [FLT]   Gene name text size (def = $gene_text_size)
   -r           Do not randomize cell order (def = randomize)
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
	if ($opt{'G'} eq "hg38") {
		$opt{'G'} = $hg38_ref; $opt{'G'} =~ s/\.fa$/\.refGene.txt/;
	} elsif ($opt{'G'} eq "hg19") {
		$opt{'G'} = $hg19_ref; $opt{'G'} =~ s/\.fa$/\.refGene.txt/;
	} elsif ($opt{'G'} eq "mm10") {
		$opt{'G'} = $mm10_ref; $opt{'G'} =~ s/\.fa$/\.refGene.txt/;
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

$newID = 0;
foreach $annot (keys %ANNOT_ids) {
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

open IN, "$opt{'O'}$region_name.data";
open OUT, ">$opt{'O'}$region_name.data.full";
while ($l = <IN>) {
	chomp $l;
	($annot,$id,$posn) = split(/\t/, $l);
	$newID = $ANNOT_CELLID_newID{$annot}{$id};
	print OUT "READ\t$annot\t$newID\t$posn\n";
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
print R "
		axis.ticks.y=element_blank(),
		axis.text.y=element_blank(),
		legend.title=element_blank(),
		panel.background=element_blank(),
		plot.background=element_blank())
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


#####################################################
#################### SUBROUTINES ####################
#####################################################

sub read_annot {
	%CELLID_annot = (); %ANNOT_count = (); $annot_count = 0;
	@ANNOT_FILES = split(/,/, $_[0]);
	foreach $annot_file (@ANNOT_FILES) {
		open ANNOT, "$annot_file";
		while ($annot_line = <ANNOT>) {
			chomp $annot_line;
			($annot_cellID,$annot) = split(/\t/, $annot_line);
			$CELLID_annot{$annot_cellID} = $annot;
			if (!defined $ANNOT_count{$annot}) {
				$annot_count++;
				$ANNOT_count{$annot}=0;
			}
		} close ANNOT;
	}
}

sub read_complexity {
	%CELLID_uniq_reads = ();
	%CELLID_raw_reads = ();
	%CELLID_complexity = ();
	%CELLID_complexity_rank = ();
	@COMPLEXITY_FILES = split(/,/, $_[0]);
	foreach $complexity_file (@COMPLEXITY_FILES) {
		open COMPL, "$complexity_file";
		while ($comp_line = <COMPL>) {
			chomp $comp_line;
			($num,$cellID,$raw,$uniq,$pct) = split(/\t/, $comp_line);
			$CELLID_complexity_rank{$cellID} = $num;
			$CELLID_uniq_reads{$cellID} = $uniq;
			$CELLID_raw_reads{$cellID} = $raw;
			$CELLID_complexity{$cellID} = $pct;
		} close COMPL;
	}
}

sub read_matrix {
	%CELLID_FEATURE_value = (); @MATRIX_COLNAMES = (); @MATRIX_ROWNAMES = ();
	%MATRIX_CellID_nonZero = (); %MATRIX_feature_nonZero = ();
	%MATRIX_CellID_signal = (); %MATRIX_feature_signal = ();
	$matrix_colNum = 0; $matrix_rowNum = 0;
	open MAT, "$_[0]";
	$column_line = <MAT>; chomp $column_line; @MATRIX_COLNAMES = split(/\t/, $column_line);
	foreach $column (@MATRIX_COLNAMES) {$matrix_colNum++};
	while ($matrix_line = <MAT>) {
		chomp $matrix_line;
		@MATRIX_ROW = split(/\t/, $matrix_line);
		$featureName = shift(@MATRIX_ROW);
		push @MATRIX_ROWNAMES, $featureName;
		for ($colNum = 0; $colNum < @MATRIX_ROW; $colNum++) {
			$CELLID_FEATURE_value{$MATRIX_COLNAMES[$colNum]}{$featureName} = $MATRIX_ROW[$colNum];
			if ($MATRIX_ROW[$colNum]>0) {
				$MATRIX_CellID_nonZero{$MATRIX_COLNAMES[$colNum]}++;
				$MATRIX_feature_nonZero{$featureName}++;
				$MATRIX_CellID_signal{$MATRIX_COLNAMES[$colNum]}+=$MATRIX_ROW[$colNum];
				$MATRIX_feature_signal{$featureName}+=$MATRIX_ROW[$colNum];
			}
		}
		$matrix_rowNum++;
	} close MAT;
}

sub read_color_string {
	%ANNOT_color = ();
	@COL_STRING = split(/,/, $_[0]);
	$color_mapping = "\"Cell\" = \"lightsteelblue4\",";
	foreach $color_assignment (@COL_STRING) {
		($annot,$color) = split(/=/, $color_assignment);
		$ANNOT_color{$annot} = $color;
		$color_mapping .= "\"$annot\" = \"$color\",";
	} $color_mapping =~ s/,$//;
}

sub read_color_file {
	%ANNOT_color = ();
	open COL_FILE, "$_[0]";
	$color_mapping = "\"Cell\" = \"lightsteelblue4\",";
	while ($color_assignment = <COL_FILE>) {
		chomp $color_assignment;
		($annot,$color) = split(/\t/, $color_assignment);
		$ANNOT_color{$annot} = $color;
		$color_mapping .= "\"$annot\" = \"$color\",";
	} $color_mapping =~ s/,$//; close COL_FILE;
}

sub read_dims {
	%CELLID_DIMS = ();
	@DIM_FILES = split(/,/, $_[0]);
	foreach $dim_file (@DIM_FILES) {
		open DIMS, "$dim_file";
		while ($dim_line = <DIMS>) {
			chomp $dim_line;
			@DIM_SET = split(/\t/, $dim_line);
			$cellID = $DIM_SET[0];
			@{$CELLID_DIMS{$cellID}} = @DIM_SET;
			$Ndims = @DIM_SET;
		} close DIMS;
	}
}

sub read_pcurve_dims {
	%CELLID_PCURVE_DIMS = ();
	@PC_DIM_FILES = split(/,/, $_[0]);
	foreach $pc_dim_file (@PC_DIM_FILES) {
		open DIMS, "$pc_dim_file";
		while ($dim_line = <DIMS>) {
			chomp $dim_line;
			@DIM_SET = split(/\t/, $dim_line);
			$cellID = $DIM_SET[0];
			@{$CELLID_PCURVE_DIMS{$cellID}} = @DIM_SET;
		} close DIMS;
	}
}

sub read_values {
	%CELLID_value = ();
	@VALUES = ();
	$value_min = "na"; $value_max = "na";
	$value_mean = 0; $value_sum = 0; $value_median = 0;
	@VAL_FILES = split(/,/, $_[0]);
	foreach $val_file (@VAL_FILES) {
		open VALS, "$val_file";
		while ($val_line = <VALS>) {
			chomp $val_line;
			($cellID,$value) = split(/\t/, $val_line);
			$CELLID_value{$cellID} = $value;
			push @VALUES, $value;
			$value_sum += $value;
			if ($value > $value_max || $value_max eq "na") {$value_max = $value};
			if ($value < $value_min || $value_min eq "na") {$value_min = $value};
		} close VALS;
	}
	$value_range = $value_max - $value_min;
	$value_mean = $value_sum/@VALUES;
	@SORTED_VALUES = sort {$a<=>$b} @VALUES;
	$value_median = $SORTED_VALUES[int(@VALUES/2)];
}

sub read_ranges {
	@RANGE_VALUES = ();
	$range_R_set = "c(";
	@RANGE_SETS = split(/,/, $_[0]);
	for ($range_pos = 0; $range_pos < @RANGE_SETS; $range_pos++) {
		$range_set = $RANGE_SETS[$range_pos];
		if ($range_set =~ /-/) {
			($range_start,$range_end) = split(/-/, $range_set);
			for ($range_value = $range_start; $range_value <= $range_end; $range_value++) {
				push @RANGE_VALUES, $range_value;
				$range_R_set .= "$range_value,";
			}
		} else {
			push @RANGE_VALUES, $range_set;
			$range_R_set .= "$range_set,";
		}
	}
	$range_R_set =~ s/,$//;
	$range_R_set .= ")";
}

sub read_indexes {
	%INDEX_POS_SEQ_id = ();
	%INDEX_POS_SEQ_well = ();
	open INDEX, "$_[0]";
	while ($index_line = <INDEX>) {
		chomp $index_line;
		($index_id,$index_pos,$index_seq) = split(/\t/, $index_line);
		($id_tier,$id_set,$id_side,$id_well) = split(/_/, $index_id);
		$INDEX_POS_SEQ_id{$index_pos}{$index_seq} = $id_set;
		$INDEX_POS_SEQ_well{$index_pos}{$index_seq} = $id_well;
	} close INDEX;
}

sub read_refgene {
	%GENENAME_coords = ();
	%GENEID_coords = ();
	%GENECOORDS_refGene = ();
	open REFGENE, "$_[0]";
	while ($refgene_line = <REFGENE>) {
		chomp $refgene_line;
		@REFGENE = split(/\t/, $refgene_line);
		$gene_coords = "$REFGENE[2]:$REFGENE[4]-$REFGENE[5]";
		$GENEID_coords{$REFGENE[1]} = $gene_coords;
		$GENENAME_coords{$REFGENE[12]} = $gene_coords;
		$GENECOORDS_refGene{$gene_coords} = $refgene_line;
	} close REFGENE;
}

sub get_gradient {
	if ($_[0] =~ /,/) {
		@GRADIENT_COLORS = split(/,/, $_[0]);
		$gradient_specification = "gradient_funct<-colorRampPalette(c(";
		for ($grad_pos = 0; $grad_pos < @GRADIENT_COLORS; $grad_pos++) {
			$gradient_specification .= "\"$GRADIENT_COLORS[$grad_pos]\",";
		} $gradient_specification =~ s/,$//;
		$gradient_specification .= "))";
	} elsif (defined $COLOR_GRADIENT{$_[0]}) {
		$gradient_specification = $COLOR_GRADIENT{$_[0]};
	} else {
		die "
ERROR: Color gradients must either be one of the pre-set scitools color gradients,
       or defined as a comma separated list of two or more colors. For more information
       run 'scitools gradient'.
";
	}
	return $gradient_specification;
}

sub load_gradient_defaults {
# COLOR GRADIENT DEFAULTS:
%COLOR_GRADIENT = ();

# diverging
$COLOR_GRADIENT{'PuOr'} = "gradient_funct<-colorRampPalette(c(\"#542788\",\"#8073ac\",\"#b2abd2\",\"#d8daeb\",\"#f7f7f7\",\"#fee0b6\",\"#fdb863\",\"#e08214\",\"#b35806\"))";
$COLOR_GRADIENT{'BuRd'} = "gradient_funct<-colorRampPalette(c(\"#2166ac\",\"#4393c3\",\"#92c5de\",\"#d1e5f0\",\"#f7f7f7\",\"#fddbc7\",\"#f4a582\",\"#d6604d\",\"#b2182b\"))";
$COLOR_GRADIENT{'PuRd'} = "gradient_funct<-colorRampPalette(c(\"#542788\",\"#8073ac\",\"#b2abd2\",\"#d8daeb\",\"#f7f7f7\",\"#fddbc7\",\"#f4a582\",\"#d6604d\",\"#b2182b\"))";
$COLOR_GRADIENT{'RdYlGn'} = "gradient_funct<-colorRampPalette(c(\"#d73027\",\"#f46d43\",\"#fdae61\",\"#fee08b\",\"#ffffbf\",\"#d9ef8b\",\"#a6d96a\",\"#66bd63\",\"#1a9850\"))";
$COLOR_GRADIENT{'BrBG'} = "gradient_funct<-colorRampPalette(c(\"#8c510a\",\"#bf812d\",\"#dfc27d\",\"#f6e8c3\",\"#f5f5f5\",\"#c7eae5\",\"#80cdc1\",\"#35978f\",\"#01665e\"))";
$COLOR_GRADIENT{'PiYG'} = "gradient_funct<-colorRampPalette(c(\"#c51b7d\",\"#de77ae\",\"#f1b6da\",\"#fde0ef\",\"#f7f7f7\",\"#e6f5d0\",\"#b8e186\",\"#7fbc41\",\"#4d9221\"))";
$COLOR_GRADIENT{'PuRdGn'} = "gradient_funct<-colorRampPalette(c(\"#762a83\",\"#9970ab\",\"#c2a5cf\",\"#e7d4e8\",\"#f7f7f7\",\"#d9f0d3\",\"#a6dba0\",\"#5aae61\",\"#1b7837\"))";
$COLOR_GRADIENT{'BuYlRd'} = "gradient_funct<-colorRampPalette(c(\"#2c7bb6\",\"#abd9e9\",\"#ffffbf\",\"#fdae61\",\"#d7191c\"))";

# sequential multi-hue
$COLOR_GRADIENT{'YlOrRd'} = "gradient_funct<-colorRampPalette(c(\"#ffffcc\",\"#ffeda0\",\"#fed976\",\"#feb24c\",\"#fd8d3c\",\"#fc4e2a\",\"#e31a1c\",\"#bd0026\",\"#800026\"))";
$COLOR_GRADIENT{'WtYlOrRd'} = "gradient_funct<-colorRampPalette(c(\"white\",\"#ffffcc\",\"#ffeda0\",\"#fed976\",\"#feb24c\",\"#fd8d3c\",\"#fc4e2a\",\"#e31a1c\",\"#bd0026\",\"#800026\"))";
$COLOR_GRADIENT{'RdPu'} = "gradient_funct<-colorRampPalette(c(\"#fff7f3\",\"#fde0dd\",\"#fcc5c0\",\"#fa9fb5\",\"#f768a1\",\"#dd3497\",\"#ae017e\",\"#7a0177\",\"#49006a\"))";
$COLOR_GRADIENT{'YlGnBu'} = "gradient_funct<-colorRampPalette(c(\"#ffffd9\",\"#edf8b1\",\"#c7e9b4\",\"#7fcdbb\",\"#41b6c4\",\"#1d91c0\",\"#225ea8\",\"#253494\",\"#081d58\"))";
$COLOR_GRADIENT{'BuGnYl'} = "gradient_funct<-colorRampPalette(c(\"#081d58\",\"#253494\",\"#225ea8\",\"#1d91c0\",\"#41b6c4\",\"#7fcdbb\",\"#c7e9b4\",\"#edf8b1\",\"#ffffd9\"))";
$COLOR_GRADIENT{'YlGn'} = "gradient_funct<-colorRampPalette(c(\"#ffffe5\",\"#f7fcb9\",\"#d9f0a3\",\"#addd8e\",\"#78c679\",\"#41ab5d\",\"#238443\",\"#006837\",\"#004529\"))";
$COLOR_GRADIENT{'Spect'} = "gradient_funct<-colorRampPalette(c(\"#3288bd\",\"#66c2a5\",\"#abdda4\",\"#e6f598\",\"#ffffbf\",\"#fee08b\",\"#fdae61\",\"#f46d43\",\"#d53e4f\"))";
$COLOR_GRADIENT{'WtSpect'} = "gradient_funct<-colorRampPalette(c(\"white\",\"#3288bd\",\"#66c2a5\",\"#abdda4\",\"#e6f598\",\"#ffffbf\",\"#fee08b\",\"#fdae61\",\"#f46d43\",\"#d53e4f\"))";


# sequential single-hue
$COLOR_GRADIENT{'WtPu'} = "gradient_funct<-colorRampPalette(c(\"white\",\"#fcfbfd\",\"#efedf5\",\"#dadaeb\",\"#bcbddc\",\"#9e9ac8\",\"#807dba\",\"#6a51a3\",\"#54278f\",\"#3f007d\"))";
$COLOR_GRADIENT{'WtRd'} = "gradient_funct<-colorRampPalette(c(\"white\",\"#fff5f0\",\"#fee0d2\",\"#fcbba1\",\"#fc9272\",\"#fb6a4a\",\"#ef3b2c\",\"#cb181d\",\"#a50f15\",\"#67000d\"))";

# non colorbrewer
$COLOR_GRADIENT{'BuGo'} = "gradient_funct<-colorRampPalette(c(\"lightblue\",\"blue\",\"red\",\"orange\",\"gold\"))";
$COLOR_GRADIENT{'BuG90Rd'} = "gradient_funct<-colorRampPalette(c(\"blue\",\"gray90\",\"red\"))";
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
The default executables can be altered at the start of the scitools
code in the # GLOBAL DEFAULTS section.

Executables:
   gzip         For gzipped fastq files. Default: $gzip & $zcat
                (these are hardoded at the beginning of the scitools
                 code and not available as options)
   samtools     Bam-related commands. Default: $samtools
   bedtools     Bed-related commands. Default: $bedtools
   Rscript      Numerous R operations. Default: $Rscript
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
      $fastq_input_directory
   Output fastq directory (for processed SCI fastq files)
      $SCI_fastq_directory
   SCI index file (should comtain all barcodes in the proper format)
      $SCI_index_file
   hg19 ([reference.fa], [reference.fa].fai, [reference.fa].[bwa_index], and [reference].refGene.txt)
      $hg19_ref
   hg38 ([reference.fa], [reference.fa].fai, [reference.fa].[bwa_index], and [reference].refGene.txt)
      $hg38_ref
   mm10 ([reference.fa], [reference.fa].fai, [reference.fa].[bwa_index], and [reference].refGene.txt)
      $mm10_ref

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

   3) Align fastq files (if reference defaults are set up use hg38, hg19, or mm10)
      scitools align -t [n threads] [reference_prefix or shortcut] [reads.sampleA.1.fq.gz] [reads.sampleA.2.fq.gz] [sampleA_prefix]

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



#################### END SCITOOLS ####################