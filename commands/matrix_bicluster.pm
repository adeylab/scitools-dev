package commands::matrix_bicluster;

use commands::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_bicluster");

sub matrix_bicluster {

# Defaults
%COLS = ();
$colfunc = "BuRd";
($width,$height) = (12,12);
$imageType = "png";
$res = 600;
$gradient_def = "BuRd";

@ARGV = @_;
getopts("", \%opt);

$die2 = "
scitools matrix-bicluster [options] [input matrix]
   or    bicluster-matrix
   or    bicluster

Produces a sorted bam file with read name and RG barcodes.

Options:
   -O   [STR]   Output (default is matrix file prefix)
   -A   [STR]   Annotation file (for column coloring)
   -C   [STR]   Color coding file (annot (tab) #hexColor)
   -c   [STR]   Color coding string
                  Annot=#hexColor,Annot2=#hexColor
   -N   [STR]   Row annotation file
   -n   [STR]   Color coding file for rows
   -t   [STR]   Color coding string for rows
   -G   [GRD]   Color gradient (def = $gradient_def)
                  For all available gradients, run 'scitools gradient'
   -R   [STR]   Rscript call (def = $Rscript)
   -V           Do not cluster columns
   -v           Do not cluster rows
   -X           Do not delete intermediate files (def = delete)
   
";

if (!defined $ARGV[0]) {die $die2};



#####################################################################


for ($i = 0; $i < @ARGV; $i++) {
	($type,$field) = split(/=/, $ARGV[$i]);
	if ($type =~ /ROWC/i) {
		$COLS{'ROW'} = $field;
	} elsif ($type =~ /COLC/i) {
		$COLS{'COL'} = $field;
	} elsif ($type =~ /GRA/i) {
		$colfunc = $field;
	} elsif ($type =~ /SIZE/i) {
		($width,$height) = split(/,/, $field);
	} elsif ($type =~ /TYPE/i) {
		$imageType = $field;
	} elsif ($type =~ /ROWV/i) {
		$ROWV = $field;
	} elsif ($type =~ /COLV/i) {
		$COLV = $field;
	} elsif ($type =~ /RES/i) {
		$res = $field;
	} elsif ($type =~ /MAT/i) {
		$matrix = $field;
	} elsif ($type =~ /OUT/i) {
		$out = $field;
	} else {
		die "$die\n";
	}
}

if (!defined $matrix) {die "Must provide matrix as MAT=[file] argument!\n$die"};

if (!defined $out) {
	$out = $matrix;
	$out =~ s/\.matrix//;
}

open R, ">$out.heatmap2.r";

print R "
library(gplots)
library(grid)
source(\"/home/users/adey/src/load_color_functions.r\")
colfunc <- $colfunc
IN <- read.table(\"$matrix\")\n";

if (defined $COLS{'ROW'}) {print R "ROW <- read.table(\"$COLS{'ROW'}\")\n"}
if (defined $COLS{'COL'}) {print R "COL <- read.table(\"$COLS{'COL'}\")\n"}

if ($imageType =~ /pdf/i) {
	print R "$imageType(\"$out.heatmap2.$imageType\",width=$width,height=$height)\n";
} else {
	print R "$imageType(\"$out.heatmap2.$imageType\",width=$width,height=$height,units=\"in\",res=$res)\n";
}

print R "HM2 <- heatmap.2(as.matrix(IN),
	trace=\"none\",
	col=colfunc(99),
	tracecol=\"black\"";

if (defined $COLS{'ROW'} && defined $COLS{'COL'}) {
	print R ",
	ColSideColors=as.vector(COL\$V1),
	RowSideColors=as.vector(ROW\$V1)";
} elsif (defined $COLS{'ROW'}) {
	print R ",
	RowSideColors=as.vector(ROW\$V1)";
} elsif (defined $COLS{'COL'}) {
	print R ",
	ColSideColors=as.vector(COL\$V1)";
}

if ($ROWV =~ /f/i) {
	print R ",
	Rowv=FALSE";
}
if ($COLV =~ /f/i) {
	print R ",
	Colv=FALSE";
}

print R ")\ndev.off()\n";

close R;

system("Rscript $out.heatmap2.r");