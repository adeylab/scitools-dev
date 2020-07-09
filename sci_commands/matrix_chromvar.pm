package sci_commands::matrix_chromvar;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_chromvar");

sub matrix_chromvar {

    @ARGV = @_;

# Defaults
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
scitools matrix-chromvar [options] [dense matrix] [peaks bed file]

This wrapper will execute chromVAR and output
deviation and deviation z-score matrix files, along with a variaiton
file for each motif or peak set assessed.

Options:
   -O   [STR]   Output prefix (default is bam prefix)
                (creates prefix.chromVAR directory)
   -g   [STR]   Genome (hg38, hg19, or mm10; def = hg38)
   -M   [STR]   Motif set (def = $motif_set)
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
    if (!defined $GENOMES{$opt{'g'}}) {die "\n\nERROR: Genome $opt{'g'} is not a proper genome reference! Exiting!\n"} else {$genome = $GENOMES{$opt{'g'}}};
    if (!defined $opt{'M'}) {$opt{'M'} = "human_pwms_v2"};
    if (!defined $MOTIF_SETS{$opt{'M'}}) {die "\n\nERROR: The motif set provided ($opt{'M'}) does not exist in chromVARmotifs\n"};

    if (-e "$opt{'O'}.chromVAR") {
	if (defined $opt{'X'}) {
	    print STDERR "\n\nWARNING: $opt{'O'}.chromVAR exists, will rewrite contents!\n";
	} else {
	    die "\n\nERROR: $opt{'O'}.chromVAR exists! Exiting!\n";
	}
    }

    system("mkdir $opt{'O'}.chromVAR");

    $matrixfile = $ARGV[0];
    $peakfile = $ARGV[1];

    open R, ">$opt{'O'}.chromVAR/chromVAR.r";
    $ts = localtime(time);
print R "# chromVAR script generated: $ts
# matrix: $ARGV[0]
# peak bed: $ARGV[1]
# genome: $opt{'g'}
# motif set: $opt{'M'}

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
data($opt{'M'})

# read in peaks and filter them
peaks <- getPeaks(\"$ARGV[1]\",sort=TRUE)
peaks <- resize(peaks, width = 500, fix = \"center\")

# make counts matrix and normalize / filter
counts_matrix <- as.matrix(read.delim(\"$ARGV[0]\",row.names=1,header=T))
counts <- SummarizedExperiment(assays = list(counts = counts_matrix), rowRanges = peaks)

counts <- addGCBias(counts, genome = $genome)

# make motif index
motif_ix <- matchMotifs($opt{'M'}, counts, genome = $genome)

# calculate & print deviations
dev <- computeDeviations(object = counts, annotations = motif_ix)
write.table(as.matrix(deviations(dev)),file = \"$opt{'O'}.chromVAR/deviations.matrix\", col.names = TRUE, row.names = TRUE, sep = \"\\t\", quote = FALSE)
write.table(as.matrix(deviationScores(dev)),file = \"$opt{'O'}.chromVAR/deviation_scores.matrix\", col.names = TRUE, row.names = TRUE, sep = \"\\t\", quote = FALSE)

# calculate & print variabilities
var <- computeVariability(dev)
write.table(as.matrix(var),file = \"$opt{'O'}.chromVAR/variability.txt\", col.names = TRUE, row.names = TRUE, sep = \"\\t\", quote = FALSE)
plot<-plotVariability(var,use_plotly=FALSE)
ggsave(plot,file = \"$opt{'O'}.chromVAR/variability.png\")

"; close R;

    system("$Rscript $opt{'O'}.chromVAR/chromVAR.r >> $opt{'O'}.chromVAR/chromVAR.log 2>> $opt{'O'}.chromVAR/chromVAR.log");

}
1;
