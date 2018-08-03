package sci_commands::cds_cicero;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("cds_cicero");

sub cds_cicero {

@ARGV = @_;
# Defaults

getopts("O:R:XG:", \%opt);

$die2 = "
scitools cds-cicero [options] [directory containing cds files]
   or    cicero-cds

Prior to running, ensure you have ran scitools matrix-makecds. 
That function will convert matrix files into the CDS format that Cicero requires.
Runs the current version of Cicero within a directory containing files for CDS format input. 


Options:
   -O    [STR]   Output Directory (default is [current working directory]/cicero_output)
   -R    [STR]   Rscript call (def = $Rscript)
   -G    [STR]   Genome to be used. Must be of list:
                  human.hg38.genome,human.hg19.genome,human.hg38.genome,mouse.mm10.genome
                  Default: human.hg38.genome
   -X           Retain intermediate files (Default = delete)
                  
";


if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = "$ARGV[0]/cicero_output"};
if (!defined $opt{'D'}) {$opt{'G'} = "human.hg38.genome"};

system("mkdir $opt{'O'}");

open R, ">$opt{'O'}/cicero.r";

print R "
suppressWarnings(library(cicero))
message(\"Loading Cicero\")
# reading in matrix, annotation, and dimension data

cds_cell_data <- read.delim(\"$ARGV[0]/cds_cell_data.txt\")
cds_dims_data <- read.delim(\"$ARGV[0]/cds_dims_data.txt\",header=F)
cds_site_data <- read.delim(\"$ARGV[0]/cds_site_data.txt\")
cds_counts_matrix <- read.table(\"$ARGV[0]/cds_counts_matrix.txt\")
message(\"Read in CDS Files.\")

# generating cell data set
feature_data <- new(\"AnnotatedDataFrame\", data = cds_site_data)
sample_data <- new(\"AnnotatedDataFrame\", data = cds_cell_data)

if(ncol(cds_dims_data)==4){
colnames(cds_dims_data)<-c(\"cellID\",\"Dimension_1\",\"Dimension_2\",\"Dimension_3\")
row.names(cds_dims_data)<-cds_dims_data\$cellID
dimension_reduction<-cds_dims_data[2:ncol(cds_dims_data)]
}else{
colnames(cds_dims_data)<-c(\"cellID\",\"Dimension_1\",\"Dimension_2\")
row.names(cds_dims_data)<-cds_dims_data\$cellID
dimension_reduction<-cds_dims_data[2:ncol(cds_dims_data)]
}

message(\"Setting up CDS matrix, binarized for Cicero.\")


cds <- suppressWarnings(newCellDataSet(as.matrix(cds_counts_matrix), phenoData = sample_data, featureData = feature_data,expressionFamily = negbinomial.size(),lowerDetectionLimit = 0))

cds\@expressionFamily <- binomialff()
cds\@expressionFamily\@vfamily <- \"binomialff\"

pData(cds)\$temp <- NULL
fData(cds)\$chr <- as.numeric(as.character(fData(cds)\$chr))
fData(cds)\$bp1 <- as.numeric(as.character(fData(cds)\$bp1))
fData(cds)\$bp2 <- as.numeric(as.character(fData(cds)\$bp2))
cds <- cds[order(fData(cds)\$chr, fData(cds)\$bp1),]

set.seed(2017)
pData(cds)\$cells <- NULL
cds <- aggregate_nearby_peaks(cds, distance = 10000)
cds <- detectGenes(cds)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

message(\"Using Supplied Dims File.\")
cds <- make_cicero_cds(cds, reduced_coordinates = dimension_reduction)

message(\"Running Cicero\")
data(\"human.hg38.genome\")
conns <- run_cicero(cds, human.hg38.genome) 
# Takes a few minutes to run

message(\"Sample Cicero Output:\")
head(conns)

message(\"Generating CCANS.\")
CCAN_assigns <- generate_ccans(conns)

message(\"Sample CCANs Output:\")
head(CCAN_assigns)

saveRDS(cds,\"$opt{'O'}/cicero.CDS.rds\")
write.table(as.data.frame(conns), file=\"$opt{'O'}/cicero.output.txt\",quote=F,sep=\"\\t\",row.names=T,col.names=T)
write.table(as.data.frame(CCAN_assigns), file=\"$opt{'O'}/cicero.CCANS.txt\",quote=F,sep=\"\\t\",row.names=T,col.names=T)
";

close R;
system("$Rscript $opt{'O'}/cicero.r");
if (!defined $opt{'X'}) {
    system("rm -f $opt{'O'}/cicero.r");
}
}
1;
