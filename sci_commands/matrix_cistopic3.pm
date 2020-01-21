package sci_commands::matrix_cistopic3;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_cistopic3");

sub matrix_cistopic3 {

    @ARGV = @_;
    getopts("O:c:r:n:T:XR:P:a:b:L:W:", \%opt);

    $thrP = 0.975;
    $cis_alpha = 50;
    $cis_beta = 0.1;

$die2 = "
scitools matrix-cistopic3 [options] [input matrix]

[input matrix] =A counts matrix generated from atac_count. 
    Can be filtered prior to this command call, or 
    can supply the filtering options during cisTopic processing.

cisTopic serves as an alternative textmining algorithm to TFIDF-LSI processing.
Version 3 uses Warp LDA implementation for topic modeling.
It is to be run on a sciATAC counts matrix. For more information see:
github.com/aertslab/cisTopic/

Outputs a matrix file similar to matrix-irlba function call. To be processed through
[matrix_tsne|matrix_umap|matrix_PCA|matrix_SWNE]

cisTopic consists of 4 main steps: 
(1) generation of a binary accessibility matrix as input for LDA; 
(2) LDA and model selection; 
(3) cell state identification using the topic-cell distributions from LDA and 
(4) exploration of the region-topic distributions.

    Options:
   -O   [STR]  Output prefix (default is [input].cistopic.dims)
   -c   [INT]  Number of nonZero sites per column (cell) to retain (def = 1)
   -r   [INT]  Number of nonZero sites per row (peak) to retain (def = 1)
   -n   [INT]  Number of cores for parallel processing. (def=1)
   -L   [STR]  Columns file for sparseMatrix (will try to auto-detect), currently unsupported
   -W   [STR]  Rows file for sparseMatrix (will try to auto-detect), currently unsupported
   -T   [INT]  User defined number of Topics to use. 
                  If unspecified: will generate 15,20,25,30,50,65,100 Topics,
                  and use log-liklihood estimators to select the best.
                  Specification can be a single number of a comma separated list.
                  Will use a core for each number supplied (DO NOT EXCEED A LIST LENGTH OF 10)
   -P   [FLT]  ThrP (def = $thrP)
   -a   [FLT]  cistopic alpha (def = $cis_alpha)
   -b   [FLT]  cistopic beta (def = $cis_beta)
   -X          Retain intermediate files (def = delete)
   -R   [STR]  Rscript call (def = $Rscript)

    Note: Requires cisTopic v3 R package

";

###

if (!defined $ARGV[0]) {die $die2};

$prefix = $ARGV[0];
$prefix =~ s/\.gz$//;
$prefix =~ s/\.matrix$//;
$prefix =~ s/\.values$//;
$prefix =~ s/\.sparseMatrix$//;
if (!defined $opt{'O'}) {
   $opt{'O'} = $prefix;
}

if ($ARGV[0] =~ /sparseMatrix/i) {
   $sparse = 1;
   if (defined $opt{'L'}) {
      $col_file = $opt{'L'};
   } else {
      if (-e "$prefix.sparseMatrix.cols") {
         $col_file = "$prefix.sparseMatrix.cols";
      } elsif (-e "$prefix.sparseMatrix.cols.gz") {
         $col_file = "$prefix.sparseMatrix.cols.gz";
   } else {
      die "ERROR: Cannot detect cols file (e.g. $prefix.sparseMatrix.cols), please provide as -C\n";
   }
}
if (defined $opt{'W'}) {
   $row_file = $opt{'W'};
} else {
   if (-e "$prefix.sparseMatrix.rows") {
      $row_file = "$prefix.sparseMatrix.rows";
   } elsif (-e "$prefix.sparseMatrix.rows.gz") {
      $row_file = "$prefix.sparseMatrix.rows.gz";
   } else {
      die "ERROR: Cannot detect rows file (e.g. $prefix.sparseMatrix.rows), please provide as -C\n";
   }
}
} else {$sparse = 0};

if (!defined $opt{'c'}) {$opt{'c'} = 1};
if (!defined $opt{'r'}) {$opt{'r'} = 1};
if (!defined $opt{'n'}) {$opt{'n'} = 1};
if (defined $opt{'P'}) {$thrP = $opt{'P'}};
if (defined $opt{'a'}) {$cis_alpha = $opt{'a'}};
if (defined $opt{'b'}) {$cis_beta = $opt{'b'}};
if (!defined $opt{'T'}) {$opt{'T'} = "15,20,25,30,50,65,100"};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};

###

open R, ">$opt{'O'}.cistopic.r";
print R "
suppressWarnings(library(cisTopic))
library(Matrix)
library(ComplexHeatmap)
";

if ($ARGV[0] =~ /sparseMatrix/i) {
print R "IN<-as.matrix(read.table(\"$ARGV[0]\"))
IN<-sparseMatrix(i=IN[,1],j=IN[,2],x=IN[,3])
COLS<-read.table(\"$col_file\")
colnames(IN)<-COLS\$V1
ROWS<-read.table(\"$row_file\")
row.names(IN)<-ROWS\$V1\n";
} else {
print R "IN<-as.matrix(read.table(\"$ARGV[0]\"))\n";
}

print R "
row.names(IN)<-sub(\"_\",\"-\",sub(\"_\",\":\",row.names(IN)))
#Set up filtered binarized counts matrix
cisTopicObject <- createcisTopicObject(IN,min.cells=$opt{'r'},min.regions=$opt{'c'}, keepCountsMatrix=FALSE)
";

print R "
cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=c($opt{'T'}), seed=2020, nCores=$opt{'n'}, iterations=500, alpha=$cis_alpha, beta=$cis_beta, tmp=\"$opt{'O'}.temp\")
saveRDS(cisTopicObject,\"$opt{'O'}.cistopicObject.rds\")

# output all matrices
for (i in 1:length(cisTopicObject\@models)) {
    model <- cisTopicObject\@models[i]
    modelMat <- scale(model[[1]]\$document_expects, center = TRUE, scale = TRUE)
    tModelmat<-as.data.frame(t(modelMat))
    Modeldf<-as.data.frame(modelMat)
    rownames(tModelmat)<-cisTopicObject\@cell.names
    colnames(Modeldf)<-cisTopicObject\@cell.names
    row.names(Modeldf)<-paste0(\"Topic_\",row.names(Modeldf))
    write.table(Modeldf,file=paste(\"$opt{'O'}.cistopic.\",nrow(Modeldf),\".matrix\",sep=\"\"),col.names=T,row.names=T,quote=F,sep=\"\t\")
}

pdf(file=\"$opt{'O'}.cistopic.modelselection.pdf\")
par(mfrow=c(3,1))
cisTopicObject <- selectModel(cisTopicObject, type='derivative')
dev.off()

cisTopicObject <- getRegionsScores(cisTopicObject, method='Z-score', scale=TRUE)
cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=$thrP, plot=FALSE)
getBedFiles(cisTopicObject, path='$opt{'O'}.topics')
saveRDS(cisTopicObject,\"$opt{'O'}.cistopicObject.rds\")
";

print R "
png(\"$opt{'O'}.Heatmap_prob_cistopic.png\",width=12,height=12,units=\"in\",res=600)
cellTopicHeatmap(cisTopicObject, method=\'Probability\')
dev.off()
pdf(\"$opt{'O'}.Heatmap_prob_cistopic.pdf\",width=12,height=12)
cellTopicHeatmap(cisTopicObject, method=\'Probability\')
dev.off()
png(\"$opt{'O'}.Heatmap_zscore_cistopic.png\",width=12,height=12,units=\"in\",res=600)
cellTopicHeatmap(cisTopicObject, method=\'Z-score\')
dev.off()
pdf(\"$opt{'O'}.Heatmap_zscore_cistopic.pdf\",width=12,height=12)
cellTopicHeatmap(cisTopicObject, method=\'Z-score\')
dev.off()
";

close R;

system("$Rscript $opt{'O'}.cistopic.r");

if (!defined $opt{'X'}) {
system("rm -f $opt{'O'}.cistopic.r");
}

}
1;
