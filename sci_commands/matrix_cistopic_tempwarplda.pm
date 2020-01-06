package sci_commands::matrix_cistopic_tempwarplda;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_cistopic_tempwarplda");

sub matrix_cistopic_tempwarplda {

    @ARGV = @_;
    getopts("O:A:c:r:n:T:t:XR:P:L:W:", \%opt);

    $thrP = 0.975;

$die2 = "
scitools matrix-cistopic-tempwarplda [options] [input matrix]
   or    cistopic

[input matrix] =A counts matrix generated from atac_count. 
    Can be filtered prior to this command call, or 
    can supply the filtering options during cisTopic processing.

cisTopic serves as an alternative textmining algorithm to TFIDF-LSI processing.
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
   -n   [INT] Number of cores for parallel processing. (def=1)
   -L   [STR]  Columns file for sparseMatrix (will try to auto-detect)
   -W   [STR]  Rows file for sparseMatrix (will try to auto-detect)
   -t   [STR]  prefix for temporary topic RDS files for each model
   -T   [INT]  User defined number of Topics to use. 
                  If unspecified: will generate 15,20,25,30,50,65,100 Topics,
                  and use log-liklihood estimators to select the best.
                  Specification can be a single number of a comma separated list.
                  Will use a core for each number supplied (DO NOT EXCEED A LIST LENGTH OF 10)
   -P   [FLT]  ThrP (def = $thrP)
   -X          Retain intermediate files (def = delete)
   -R   [STR]  Rscript call (def = $Rscript)

    Note: Requires cisTopic R package, use: /home/groups/oroaklab/src/R/R-3.5.1/library_arsn

";


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
if (!defined $opt{'T'}) {$opt{'T'} = "15,20,25,30,50,65,100"};
if (!defined $opt{'t'}) {$opt{'t'} = "temp_cistopics_warplda/"};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};

###

system("mkdir $opt{'t'}");

open R, ">$opt{'O'}.cistopic.r";
print R "
library(plyr)
library(cisTopic)
library(Matrix)
source(\"/home/groups/oroaklab/nishida/scripts/cistopics_runWarpLDA.R\")
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
cisTopicObject <- createcisTopicObject(IN,min.cells=$opt{'r'},min.regions=$opt{'c'}, keepCountsMatrix=FALSE)
cisTopicObject <- runWarpLDA(cisTopicObject, topic=c($opt{'T'}), seed=2018, nCores=$opt{'n'}, iterations = 300, tmp=\"$opt{'t'}\")

if (length(c($opt{'T'}))>1){
cisTopicObject <- selectModel(cisTopicObject)
modelMat <- scale(cisTopicObject\@selected.model\$document_expects, center = TRUE, scale = TRUE)
} else {
modelMat<-scale(cisTopicObject\@models\$document_expects,center=TRUE,scale=TRUE)
}

#Print out cisTopic Matrix#
tModelmat<-as.data.frame(t(modelMat))
Modeldf<-as.data.frame(modelMat)
rownames(tModelmat)<-cisTopicObject\@cell.names
colnames(Modeldf)<-cisTopicObject\@cell.names
row.names(Modeldf)<-paste0(\"Topic_\",row.names(Modeldf))

write.table(Modeldf,file=\"$opt{'O'}.cistopic_warplda.matrix\",col.names=T,row.names=T,quote=F,sep=\"\\t\")
saveRDS(cisTopicObject,\"$opt{'O'}.cistopicObject_warplda.rds\")
";

close R;

system("$Rscript $opt{'O'}.cistopic.r");

if (!defined $opt{'X'}) {
    system("rm -f $opt{'O'}.cistopic.r");
}

}
1;
