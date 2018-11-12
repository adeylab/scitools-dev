package sci_commands::matrix_cistopic;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_cistopic");

sub matrix_cistopic {

@ARGV = @_;
getopts("O:A:c:r:n:T:SXR:", \%opt);

$die2 = "
scitools matrix-cistopic [options] [input matrix]
   or    cistopic

[input matrix] =	A counts matrix generated from atac_count. 
					Can be filtered prior to this command call, or 
					can supply the filtering options during cisTopic processing.

cisTopic serves as an alternative textmining algorithm to TFIDF-LSI processing.
It is to be run on a sciATAC counts matrix. For more information see:
https://github.com/aertslab/cisTopic/

Outputs a matrix file similar to matrix-irlba function call. To be processed through
[matrix_tsne|matrix_umap|matrix_PCA|matrix_SWNE]

cisTopic consists of 4 main steps: 
(1) generation of a binary accessibility matrix as input for LDA; 
(2) LDA and model selection; 
(3) cell state identification using the topic-cell distributions from LDA and 
(4) exploration of the region-topic distributions.

Options:
   -O   [STR]   Output prefix (default is [input].cistopic.dims)
   -c   [INT]   Number of nonZero sites per column (cell) to retain (def = 1)
   -r   [INT]   Number of nonZero sites per row (peak) to retain (def = 1)
   -n 	[INT] 	Number of cores for parallel processing. (def=1)
   -A   [STR]   Annotation file is useful. cisTopic will provide a PCA with the influence of Topics
   -T    [INT]    User defined number of Topics to use. 
                  If unspecified: will generate 15, 20, 25, 30, 50, 65, 100 Topics, 
                  and use log-liklihood estimators to select the best.
                  Specification can be a single number of a comma separated list.
                  Will use a core for each number supplied (DO NOT EXCEED A LIST LENGTH OF 10)
   -S			If defined CDS will be retained in RDS format for further analysis			
   -X           Retain intermediate files (def = delete)
   -R   [STR]   Rscript call (def = $Rscript)

Note: Requires cisTopic R package

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};
if (!defined $opt{'c'}) {$opt{'c'} = 1};
if (!defined $opt{'r'}) {$opt{'r'} = 1};
if (!defined $opt{'n'}) {$opt{'n'} = 1};
if (!defined $opt{'T'}) {$opt{'T'} = "15,20,25,30,50,65,100"};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};

open R, ">$opt{'O'}.cistopic.r";
print R "
library(plyr)
library(cisTopic)
IN<-as.matrix(read.table(\"$ARGV[0]\"))
row.names(IN)<-sub(\"_\",\"-\",sub(\"_\",\":\",row.names(IN)))

#Set up filtered binarized counts matrix
cisTopicObject <- createcisTopicObject(IN,min.cells=$opt{'r'},min.regions=$opt{'c'}, keepCountsMatrix=FALSE)
";
if (defined $opt{'A'})
{
print R "
annot<-read.table(file=\"$opt{'A'}\",header=F)
#match rows of annot
annot<-annot[match(colnames(IN), annot\$V1),]
row.names(annot)<-anno\$V1
names(annot)<-c(\"cellname\",\"LineType\")
cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = annot)
";

}



print R "
cisTopicObject <- runModels(cisTopicObject, topic=c($opt{'T'}), seed=2018, nCores=$opt{'n'}, burnin = 250, iterations = 300)

if (length(c($opt{'T'}))>1){
#Future update:Plot model log likelihood (P(D|T)) at the last iteration
pdf(file=\"$opt{'O'}.cistopic.modelselection.pdf\")
cisTopicObject <- selectModel(cisTopicObject)
logLikelihoodByIter(cisTopicObject)
dev.off()  
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
write.table(Modeldf,file=\"$opt{'O'}.cistopic.matrix\",col.names=T,row.names=T,quote=F,sep=\"\\t\")
#adding part where the contribution matrix is calculated, we use a binarization method to select for peaks that contribute to each topic
cisTopicObject <- getRegionsScores(cisTopicObject, method='Zscore', scale=TRUE)
cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=0.975, plot=FALSE)
getBedFiles(cisTopicObject, path='$opt{'O'}')
saveRDS(cisTopicObject,\"$opt{'O'}.cistopicObject.rds\")
";

if (defined $opt{'A'}){
print R "
cisTopicObject <- runPCA(cisTopicObject)
coordinates <- object\@dr[[\'PCA\']]\$ind.coord
write.table(coordinates,file=\"$opt{'O'}.PCA.internal.dims\",col.names=T,row.names=T,quote=F,sep=\"\\t\")
png(file=\"PCA_cistopic.png\",width=12,height=12,units=\"in\",res=300)
plotCellStates(cisTopicObject, method=\'Biplot\', topic_contr=\'Zscore\',topics=\'all\', colorBy=c(\'LineType\'))
dev.off()
plotCellStates(cisTopicObject, method=\'Biplot\', topic_contr=\'Zscore\', topics=\'all\', colorBy=c(\'LineType\'))
";
}

close R;

system("$Rscript $opt{'O'}.cistopic.r");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.cistopic.r");
}

}
1;
