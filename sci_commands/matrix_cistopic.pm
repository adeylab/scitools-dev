package sci_commands::matrix_cistopic;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_cistopic");

sub matrix_cistopic {

@ARGV = @_;
getopts("O:c:r:n:T:XR:", \%opt);

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
   -T    [INT]    User defined number of Topics to use. 
                  If unspecified: will generate 15, 20, 25, 30, 50, and 65 Topics, 
                  and use log-liklihood estimators to select the best.
   -D   [INT]   Max dimensions to compute (def = $dims)
   -X           Retain intermediate files (def = delete)
   -R   [STR]   Rscript call (def = $Rscript)

Note: Requires cisTopic R package

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};
if (!defined $opt{'c'}) {$opt{'c'} = 1};
if (!defined $opt{'r'}) {$opt{'r'} = 1};
if (!defined $opt{'n'}) {$opt{'n'} = 1};

if (defined $opt{'R'}) {$Rscript = $opt{'R'}};

open R, ">$opt{'O'}.cistopic.r";
print R "
library(cisTopic)
IN<-as.matrix(read.table(\"$ARGV[0]\"))
row.names(IN)<-sub(\"_\",\"-\",sub(\"_\",\":\",row.names(IN)))

#Set up filtered binarized counts matrix
cisTopicObject <- createcisTopicObject(IN,min.cells=$opt{'r'},min.regions=$opt{'c'}, keepCountsMatrix=FALSE)

";

if (!defined $opt{'T'}) {
print R "
#Run Multiple models with different topic numbers. Optimal topic number is generally slightly bigger than the potential cell states in the data set
cisTopicObject <- runModels(cisTopicObject, topic=c(15, 20, 25, 30, 50, 65, 100), seed=2018, nCores=$opt{'n'}, burnin = 250, iterations = 300)
saveRDS(cisTopicObject,\"$opt{'O'}.cistopicObject.rds\")

#Future update:Plot model log likelihood (P(D|T)) at the last iteration
pdf(file=\"$opt{'O'}.cistopic.modelselection.pdf\")
cisTopicObject <- selectModel(cisTopicObject)
logLikelihoodByIter(cisTopicObject)
dev.off()
";
} else {
print R "
cisTopicObject <- runModels(cisTopicObject, topic=$opt{'T'}, seed=2018, nCores=$opt{'n'}, burnin = 250, iterations = 300)
cisTopicObject <- selectModel(cisTopicObject)
";
}

print R "
#Print out cisTopic Matrix#
modelMat <- scale(cisTopicObject\@selected.model\$document_expects, center = TRUE, scale = TRUE)
tModelmat<-as.data.frame(t(modelMat))
Modeldf<-as.data.frame(modelMat)
rownames(tModelmat)<-cisTopicObject\@cell.names
colnames(Modeldf)<-cisTopicObject\@cell.names
row.names(Modeldf)<-paste0(\"Topic_\",row.names(Modeldf))
write.table(Modeldf,file=\"$opt{'O'}.cistopic.matrix\",col.names=T,row.names=T,quote=F,sep=\"\\t\")
";
close R;

system("$Rscript $opt{'O'}.cistopic.r");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.cistopic.r");
}

}
1;
