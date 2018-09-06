package sci_commands::matrix_cistopic;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_cistopic");

sub matrix_cistopic {

$dims = 50;

@ARGV = @_;
getopts("O:XR:D:", \%opt);

$die2 = "
scitools matrix-cistopic [options] [input matrix]
   or    cistopic

[input matrix] =	A counts matrix generated from atac_count. 
					Can be filtered prior to this command call, or 
					can supply the filtering options during cisTopic processing.

cisTopic serves as an alternative textmining algorithm to TFIDF-LSI processing.
It is to be run on a sciATAC counts matrix. For more information see:
https://github.com/aertslab/cisTopic/

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
   -D   [INT]   Max dimensions to compute (def = $dims)
   -X           Retain intermediate files (def = delete)
   -R   [STR]   Rscript call (def = $Rscript)

Note: Requires cisTopic R package

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};
if (!defined $opt{'c'}) {$opt{'c'} = 1;
if (!defined $opt{'r'}) {$opt{'r'} = 1;
if (!defined $opt{'n'}) {$opt{'n'} = 1;

if (defined $opt{'R'}) {$Rscript = $opt{'R'}};

open R, ">$opt{'O'}.cistopic.r";
print R "
library(cisTopic)
library(Rtsne)
IN<-as.matrix(read.table(\"$ARGV[0]\"))
row.names(IN)<-sub(\"_\",\"-\",sub(\"_\",\":\",row.names(IN)))

#Set up filtered binarized counts matrix
cisTopicObject <- createcisTopicObject(IN,min.cells=$opt{'r'},min.regions=$opt{'c'}, keepCountsMatrix=FALSE)

#Run Multiple models with different topic numbers. Optimal topic number is generally slightly bigger than the potential cell states in the data set
cisTopicObject <- runModels(cisTopicObject, topic=c(2, 5, 10, 15, 20, 25, 30), seed=2018, nCores=$opt{'n'}, burnin = 250, iterations = 500)

#Plot and select model based on the highest log likelihood (P(D|T)) at the last iteration
cisTopicObject <- selectModel(cisTopicObject)
logLikelihoodByIter(cisTopicObject)

#Run Native tSNE
cisTopicObject<-runtSNE(cisTopicObject,seed=2018)
write.table(cisTopicObject@dr,file=\"$opt{'O'}.cistopic.tsne.dims\",col.names=T,row.names=T,quote=F,sep=\"\\t\")


#Print out cisTopic Matrix#
modelMat <- scale(cisTopicObject@selected.model$document_expects, center = TRUE, scale = TRUE)
tModelmat<-as.data.frame(t(modelMat))
Modeldf<-as.data.frame(modelMat)
rownames(tModelmat)<-cisTopicObject@cell.names
colnames(Modeldf)<-cisTopicObject@cell.names
write.table(tModelmat,file=\"$opt{'O'}.cistopic.dim\",col.names=T,row.names=T,quote=F,sep=\"\\t\")
row.names(Modeldf)<-paste0(\"Dim_\",row.names(Modeldf))
write.table(Modeldf,file=\"$opt{'O'}.cistopic.matrix\",col.names=T,row.names=T,quote=F,sep=\"\\t\")
";
close R;

system("$Rscript $opt{'O'}.cistopic.r");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.cistopic.r");
}

}
1;
