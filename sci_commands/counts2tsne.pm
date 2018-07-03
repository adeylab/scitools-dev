package sci_commands::counts2tsne;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("counts2tsne");

sub counts2tsne {

@ARGV = @_;
# Deafults
$irlba_dims = 50;
$tsne_dims = 2;
$perp = 30;

getopts("O:XD:R:P:d:", \%opt);

$die2 = "
scitools counts2tsne [options] [input counts matrix]

Takes in a counts matrix (filtered) and will run tfidf,
irlba, and then rtsne.

Options:
   -O   [STR]   Output prefix (default is [input].tSNE)
   -D   [INT]   Dimensions to use from irlba (def = $irlba_dims)
   -d   [INT]   Dimensions to embed tSNE (def = $tsne_dims)
   -P   [INT]   Perplexity (def = $perp)
   -X           Retain intermediate files (def = delete)
   -R   [STR]   Rscript call (def = $Rscript)

Note: Requires irlba and Rtsne R packages

";

if (!defined $ARGV[0]) {die $die2};
if ($ARGV[0] =~ /\.h5$/) {$opt{'H'} = 1};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};
if (defined $opt{'D'}) {$irlba_dims = $opt{'D'}};
if (defined $opt{'d'}) {$tsne_dims = $opt{'d'}};
if (defined $opt{'P'}) {$perp = $opt{'P'}};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};

open R, ">$opt{'O'}.counts2tsne.r";

print R "
library(irlba)
library(Rtsne)

# load in counts matrix
COUNTS<-as.matrix(read.table(\"$ARGV[0]\"))

# tfidf
TFIDF<-(COUNTS/colSums(COUNTS))*(log(1+(ncol(COUNTS)/(rowSums(COUNTS)+1))))

# irlba
IRLBA<-(TFIDF,$irlba_dims)

# tSNE
TSNE(IRLBA\$v,dims=$tsne_dims,perplexity=$perp,check_duplicates=FALSE,pca=FALSE)

# output
rownames(TSNE\$Y)<-colnames(COUNTS)
write.table(TSNE\$Y,file=\"$opt{'O'}.counts2tsne.dims\",sep=\"\\t\",col.names=FALSE,row.names=TRUE,quote=FALSE)

";

close R;

system("$Rscript $opt{'O'}.counts2tsne.r 2>/dev/null");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.counts2tsne.r");
}

}
1;