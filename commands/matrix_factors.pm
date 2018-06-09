package commands::matrix_factors;

use commands::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_factors");

sub matrix_factors {

@ARGV = @_;
getopts("O:XR:", \%opt);

$die2 = "
scitools matrix-approx_factors [options] [input tf-idf matrix]
   or    factors

Produces a k-range vs error txt of NMF

Options:
   -O   [STR]   Output prefix (default is [input].factors.txt)
   -X           Retain intermediate files (def = delete)
   -R   [STR]   Rscript call (def = $Rscript)

Note: Requires Seurat and swne R packages so dependencies need to look for those

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};

open R, ">$opt{'O'}.factors.r";
print R "
library(Seurat)
library(swne)
norm<-as.matrix(read.table(file=\"$ARGV[0]\",row.names=1), \"dgCMatrix\")

#read in tf-idf matrix into dgc format, might want to do original count matrix
norm.counts<-as.matrix(norm,\"dgCMatrix\")

#or do original matrix and maybe frequency transform
#norm.counts <- ScaleCounts(counts, batch = NULL, method = \"ft\", adj.var = T)

## Unguided NMF
loss <- \"mse\" ## Loss function
#loss <- \"mkl\" ## Loss function
k.range <- seq(1,100,1) ## Range of factors to iterate over
n.cores <- 25 ## Number of cores to use
seed <- 32566 ## Set seed for 

## Identify optimal number of factors

n.comp.res <- FindNumFactors(norm.counts, k.range = k.range, n.cores = n.cores, do.plot = F, loss = loss,max.iter = 10000)
output=data.frame(\"k\"=k.range,\"recon.err\"=n.comp.res\$err[loss,])
output2=data.frame(\"k\"=k.range,n.comp.res\$err)

#output k vs error txt to plot with plot-k
write.table(output,file=\"$opt{'O'}.factors.txt\",sep=\"\\t\",row.names=FALSE,col.names=TRUE,quote=FALSE)

#output k vs error txt to plot with plot-k
write.table(output2,file=\"$opt{'O'}.allerr.factors.txt\",sep=\"\\t\",row.names=FALSE,col.names=TRUE,quote=FALSE)


";
close R;

system("$Rscript $opt{'O'}.factors.r");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.factors.r");
}

}
1;
