package sci_commands::matrix_irlba;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_irlba");

sub matrix_irlba {

$dims = 50;

@ARGV = @_;
getopts("O:XR:D:", \%opt);

$die2 = "
scitools matrix-irlba [options] [input matrix]
   or    irlba

Produces a dims file with irlba coordinates

Options:
   -O   [STR]   Output prefix (default is [input].irlba.dims)
   -D   [INT]   Max dimensions to compute (def = $dims)
   -X           Retain intermediate files (def = delete)
   -R   [STR]   Rscript call (def = $Rscript)

Note: Requires irlba R package

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};

open R, ">$opt{'O'}.irlba.r";
print R "
library(irlba)
library(ggplot2)
IN<-as.matrix(read.table(\"$ARGV[0]\"))
I<-irlba(IN, $dims)
rownames(I\$v)<-colnames(IN)

var<-as.data.frame(I\$d,ncol=1)
colnames(var)<-c(\"d\")
var\$var_names<-paste0(\"PC_\",row.names(var))
var\$var_names<-factor(var\$var_names,levels=var\$var_names)

var\$variance<-prop.table(var\$d^2)
ggplot()+geom_histogram(aes(x=var\$var_names,y=var\$variance),stat=\"identity\")+xlab(\"PC\")+ylab(\"Variance Explained\")+theme_minimal()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(\"$opt{'O'}.irlba_variance.png\")

df<-as.data.frame(I\$v)
PCA_colnames<-c()
for (i in 1:$dims) {
  p<-paste(\"PCA\",i, sep = \"\")
  PCA_colnames <- append(PCA_colnames, p)
}
colnames(df)<-PCA_colnames
dft<-t(df)

write.table(dft, file=\"$opt{'O'}.irlba_$dims.dims\", quote=FALSE, sep=\"\\t\")

";
close R;

system("$Rscript $opt{'O'}.irlba.r");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.irlba.r");
}

}
1;
