package sci_commands::matrix_pg;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_pg");

sub matrix_pg {

@ARGV = @_;


getopts("O:X:P:", \%opt);

$die2 = "
scitools matrix-pg [options] [input matrix or dims] 

This scripts uses Louvain method to find clusters of cells based on provided dims or matrix

Options:
   -O   [STR]   Output prefix (default is [input].pg.annot)
   -X           Retain intermediate files (def = delete)
   -P   [STR]   python script call (def = $Pscript)

Note: Requires numpy and phenograph python packages

";


if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.matrix$//;
$opt{'O'} =~ s/\.dims$//;
if (defined $opt{'P'}) {$Pscript = $opt{'P'}};
read_matrix($ARGV[0]);

open OUT, ">$opt{'O'}.pg.py";
print OUT"
import phenograph
import numpy
";
if ($ARGV[0] =~ /\.dims/) {
print OUT "
data_matrix=numpy.loadtxt(\"$ARGV[0]\",skiprows=1,usecols=range(1,$matrix_colNum))\n";
} else {
print OUT "ori_matrix=numpy.loadtxt(\"$ARGV[0]\",skiprows=1,usecols=range(1,$matrix_colNum))
data_matrix=ori_matrix.transpose()\n";
}


#detext communities. Graph, later
print OUT "

communities, graph, Q = phenograph.cluster(data_matrix)
numpy.savetxt(\"$opt{'O'}.temp.pg.annot\",communities,delimiter=\"\\t\")
";
close OUT;
system("$Pscript $opt{'O'}.pg.py");

open OUT, ">$opt{'O'}.UMAP.dims";
$counter=0;
open IN, "$opt{'O'}.temp.pg.annot"; 
while($l=<IN>)
{
$l =~ s/"//g;   
print OUT $MATRIX_COLNAMES[$counter]."\t".$l;
$counter++;
}
close(IN);
close OUT;


if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.temp.pg.annot $opt{'O'}.pg.py");
}

}
1;