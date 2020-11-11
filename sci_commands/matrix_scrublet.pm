package sci_commands::matrix_scrublet;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_scrublet");

sub matrix_scrublet {

@ARGV = @_;

getopts("O:P:", \%opt);

$die2 = "
scitools scrublet [options] [matrix (dense only, eg topic matrix)]

Options:
   -O   [STR]   Output prefix (default is matrix file prefix)
   -P   [STR]   python script call (def = $Pscript)

   
Note: Requires python, numpy, and scrublet to be installed and callable
      This command is a wrapper for executing the python code.

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.matrix$//; $opt{'O'} =~ s/\.dims$//;
if (defined $opt{'P'}) {$Pscript = $opt{'P'}};
read_matrix($ARGV[0]);

open OUT, ">$opt{'O'}.scrublet.py";
print OUT "
import numpy
import scrublet as scr
data_matrix=numpy.loadtxt(\"$ARGV[0]\",skiprows=1,usecols=range(1,$matrix_colNum))
scrub = scr.Scrublet(data_matrix.T)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
numpy.savetxt(\"$opt{'O'}.temp.doublet_scores\",doublet_scores,delimiter=\"\\t\")
numpy.savetxt(\"$opt{'O'}.temp.predicted_doublets\",predicted_doublets,delimiter=\"\\t\")
";

system("$Pscript $opt{'O'}.scrublet.py");

}
1;
