package sci_commands::fastq_barcode_collapse;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("fastq_barcode_collapse");

sub fastq_barcode_collapse {

# defaults
$hdist = 1;

@ARGV = @_;
getopts("O:", \%opt);

$die2 = "
scitools fastq-barcode-collapse [options] [index_read_fastq] [read1_fastq] (read2_fastq)

Options:
   -O   [STR]   Output prefix
   -H   [STR]   Hamming distance (def = $hdist)

";

if (!defined $ARGV[0]) {die $die2};
if (defined $opt{'H'}) {$hdist = $opt{'H'}};


}
1;
