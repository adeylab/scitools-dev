package sci_commands::fastq_simulate;

use sci_utils::general;
use sci_utils::modes;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("fastq_simulate");

sub fastq_simulate {

@ARGV = @_;
getopts("O:M:N:m:I:", \%opt);

# Defaults
$mode_name = "sci";
$num_reads = 100000;

$die2 = "
scitools fastq-simulate [options] [output prefix]
   or    simulate-fastq

Takes sequencer fastq files (from bcl2fastq) and will format
them into fastq files with matched barcodes.

Options:
   -M   [STR]   Mode (def = $mode_name)
   -N   [INT]   Number of reads to simulate (def = $num_reads)
   -m   [STR]   Run format specification file
         (def = $VAR{'sci_modes'})
   -I   [STR]   Index files/directory - comma separated
         Def: DIR=$VAR{'index_directory'}

";

if (!defined $ARGV[0] && !defined $opt{'O'}) {die $die2};
if (defined $opt{'O'}) {$ARGV[0] = $opt{'O'}};
if (defined $opt{'M'}) {$mode_name = $opt{'M'}};
if (defined $opt{'N'}) {$num_reads = $opt{'N'}};

}
1;