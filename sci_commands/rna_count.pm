package sci_commands::rna_count;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("rna_count");

sub rna_count {

@ARGV = @_;

getopts("s:O:m:r:", \%opt);

$die2 = "
scitools rna-count [options] [GTF file / shortcut] [rna-align sorted bam]
   or    count-rna

Options:
   -r   [INT]   Threads for sorting (def = $sort_threads)
   -s   [STR]   Samtools call (def = $samtools)
   -m   [MEM]   Samtools sort mex memory per thread, K/M/G (def = $memory)

Reference shortcuts:
$ref_shortcuts

";

if (defined $opt{'O'}) {
	$out_prefix = $opt{'O'};
} else {
	$out_prefix = $ARGV[1];	
}
$out_prefix =~ s/\.bam$//;
$out_prefix =~ s/\.$//;

if (defined $opt{'s'}) {$samtools = $opt{'s'}};
if (!defined $opt{'m'}) {$opt{'m'} = $memory};
if (defined $opt{'r'}) {$sort_threads = $opt{'r'}};


}
1;