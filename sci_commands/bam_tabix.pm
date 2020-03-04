package sci_commands::bam_tabix;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("bam_tabix");

sub bam_tabix {

    @ARGV = @_;

# Defaults:
    $threads = 1;

    getopts("s:t:O", \%opt);

$die2 = "
scitools bam-tabix [options] [duplicate removed and filtered bam file]

Requires tabix and bgzip from htslib.

Options:
   -O   [STR]   Output prefix (default is bam file prefix)
   -t   [INT]   number of threads to use
   -s   [STR]   Samtools call (def = $samtools)
";

    if (!defined $ARGV[0]) {die $die2};
    if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.bam$//};
    if (defined $opt{'t'}) {$threads = $opt{'t'}};
    if (defined $opt{'s'}) {$samtools = $opt{'s'}};

    system("samtools view --threads $threads $ARGV[0] | awk 'OFS=\"\\t\" {split(\$1,a,\":\"); print \$3,\$4,\$8,a[1],1}' | sort -S 2G -T . --parallel=$threads -k1,1 -k2,2n -k3,3n | bgzip > $opt{'O'}.fragments.tsv.gz");
    system("tabix -p bed $opt{'O'}.fragments.tsv.gz");
}
1;
