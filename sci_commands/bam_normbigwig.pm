package sci_commands::bam_normbigwig;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("bam_normbigwig");

sub bam_normbigwig {

    @ARGV = @_;

# Defaults:
    $scale = 1000000;

    getopts("s:b:O", \%opt);

$die2 = "
scitools bam-normbigwig [options] [duplicate removed and filtered bam file] [genome]

Requires bedtools, samtools, and ucsc-kent-toolkit.
Genome must be: hg19, hg38, mm10, dm6.

Options:
   -O   [STR]   Output prefix (default is bam file prefix)
   -S   [INT]   scaling factor to use, default 1,000,000
   -s   [STR]   Samtools call (def = $samtools)
   -b   [STR]   Bedtools call (def = $bedtools)
";

    if (!defined $ARGV[0]) {die $die2};
    if (!defined $ARGV[1]) {
	die $die2;
    } else {
	if ($ARGV[1] eq "hg19") {
	    $genomefile = "/home/groups/oroaklab/refs/hg19/hg19.chromInfo";
	} elsif ($ARGV[1] eq "hg38") {
	    $genomefile = "/home/groups/oroaklab/refs/hg38/hg38.chromInfo";
	} elsif ($ARGV[1] eq "mm10") {
	    $genomefile = "/home/groups/oroaklab/refs/mm10/mm10.chromInfo";
	} elsif ($ARGV[1] eq "dm6") {
	    $genomefile = "/home/groups/oroaklab/refs/dm6/dm6.chromInfo";
	} else {
	    die $die2;
	}
    }
    if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.bam$//};
    if (defined $opt{'S'}) {$scale = $opt{'S'}};
    if (defined $opt{'s'}) {$samtools = $opt{'s'}};
    if (defined $opt{'b'}) {$bedtools = $opt{'b'}};

    $bamcount = qx(samtools view -c $ARGV[0]);
    $scalingfactor = $scale/$bamcount;

    system("genomeCoverageBed -ibam $ARGV[0] -bg -scale $scalingfactor -g $genomefile > temp.$opt{'O'}.norm.bedgraph");
    system("sort-bed temp.$opt{'O'}.norm.bedgraph > temp.$opt{'O'}.norm.sorted.bedgraph");
    system("bedGraphToBigWig temp.$opt{'O'}.norm.sorted.bedgraph $genomefile $opt{'O'}.norm.bigwig");
    system("rm temp.$opt{'O'}.norm.bedgraph");
    system("rm temp.$opt{'O'}.norm.sorted.bedgraph");

}
1;
