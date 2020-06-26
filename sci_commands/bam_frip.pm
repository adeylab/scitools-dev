package sci_commands::bam_frip;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("bam_frip");

sub bam_frip {

    @ARGV = @_;

# Defaults:
    $bedtools = "bedtools";

    getopts("O:b:X", \%opt);

$die2 = "
scitools bam_frip [options] [duplicate removed and filtered bam file] [bed file]

Options:
   -O   [STR]   Output prefix (default is bam file prefix)
   -b   [STR]   Bedtools call (def = $bedtools)
   -X           Retain intermediate files, this includes individual value files for the TSS and background (def = remove)

";

    if (!defined $ARGV[0]) {die $die2};
    if (!defined $ARGV[1]) {die $die2};
    if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.bam$//};
    if (defined $opt{'b'}) {$bedtools = $opt{'b'}};

    system("$bedtools intersect -bed -u -a $ARGV[0] -b $ARGV[1] | cut -f 4 | sed -e 's/:.*//g' | sort | uniq -c | sed -e 's/^[ ]*//g' | awk '{print \$2,\$1}' | tr ' ' '\t' > $opt{'O'}.num_reads_in_bed.value");
    system("samtools view $ARGV[0] | cut -f 1 | sed -e 's/:.*//g' | sort | uniq -c | sed -e 's/^[ ]*//g' | awk '{print \$2,\$1}' | tr ' ' '\t' > $opt{'O'}.num_reads_in_bam.value");
    system("join -e \"0\" -a1 -a2 -o \"0,1.2,2.2\" $opt{'O'}.num_reads_in_bed.value $opt{'O'}.num_reads_in_bam.value | awk '{print \$1,\$2/\$3}' | tr ' ' '\t' > $opt{'O'}.frip.value");

    if (!defined $opt{'X'}) {
        system("rm -f $opt{'O'}.xoxo");
    }
    
}

1;
