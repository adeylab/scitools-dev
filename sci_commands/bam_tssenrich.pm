package sci_commands::bam_tssenrich;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("bam_tssenrich");

sub bam_tssenrich {

    @ARGV = @_;

# Defaults:
    $bedtools = "bedtools";

    getopts("O:T:B:b:X", \%opt);

$die2 = "
scitools bam_tssenrich [options] [duplicate removed and filtered bam file] [genome -- either 'hg38' or 'mm10']

Options:
   -O   [STR]   Output prefix (default is bam file prefix)
   -T   [STR]   An alternative BED file for TSS signal
   -B   [STR]   An alternative BED file for background signal
   -b   [STR]   Bedtools call (def = $bedtools)
   -X           Retain intermediate files, this includes individual value files for the TSS and background (def = remove)

";

    if (!defined $ARGV[0]) {die $die2};
    if (!defined $ARGV[1]) {die $die2};
    if ($ARGV[1] eq "hg38") {
	$tss_signal = "/home/groups/oroaklab/refs/hg38/ensembl_tss/ensembl.hg38.tss.chr100bpWINDOW.bed";
	$bg_signal = "/home/groups/oroaklab/refs/hg38/ensembl_tss/ensembl.hg38.tss.chr.b1_b2.bed";
    } elsif ($ARGV[1] eq "mm10") {
	$tss_signal = "/home/groups/oroaklab/refs/mm10/ensembl_tss/ensembl.mm10.tss.chr100bpWINDOW.bed";
	$bg_signal = "/home/groups/oroaklab/refs/mm10/ensembl_tss/ensembl.mm10.tss.chr.b1_b2.bed";
    } else {
	die $die2;
    }
    if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.bam$//};
    if (defined $opt{'b'}) {$bedtools = $opt{'b'}};
    if (defined $opt{'T'}) {$tss_signal = $opt{'T'}};
    if (defined $opt{'B'}) {$bg_signal = $opt{'B'}};

    system("$bedtools intersect -bed -u -a $ARGV[0] -b $tss_signal | cut -f 4 | sed -e 's/:.*//g' | sort | uniq -c | sed -e 's/^[ ]*//g' | awk '{print \$2,\$1}' | tr ' ' '\t' > $opt{'O'}.tss_reads.value");
    system("$bedtools intersect -bed -u -a $ARGV[0] -b $bg_signal | cut -f 4 | sed -e 's/:.*//g' | sort | uniq -c | sed -e 's/^[ ]*//g' | awk '{print \$2,\$1}' | tr ' ' '\t' > $opt{'O'}.bg_reads.value");
    system("join -e \"0\" -a1 -a2 -o \"0,1.2,2.2\" $opt{'O'}.tss_reads.value $opt{'O'}.bg_reads.value | awk '{print \$1,\$2/(\$3+1)}' | tr ' ' '\t' > $opt{'O'}.TSSenrich.value");
    system("cat $opt{'O'}.TSSenrich.value | awk '{sum+=\$2}END{print \"Average TSS enrichment: \" sum/NR}' > $opt{'O'}.TSSenrich.log");

    if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.tss_reads.value $opt{'O'}.bg_reads.value");
    }

}

1;
