package sci_commands::rna_align;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("rna_align");

sub rna_align {

@ARGV = @_;

# Defaults
$threads = 1;
$memory = "2G";
$sort_threads = 1;

getopts("t:S:s:O:m:r:p", \%opt);

$die2 = "
scitools rna-align [options] [SNAP genomeDir or shortcut] [output_prefix] [rna_read.fq.gz]
   or    align-rna

or scitools rna-align [options] -p [SNAP.sam Output]; just process after alignment

Takes in the single fastq file in sci format and produces and aligned bam with
UMI and cell barcode fields.

If genomeDir shortcut is used, it must be the path to the fastq file, as with
bwa aligner. The fastq file will be removed and the directory used.

Options:
   -t   [INT]   Threads for alignment (def = $threads)
   -r   [INT]   Threads for sorting (def = $sort_threads)
   -S   [STR]   STAR call (def = $STAR)
   -s   [STR]   Samtools call (def = $samtools)
   -m   [MEM]   Samtools sort mex memory per thread, K/M/G (def = $memory)

Reference shortcuts:
$ref_shortcuts

";

if (!defined $opt{'p'}) {
	if (!defined $ARGV[2]) {die $die2};
	if (defined $REF{$ARGV[0]}) {
		$genomeDir = $REF{$ARGV[0]};
	} else {
		$genomeDir = $ARGV[0];
	}
} else {
	if (!defined $ARGV[0]) {die $die2};
}

if (defined $opt{'O'}) {
	$out_prefix = $opt{'O'};
} else {
	$out_prefix = $ARGV[1];	
}
$out_prefix =~ s/\.bam$//;
$out_prefix =~ s/\.$//;

if (defined $opt{'S'}) {$STAR = $opt{'S'}};
if (defined $opt{'s'}) {$samtools = $opt{'s'}};
if (!defined $opt{'t'}) {$opt{'t'} = $threads};
if (!defined $opt{'m'}) {$opt{'m'} = $memory};
if (defined $opt{'r'}) {$sort_threads = $opt{'r'}};

$align_command = "$STAR mem --runThreadN $opt{'t'} --genomeDir $genomeDir --outFileNamePrefix $out_prefix\. --readFilesIn $ARGV[1] --readFilesCommand $zcat 2> $out_prefix.STAR.log";

if (!defined $opt{'p'}) {
	#print "Running: $align_command\n";
	system($align_command);
	open SAM, "$out_prefix.Aligned.sam";
} else {
	open SAM, "$ARGV[0]";
}
open OUT, "| $samtools sort -@ $sort_threads -m $opt{'m'} -T $out_prefix.TMP - 2>> $out_prefix.bam";
while ($l = <SAM>) {
	chomp $l;
	if ($l ~ /^\@/) {
		print OUT "$l\n";
	} else {
		@P = split(/\t/, $l) {
			$barc = $P[0]; $barc =~ s/:.+$//;
			$umi = $P[0]; $umi =~ s/^.+:UMI=//; $umi =~ s/\..+$//;
			$l .= "\tCB:Z:$barc\tUB:Z:$umi";
			print OUT "$l\n";
		}
	}
} close SAM; close OUT;

}
1;
