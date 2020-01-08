package sci_commands::fastq_align;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("fastq_align");

sub fastq_align {

@ARGV = @_;

# Defaults
$threads = 1;
$memory = "2G";
$sort_threads = 1;

getopts("t:b:s:A:L:O:m:nr:Sp:", \%opt);

$die2 = "
scitools fastq-align [options] [aligner reference] [output_prefix] [read1.fq] (read2.fq)
   or    align-fastq
   or    align

Produces a sorted bam file. Read 2 is optional.

Options:
   -t   [INT]   Threads for alignment (def = $threads)
   -r   [INT]   Threads for sorting (def = $sort_threads)
   -b   [STR]   Bwa call (def = $bwa)
   -s   [STR]   Samtools call (def = $samtools)
   -S           Instead use the snap-aligner
                (reference shortcuts must have a snap folder with
                 ref files in the same directory as the fasta)
   -p   [STR]   Snap-aligner call (def = $snap_aligner)
   -n           Sort reads by name / index (def = coordinate)
   -m   [MEM]   Samtools sort mex memory per thread, K/M/G (def = $memory)

Reference shortcuts:
$ref_shortcuts

";

if (!defined $ARGV[2]) {die $die2};
if (defined $REF{$ARGV[0]}) {
	$ref_file = $REF{$ARGV[0]};
} else {
	$ref_file = $ARGV[0];
}

if (defined $opt{'S'}) {
	$ref_file =~ s/[^\/]+$/snap/;
	die "DEBUG: $ref_file is ref file. snap-aligner is $snap_aligner\n";
}

if (defined $opt{'O'}) {
	$out_prefix = $opt{'O'};
} else {
	$out_prefix = $ARGV[1];	
}
$out_prefix =~ s/\.bam$//;
$out_prefix =~ s/\.$//;

if (defined $opt{'b'}) {$bwa = $opt{'b'}};
if (defined $opt{'s'}) {$samtools = $opt{'s'}};
if (!defined $opt{'t'}) {$opt{'t'} = $threads};
if (!defined $opt{'m'}) {$opt{'m'} = $memory};
if (defined $opt{'p'}) {$snap_aligner = $opt{'p'}};
if (defined $opt{'r'}) {$sort_threads = $opt{'r'}};

if (defined $ARGV[3]) {
	if (defined $opt{'n'}) {
		$align_command = "$bwa mem -t $opt{'t'} $ref_file $ARGV[2] $ARGV[3] 2>> $out_prefix.align.log | $samtools view -bSu - 2>> $out_prefix.align.log | $samtools sort -@ $sort_threads -m $opt{'m'} -T $out_prefix.TMP -n - > $out_prefix.nsrt.bam 2>> $out_prefix.align.log";
	} else {
		$align_command = "$bwa mem -t $opt{'t'} $ref_file $ARGV[2] $ARGV[3] 2>> $out_prefix.align.log | $samtools view -bSu - 2>> $out_prefix.align.log | $samtools sort -@ $sort_threads -m $opt{'m'} -T $out_prefix.TMP - > $out_prefix.bam 2>> $out_prefix.align.log";
	}
} else { # single ended
	if (defined $opt{'n'}) {
		$align_command = "$bwa mem -t $opt{'t'} $ref_file $ARGV[2] 2>> $out_prefix.align.log | $samtools view -bSu - 2>> $out_prefix.align.log | $samtools sort -@ $sort_threads -m $opt{'m'} -T $out_prefix.TMP -n - > $out_prefix.nsrt.bam 2>> $out_prefix.align.log";
	} else {
		$align_command = "$bwa mem -t $opt{'t'} $ref_file $ARGV[2] 2>> $out_prefix.align.log | $samtools view -bSu - 2>> $out_prefix.align.log | $samtools sort -@ $sort_threads -m $opt{'m'} -T $out_prefix.TMP - > $out_prefix.bam 2>> $out_prefix.align.log";
	}
}

$ts = localtime(time);
open LOG, ">$out_prefix.align.log";
print LOG "$ts Alignment started.
Alignment command:
$align_command

"; close LOG;

#print "Running: $align_command\n";
system($align_command);

}
1;
