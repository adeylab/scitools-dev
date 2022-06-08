package sci_commands::triplet2ascii;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("triplet2ascii");

sub triplet2ascii {

@ARGV = @_;

# defaults

getopts("D:c:", \%opt);

$die2 = "
scitools triplet2ascii (options) [input_file] [output_file]

Collapses cell barcodes to ascii code. Will try to auto-detect file type.
Will only output a file of the same type. (e.g. bam->bam, fq.gz->fq.gz)

Options:
   -D   [STR]   Print all detected barcodes to provided annotation file
   -c   [INT]   1-based column of a flatfile for the cellID

";

if (!defined $ARGV[1]) {die $die2};

# load key
load_triplet2ascii();

# determine input file type & process it

if ($ARGV[0] =~ /bam$/) { # bam
	open IN, "$samtools view -h $ARGV[0] |";
	open OUT, "| $samtools view -b - > $ARGV[1]";
	while ($l = <IN>) {
		chomp $l;
		if ($l =~ /^\@/) {
			print OUT "$l\n";
		} else {
			@P = split(/\t/, $l);
			@B = split(/:/, $P[0]);
			$collapsed = collapse_barcode($B[0]);
			$BARCODE_ascii{$B[0]} = $collapsed;
			$B[0] = $collapsed;
			$P[0] = join(":", @B);
			$l = join("\t", @P);
			print OUT "\n";
		}
	} close IN; close OUT;
} elsif ($ARGV[0] =~ /fq.gz$/ || $ARGV[0] =~ /fastq.gz$/) { # fastq
	open IN, "$zcat $ARGV[0] |";
	open OUT, "| $gzip > $ARGV[1]";
	while ($tag = <IN>) {
		chomp $tag;
		@B = split(/:/, $tag); $B[0] =~ s/^\@//;
		$collapsed = collapse_barcode($B[0]);
		$BARCODE_ascii{$B[0]} = $collapsed;
		$B[0] = $collapsed;
		$tag = join(":", @B);
		print OUT "\@$tag\n";
		$seq = <IN>; print OUT "$seq";
		$null = <IN>; print OUT "$null";
		$qual = <IN>; print OUT "$qual";
	} close IN; close OUT;
} else { # assume a text file, assume 1st column is index unless specified.
	if (!defined $opt{'c'}) {$opt{'c'} = 1};
	$col = $opt{'c'};
	if ($ARGV[0] =~ /complexity.txt$/) {$col = 2};
	$col--;
	open IN, "$ARGV[0]";
	open OUT, ">$ARGV[1]";
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		$collapsed = collapse_barcode($P[$col]);
		$BARCODE_ascii{$P[$col]} = $collapsed;
		$P[$col] = $collapsed;
		$l = join("\t", @P);
		print OUT "$l\n";
	} close IN; close OUT;
}

if (defined $opt{'D'}) {
	open KEY, ">$opt{'D'}";
	foreach $barc (keys %BARCODE_ascii) {
		print KEY "$barc\t$BARCODE_ascii{$barc}\n";
	} close KEY;
}

}
1;
