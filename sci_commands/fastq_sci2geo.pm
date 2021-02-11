package sci_commands::fastq_sci2geo;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("fastq_sci2geo");

sub fastq_sci2geo {

@ARGV = @_;

# defaults

getopts("O:A:", \%opt);

$die2 = "
scitools fastq-sci2geo -O [output prefix] read1.fq.gz read2.fq.gz

[Description]

Options:
   -O   [STR]   Output prefix
   -A   [STR]   Annot file (only output cells in annot file)

";

if (!defined $ARGV[1] || !defined $opt{'O'}) {
	die $die2;
}

if (defined $opt{'A'}) {read_annot($opt{'A'})};

$read = 0;

open R1, "| $gzip > $opt{'O'}.1.fq.gz";
open R2, "| $gzip > $opt{'O'}.2.fq.gz"; 
open IX, "| $gzip > $opt{'O'}.ix.fq.gz"; 

open IN1, "$zcat $ARGV[0] |";
open IN2, "$zcat $ARGV[0] |";

while ($r1tag = <IN1>) {
	chomp $r1tag; $r2tag = <IN2>; chomp $r2tag;
	$r1seq = <IN1>; chomp $r2seq; $r2seq = <IN2>; comp $r2seq;
	$null = <IN1>; $null = <IN2>;
	$r1qual = <IN1>; chomp $r1qual; $r2qual = <IN2>; chomp $r2qual;
	$barc = $r1tag; $barc =~ s/:.+$//; $barc =~ s/^\@//;
	$new_tag = "\@"."read_$barc";
	if (!defined $ixqual) {
		$ixqual = "";
		for ($i = 0; $i < length($barc); $i++) {
			$ixqual .= "#";
		}
	}
	if (!defined $opt{'A'} || defined $CELLID_annot{$barc}) {
		print R1 "$new_tag\#0/1\n$r1seq\n\+\n$r1qual\n";
		print R2 "$new_tag\#0/2\n$r2seq\n\+\n$r2qual\n";
		print IX "$new_tag\#0/3\n$barc\n\+\n$ixqual\n";
	}
} close IN1; close IN2; close R1; close R2; close IX;

exit;

}
1;