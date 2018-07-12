package sci_commands::atac_distal;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("atac_distal");

sub atac_distal {

@ARGV = @_;

# Defaults
$dist = 20000;
getopts("O:", \%opt);

$die2 = "
scitools atac-distal [options] [input matrix/bed]

Filters a matrix or bed file to only include sites that are
considered distal based on a threshold distance from the
TSS of genes.

Options: Either -B OR -G must be specified.
   -O   [STR]   Output prefix (default is input prefix)
   -B   [BED]   Bed file of TSSs / gene bodies
                 chr start end (name) (strand)
                 if no strand will assume start as TSS
   -G   [STR]   Gene info (refGene.txt formats)
                 Shortcut eg: hg38, hg19, mm10
   -D   [BP]    Minimum distance (in bp) from a TSS for a
                 site to be considered distal (def = $dist)

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//; $opt{'O'} =~ s/\.bed$//;};
if (defined $opt{'D'}) {$dist = $opt{'D'}};
if (!defined $opt{'G'} && !defined $opt{'B'}) {die "ERROR: A TSS file MUST be specified! -G or -B\n$die2"};

if (defined $opt{'G'}) {
	if (defined $REF{$opt{'G'}}) {
		$ref_file = $REF{$opt{'G'}};
		$opt{'G'} = $ref_file;
		$opt{'G'} =~ s/\.fa$/\.refGene.txt/;
	}
	read_refgene($opt{'G'});
}

##############################


}
1;
