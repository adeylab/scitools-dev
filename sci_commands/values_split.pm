package sci_commands::values_split;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("values_split");

sub values_split {

@ARGV = @_;

getopts("A:O:c:", \%opt);

$die2 = "
scitools values-split [options] -A [annot file] [values file]
   or    values-fastq

Will split your values or other files by annotation.

Options:
   -A   [STR]   Annotation file (comma separated for more than one)
                (If multiple, must be non-conflicting)
				(Required)
   -O   [STR]   Output prefix (def = values file prefix)
   -c   [INT]   Column for cellID in values file (def = 1)

";

if (!defined $ARGV[0] || !defined $opt{'A'}) {die $die2};
if (!defined $opt{'O'}) {
	$opt{'O'} = $ARGV[0];
	$opt{'O'} =~ s/\.vals$//;
	$opt{'O'} =~ s/\.values$//;
}
if (!defined $opt{'c'}) {$opt{'c'} = 1};

read_annot($opt{'A'});
print STDERR "$annot_count total annotations found.\n";

# setup output files
foreach $annot (keys %ANNOT_count) {
	$out_handle = "$annot";
	open $out_handle, "> $opt{'O'}.$annot.values";
	$HANDLE{$annot} = $out_handle;
}

open IN, $ARGV[0];
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$cellID = $P[$opt{'c'}-1];
	if (defined $CELLID_annot{$cellID}) {
		$annot = $CELLID_annot{$cellID};
		$ANNOT_count{$annot}++;
		$out_handle = $HANDLE{$annot};
		print $out_handle "$l\n";
	}
} close IN;

foreach $annot (keys %ANNOT_count) {
	$out_handle = $HANDLE{$annot};
	close $out_handle;
}

}
1;
