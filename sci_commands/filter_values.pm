package sci_commands::filter_values;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("filter_values");

sub filter_values {

@ARGV = @_;

$cellID_col = 1;

getopts("O:A:a:V:L:C:v", \%opt);

$die2 = "
scitools filter-values [options] [file with a column containing cellID]

Options:
   -O   [STR]   Output file name (default is [input file].filt.txt)
   -C   [INT]   Column of cellID (def = $cellID_col)
   -A   [STR]   Annot file
   -a   [STR]   Comma separated list of annots to include
   -L   [STR]   List of cellIDs to include (or other file with 1st col as cellIDs)
   -V   [STR]   Values to filter on, comma separated
                Format: [column]:[min value]:[max value],[coulmn2]:[min value]:[max value],etc...
                Min or max values can be '-inf' or 'inf'
   -Q   [STR]   Qualitative values to filter, comma separated
                Format: [column]:[included value1]:[included value2]:etc,[column2]:[included value1]:etc...
   -v           Only include the inverse of specified filters

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.txt$//; $opt{'O'} .= ".filt.txt"};
if (defined $opt{'C'}) {$cellID_col = $opt{'C'}};
if (defined $opt{'a'} && !defined $opt{'A'}) {die "\nMust provide an annotation file (-A) if specifying annotations to plot (-a)!\n$die2"};
if (defined $opt{'A'}) {read_annot($opt{'A'})};
if (defined $opt{'a'}) {
	@ANNOT_LIST = split(/,/, $opt{'a'});
	foreach $annot (@ANNOT_LIST) {
		$ANNOT_include{$annot} = 1;
	}
}
if (defined $opt{'L'}) {
	open IN, "$opt{'L'}";
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\s/, $l);
		$CELLID_include{$P[0]} = 1;
	} close IN;
}
if (defined $opt{'V'}) {
	@V = split(/,/, $opt{'V'});
	foreach $val_filt (@V) {
		($col,$min,$max) = split(/:/, $val_filt);
		if (defined $VAL_col_min{$col}) {die "ERROR: Cannot provide two filters for the same column! (problem column = $col)\n"};
		$VAL_col_min{$col} = $min;
		$VAL_col_min{$col} = $max;
	}
}
if (defined $opt{'Q'}) {
	@Q = split(/,/, $opt{'Q'});
	foreach $qual_filt (@Q) {
		@R = split(/:/, $qual_filt);
		if (defined $QUAL_col_include{$R[0]}) {die "ERROR: Cannot provide two filters for the same column! (problem column = $R[0])\n"};
		for ($i = 1; $i < @R; $i++) {
			$QUAL_col_include{$R[0]}{$R[$i]} = 1;
		}
	}
}

open OUT, ">$opt{'O'}";
open IN, "$ARGV[0]";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$pass = 1;
	$cellID = $P[$cellID_col-1];
	if (defined $opt{'A'}) {
		if (!defined $CELLID_annot{$cellID}) {
			$pass = 0;
		} elsif (defined $opt{'a'}) {
			if (!defined $ANNOT_include{$CELLID_annot{$cellID}}) {
				$pass = 0;
			}
		}
	}
	if (defined $opt{'L'} && !defined $CELLID_include{$cellID}) {
		$pass = 0;
	}
	if (defined $opt{'V'}) {
		foreach $col (keys %VAL_col_min) {
			if ($P[$col-1] < $VAL_col_min{$col} && $VAL_col_min{$col} !~ /-inf/) {$pass = 0};
			if ($P[$col-1] > $VAL_col_max{$col} && $VAL_col_max{$col} !~ /inf/) {$pass = 0};
		}
	}
	if (defined $opt{'Q'}) {
		foreach $col (keys %QUAL_col_include) {
			if (!defined $QUAL_col_include{$col}{$P[$col-1]}) {$pass = 0};
		}
	}
	if (($pass > 0 && !defined $opt{'v'}) || (defined $opt{'v'} && $pass < 1)) {
		print OUT "$l\n";
	}
} close IN; close OUT;

}
1;