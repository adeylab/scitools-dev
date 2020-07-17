package sci_commands::matrix_tf;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_tf");

sub matrix_tf {

@ARGV = @_;

getopts("O:N:", \%opt);

$die2 = "
scitools matrix-tf [options] [counts matrix]
   or    tf

Term frequency normalization of matrix.

If a regular or a sparseMatrix is provided, the output will be in the same
format as the input.

Options:
   -O   [STR]   Output prefix (default is [input].tf)
   -N   [VAL]   Normalization constant (def = row number)

";

if (!defined $ARGV[0]) {die $die2};

if ($ARGV[0] =~ /\.matrix/ || $ARGV[0] =~ /\.matrix\.gz/) {

	if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};

	read_matrix_stats($ARGV[0]);
	if (!defined $opt{'N'}) {$opt{'N'} = $matrix_rowNum};

	open IN, "$ARGV[0]";
	open OUT, ">$opt{'O'}.tf";
	$h = <IN>; print OUT "$h";
	while ($l = <IN>) {
		chomp $l;
		@P = split (/\t/, $l);
		$rowID = shift(@P);
		print OUT "$rowID";
		for ($i = 0; $i < @P; $i++) {
			if ($P[$i]>0) {
				$tf = ($P[$i]/$MATRIX_CellID_signal{$MATRIX_COLNAMES[$i]});
				$score = sprintf("%.6f", $tf*$opt{'N'});
			} else {
				$score = "0";
			}
			print OUT "\t$score";
		}
		print OUT "\n";
	} close IN; close OUT;

} elsif ($ARGV[0] =~ /(sparseMatrix|values)/) {

	if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.gz$//; $opt{'O'} =~ s/\.values$//; $opt{'O'} =~ s/\.sparseMatrix$//};

	# open matrix
	if ($ARGV[0] =~ /\.gz$/) {
		open IN, "$zcat $ARGV[0] |";
	} else {
		open IN, "$ARGV[0]";
	}
	
	# read matrix once to get col counts
	while ($l = <IN>) {
		chomp $l;
		($row,$col,$val) = split(/\s/, $l);
		$COLID_sum{$col} += $val;
		$lastrow = $row;
	}
	close IN;

	# define normalization factor, default num rows
	if (!defined $opt{'N'}) {$opt{'N'} = $lastrow};

	# open matrix
	if ($ARGV[0] =~ /\.gz$/) {
		open IN, "$zcat $ARGV[0] |";
	} else {
		open IN, "$ARGV[0]";
	}
	
	# open output
	if (!defined $opt{'z'}) {
		open OUT, "| $gzip > $opt{'O'}.sparseMatrix.tf.gz";
	} else {
		open OUT, ">$opt{'O'}.sparseMatrix.tf";
	}

	print "test info\n";
	foreach my $key (keys %COLID_sum) {
		print "$key\t$COLID_sum{$key}\n";
	}
	print "$lastrow\n";
		
	# read matrix values and normalize
	while ($l = <IN>) {
		chomp $l;
		($row,$col,$val) = split(/\s/, $l);
		$tf = $val/$COLID_sum{$col};
		$score = sprintf("%.6f", $tf*$opt{'N'});
		print OUT "$row\t$col\t$score\n";
	}
	close IN; close OUT;

}

}
1;