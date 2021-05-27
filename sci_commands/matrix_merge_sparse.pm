package sci_commands::matrix_merge_sparse;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_merge_sparse");

sub matrix_merge_sparse {

    @ARGV = @_;
    getopts("O:z", \%opt);

$die2 = "

scitools matrix-merge-sparse (options) [matrix1.sparseMatrix.values.gz] [matrix2.sparseMatrix.values.gz]

Options:
   -O   [STR]   Output prefix (def = combination of two prefixes)
   -z           Turn off gzip output (default: output is gzipped)

";

    if (!defined $ARGV[0]) {die $die2};
    if (!defined $ARGV[1]) {die $die2};
    my $name1 = $ARGV[0];
    my $name2 = $ARGV[1];
    $name1 =~ s/\.sparseMatrix\.values\.gz//;
    $name2 =~ s/\.sparseMatrix\.values\.gz//;
    if (!defined $opt{'O'}) { $opt{'O'} = $name1 . "_" . $name2; }

# infer file names from input
    my $mat1_row = $name1 . ".sparseMatrix.rows.gz";
    my $mat1_col = $name1 . ".sparseMatrix.cols.gz";
    my $mat1_val = $name1 . ".sparseMatrix.values.gz";
    my $mat2_row = $name2 . ".sparseMatrix.rows.gz";
    my $mat2_col = $name2 . ".sparseMatrix.cols.gz";
    my $mat2_val = $name2 . ".sparseMatrix.values.gz";

# initalize counts and dictionaries
    my %row1; my %row2; my %mergedrow;
    my %col1; my %col2; my %mergedcol;
    my %val;
    my $row1_count = 1; my $row2_count = 1; my $mergedrow_count = 1;
    my $col1_count = 1; my $col2_count = 1; my $mergedcol_count = 1;

# read matrix1 rows
    open(IN_ROW, "zcat $mat1_row |") or die "gzip $mat1_row: $!";
    while(<IN_ROW>) {
	chomp $_;
	$row1{$row1_count} = $_;
	$row1_count++;
	if (! defined $mergedrow{$_}) {
	    $mergedrow{$_} = 0;
	}
    }
    close IN_ROW;

# read matrix2 rows
    open(IN_ROW, "zcat $mat2_row |") or die "gzip $mat2_row: $!";
    while(<IN_ROW>) {
	chomp $_;
	$row2{$row2_count} = $_;
	$row2_count++;
	if (! defined$mergedrow{$_}) {
	    $mergedrow{$_} = 0;
	}
    }
    close IN_ROW;

# read matrix1 cols
    open(IN_COL, "zcat $mat1_col |") or die "gzip $mat1_col: $!";
    while(<IN_COL>) {
	chomp $_;
	$col1{$col1_count} = $_;
	$col1_count++;
	if (! defined$mergedcol{$_}) {
	    $mergedcol{$_} = 0;
	}
    }
    close IN_COL;

# read matrix2 cols
    open(IN_COL, "zcat $mat2_col |") or die "gzip $mat2_col: $!";
    while(<IN_COL>) {
	chomp $_;
	$col2{$col2_count} = $_;
	$col2_count++;
	if (! defined$mergedcol{$_}) {
	    $mergedcol{$_} = 0;
	}
    }
    close IN_COL;

# write new rows
    if (! defined $opt{'z'}) {
	open OUTROW, "| $gzip > $opt{'O'}.rows.gz";
    } else {
	open OUTROW, ">$opt{'O'}.rows";
    }
    foreach my $key (sort {$a cmp $b} keys(%mergedrow)) {
	$mergedrow{$key} = $mergedrow_count;
	$mergedrow_count++;
	print OUTROW "$key\n";
    }
    close OUTROW;

# write new cols
    if (! defined $opt{'z'}) {
	open OUTCOL, "| $gzip > $opt{'O'}.cols.gz";
    } else {
	open OUTCOL, ">$opt{'O'}.cols";
    }
    foreach my $key (sort {$a cmp $b} keys(%mergedcol)) {
	$mergedcol{$key} = $mergedcol_count;
	$mergedcol_count++;
	print OUTCOL "$key\n";
    }
    close OUTCOL;

# read values1 and translate
    open(IN_VAL, "zcat $mat1_val |") or die "gzip $mat1_val: $!";
    while(<IN_VAL>) {
	chomp $_;
	my ($a, $b, $c) = split("\t",$_);
	my $trans_a = $mergedrow{$row1{$a}};
	my $trans_b = $mergedcol{$col1{$b}};
	my $coords = $trans_a . "\t" . $trans_b;
	if (defined $val{$coords}) {
	    $val{$coords} = $val{$coords} + $c;
	} else {
	    $val{$coords} = $c;
	}
    }

# read values2 and translate
    open(IN_VAL, "zcat $mat2_val |") or die "gzip $mat2_val: $!";
    while(<IN_VAL>) {
	chomp $_;
	my ($a, $b, $c) = split("\t",$_);
	my $trans_a = $mergedrow{$row2{$a}};
	my $trans_b= $mergedcol{$col2{$b}};
	my $coords = $trans_a . "\t" . $trans_b;
	if (defined$val{$coords}) {
	    $val{$coords} =$val{$coords} +$c;
	} else {
	    $val{$coords} =$c;
	}
    }

# output values
    if (! defined $opt{'z'}) {
	open OUTVAL, "| $gzip > $opt{'O'}.values.gz";
    } else {
	open OUTVAL, ">$opt{'O'}.values";
    }
    foreach my $key (keys %val) {
	print OUTVAL "$key\t$val{$key}\n";
    }
    close OUTVAL;

}
1;
