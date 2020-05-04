package sci_commands::matrix_makeMM;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_makeMM");

sub matrix_makeMM {

    @ARGV = @_;
    getopts("O:R:L:W:", \%opt);

    $thrP = 0.975;
    $cis_alpha = 50;
    $cis_beta = 0.1;

$die2 = "
scitools matrix-makeMM [options] [input matrix]

    Options:
   -O   [STR]  Output prefix (default is [input].cistopic.dims)
   -L   [STR]  Columns file for sparseMatrix (will try to auto-detect)
   -W   [STR]  Rows file for sparseMatrix (will try to auto-detect)
   -R   [STR]  Rscript call (def = $Rscript)

";

###

    if (!defined $ARGV[0]) {die $die2};

    $prefix = $ARGV[0];
    $prefix =~ s/\.gz$//;
    $prefix =~ s/\.matrix$//;
    $prefix =~ s/\.values$//;
    $prefix =~ s/\.sparseMatrix$//;
    if (!defined $opt{'O'}) {
	$opt{'O'} = $prefix;
    }

    if ($ARGV[0] =~ /sparseMatrix/i) {
	$sparse = 1;
	if (defined $opt{'L'}) {
	    $col_file = $opt{'L'};
	} else {
	    if (-e "$prefix.sparseMatrix.cols") {
		$col_file = "$prefix.sparseMatrix.cols";
	    } elsif (-e "$prefix.sparseMatrix.cols.gz") {
		$col_file = "$prefix.sparseMatrix.cols.gz";
	    } else {
		die "ERROR: Cannot detect cols file (e.g. $prefix.sparseMatrix.cols), please provide as -C\n";
	    }
	}
	if (defined $opt{'W'}) {
	    $row_file = $opt{'W'};
	} else {
	    if (-e "$prefix.sparseMatrix.rows") {
		$row_file = "$prefix.sparseMatrix.rows";
	    } elsif (-e "$prefix.sparseMatrix.rows.gz") {
		$row_file = "$prefix.sparseMatrix.rows.gz";
	    } else {
		die "ERROR: Cannot detect rows file (e.g. $prefix.sparseMatrix.rows), please provide as -C\n";
	    }
	}
    } else {die "Not a sparse matrix."};

    if (defined $opt{'R'}) {$Rscript = $opt{'R'}};

###

    open R, ">$opt{'O'}.makeMM.r";

    print R "
library(Matrix)
IN<-as.matrix(read.table(\"$ARGV[0]\"))
IN<-sparseMatrix(i=IN[,1],j=IN[,2],x=IN[,3])
COLS<-read.table(\"$col_file\")
colnames(IN)<-COLS\$V1
ROWS<-read.table(\"$row_file\")
row.names(IN)<-ROWS\$V1
writeMM(IN,\"$opt{'O'}.mtx\")\n
";
    close R;

    system("$Rscript $opt{'O'}.makeMM.r");
    system("rm -f $opt{'O'}.makeMM.r");

}
1;
