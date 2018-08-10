package sci_commands::matrix_match;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_match");

sub matrix_match {

$| = 1;
@ARGV = @_;

# Defaults
getopts("O:V", \%opt);

$die2 = "
scitools matrix-match [options] [counts matrix 1] [counts matrix 2]
   or    match-matrix

Calculates the distance between counts matrix cells and makes an
annotation file for the top hit match based on the shared nonzero
sites.

Note: rows must be identical

Options:
   -O   [STR]   Output prefix (default is [mat1_vs_mat2])
   -V           Verbose (prints a log file)

";

if (!defined $ARGV[1]) {die $die2};
if (!defined $opt{'O'}) {
	$opt{'O'} = $ARGV[0];
	$opt{'O'} =~ s/\.matrix$//;
	$opt{'O'} .= "_vs_".$ARGV[1];
	$opt{'O'} =~ s/\.matrix$//;
}

if (defined $opt{'V'}) {
	open LOG, ">$opt{'O'}.matching.log";
	$ts = localtime(time);
	print LOG "$ts matrix-match called
\tMatrix 1 = $ARGV[0]
\tMatrix 2 = $ARGV[1]
\tChecking row count ... ";
	# get line count for progress meter:
	open M1, "$ARGV[0]";
	while ($l = <M1>) {$lineCT++};
	close M1;
	$lineCT--; # header line
	print LOG "$lineCT rows.\n";
	$ts = localtime(time);
	print LOG "$ts Building empty matrix ...\n";
}

open M1, "$ARGV[0]";
$h1 = <M1>; chomp $h1; @H1 = split(/\t/, $h1);
open M2, "$ARGV[1]";
$h2 = <M2>; chomp $h2; @H2 = split(/\t/, $h2);

for ($m1 = 0; $m1 < @H1; $m1++) {
	for ($m2 = 0; $m2 < @H2; $m2++) {
		$M1_M2_nonzero{$m1}{$m2} = 0;
	}
}

if (defined $opt{'V'}) {
	$ts = localtime(time);
	print LOG "$ts Starting matching ...\n";
	$report_increment = 0.1; $report = $report_increment;
}

$lineID = 0;
while ($l1 = <M1>) {
	chomp $l1; $l2 = <M2>; chomp $l2;
	@P1 = split(/\t/, $l1); @P2 = split(/\t/, $l2);
	$siteID1 = shift(@P1); $siteID2 = shift(@P2);
	$lineID++;
	if ($siteID1 ne $siteID2) {
		die "ERROR: Rows must be identical between the two matrixes! $siteID1 ne $siteID2! (line $lineID)\n";
	}
	for ($m1 = 0; $m1 < @H1; $m1++) {
		for ($m2 = 0; $m2 < @H2; $m2++) {
			if ($P1[$m1] != 0 && $P2[$m2] != 0) {
				$M1_M2_nonzero{$m1}{$m2}++;
			}
		}
	}
	if (defined $opt{'V'}) {
		if (($lineID/$lineCT)>=$report) {
			$ts = localtime(time);
			print LOG "\t$report fraction complete, $lineID lines processed. ($ts)\n";
			$report += $report_increment;
		}
	}
}
close M1; close M2;

if (defined $opt{'V'}) {
	$ts = localtime(time);
	print LOG "$ts Matching complete - reporting best hits and matching matrix.\n";
	$report_increment = 0.1; $report = $report_increment;
}

open D, ">$opt{'O'}.all_dist.txt";
open A, ">$opt{'O'}.best_match.annot";
open R, ">$opt{'O'}.best_match_rev.annot";

for ($m1 = 0; $m1 < @H1; $m1++) {
	$cellID1 = $H1[$m1];
	print D "$cellID1";
	$win = "null";
	# sort by the one with the most shared non-zero sites
	foreach $m2 (sort {$M1_M2_nonzero{$m1}{$b}<=>$M1_M2_nonzero{$m1}{$a}} keys %{$M1_M2_nonzero{$m1}}) {
		if ($win eq "null") {$win = $H2[$m2]; print D "\t$win"};
		$av_d = $M1_M2_nonzero{$m1}{$m2};
		print D "\t$av_d";
	}
	print D "\n";
	print A "$cellID1\t$win\n";
	print R "$win\t$cellID1\n";
	if (defined $opt{'V'}) {
		if (($m1/@H1)>=$report) {
			$ts = localtime(time);
			print LOG "\t$report fraction complete, $m1 cellIDs from Matrix 1 processed. ($ts)\n";
			$report += $report_increment;
		}
	}
}

close D; close A; close R;

if (defined $opt{'V'}) {
	$ts = localtime(time);
	print LOG "$ts Complete!\n";
	close LOG;
}

}
1;