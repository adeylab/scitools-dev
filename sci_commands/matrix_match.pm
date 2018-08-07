package sci_commands::matrix_match;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_match");

sub matrix_match {

@ARGV = @_;

# Defaults
getopts("O:", \%opt);

$die2 = "
scitools matrix-match [options] [counts matrix 1] [counts matrix 2]
   or    match-matrix

Calculates the distance between counts matrix cells and makes an
annotation file for the top hit match based on the shared nonzero sites.

Options:
   -O   [STR]   Output prefix (default is [mat1_vs_mat2])

";

if (!defined $ARGV[1]) {die $die2};
if (!defined $opt{'O'}) {
	$opt{'O'} = $ARGV[0];
	$opt{'O'} =~ s/\.matrix$//;
	$opt{'O'} .= "_vs_".$ARGV[1];
	$opt{'O'} =~ s/\.matrix$//;
}

open M1, "$ARGV[0]";
$h1 = <M1>; chomp $h1; @H1 = split(/\t/, $h1);
while ($l = <M1>) {
	chomp $l;
	@P = split(/\t/, $l);
	$siteID = shift(@P);
	@{$SITEID_row{$siteID}} = @P;
} close M1;

open M2, "$ARGV[1]";
$h2 = <M2>; chomp $h2; @H2 = split(/\t/, $h2);
for ($m1 = 0; $m1 < @H1; $m1++) {
	for ($m2 = 0; $m2 < @H2; $m2++) {
		$M1_M2_dist{$m1}{$m2} = 0;
		$M1_M2_nonzero{$m1}{$m2} = 0;
	}
}
while ($l = <M2>) {
	chomp $l;
	@P = split(/\t/, $l);
	$siteID = shift(@P);
	if (defined $SITEID_row{$siteID}[0]) {
		for ($m1 = 0; $m1 < @H1; $m1++) {
			for ($m2 = 0; $m2 < @H2; $m2++) {
				$M1_M2_dist{$m1}{$m2} += abs($P[$m2]-$SITEID_row{$siteID}[$m1]);
				if ($P[$m2] != 0 && $SITEID_row{$siteID}[$m1] != 0) {
					$M1_M2_nonzero{$m1}{$m2}++;
				}
			}
		}
		$siteCT++;
	}
} close M2;

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
}

close D; close A; close R;

}
1;