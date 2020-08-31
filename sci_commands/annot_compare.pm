package sci_commands::annot_compare;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("annot_compare");

sub annot_compare {

@ARGV = @_;
getopts("O:S", \%opt);

$die2 = "
scitools annot-compare [options] [annot1] [annot2]
   or    compare-annot

Will compare the proportions in each of the annotaitons across files.
CellIDs must be the same (will take the intersect only)

Options:
   -O   [STR]   Output prefix (def = annot1.annot2)
                Adds .compare.txt suffix
   -S           Only print proportions of Annot2 within Annot1 (instead of all)
   -M           Print a matrix of the comparisons (as in annot-proportions)

";

if (!defined $ARGV[1]) {die $die2};

if (!defined $opt{'O'}) {
	$pfx1 = $ARGV[0]; $pfx1 =~ s/\.annot$//;
	$pfx2 = $ARGV[1]; $pfx2 =~ s/\.annot$//;
	$opt{'O'} = "$pfx1.$pfx2";
}

if (defined $opt{'M'}) { # matric (annot-proportions)

$out_header = "#"; @ANNOT1 = ();
open IN, "$ARGV[0]";
while ($l = <IN>) {
	chomp $l;
	($cellID,$annot) = split(/\t/, $l);
	$CELLID_annot1{$cellID} = $annot;
	if (!defined $ANNOT1_check{$annot}) {
		$out_header .= "\t$annot";
		push @ANNOT1, $annot;
		$ANNOT1_check{$annot} = 1;
	}
} close IN;
$cell_match_ct = 0; @ANNOT2 = ();
open IN, "$ARGV[1]";
while ($l = <IN>) {
	chomp $l;
	($cellID,$annot2) = split(/\t/, $l);
	if (!defined $ANNOT2_check{$annot2}) {
		push @ANNOT2, $annot2;
		$ANNOT2_check{$annot2} = 1;
	}
	if (defined $CELLID_annot1{$cellID}) {
		$annot1 = $CELLID_annot1{$cellID};
		$cell_match_ct++;
		$ANNOT1_ANNOT2_ct{$annot1}{$annot2}++;
	}
} close IN;
open OUT, ">$opt{'O'}.compare.txt";
print OUT "$out_header\n";
foreach $annot2 (@ANNOT2) {
	print OUT "$annot2";
	for ($i = 0; $i < @ANNOT1; $i++) {
		if (!defined $ANNOT1_ANNOT2_ct{$ANNOT1[$i]}{$annot2}) {$ANNOT1_ANNOT2_ct{$ANNOT1[$i]}{$annot2}=0};
		print OUT "\t$ANNOT1_ANNOT2_ct{$ANNOT1[$i]}{$annot2}";
	}
	print OUT "\n";
}
close OUT;

} else { # regular

$opt{'O'} =~ s/\.txt$//; $opt{'O'} =~ s/\.compare$//;

open IN, "$ARGV[0]";
while ($l = <IN>) {
	chomp $l;
	($cellID,$annot) = split(/\t/, $l);
	$CELLID_annot1{$cellID} = $annot;
} close IN;

$cell_match_ct = 0;
open IN, "$ARGV[1]";
while ($l = <IN>) {
	chomp $l;
	($cellID,$annot2) = split(/\t/, $l);
	if (defined $CELLID_annot1{$cellID}) {
		$annot1 = $CELLID_annot1{$cellID};
		$cell_match_ct++;
		$ANNOT1_ANNOT2_ct{$annot1}{$annot2}++;
		if (!defined $opt{'S'}) {
			$ANNOT1_ANNOT2_ct{$annot2}{$annot1}++;
		}
	}
} close IN;

open OUT, ">$opt{'O'}.compare.txt";
print OUT "Annot1\tAnnot2\tN_A2_in_A1\tPct_A2_in_A1\tPct_total\n";
foreach $annot1 (keys %ANNOT1_ANNOT2_ct) {
	$annot_ct = 0;
	foreach $annot2 (keys %{$ANNOT1_ANNOT2_ct{$annot1}}) {
		$annot_ct += $ANNOT1_ANNOT2_ct{$annot1}{$annot2};
	}
	foreach $annot2 (keys %{$ANNOT1_ANNOT2_ct{$annot1}}) {
		$A1_pct = sprintf("%.3f", ($ANNOT1_ANNOT2_ct{$annot1}{$annot2}/$annot_ct)*100);
		$G_pct = sprintf("%.3f", ($ANNOT1_ANNOT2_ct{$annot1}{$annot2}/$cell_match_ct)*100);
		print OUT "$annot1\t$annot2\t$ANNOT1_ANNOT2_ct{$annot1}{$annot2}\t$A1_pct\t$G_pct\n";
	}
} close OUT;

# plotting?

open R, ">$opt{'O'}.compare.plot.r";
print R "
library(ggplot2)
IN<-read.table(header=T,\"$opt{'O'}.compare.txt\")
PLT<-ggplot() + theme_bw() +
	geom_col(aes(x=IN\$Annot1,y=IN\$Pct_A2_in_A1,fill=IN\$Annot2)) +
	labs(fill=\"Annot 2\") + xlab(\"Annot 1\") + ylab(\"Percent Annot 2\")
ggsave(plot=PLT,filename=\"$opt{'O'}.compare.png\")
ggsave(plot=PLT,filename=\"$opt{'O'}.compare.pdf\")
";
close R;
system("Rscript $opt{'O'}.compare.plot.r && rm -f $opt{'O'}.compare.plot.r");

} # end regular

}
1;