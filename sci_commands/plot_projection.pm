package sci_commands::plot_projection;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("plot_projection");

sub plot_projection {

@ARGV = @_;

# Defaults
$gradient_def = "BuYlRd";
$point_size = 3;
$ribbon_color = "lightsteelblue4";
$width = 6;
$height = 5;
$median_width = 1;
$line_color = "black";
$ribbon_alpha = 0.5;
$ymax = 6;
$cell_width = 0.15;
$quantile_width = 0.5;

getopts("O:G:ep:r:a:c:y:q:m:l:w:h:R:X", \%opt);

$die2 = "
scitools plot-projection [options] [output directory from bam-project]
   or    plot-projections
   or    plot-project

Plots the projection of reads from bam-project output folder.

Options:
   -O   [STR]   Output prefix (def = bam-projection folder)
   -G   [GRD]   Color gradient in plot (def = $gradient_def)
                For all available gradients, run 'scitools gradient'
   -e           Plot projections for every cell (def = no)
   -p   [FLT]   Point size (def = $point_size)
   -r   [STR]   Ribbon/ind. cell color (def = $ribbon_color)
   -a   [FLT]   Ribbon/ind cell alpha (def = $ribbon_alpha)
   -c   [FLT]   Individual cell line width (if -e, def = $cell_width)
   -m   [FLT]   Median cell line width (def = $median_width)
   -l   [STR]   Line color (def = $line_color)
   -y   [FLT]   Y-axis max (def = $ymax)
   -q   [FLT]   Quantile line width (if no -e, def = $quantile_width)
   -w   [FLT]   Plot width (def = $width)
   -h   [FLT]   Plot height (def = $height)
   -R   [STR]   Rscript call (def = $Rscript)
   -X           Keep intermediate files (def = delete)

  ";

if (!defined $ARGV[0]) {die $die2};
$ARGV[0] =~ s/\/$//;
if (!defined $opt{'O'}) {$opt{'O'} = "$ARGV[0]/projected_complexity_plot"};
if (!defined $opt{'G'}) {$opt{'G'} = $gradient_def};
$gradient_function = get_gradient($opt{'G'});
if (defined $opt{'p'}) {$point_size = $opt{'p'}};
if ($point_size>1) {
	$inner_size = $point_size-1;
} else {
	$inner_size = $point_size*0.5;
}
if (defined $opt{'r'}) {$ribbon_color = $opt{'r'}};
if (defined $opt{'a'}) {$ribbon_alpha = $opt{'a'}};
if (defined $opt{'c'}) {$cell_width = $opt{'c'}};
if (defined $opt{'y'}) {$ymax = $opt{'y'}};
if (defined $opt{'q'}) {$quantile_width = $opt{'q'}};
if (defined $opt{'m'}) {$median_width = $opt{'m'}};
if (defined $opt{'l'}) {$line_color = $opt{'l'}};
if (defined $opt{'w'}) {$width = $opt{'w'}};
if (defined $opt{'h'}) {$height = $opt{'h'}};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};

if (-e "$ARGV[0]/cell_summaries.txt") {print STDERR "$ARGV[0]/cell_summaries.txt found!\n"} else {die "ERROR: Cannot find $ARGV[0]/cell_summaries.txt!\n"};
if (-e "$ARGV[0]/summary_projections.txt") {print STDERR "$ARGV[0]/summary_projections.txt found!\n"} else {die "ERROR: Cannot find $ARGV[0]/summary_projections.txt!\n"};
if (-e "$ARGV[0]/cell_projections.txt") {print STDERR "$ARGV[0]/cell_projections.txt found!\n"} else {die "ERROR: Cannot find $ARGV[0]/cell_projections.txt!\n"};

open R, ">$opt{'O'}.r";

print R "library(ggplot2)
$gradient_function
";

if (!defined $opt{'e'}) {
print R "
RANGE<-read.table(\"$ARGV[0]/cell_summaries.txt\")
SUMMARY<-read.table(\"$ARGV[0]/summary_projections.txt\")
PLT<-ggplot() + theme_bw() +
	geom_ribbon(aes(x=(RANGE\$V1*100),ymin=log10(RANGE\$V6),ymax=log10(RANGE\$V5)),fill=\"$ribbon_color\",alpha=$ribbon_alpha) +
	geom_line(aes((RANGE\$V1*100),log10(RANGE\$V5)),color=\"$line_color\",size=$quantile_width,linetype=\"dashed\") +
	geom_line(aes((RANGE\$V1*100),log10(RANGE\$V6)),color=\"$line_color\",size=$quantile_width,linetype=\"dashed\") +
	geom_line(aes((RANGE\$V1*100),log10(RANGE\$V7)),color=\"$line_color\",size=$quantile_width,linetype=\"dashed\") +
	geom_line(aes((RANGE\$V1*100),log10(RANGE\$V8)),color=\"$line_color\",size=$quantile_width,linetype=\"dashed\") +";
} else {
print R "
CELLS<-read.table(\"$ARGV[0]/cell_projections.txt\")
SUMMARY<-read.table(\"$ARGV[0]/summary_projections.txt\")
PLT<-ggplot() + theme_bw() +
	geom_line(aes((CELLS\$V5*100),log10(CELLS\$V4),group=CELLS\$V1),alpha=$ribbon_alpha,size=$cell_width,color=\"$ribbon_color\") +";
}

print R "
	geom_line(aes((SUMMARY\$V1*100),log10(SUMMARY\$V3)),color=\"$line_color\",size=$median_width) +
	geom_point(aes((SUMMARY\$V1*100),log10(SUMMARY\$V3)),color=\"$line_color\",size=$point_size) +
	geom_point(aes((SUMMARY\$V1*100),log10(SUMMARY\$V3),color=log10(SUMMARY\$V2)),size=$inner_size) +
	scale_color_gradientn(colours=gradient_funct(99)) +
	scale_x_continuous(limits=c(0,100)) +
	scale_y_continuous(limits=c(2,$ymax)) +
	xlab(\"Complexity\") +
	ylab(\"log10 Passing Reads\") +
	labs(color=\"Log10 Total\\nReads\")
ggsave(plot=PLT,filename=\"$opt{'O'}.png\",width=$width,height=$height)
ggsave(plot=PLT,filename=\"$opt{'O'}.pdf\",width=$width,height=$height)
";

close R;

system("$Rscript $opt{'O'}.r 2>/dev/null");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.r");
}

}
1;