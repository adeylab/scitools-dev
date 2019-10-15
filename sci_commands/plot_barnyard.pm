package sci_commands::plot_barnyard;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("plot_barnyard");

sub plot_barnyard {

@ARGV = @_;
getopts("O:f:p:a:C:w:h:R:", \%opt);

$maxF = 0.1;
$ptSize = 0.75;
$alpha = 1;
$colors = "red3,blue3,purple3";
$width = 4;
$height = 4;

$die2 = "
scitools plot-barnyard [options] [barnyard_cells.txt file]

Options:
   -O   [STR]   Output prefix (default is input file prefix)
   -f   [FLT]   Max fraction of other species to be considered pure (def = $maxF)
   -p   [FLT]   Point size (def = $ptSize)
   -a   [FLT]   Alpha for plotting points (def = $alpha)
   -C   [STR]   Colors: human,mouse,mix (def = $colors)
   -w   [FLT]   Plot width (def = $width)
   -h   [FLT]   Plot height (def = $height)
   -R   [STR]   Rscript call (def = $Rscript)

";

if (!defined $ARGV[0]) {die $die2};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (defined $opt{'C'}) {$colors = $opt{'C'}};
($human_color,$mouse_color,$mix_color) = split(/,/, $colors);
if (defined $opt{'p'}) {$ptSize = $opt{'p'}};
if (defined $opt{'a'}) {$alpha = $opt{'a'}};
if (defined $opt{'h'}) {$height = $opt{'h'}};
if (defined $opt{'w'}) {$width = $opt{'w'}};
if (defined $opt{'f'}) {$maxF = $opt{'f'}};
$minF = 1-$maxF;
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.txt$//;

open R, ">$opt{'O'}.plot.r";
print R "
library(ggplot2)

IN<-read.table(\"$ARGV[0]\")
HUM<-subset(IN,\$V5>=$minF)
MUS<-subset(IN,\$V5<=$maxF)
MIX<-subset(IN,\$V5<$minF&\$V5>$maxF)

PLT<-ggplot() +
	geom_point(aes(MIX\$V3,MIX\$V4),color=\"$mix_color\",alpha=$alpha,size=$ptSize) +
	geom_point(aes(HUM\$V3,HUM\$V4),color=\"$human_color\",alpha=$alpha,size=$ptSize) +
	geom_point(aes(MUS\$V3,MUS\$V4),color=\"$mouse_color\",alpha=$alpha,size=$ptSize) +
	xlab(\"Human Passing Reads\") +
	ylab(\"Mouse Passing Reads\")

ggsave(plot=PLT,filename=\"$opt{'O'}.plot.png\",width=$width,height=$height)
ggsave(plot=PLT,filename=\"$opt{'O'}.plot.pdf\",width=$width,height=$height)
";

close R;

system("$Rscript $opt{'O'}.plot.r");

}
1;