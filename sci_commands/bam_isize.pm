package sci_commands::bam_isize;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("bam_isize");

sub bam_isize {

@ARGV = @_;
getopts("s:A:a:O:R:M:", \%opt);

$max = "1000";
$width = 4.5;
$height = 4;

$die2 = "
scitools bam-isize [options] [input bam]
   or    isize

Produces a file & plot of insert size distribution.

Options:
   -O   [STR]   Output (default is bam file prefix)
   -M   [INT]   Max isize to include (def = $max)
   -A   [STR]   Annotation file (will split by annot)
   -a   [STR]   Comma separated list of annoations to include in plot
                 (requires -A to be specified)
   -C   [STR]   Color coding file (annot (tab) #hexColor)
   -c   [STR]   Color coding string
                 Annot=#hexColor,Annot2=#hexColor
   -w   [FLT]   Plot width (def = $width)
   -h   [FLT]   Plot height (def = $height)
   -s   [STR]   Samtools call (def = $samtools)
   -R   [STR]   Rscript call (def = $Rscript)

";

if (!defined $ARGV[0]) {die $die2};
if (defined $opt{'s'}) {$samtools = $opt{'s'}};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.bam$//;
if (defined $opt{'C'}) {read_color_file($opt{'C'})};
if (defined $opt{'c'}) {read_color_string($opt{'c'})};
if (defined $opt{'M'}) {$max = $opt{'M'}};
if (defined $opt{'h'}) {$height = $opt{'h'}};
if (defined $opt{'w'}) {$width = $opt{'w'}};

if (defined $opt{'A'}) {read_annot($opt{'A'})};
if (defined $opt{'a'}) {
	@ANNOT_LIST = split(/,/, $opt{'a'});
	foreach $annot (@ANNOT_LIST) {
		$ANNOT_include{$annot} = 1;
	}
}

open OUT, ">$opt{'O'}.isize.values";

open IN, "$samtools view $ARGV[0] 2>/dev/null |";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$cellID = $P[0]; $cellID =~ s/:.+$//;
	if ((defined $ANNOT_include{$CELLID_annot{$cellID}} || !defined $opt{'a'}) && $P[8] > 0 && $P[8] <= $max) {
		if (defined $opt{'A'}) {
			print OUT "$CELLID_annot{$cellID}\t$P[8]\n";
		} else {
			print OUT "$P[8]\n";
		}
	}
} close IN; close OUT;

open R, ">$opt{'O'}.isize.plot.r";
if (!defined $opt{'A'}) {
print R "
library(ggplot2)
IN<-read.table(\"$opt{'O'}.isize.values\")
PLT<-ggplot() + theme_bw() +
	geom_density(aes(IN\$V1)) +
	scale_x_continuous(limits=c(0,$max)) +
	xlab(\"Insert Size\") +
	ylab(\"Density\") +
	ggtitle(\"Insert Size Distribution\")
ggsave(plot=PLT,filename=\"$opt{'O'}.isize.plot.pdf\",width=$width,height=$height)
ggsave(plot=PLT,filename=\"$opt{'O'}.isize.plot.png\",width=$width,height=$height)
";
} else {
print R "
library(ggplot2)
IN<-read.table(\"$opt{'O'}.isize.values\")
PLT<-ggplot() + theme_bw() +
	geom_density(aes(IN\$V2,colour=IN\$V1)) +";
if (defined $opt{'C'} || defined $opt{'c'}) {
	print R "
	scale_colour_manual(values = c($color_mapping)) +
	guides(colour = guide_legend(override.aes = list(size=4))) +";
}
print R "
	scale_x_continuous(limits=c(0,$max)) +
	xlab(\"Insert Size\") +
	ylab(\"Density\") +
	ggtitle(\"Insert Size Distribution\") +
	theme(legend.background=element_blank(),legend.title=element_blank())
ggsave(plot=PLT,filename=\"$opt{'O'}.isize.plot.pdf\",width=$width,height=$height)
ggsave(plot=PLT,filename=\"$opt{'O'}.isize.plot.png\",width=$width,height=$height)
";	
}
close R;

system("$Rscript $opt{'O'}.isize.plot.r && rm -f $opt{'O'}.isize.plot.r");

}
1;
