package sci_commands::plot_values;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("plot_values");

sub plot_values {

$type = "boxplot";
$width = 5;
$height = 4;
$ptSize = 0.5;
$color_mapping = "none";
$alpha = 0.5;
$ylab = "Feature value";
$val_col = 2;
$cell_col = 1;

@ARGV = @_;
getopts("O:A:a:C:c:R:T:M:Xs:p:f:w:h:Wr:L:l:m:P:n:i:", \%opt);

$die2 = "
scitools plot-values [options] [values file / null if -M]

Options - general:
   -O   [STR]   Output prefix (def: vals prefix; adds type suffix)
   -T   [STR]   Type, comma separated list: (def = $type)
                  violin/vln/v
                  boxplot/box/b
                  histogram/hist/h
                  density/dens/d
   -l   [INT]   Column with values
                (def = $val_col, i.e. cellID (tab) value)
   -i   [INT]   CellID column (def = $cell_col)
   -m   [FLT]   Maximum value (def = none)
   -n   [FLT]   Minimum value (def = none)
   -P   [STR]   Operation (def = none)
                  log2, log10, inv/inverse, ln
   -p   [FLT]   Point size (def = $ptSize)
   -f   [FLT]   Alpha for fill (def = $alpha)
   -L   [STR]   Name for values axis (def = $ylab)
   -w   [FLT]   Plot width (inches, def = $width)
   -h   [FLT]   Plot height (inches, def = $height)
   -W           If multiple types, plot in same file.
                 Note: width is for each individual panel
   -r   [INT]   Panel row number (def = sqrt(types))

Plotting by annotations:
   -A   [STR]   Annotation file (to color code points)
   -a   [STR]   Comma separated list of annotations to include in plot
                  (requires -A to be specified)
   -C   [STR]   Color coding file (annot (tab) #hexColor)
   -c   [STR]   Color coding string
                  Annot=#hexColor,Annot2=#hexColor

To plot multiple values specified in a matrix file:
   -M   [STR]   Matrix file (e.g. deviation z-scores from chromVAR)
                Will create a folder as -O option and produce
                a plot for each row of the matrix.
                Only includes rows with at least one non-zero value.

Other options:
   -R   [STR]   Rscript call (def = $Rscript)
   -s   [STR]   scitools call (def = $scitools)
   -X           Do not delete intermediate files (def = delete)

Note: Requires ggplot2 R package

";

if (!defined $ARGV[0] && !defined $opt{'M'}) {die $die2};
if (defined $opt{'a'} && !defined $opt{'A'}) {die "\nMust provide an annotation file (-A) if specifying annotations to plot (-a)!\n$die2"};
if (defined $opt{'C'} && defined $opt{'c'}) {die "\nSpecify either a color string (-c) or a color coding file (-C), not both!\n$die2"};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.txt$//; $opt{'O'} =~ s/\.values$//; $opt{'O'} =~ s/\.matrix$//;
if (defined $opt{'p'}) {$ptSize = $opt{'p'}};
if (defined $opt{'f'}) {$alpha = $opt{'f'}};
if (defined $opt{'h'}) {$height = $opt{'h'}};
if (defined $opt{'w'}) {$width = $opt{'w'}};
if (defined $opt{'r'}) {$panel_nrow = $opt{'r'}};
if (defined $opt{'L'}) {$ylab = $opt{'L'}};
if (defined $opt{'l'}) {$val_col = $opt{'l'}};
if (defined $opt{'i'}) {$cell_col = $opt{'i'}};

if (defined $opt{'M'}) {
	print STDERR "SCITOOLS: Matrix file plotting detected! Will plot a separate value-based plot for each row entry of the matrix.\n";
	
	$matrix_out = "matrix_plots";
	
	if (-e "$opt{'O'}.$matrix_out") {
		die "\nFATAL: $opt{'O'}.$matrix_out directory already exists! Exiting!\n$die2";
	}
	system("mkdir $opt{'O'}.$matrix_out");
	
	read_matrix($opt{'M'});
	
	$common_opts = "";
	if (defined $opt{'A'}) {$common_opts .= "-A $opt{'A'} "};
	if (defined $opt{'C'}) {$common_opts .= "-C $opt{'C'} "};
	if (defined $opt{'c'}) {$common_opts .= "-c $opt{'c'} "};
	if (defined $opt{'a'}) {$common_opts .= "-a $opt{'a'} "};
	if (defined $opt{'T'}) {$common_opts .= "-T $opt{'T'} "};
	if (defined $opt{'R'}) {$common_opts .= "-R $opt{'R'} "};
	if (defined $opt{'p'}) {$common_opts .= "-p $opt{'p'} "};
	if (defined $opt{'f'}) {$common_opts .= "-f $opt{'f'} "};
	if (defined $opt{'w'}) {$common_opts .= "-w $opt{'w'} "};
	if (defined $opt{'h'}) {$common_opts .= "-h $opt{'h'} "};
	if (defined $opt{'r'}) {$common_opts .= "-r $opt{'r'} "};
	if (defined $opt{'X'}) {$common_opts .= "-X "};
	if (defined $opt{'W'}) {$common_opts .= "-W "};
	
	$common_opts =~ s/\s$//;
	
	foreach $feature (keys %MATRIX_feature_nonZero) {
		$feature_polished = $feature;
		$feature_polished =~ s/(:|;|'|\.|\(|\)|\{|\}|\[|\]|\s|\|)/_/;
		open VALS, ">$opt{'O'}.$matrix_out/$feature_polished.values";
		foreach $cellID (keys %CELLID_FEATURE_value) {
			if (defined $CELLID_FEATURE_value{$cellID}{$feature}) {
				print VALS "$cellID\t$CELLID_FEATURE_value{$cellID}{$feature}\n";
			}
		} close VALS;
		system("$scitools plot-values $common_opts -O $opt{'O'}.$matrix_out/$feature_polished $opt{'O'}.$matrix_out/$feature_polished.values");
		if (!defined $opt{'X'}) {system("rm -f $opt{'O'}.$matrix_out/$feature_polished.values")};
	}
	
	exit;
	
}

if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (defined $opt{'s'}) {$scitools = $opt{'s'}};
if (defined $opt{'A'}) {read_annot($opt{'A'})};
if (defined $opt{'a'}) {
	@ANNOT_LIST = split(/,/, $opt{'a'});
	foreach $annot (@ANNOT_LIST) {
		$ANNOT_include{$annot} = 1;
	}
}
if (defined $opt{'C'}) {read_color_file($opt{'C'})};
if (defined $opt{'c'}) {read_color_string($opt{'c'})};
if (defined $opt{'T'}) {$type = $opt{'T'}};
@TYPES = split(/,/, $type);

$cell_col--;
$val_col--;
%CELLID_value = ();
open IN, "$ARGV[0]";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$cellID = $P[$cell_col];
	$value = $P[$val_col];
	$pass = 1;
	if (defined $opt{'m'} && $value > $opt{'M'}) {$pass = 0};
	if (defined $opt{'n'} && $value < $opt{'n'}) {$pass = 0};
	if ($pass>0) {
		if (defined $opt{'P'}) {
			$oval = $value;
			if ($opt{'P'} =~ /log2/) {$value = log($oval)/log(2)};
			if ($opt{'P'} =~ /log10/) {$value = log($oval)/log(10)};
			if ($opt{'P'} =~ /ln/) {$value = log($oval)};
			if ($opt{'P'} =~ /^inv/) {$value = 1/$oval};
		}
		$CELLID_value{$cellID} = $value;
	}
} close IN;

open OUT, ">$opt{'O'}.plot.txt";
foreach $cellID (keys %CELLID_value) {
	if (!defined $opt{'A'}) {$CELLID_annot{$cellID} = "Cell"};
	if (defined $CELLID_annot{$cellID}) {
		if (!defined $opt{'a'} || (defined $opt{'a'} && defined $ANNOT_include{$CELLID_annot{$cellID}})) {
			print OUT "$cellID\t$CELLID_annot{$cellID}\t$CELLID_value{$cellID}\n";
		}
	}
} close OUT;

open R, ">$opt{'O'}.plot.r";
print R "
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

IN<-read.table(\"$opt{'O'}.plot.txt\")

";

if (defined $opt{'W'}) {$grid_list = ""};

for ($typeID = 0; $typeID < @TYPES; $typeID++) {
	
	# VIOLIN
	if ($TYPES[$typeID] =~ /^v/) {
print R "VLN<-ggplot() + theme_bw() +
	geom_violin(aes(IN\$V2,IN\$V3,fill=IN\$V2),alpha=$alpha,color=\"gray30\",size=0.5) +
	geom_jitter(aes(IN\$V2,IN\$V3),color=\"gray30\",size=0.15) +";
if ($color_mapping !~ /none/i) {
	print R "
	scale_fill_manual(values = c($color_mapping)) +";
}
print R "
	guides(fill=FALSE,colour=FALSE) +
	xlab(\"Annot\") +
	ylab(\"$ylab\") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
";
if (!defined $opt{'W'}) {
print R "ggsave(plot=VLN,filename=\"$opt{'O'}.violin.png\",width=$width,height=$height,dpi=900)
ggsave(plot=VLN,filename=\"$opt{'O'}.violin.pdf\",width=$width,height=$height)
";
} else {
	$grid_list .= "VLN,";
}

	# BOXPLOT
	} elsif ($TYPES[$typeID] =~ /^b/) {
print R "BOX<-ggplot() + theme_bw() +
	geom_boxplot(aes(IN\$V2,IN\$V3,fill=IN\$V2),size=0.5) +";
if ($color_mapping !~ /none/i) {
	print R "
	scale_fill_manual(values = c($color_mapping)) +";
}
print R "
	guides(fill=FALSE,colour=FALSE) +
	xlab(\"Annot\") +
	ylab(\"$ylab\") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
";
if (!defined $opt{'W'}) {
print R "ggsave(plot=BOX,filename=\"$opt{'O'}.boxplot.png\",width=$width,height=$height,dpi=900)
ggsave(plot=BOX,filename=\"$opt{'O'}.boxplot.pdf\",width=$width,height=$height)
";
} else {
	$grid_list .= "BOX,";
}
	
	# HISTOGRAM
	} elsif ($TYPES[$typeID] =~ /^h/) {
print R "HIST<-ggplot() + theme_bw() +
	geom_histogram(aes(IN\$V3,fill=IN\$V2),size=0.5,alpha=1) +";
if ($color_mapping !~ /none/i) {
	print R "
	scale_fill_manual(values = c($color_mapping)) +";
}
print R "
	guides(fill=FALSE,colour=FALSE) +
	xlab(\"$ylab\") +
	ylab(\"Count\")
";
if (!defined $opt{'W'}) {
print R "ggsave(plot=HIST,filename=\"$opt{'O'}.histogram.png\",width=$width,height=$height,dpi=900)
ggsave(plot=HIST,filename=\"$opt{'O'}.histogram.pdf\",width=$width,height=$height)
";
} else {
	$grid_list .= "HIST,";
}
	
	# DENSITY
	} elsif ($TYPES[$typeID] =~ /^d/) {
print R "DENS<-ggplot() + theme_bw() +
	geom_density(aes(IN\$V3,fill=IN\$V2,colour=IN\$V2),alpha=$alpha,size=0.5) +";
if ($color_mapping !~ /none/i) {
	print R "
	scale_fill_manual(values = c($color_mapping)) +
	scale_colour_manual(values = c($color_mapping)) +";
}
print R "
	guides(fill=FALSE,colour=FALSE) +
	xlab(\"$ylab\") +
	ylab(\"Density\")
";
if (!defined $opt{'W'}) {
print R "ggsave(plot=DENS,filename=\"$opt{'O'}.density.png\",width=$width,height=$height,dpi=900)
ggsave(plot=DENS,filename=\"$opt{'O'}.density.pdf\",width=$width,height=$height)
";
} else {
	$grid_list .= "DENS,";
}

	# NOT A TYPE
	} else {
		print STDERR "ERROR: Cannot determine type \"$TYPES[$typeID]\"\n";
	}
}

# PLOT PANELS
if (defined $opt{'W'}) {
	
	if (!defined $opt{'r'}) {
		$panel_nrow = int(sqrt(@TYPES));
	}
	
	$ncol_factor = (@TYPES+1)/$panel_nrow;
	$panel_ncol = int($ncol_factor);
	if ($ncol_factor>$panel_ncol) {$panel_ncol++};
	$grid_width = $width*$panel_ncol;
	$grid_height = $height*$panel_nrow;
	
	$grid_list =~ s/,$//;
	
	print R "\nPLT_grid<-grid.arrange($grid_list,nrow=$panel_nrow)
ggsave(plot=PLT_grid,filename=\"$opt{'O'}.plot.png\",width=$grid_width,height=$grid_height,dpi=900)
ggsave(plot=PLT_grid,filename=\"$opt{'O'}.plot.pdf\",width=$grid_width,height=$grid_height)
";
}

close R;

system("$Rscript $opt{'O'}.plot.r");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.plot.r $opt{'O'}.plot.txt");
}

}
1;
