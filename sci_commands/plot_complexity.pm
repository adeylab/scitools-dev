package sci_commands::plot_complexity;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("plot_complexity");

sub plot_complexity {

@ARGV = @_;
getopts("O:A:a:C:c:R:Mf:p:k:h:w:T:K:m:y:t:D:HN", \%opt);

#defaults
$kmean_centers = 3;
$alpha = 0.3;
$ptSize = 1;
$height = 6;
$width = 7;
$max = 6;
$title = "Library Complexity";
$contourCT = 50;
$die2 = "
scitools plot-complexity [options] [complexity file(s) can be comma separated]

Options:
   -O   [STR]   Output prefix (default is complexity file 1 prefix)
   -M           Run mixed model to determine read count cutoff for cells (def = no)
   -m   [INT]   Number of k-means clusters to use. (def=3)
   -K 	[INT] 	Force a minimum number of cells to be called by the knee analysis. (def=500)
   -N           Generate knee plot (def = no)
   -A   [STR]   Annotation file (to color code points)
   -a   [STR]   Comma separated list of annoations to include in plot
                 (requires -A to be specified)
   -C   [STR]   Color coding file (annot (tab) #hexColor)
   -c   [STR]   Color coding string
                 Annot=#hexColor,Annot2=#hexColor
   -p   [FLT]   Point size (def = $ptSize)
   -f   [FLT]   Alpha for plotting points (def = $alpha)
   -k   [STR]   If defined will color density lines the specified color (def = same as points)
                 either #hexcolor, or colorName
   -t   [INT]   Number of contours for 2d density (def = $contourCT)
   -w   [FLT]   Plot width (def = $width)
   -h   [FLT]   Plot height (def = $height)
   -y   [INT]   Max scale for plot in log10 unique reads (def = $max)
   -T   [STR]   Title (def = $title)
   -R   [STR]   Rscript call (def = $Rscript)

Note: Requires ggplot2 R package

";

if (!defined $ARGV[0]) {die $die2};
if (defined $opt{'a'} && !defined $opt{'A'}) {die "\nMust provide an annotaiton file (-A) if specifying annotations to plot (-a)!\n$die2"};
if (defined $opt{'C'} && defined $opt{'c'}) {die "\nSpecify either a color string (-c) or a color coding file (-C), not both!\n$die2"};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (defined $opt{'A'}) {read_annot($opt{'A'})};
if (defined $opt{'a'}) {
	@ANNOT_LIST = split(/,/, $opt{'a'});
	foreach $annot (@ANNOT_LIST) {
		$ANNOT_include{$annot} = 1;
	}
}
if (defined $opt{'C'}) {read_color_file($opt{'C'})};
if (defined $opt{'c'}) {read_color_string($opt{'c'})};
if (defined $opt{'m'}) {$kmean_centers = $opt{'m'}};
if (defined $opt{'p'}) {$ptSize = $opt{'p'}};
if (defined $opt{'f'}) {$alpha = $opt{'f'}};
if (defined $opt{'k'}) {
	if ($opt{'k'} =~ /^#/) {$cont_col = $opt{'k'}}
	else {$cont_col = "\"$opt{'k'}\""};
}
if (defined $opt{'h'}) {$height = $opt{'h'}};
if (defined $opt{'w'}) {$width = $opt{'w'}};
if (defined $opt{'T'}) {$title = $opt{'T'}; $title =~ s/"//g};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
if (!defined $opt{'K'}) {$opt{'K'} = 500};
if (defined $opt{'y'}) {$max = $opt{'y'}};
if (defined $opt{'t'}) {$contourCT = $opt{'t'}};
if (defined $opt{'D'}) {$hexBins = $opt{'D'}; $opt{'H'} = 1};

if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.txt$//;

read_complexity($ARGV[0]);

open OUT, ">$opt{'O'}.plot.txt";
foreach $cellID (keys %CELLID_complexity) {
	if (defined $opt{'a'}) {
		$annot = $CELLID_annot{$cellID};
		if (defined $ANNOT_include{$annot} && defined $CELLID_annot{$cellID}) {
			print OUT "$cellID\t$CELLID_annot{$cellID}\t$CELLID_uniq_reads{$cellID}\t$CELLID_complexity{$cellID}\n";
		}
	} elsif (defined $opt{'A'} && defined $CELLID_annot{$cellID}) {
		print OUT "$cellID\t$CELLID_annot{$cellID}\t$CELLID_uniq_reads{$cellID}\t$CELLID_complexity{$cellID}\n";
	} else {
		print OUT "$cellID\tCell\t$CELLID_uniq_reads{$cellID}\t$CELLID_complexity{$cellID}\n";
	}
} close OUT;

open R, ">$opt{'O'}.plot.r";
print R "
library(MASS)
library(mixtools)
library(ggplot2)
IN<-read.table(\"$opt{'O'}.plot.txt\")
PLT<-ggplot(data=subset(IN,V4<100&V4>0)) + theme_bw() +
";
if (!defined $opt{'H'}) {
	print R "   geom_point(aes(V4,log10(V3),color=V2),size=$ptSize,alpha=$alpha,shape=16) +\n";
	if (defined $opt{'k'}) {
		print R "   geom_density2d(aes(V4,log10(V3),color=$cont_col,bins=$contourCT),size=0.3) +\n";
		} else {
		print R "   geom_density2d(aes(V4,log10(V3),color=V2,bins=$contourCT),size=0.3) +\n";
	}
	if (defined $opt{'C'} || defined $opt{'c'}) {
	print R "   scale_colour_manual(values = c($color_mapping)) +
guides(colour = guide_legend(override.aes = list(size=4))) +\n";
	}
} else {
	if (!defined $opt{'c'} && !defined $opt{'C'} && !defined $opt{'A'}) {
		print R "   geom_hex(aes(V4,log10(V3)),color=\"lightsteelblue4\"),bins=$hexBins) +\n";
	} else {
		print R "	geom_hex(aes(V4,log10(V3),fill=annot),bins=$hexBins) +
	guides(colour = guide_legend(override.aes = list(size=4))) +\n";
		if (defined $opt{'C'} || defined $opt{'c'}) {
				print R "	scale_fill_manual(values = c($color_mapping)) +\n";
		}
	}
}


print R "
   scale_x_continuous(limits=c(0,100)) +
   scale_y_continuous(limits=c(0,$max)) +
   xlab(\"Percent Passing Reads\") +
   ylab(\"log10 Passing Reads\") +
   ggtitle(\"$title\") +";
if (defined $opt{'A'}) {
print R "
	theme(legend.background=element_blank(),legend.title=element_blank())";
} else {
print R "
	theme(legend.position=\"none\")";
}
print R "
ggsave(plot=PLT,filename=\"$opt{'O'}.png\",width=$width,height=$height)
ggsave(plot=PLT,filename=\"$opt{'O'}.pdf\",width=$width,height=$height)
";

if (defined $opt{'M'}) {
	print R "
IN_sub=subset(IN,V4<100&V4>0)

#take unique aligned 
IN_sub\$unique_aligned<-IN_sub\$V3 
x1 <- IN_sub\$unique_aligned[IN_sub\$unique_aligned != 0]

# trimodal fit
km <- kmeans(log10(x1),centers=$kmean_centers)
clustr <- as.factor(km\$cluster)
data_1<-data.frame(\"val\"=log10(x1),\"km\"=clustr)
highest_cluster<-which.max(km\$centers)
second<-km\$center[-highest_cluster,1]
second_highest_cluster<-as.numeric(names(which.max(second)))
cluster_top<-data_1\$val[which(data_1\$km==highest_cluster)]
border1=mean(max(data_1\$val[which(data_1\$km==second_highest_cluster)]),min(data_1\$val[which(data_1\$km==highest_cluster)]))

#in case the groups seperate we do normal fit
sink(\"$opt{'O'}.log\")
fit <- fitdistr(cluster_top, \"normal\")
para <- fit\$estimate
threshold_final_1=para[1]-1.96*para[2]

#in case 3 groups do not sepearte well do a mixed model instead of normal
MMS<-normalmixEM(log10(cluster_top),maxrestarts=1000,maxit = 10000)
threshold_final_2 =10**(MMS\$mu[2]-1.96*MMS\$sigma[2])

p1<-ggplot(data_1, aes(x=data_1\$val)) + geom_histogram(aes(fill=data_1\$km,y=..count../sum(..count..)),color=\"grey50\")+ylab(\"Density\")+stat_density(geom=\"line\",color=\"red\")+geom_vline(xintercept = threshold_final_1)+xlab(\"Log10 Unique Reads\")+theme_bw()+theme(legend.background=element_blank(),text = element_text(size=10),legend.title=element_blank())
ggsave(p1,filename = \"$opt{'O'}.dist.threshold_normal.pdf\")
ggsave(p1,filename = \"$opt{'O'}.dist.threshold_normal.png\")

p2<-ggplot(data_1, aes(x=data_1\$val)) + geom_histogram(aes(fill=data_1\$km,y=..count../sum(..count..)),color=\"grey50\")+ylab(\"Density\")+stat_density(geom=\"line\",color=\"red\")+geom_vline(xintercept = threshold_final_2)+xlab(\"Log10 Unique Reads\")+theme_bw()+theme(legend.background=element_blank(),text = element_text(size=10),legend.title=element_blank())
ggsave(p2,filename = \"$opt{'O'}.dist.threshold_mixed.pdf\")
ggsave(p2,filename = \"$opt{'O'}.dist.threshold_mixed.png\")

print(paste(\"the number of total cells: \",length(IN_sub\$unique_aligned)))
print(paste(\"threshold 1: \",10**threshold_final_1))
print(paste(\"the number of cells above this threshold 1: \",length(which(IN_sub\$unique_aligned>round(10**threshold_final_1)))))

print(paste(\"threshold 2: \",10**threshold_final_2))
print(paste(\"the number of cells above this threshold 2: \",length(which(IN_sub\$unique_aligned>round(10**threshold_final_2)))))
sink()
";
}

print R "
PLT<-ggplot(data=subset(IN,V4<100&V4>0)) + theme_bw() +
	geom_histogram(aes(log10(V3),fill=V2)) +";
if (defined $opt{'C'} || defined $opt{'c'}) {
	print R "   scale_fill_manual(values = c($color_mapping)) +";
}
print R "
	xlab(\"log10 Passing Reads\") +
	ylab(\"Counts\") +
	ggtitle(\"$title\") +
	scale_x_continuous(limits=c(0,$max)) +";
if (defined $opt{'A'}) {
print R "
	theme(legend.background=element_blank(),legend.title=element_blank())";
} else {
print R "
	theme(legend.position=\"none\")";
}
print R "
ggsave(plot=PLT,filename=\"$opt{'O'}.hist.png\",width=$width,height=$height)
ggsave(plot=PLT,filename=\"$opt{'O'}.hist.pdf\",width=$width,height=$height)
";

if (defined $opt{'N'}) {

	# knee plotting
	print R "
library(inflection)
IN_sub<-subset(IN,V4<100&V4>0)
IN_sub\$cell_order<-NA
IN_sub[order(IN_sub\$V3,decreasing=T),]\$cell_order<-c(1:nrow(IN_sub))
IN_sub_forcedcells<-subset(IN_sub,cell_order>$opt{'K'})
#first order then call knee
order_of_cells<-order(IN_sub_forcedcells\$cell_order)
IN_sub_forcedcells<-IN_sub_forcedcells[order_of_cells,]

kneecalling_xintercept<-uik(x=log10(IN_sub_forcedcells\$cell_order),y=log10(IN_sub_forcedcells\$V3))

cell_count_cutoff<-10^kneecalling_xintercept
read_cutoff<-IN_sub[IN_sub\$cell_order==as.integer(cell_count_cutoff),]\$V3
PLT<-ggplot(data=IN_sub) + theme_bw() +
	geom_point(aes(log10(cell_order),log10(V3))) + 
	geom_vline(xintercept=kneecalling_xintercept,color=\"red\") +
	geom_text(aes(label=paste(\"Cells:\",cell_count_cutoff,\"\\n Read Cutoff:\",read_cutoff),x=max(log10(IN_sub\$cell_order))-0.5,y=max(log10(IN_sub\$V3))-0.5)) +
";

	if (defined $opt{'C'} || defined $opt{'c'}) {
		print R "   scale_colour_manual(values = c($color_mapping)) +";
	}
	print R "
	xlab(\"log10 Cell Order\") +
	ylab(\"log10 Passing Reads\") +
	ggtitle(\"$title\") +";
	if (defined $opt{'A'}) {
	print R "
	theme(legend.background=element_blank(),legend.title=element_blank())";
	} else {
	print R "
	theme(legend.position=\"none\")";
	}
	print R "
ggsave(plot=PLT,filename=\"$opt{'O'}.knee.png\",width=$width,height=$height)
ggsave(plot=PLT,filename=\"$opt{'O'}.knee.pdf\",width=$width,height=$height)
	";
}

close R;

system("$Rscript $opt{'O'}.plot.r");

}
1;
