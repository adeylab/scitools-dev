$Rscript = "Rscript"; #DEFAULT=Rscript
$Bedtools = "bedtools"; #DEFAULT=Bedtools

$hg19_ref_bedtools = "/home/groups/oroaklab/refs/hg19/hg19.fa.genome"; #DEFAULT=hg19_ref_bismark
$hg38_ref_bedtools = "/home/groups/oroaklab/refs/hg38/hg38.bedtools.genome"; #DEFAULT=hg38_ref_bismark
$mm10_ref_bedtools = "/home/groups/oroaklab/refs/mm10/mm10_bedtools_genome"; #DEFAULT=mm10_ref_bismark

use Getopt::Std; %opt = ();
getopts("R:B:A:F:f:c:p:", \%opt);
$die2 = "
scitools probe-view [options] [Feature_Probe_Set] [Bedtools Genome File] [NOMe.sorted.bed.gz File]

sciMET tool for looking at methylation percentage over feature probe space.Multiple strategies for clutering of cellIDs by WindowMatrix set. Will output Feature Summary statistics as well as clustering attempts.

[Feature_Probe_Set]             =               Bed file of feature set for methylation probe window generation.
[Bedtools Genome File]          =               Bedtools prepared genome file.
[NOMe.sorted.bed.gz File]   =               [CG/CH/GC].sorted.bed file generated from anyC_NOMe_processing script.

Options:
   -R   [STR]   Rscript call (Default: $Rscript)
   -B   [STR]   Bedtools call (Default: $Bedtools)
   -A   [STR]   Annotation File describing coloring for plots.
                                (Currently Required)
   -f   [INT]   Flank size for either side of the probe.
                (Default: 100)
   -F   [STR]   Feature name for output.
                (Default: all text within last \"/\" and first \"\.\"\ of [Feature_Probe_Set])
   -c   [STR]   Cytosine context for output.
                (Default: all text within last \"/\" and second \"\.\"\ of [NOMe.sorted.bed.gz File])
   -p   [INT]   Number of bins to split feature probe set into. This is useful for features of uneven size.
                If not specified, will generate probe window on base-pair resolution

Bedtools reference shortcuts:
   hg19  $hg19_ref_bedtools
   hg38  $hg38_ref_bedtools
   mm10  $mm10_ref_bedtools
";

%GENOMES = (
   "hg19" => "$hg19_ref_bedtools",
   "hg38" => "$hg38_ref_bedtools",
   "mm10" => "$mm10_ref_bedtools"
);

if (!defined $ARGV[0]) {die $die2};
if (!defined $ARGV[1]) {die $die2};
if (!defined $ARGV[2]) {die $die2};

if (!defined $opt{'A'}) {die $die2};

if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (defined $opt{'B'}) {$Bedtools = $opt{'B'}};
if (!defined $opt{'f'}) {$opt{'f'} = 100};
if (!defined $GENOMES{$ARGV[1]}) {die "\n\nERROR: Genome $ARGV[1] is not a proper genome reference! Exiting!\n"} else {$genome = $GENOMES{$ARGV[1]}};
if (!defined $opt{'F'}) {$opt{'F'}=$ARGV[0]; my @F = split(/\//,$opt{'F'}); $opt{'F'}=$F[-1]; @F = split(/\./,$opt{'F'}); $opt{'F'}=$F[0]};
if (!defined $opt{'c'}) {$opt{'c'}=$ARGV[2]; my @C = split(/\//,$opt{'c'}); $opt{'c'}=$C[-1]; @C = split(/\./,$opt{'c'}); $opt{'c'}=$C[1]};

$upstream="
$Bedtools flank -i $ARGV[0] -g $genome -l $opt{'f'} -r 0 | sort -k1,1 -k2,2n -T . - | $Bedtools intersect -sorted -wa -wb -a stdin -b $ARGV[2] | awk 'OFS=\"\\t\" {a=\$3-\$6; print \$1\"_\"\$2\"\_\"\$3, a, \$8, \$9, \$10, \"upstream\"}' - > $opt{'c'}.$opt{'F'}.bed";

$downstream="
$Bedtools flank -i $ARGV[0] -g $genome -r $opt{'f'} -l 0 | sort -k1,1 -k2,2n -T . - | $Bedtools intersect -sorted -wa -wb -a stdin -b $ARGV[2] | awk 'OFS=\"\\t\" {a=\$6-\$2; print \$1\"_\"\$2\"\_\"\$3, a, \$8, \$9, \$10, \"downstream\"}' - >> $opt{'c'}.$opt{'F'}.bed";

if (!defined $opt{'p'}){
$probe="sort -k1,1 -k2,2n -T $ARGV[1] - | $Bedtools intersect -sorted -wa -wb -a stdin -b $ARGV[2] | awk 'OFS=\"\\t\" {a=\$6-(\$2+(\$3-\$2)/2); print \$1\"_\"\$2\"\_\"\$3, a, \$8, \$9, \$10, \"probe\"}' - >> $opt{'c'}.$opt{'F'}.bed";
} else {
$probe="$Bedtools makewindows -b $ARGV[0] -n $opt{'p'} -i winnum -s 1 | sort -k1,1 -k2,2n -T . - | $Bedtools intersect -sorted -wa -wb -a stdin -b $ARGV[2] | awk 'OFS=\"\\t\" {print \$1\"_\"\$2\"\_\"\$3, \$4, \$8, \$9, \$10, \"probe\"}' - >> $opt{'c'}.$opt{'F'}.bed"; }

system($upstream);
system($downstream);
system($probe);

$probe_plot="

# probe view plot for:
#    annotation file = $opt{'A'}
#    probe = $ARGV[0]
#    methylation context = $opt{'c'}

library(ggplot2)
library(reshape2)
library(plyr)
library(gridExtra)
signal<-read.table(\"$opt{'c'}.$opt{'F'}.bed\")
colnames(signal)<-c(\"featID\",\"loc\",\"nonmet\",\"met\",\"cellID\",\"window\")
annot<-read.table(\"$opt{'A'}\");
colnames(annot)<-c(\"cellID\",\"x\",\"y\",\"clusterID\")
probe<-merge(signal,annot,by=\"cellID\")
probe_m<-ddply(probe,.(window,loc,clusterID),summarize,metm=sum(met),nonmetm=sum(nonmet))
rm(probe)

probe_split<-split(probe_m,probe_m\$window)
rm(probe_m)

probe_plot<-ggplot()+
  geom_line(aes(probe_split\$probe\$loc,probe_split\$probe\$metm/(probe_split\$probe\$metm+probe_split\$probe\$nonmetm)*100,group=probe_split\$probe\$clusterID,color=as.factor(probe_split\$probe\$clusterID)),alpha=1/10)+
  geom_smooth(aes(probe_split\$probe\$loc,probe_split\$probe\$metm/(probe_split\$probe\$metm+probe_split\$probe\$nonmetm)*100,group=probe_split\$probe\$clusterID,color=as.factor(probe_split\$probe\$clusterID)),span=0.3,na.rm=T,method=\"loess\")+
  xlab(\"Across Feature\")+
  ylab(\"\")+
  theme_bw()+
  theme(plot.background=element_blank(),legend.title=element_blank())+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0),limits=c(0,100))+
  theme(legend.position=\"none\",axis.text.y = element_blank(),plot.margin = unit(c(0, 0, 0, 0), \"cm\"))


upstream_plot<-ggplot()+
  geom_line(aes(rev(probe_split\$upstream\$loc),probe_split\$upstream\$metm/(probe_split\$upstream\$metm+probe_split\$upstream\$nonmetm)*100,group=probe_split\$upstream\$clusterID,color=as.factor(probe_split\$upstream\$clusterID)),alpha=1/10)+
  geom_smooth(aes(rev(probe_split\$upstream\$loc),probe_split\$upstream\$metm/(probe_split\$upstream\$metm+probe_split\$upstream\$nonmetm)*100,group=probe_split\$upstream\$clusterID,color=as.factor(probe_split\$upstream\$clusterID)),span=0.3,na.rm=T,method=\"loess\")+
  xlab(\"Upstream (bp)\")+
  ylab(\"Methylation Percent\")+
  theme_bw()+
  theme(plot.background=element_blank(),legend.title=element_blank())+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0),limits=c(0,100))+
  theme(legend.position=\"none\",plot.margin = unit(c(0, 0, 0, 0), \"cm\"))


downstream_plot<-ggplot()+
  geom_line(aes(probe_split\$downstream\$loc,probe_split\$downstream\$metm/(probe_split\$downstream\$metm+probe_split\$downstream\$nonmetm)*100,group=probe_split\$downstream\$clusterID,color=as.factor(probe_split\$downstream\$clusterID)),alpha=1/10)+
  geom_smooth(aes(probe_split\$downstream\$loc,probe_split\$downstream\$metm/(probe_split\$downstream\$metm+probe_split\$downstream\$nonmetm)*100,group=probe_split\$downstream\$clusterID,color=as.factor(probe_split\$downstream\$clusterID)),span=0.3,na.rm=T,method=\"loess\")+
  xlab(\"Downstream (bp)\")+
  ylab(\"\")+
  theme_bw()+
  theme(plot.background=element_blank(),legend.title=element_blank())+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0),limits=c(0,100))+
  theme(legend.position=\"none\",axis.text.y = element_blank(),plot.margin = unit(c(0, 0, 0, 0), \"cm\"))
out_plot<-grid.arrange(upstream_plot,probe_plot,downstream_plot,nrow=1)
ggsave(plot=out_plot,filename=\"$opt{'c'}.$opt{'F'}.probeview.png\",width=18,height=5)
";

open R, ">$opt{'c'}.$opt{'F'}.probeview.r";
print R $probe_plot;
close R;
system("$Rscript $opt{'c'}.$opt{'F'}.probeview.r");


#perl probeviewWIP.pl -A HCG.GENE.tab.bedHCH.100KB.hippoonly.UMAP.db.dims -f 500 /home/groups/oroaklab/adey_lab/projects/sciNOMe/refs/mm10_annotations/Homer_Individual_TF/mm10_EnsemblRegBuild_CTCF_Motif.bed.sorted.bed mm10 180428_allCells.HCG.NOMe.sorted.bed
#perl probeviewWIP.pl -A HCG.GENE.tab.bedHCH.100KB.hippoonly.UMAP.db.dims -f 5000 -p 100 /home/groups/oroaklab/adey_lab/projects/sciNOMe/refs/mm10_annotations/mm10.RefGene.names.bed.sorted.bed mm10 180428_allCells.HCG.NOMe.sorted.bed

[mulqueen@ma