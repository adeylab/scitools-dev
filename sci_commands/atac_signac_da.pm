package sci_commands::atac_signac_da;


use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("atac_signac_da");

sub atac_signac_da {

@ARGV = @_;
use Getopt::Std; %opt = ();
#190425 RM Correction: Some option flags listed by unused.
getopts("O:U:N:", \%opt);

$die2 = "
scitools atac-signac-da [options] [counts matrix] [reference genome] [annotation file for grouping cells]
   or    signac_da,da_signac

This script will perform differential accessibility analysis on:
[counts matrix]:                       a counts matrix generated through atac-count command call.
[reference genome]:                    Shorthand reference genome of Signac processing. Accepts one of the following: [hg38, hg19, mm10] 
[annotation file for grouping cells]:  an annotation file for grouping cells (example: matrix-pg output). Cell names must be consistent with the counts matrix.

This scitools command wraps the new \"Signac\" package relased from the Satija lab. https://satijalab.org/signac/index.html
It will use the counts matrix to generate a Seurat object.
Append the annotation file provided for grouping cells.
Use LR (Logistic Regression) for calling differential peaks between cell groups in a one-by-all others manner.

It outputs:
  1. a tab-separated data frame with the following columns:

row.name p_val avg_logFC pct.1 pct.2 p_val_adj da_region gene_id gene_name gene_biotype seq_coord_system symbol entrezid region distance enriched_cells
See signac functions FindMarkers and ClosestFeature for a full description of columns.
The annotation group of interest (that being compared to all other cells) is listed in the \"enriched cells\" column.

  2.Test PDFs showing the top N (default: 5) DA-defined peaks per annotation group as a sanity check.
  3. A SeuratObject to further processing of your data.

WIP: Also going to add an integrated rGREAT analysis for the DA peaks per annotation group (Not yet implemented).

Options:
   -O   [STR]     Output prefix (default is [input annot].da_peaks.txt etc.)
   -U   [STR]     A 2-D umap coordinate .dims file (output from scitools matrix-umap). If given, will plot accessibity per cell over UMAP coordinates.
   -N   [INT]     Integer number of example peaks per annotation grouping to show. (Default: 5)
                  This does not affect DA tabulation, just the sanity check PDFs.

***NOTE: THIS SCRIPT CHANGES USER R LIBRARY LOOKUP TO /home/groups/oroaklab/src/R/R-3.5.1/library_arsn. TO RETURN TO NORMAL USE A NEW SHELL SESSION.***
 
";

#name output and create folder 
if (!defined $ARGV[1]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[2]; $opt{'O'} =~ s/\.annot$//};
if (!defined $opt{'N'}) {$opt{'N'} = 5};

if ($ARGV[1]=="hg38") {$ref="EnsDb.Hsapiens.v86"} if ($ARGV[1]=="hg19") {$ref="EnsDb.Hsapiens.v75"}if ($ARGV[1]=="mm10") {$ref="EnsDb.Mmusculus.v79"} else {die $die2};
system("export R_LIBS_USER=\'/home/groups/oroaklab/src/R/R-3.5.1/library_arsn\'");


print "Generating R script for analysis."."\n";        
open R, ">$opt{'O'}.da_peaks.r";
print R "
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
set.seed(1234)
library($ref)

write(\"Read in R Libraries.\", stderr())
write(\"Using Reference Genome: $ref\", stderr())

#read in counts matrix
counts <- read.table(\"$ARGV[0]\",header=T,row.names=1)
write(\"Read in counts matrix: $ARGV[0]\", stderr())

#read in annotation file for grouping
metadata <- read.table(\"$ARGV[2]\",header=F)
colnames(metadata)<-c(\"cellID\",\"pg_clus\")
row.names(metadata)<-metadata\$cellID

write(\"Read in Annot File for DA Groupings: $ARGV[2]\", stderr())

#Generate Seurat Object
write(\"Generating Seurat Object\", stderr())

dat <- CreateSeuratObject(
  counts = counts,
  assay = \'peaks\',
  project = \'ATAC\',
  meta.data = metadata,
  min.cells = 1
)";

if ($opt{'U'})
{
print R "

#Add UMAP Dims to Seurat Object
write(\"Read in UMAP Dims\", stderr())

umap<-read.table(\"$opt{'U'}\",row.names=1)
umap <- as.matrix(umap[,1:2])
colnames(umap) <- paste0(\"UMAP_\", 1:2)
dat\@reductions[[\"umap\"]] <- CreateDimReducObject(embeddings = umap, key = \"UMAP_\", assay = DefaultAssay(dat))

";
};

print R "
#Generate Gene ranges to downstream analysis
write(\"Generating Gene Ranges for reference genome supplied.\", stderr())
gene.ranges <- genes($ref)

tss.ranges <- GRanges(
  seqnames = seqnames(gene.ranges),
  ranges = IRanges(start = start(gene.ranges), width = 2),
  strand = strand(gene.ranges)
)

#Perform One vs. rest DA enrichment
write(\"Performing one vs. rest DA enrichment per annotation grouping supplied.\", stderr())

da_peaks<-list() 
for (i in unique(dat\@meta.data\$pg_clus)){
  write(paste0(\"Performing one vs. rest DA enrichment for:\",as.character(i)), stderr())
  da_peaks_tmp <- FindMarkers(object = dat,ident.1=i, group.by = \"pg_clus\",min.pct = 0.2,test.use = \'LR\',only.pos=T)
  da_peaks_tmp\$da_region<-row.names(da_peaks_tmp)
  closest_genes <- ClosestFeature(regions = row.names(da_peaks_tmp),annotation = $ref,sep = c(\':\', \'-\'))
  da_peaks_tmp<-cbind(da_peaks_tmp,closest_genes)
  da_peaks[[i]]<-da_peaks_tmp
  }

write(\"Formatting DA Table for output.\", stderr())

da_peaks<-do.call(\"rbind\",da_peaks)
da_peaks\$enriched_cells<-unlist(lapply(strsplit(row.names(da_peaks),\"[.]\"),\"[\",1))

write(\"Outputting DA Table.\", stderr())
write.table(da_peaks,file=\"$opt{'O'}.da_peaks.txt\",sep=\"\\t\",col.names=T,row.names=T,quote=F)


write(\"Outputting Seurat Object.\", stderr())
saveRDS(dat,file=\"$opt{'O'}.SeuratObject.Rds\")

write(\"Generating a sanity check PDF for enrichment. Using Top $opt{'N'} sites per annotation.\", stderr())

#grab top N peaks per cell grouping
N<-$opt{'N'}
da_peaks_topN<-do.call(\"rbind\",lapply(unique(da_peaks\$enriched_cells),function(x) head(da_peaks[da_peaks\$enriched_cells==x,][order(da_peaks\$p_val),],N)))
da_peaks_topN<-unlist(lapply(strsplit(row.names(da_peaks_topN),\"[.]\"),\"[\",2))


pdf(\"$opt{'O'}.da_peaks.violinplot.pdf\",height=100,width=100)
VlnPlot(
  object = dat,
  features = da_peaks_topN,
  ncol = N,
  group.by=\"pg_clus\"
)
dev.off()
";
if ($opt{'U'}){
print R "
pdf(\"$opt{'O'}.da_peaks.featureplot.pdf\",height=100,width=100)
FeaturePlot(
  object = dat,
  features = da_peaks_topN,
  ncol = N,
  pt.size = 4,
  max.cutoff=1
)
dev.off()
";
};

close(R);

system("Rscript $opt{'O'}.da_peaks.r");


}
1;
