package sci_commands::cds_monocle;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("cds_monocle");

sub cds_monocle {

@ARGV = @_;
# Defaults

getopts("O:R:X:P:D:", \%opt);

$die2 = "
scitools cds_monocle [options] [directory containing cds files]
   or    monocle_cds

Prior to running, ensure you have ran scitools matrix-makecds. That function will convert matrix files into the CDS format that Monocle3 requires.
Runs the current version of Monocle3 within a directory containing files for CDS format input. 
WIP: Current version performs re-clustering using the monocle pipeline, will adjust to read in scitools generated dims file in the future.


Options:
   -O   [STR]   Output Directory (default is [current working directory]/cds_files)
   -R   [STR]   Rscript call (def = $Rscript)
   -D 	[2|3] 	Dimensions to be used for final plotting (2D or 3D plotting)
   				Default: 2 Dimensions
   -P 	[INT]	Number of components to be used for dimensionality reduction by PCA to denoise.
   				Default: 50 components
   -X           Retain intermediate files (Default = delete)
                  
";


if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
if (!defined $opt{'D'}) {$opt{'D'} = 2};
if (!defined $opt{'P'}) {$opt{'P'} = 50};


open R, ">$opt{'O'}.monocle.r";

print R "
library(monocle)
# reading in matrix, annotation, and dimension data

cds_cell_data <- read.delim(\"$ARGV[0]/cds_cell_data.txt\")
#cds_dims_data <- read.delim(\"$ARGV[0]/cds_dims_data.txt\")
cds_site_data <- read.delim(\"$ARGV[0]/cds_site_data.txt\")
cds_counts_matrix <- read.table(\"$ARGV[0]/cds_counts_matrix.txt\")

# To be added with dims data update
#rownames(cds_cell_data)==rownames(cds_dims_data)
#cds_cell_data<-cbind(cds_cell_data,cds_dims_data)

# generating cell data set
feature_data <- new(\"AnnotatedDataFrame\", data = cds_site_data)
sample_data <- new(\"AnnotatedDataFrame\", data = cds_cell_data)

# To be added with dims data update
#dimensions_data <- new(\"AnnotatedDataFrame\", data = cds_dims_data)

cds <- newCellDataSet(as.matrix(cds_counts_matrix), phenoData = sample_data, featureData = feature_data)
set.seed(2017)
pData(cds)$cells <- NULL
cds <- detectGenes(cds)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds<-preprocessCDS(cds,num_dim=$opt{'P'},use_tf_idf=TRUE,verbose=T)
#Write out PCA components data frame
write.table(as.data.frame(cds\@normalized_data_projection), file=\"$opt{'O'}/$opt{'P'}_pca.dims\",quote=F,sep=\"\\t\",row.names=T,col.names=T)

cds <- reduceDimension(cds, max_components = $opt{'D'},reduction_method = \'UMAP\',metric=\"cosine\",verbose = T)
write.table(as.data.frame(t(reducedDimS(cds))), file=\"$opt{'O'}/trajectory_$opt{'D'}D.dims\",quote=F,sep=\"\\t\",row.names=T,col.names=F)
cds <- partitionCells(cds)
";
if ($opt{'D'}==3) {
print R "
#3D Plotting
cds <- learnGraph(cds, max_components = 3, RGE_method = 'SimplePPT', partition_component = T,verbose = T)
plot_3d_cell_trajectory(cds,color_by=\"annot\",webGL_filename=\"$opt{'O'}/trajectory_3D.html\",image_filename=\"$opt{'O'}/trajectory_3D.gif\",show_backbone=TRUE,useNULL_GLdev=TRUE,backbone_segment_color=\"#000000\")
";
} else {
print R "
cds <- learnGraph(cds, max_components = 2, RGE_method = 'SimplePPT', partition_component = T,verbose = T)
p<-plot_cell_trajectory(cds, color_by = \"annot\")
ggsave(plot=p,filename=\"$opt{'O'}.monocle3.timepoint_plot.png\",width=5,height=4,dpi=900)
ggsave(plot=p,filename=\"$opt{'O'}.monocle3.timepoint_plot.pdf\",width=5,height=4);
p<-plot_cell_trajectory(cds, color_by = \"State\")
ggsave(plot=p,filename=\"$opt{'O'}.monocle3.state_plot.png\",width=5,height=4,dpi=900)
ggsave(plot=p,filename=\"$opt{'O'}.monocle3.state_plot.pdf\",width=5,height=4);
#add stuff here if you want to reroot agg_cds <- orderCells(agg_cds, root_state = \"D\")
p<-plot_cell_trajectory(cds, color_by = \"Pseudotime\")
ggsave(plot=p,filename=\"$opt{'O'}.monocle3.lambda_plot.png\",width=5,height=4,dpi=900)
ggsave(plot=p,filename=\"$opt{'O'}.monocle3.lambda_plot.pdf\",width=5,height=4);
pData(input_cds)\$Pseudotime <- pData(cds)[colnames(),]$Pseudotime
pData(input_cds)\$State <- pData(cds)[colnames(cds),]$State
# writing out pseudotime etc so you can recreate everything
write.table(as.matrix(pData(agg_cds)),file=\"$opt{'O'}/monocle3_cells.txt\",col.names=TRUE,row.names=FALSE,sep=\"\\t\",quote=FALSE)
write.table(as.matrix(Biobase::exprs(agg_cds)),file=\"$opt{'O'}/monocle3_aggragated_cells_count.txt\",col.names=TRUE,row.names=TRUE,sep=\"\\t\",quote=FALSE)
write.table(as.matrix(fData(agg_cds)),file=\"$opt{'O'}/monocle3_features.txt\",col.names=TRUE,row.names=FALSE,sep=\"\\t\",quote=FALSE)
";
}
close R;
system("$Rscript $opt{'O'}.monocle.r");
if (!defined $opt{'X'}) {
    system("rm -f $opt{'O'}.monocle.r");
}
}
1;

