package sci_commands::cds_monocle;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("cds_monocle");

sub cds_monocle {

@ARGV = @_;
# Defaults

getopts("O:R:XP:D:L:", \%opt);

$die2 = "
scitools cds_monocle [options] [directory containing cds files]
   or    monocle_cds

Prior to running, ensure you have ran scitools matrix-makecds. That function will convert matrix files into the CDS format that Monocle3 requires.
Runs the current version of Monocle3 within a directory containing files for CDS format input. 
WIP: Current version performs re-clustering using the monocle pipeline, will adjust to read in scitools generated dims file in the future.


Options:
   -O   [STR]   Output Directory (default is [current working directory]/monocle_output)
   -R   [STR]   Rscript call (def = $Rscript)
   -D 	[2|3] 	Dimensions to be used for final plotting (2D or 3D plotting)
   				Default: 2 Dimensions
   -P 	[INT]	Number of components to be used for dimensionality reduction by PCA to denoise.
   				Default: 50 components 
   -L 	[STR] 	RGE method for branch analysis in monocle3. Must member of the following list:
   				[SimplePPT,L1graph,DDRTree] Default=SimplePPT
   -X           Retain intermediate files (Default = delete)
                  
";


if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = "$ARGV[0]/monocle_output"};
if (!defined $opt{'D'}) {$opt{'D'} = 2};
if (!defined $opt{'P'}) {$opt{'P'} = 50};
if (!defined $opt{'L'}) {$opt{'L'} = "SimplePPT"};

system("mkdir $opt{'O'}");

open R, ">$opt{'O'}/monocle.r";

print R "
suppressWarnings(library(monocle))
message(\"Loading Monocle3\")
# reading in matrix, annotation, and dimension data

cds_cell_data <- read.delim(\"$ARGV[0]/cds_cell_data.txt\")
#cds_dims_data <- read.delim(\"$ARGV[0]/cds_dims_data.txt\")
cds_site_data <- read.delim(\"$ARGV[0]/cds_site_data.txt\")
cds_counts_matrix <- read.table(\"$ARGV[0]/cds_counts_matrix.txt\")
message(\"Read in CDS Files.\")

# To be added with dims data update
#rownames(cds_cell_data)==rownames(cds_dims_data)
#cds_cell_data<-cbind(cds_cell_data,cds_dims_data)

# generating cell data set
feature_data <- new(\"AnnotatedDataFrame\", data = cds_site_data)
sample_data <- new(\"AnnotatedDataFrame\", data = cds_cell_data)

# To be added with dims data update
#dimensions_data <- new(\"AnnotatedDataFrame\", data = cds_dims_data)

message(\"No dims file given, using Monocle3 for dimensionality reduction\")
cds <- suppressWarnings(newCellDataSet(as.matrix(cds_counts_matrix), phenoData = sample_data, featureData = feature_data))
set.seed(2017)
pData(cds)\$cells <- NULL 
cds <- detectGenes(cds)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds<-preprocessCDS(cds,num_dim=$opt{'P'},use_tf_idf=TRUE,verbose=T)
message(\"Writing out PCA components data frame\")
write.table(as.data.frame(cds\@normalized_data_projection), file=\"$opt{'O'}/$opt{'P'}_pca.dims\",quote=F,sep=\"\\t\",row.names=T,col.names=T)

cds <- reduceDimension(cds, max_components = $opt{'D'},reduction_method = \'UMAP\',metric=\"cosine\",verbose = T)
write.table(as.data.frame(t(reducedDimS(cds))), file=\"$opt{'O'}/trajectory_$opt{'D'}D.dims\",quote=F,sep=\"\\t\",row.names=T,col.names=F)
cds <- partitionCells(cds)
";
if ($opt{'D'}==3) {
print R "
#3D Plotting
message(\"Generating 3D Plots\")
cds <- learnGraph(cds, max_components = 3, RGE_method = \"$opt{'L'}\", partition_component = T,verbose = T)
#Writing out full CDS file
saveRDS(cds,file=\"$opt{'O'}/monocle.CDS.rds\")




#Save branch point coordinates
dp_mst <- minSpanningTree(cds)
reduced_dim_coords <- reducedDimK(cds)
ica_space_df <- data.frame(Matrix::t(reduced_dim_coords[1:3,]))
colnames(ica_space_df) <- c(\"prin_graph_dim_1\", \"prin_graph_dim_2\",\"prin_graph_dim_3\")
ica_space_df\$sample_name <- row.names(ica_space_df)
ica_space_df\$sample_state <- row.names(ica_space_df)
edge_list <- as.data.frame(get.edgelist(dp_mst))
colnames(edge_list) <- c(\"source\", \"target\")
edge_df <- merge(ica_space_df, edge_list, by.x = \"sample_name\",by.y = \"source\", all = TRUE)
edge_df <- plyr::rename(edge_df, c(prin_graph_dim_1 = \"source_prin_graph_dim_1\",prin_graph_dim_2 = \"source_prin_graph_dim_2\", prin_graph_dim_3 = \"source_prin_graph_dim_3\"))
edge_df <- merge(edge_df, ica_space_df[, c(\"sample_name\",\"prin_graph_dim_1\", \"prin_graph_dim_2\", \"prin_graph_dim_3\")],by.x = \"target\", by.y = \"sample_name\", all = TRUE)
edge_df <- plyr::rename(edge_df, c(prin_graph_dim_1 = \"target_prin_graph_dim_1\",prin_graph_dim_2 = \"target_prin_graph_dim_2\", prin_graph_dim_3 = \"target_prin_graph_dim_3\"))
write.table(as.matrix(edge_df),file=\"$opt{'O'}/monocle3_branchpoints.txt\",col.names=TRUE,row.names=FALSE,sep=\"\\t\",quote=FALSE)

#Plot out 3D version
plot_3d_cell_trajectory(cds,color_by=paste(colnames(pData(cds))[1]),webGL_filename=\"$opt{'O'}/trajectory_3D.html\",image_filename=\"$opt{'O'}/trajectory_3D.gif\",show_backbone=TRUE,backbone_segment_color=\"#000000\")

";

} else {
print R "
message(\"Generating Plots\")
cds <- learnGraph(cds, max_components = 2, RGE_method = \"$opt{'L'}\", partition_component = T,verbose = T)
#Writing out full CDS file
saveRDS(cds,file=\"$opt{'O'}/monocle.CDS.rds\")

p<-plot_cell_trajectory(cds, color_by = paste(colnames(pData(cds))[1]),backbone_color=\"#000000\")
ggsave(plot=p,filename=paste0(\"$opt{'O'}/monocle3_\",paste(colnames(pData(cds))[1]),\".timepoint_plot.png\"),width=5,height=4,dpi=900)
ggsave(plot=p,filename=\"$opt{'O'}.monocle3.timepoint_plot.pdf\",width=5,height=4);

p<-plot_cell_trajectory(cds, color_by = \"State\",backbone_color=\"#000000\")
ggsave(plot=p,filename=\"$opt{'O'}/monocle3_state_plot.png\",width=5,height=4,dpi=900)
ggsave(plot=p,filename=\"$opt{'O'}/monocle3_state_plot.pdf\",width=5,height=4);

p<-plot_cell_trajectory(cds, color_by = \"Pseudotime\",backbone_color=\"#000000\")
";
}

print R "
#Determine the root state.
orderCells(cds)
pr_graph_test <- principalGraphTest(cds, k=3, cores=10)
diff_access <- dplyr::add_rownames(pr_graph_test) %>% dplyr::arrange(plyr::desc(morans_test_statistic), plyr::desc(-qval))
write.table(as.matrix(diff_access),file=\"$opt{'O'}/monocle3_diffaccess.txt\",col.names=TRUE,row.names=FALSE,sep=\"\\t\",quote=FALSE)

#Overwrite monocle.CDS file with final analysis
saveRDS(cds,file=\"$opt{'O'}/monocle.CDS.rds\")

pData(cds)\$Pseudotime <- pData(cds)[colnames(),]\$Pseudotime
pData(cds)\$State <- pData(cds)[colnames(cds),]\$State

# writing out pseudotime etc so you can recreate everything
write.table(as.matrix(pData(cds)),file=\"$opt{'O'}/monocle3_cells.txt\",col.names=TRUE,row.names=FALSE,sep=\"\\t\",quote=FALSE)
write.table(as.matrix(Biobase::exprs(cds)),file=\"$opt{'O'}/monocle3_aggragated_cells_count.txt\",col.names=TRUE,row.names=TRUE,sep=\"\\t\",quote=FALSE)
write.table(as.matrix(fData(cds)),file=\"$opt{'O'}/monocle3_features.txt\",col.names=TRUE,row.names=FALSE,sep=\"\\t\",quote=FALSE)

";
close R;
system("$Rscript $opt{'O'}/monocle.r");
if (!defined $opt{'X'}) {
    system("rm -f $opt{'O'}/monocle.r");
}
}
1;
