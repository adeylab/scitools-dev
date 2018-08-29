package sci_commands::matrix_ddrtree;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_ddrtree");

sub matrix_ddrtree {

@ARGV = @_;
# Defaults

getopts("O:X:R:C:c:", \%opt);

$die2 = "
scitools matrix-ddrtree [options] [input matrix] [annotation file] [dims file]
   or    ddrtree-matrix
 Applies ddrt to provided matrix 

Options:
   -O   [STR]   Output prefix (default is [input].ddrt.dims)
   -X           Retain intermediate files (def = delete)
   -R   [STR]   Rscript call (def = $Rscript)
   -C   [STR]   Color coding file (annot (tab) #hexColor)
   -c   [STR]   Color coding string
                  Annot=#hexColor,Annot2=#hexColor
                  
Note: Requires monocle2 and cicero R packages so dependencies need to look for those
      This works specifically with monocle2. Will be upgraded once monocle 3 is more stable. Test
";


if (!defined $ARGV[0]) {die $die2};
if (!defined $ARGV[1]) {die $die2};
if (!defined $ARGV[2]) {die $die2};

if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};
if (defined $opt{'C'} && defined $opt{'c'}) {die "\nSpecify either a color string (-c) or a color coding file (-C), not both!\n$die2"};
if (defined $opt{'C'}) {read_color_file($opt{'C'})};
if (defined $opt{'c'}) {read_color_string($opt{'c'})};


open SITE_DATA, ">./cds_site_data.txt";
open CELL_DATA, ">./cds_cell_data.txt";
open DIMS_DATA, ">./cds_dims_data.txt";
open COUNTS, ">./cds_counts_matrix.txt";
open LAMBDA_DIST, ">./Lambda_dist.txt";

read_annot($ARGV[1]);
read_dims($ARGV[2]);



#make sure annot matches matrix
open IN, "$ARGV[0]";
$h = <IN>; chomp $h; @H = split(/\t/, $h);
$h_out = "";
for ($i = 0; $i < @H; $i++) {
	if ((defined $CELLID_annot{$H[$i]}) && defined ($CELLID_DIMS{$H[$i]})) {
		$h_out .= "$H[$i]\t";
	}
}
#finish printing out header by removing last \t
$h_out =~ s/\t$//;
print COUNTS "$h_out\n";
#print out header

print SITE_DATA "site_name\tchr\tbp1\tbp2\tnum_cells_expressed\tsite_length\n";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$site = shift(@P);
	($chr,$bp1,$bp2) = split(/[_]/, $site);
	$siteName = "$chr\_$bp1\_$bp2";
	$siteOut = "$siteName";
	$siteOut2 = "$siteName";
	$chrID = $chr; $chrID =~ s/chr//;
	my $SITENAME_maxSignal=0;
	for ($i = 0; $i < @P; $i++) {
		if ((defined $CELLID_annot{$H[$i]}) && defined ($CELLID_DIMS{$H[$i]})) {
			$SITENAME_totalSignal{$siteName}+=$P[$i];
			
			if ($P[$i] > $SITENAME_maxSignal{$siteName}){
				$SITENAME_maxSignal{$siteName}=$P[$i];
				#print $siteName."\t$P[$i]\t".$SITENAME_maxSignal{$siteName}."\n";
			}
			
			if ($P[$i]>0) {
				$SITENAME_expressed{$siteName}++;
			}
			
			
			
			if ($P[$i]>0) {
				$CELLID_expressed{$H[$i]}++;
			}
			$siteOut .= "\t$P[$i]";
		}
	}
	print COUNTS "$siteOut\n";
	
	
	for ($i = 0; $i < @P; $i++) 
		{
			if ((defined $CELLID_annot{$H[$i]}) && defined ($CELLID_DIMS{$H[$i]})) 
			{
		
			if ($P[$i]>0) {
				$tempval = $P[$i]/$SITENAME_maxSignal{$siteName};
				#print $siteName."\t$P[$i]\t$SITENAME_maxSignal{$siteName}\t".$tempval."\n";
				$siteOut2 .= "\t$tempval";
			}
			else
			{
				$siteOut2 .= "\t0";
			}
			}
		
		}
	
	#I THINK YOU NEED THE NUMBER OF CELLS EXPRESSED PER SITE HERE NOT THE TOTAL SIGNAL
	#OLD: Have double checked! it works!
	#print SITE_DATA "$siteName\t$siteName\t$chrID\t$bp1\t$bp2\t$SITENAME_totalSignal{$siteName}\t".($bp2-$bp1)."\n";
	print SITE_DATA "$siteName\t$siteName\t$chrID\t$bp1\t$bp2\t$SITENAME_expressed{$siteName}\t".($bp2-$bp1)."\n";
	
} close IN;
close COUNTS; close SITE_DATA; 

print CELL_DATA "cells\ttimepoint\tnum_genes_expressed\n";
print DIMS_DATA "dimension_1\tdimension_2\n";

for ($i = 0; $i < @H; $i++) {
	if ((defined $CELLID_annot{$H[$i]}) && defined ($CELLID_DIMS{$H[$i]})) {
		print DIMS_DATA  join("\t",@{$CELLID_DIMS{$H[$i]}})."\n";
		print CELL_DATA "$H[$i]\t$H[$i]\t$CELLID_annot{$H[$i]}\t$CELLID_expressed{$H[$i]}\n";
	}
} close CELL_DATA; close DIMS_DATA;



open R, ">$opt{'O'}.ddrt.r";

print R "
library(monocle)
library(cicero)
# reading in matrix, annotation, and dimension data
cds_cell_data <- read.delim(\"./cds_cell_data.txt\")
cds_dims_data <- read.delim(\"./cds_dims_data.txt\")
cds_site_data <- read.delim(\"./cds_site_data.txt\")
cds_counts_matrix <- read.table(\"./cds_counts_matrix.txt\")

#make it binary this is just for cicero
#threshold <- 1
#bcdata <- ifelse(cds_counts_matrix < threshold, 0, 1)
#cds_counts_matrix<-bcdata

#check if these two are in the same order and combine for mean

rownames(cds_cell_data)==rownames(cds_dims_data)

cds_cell_data<-cbind(cds_cell_data,cds_dims_data)



# generating cell data set
feature_data <- new(\"AnnotatedDataFrame\", data = cds_site_data)
sample_data <- new(\"AnnotatedDataFrame\", data = cds_cell_data)
dimensions_data <- new(\"AnnotatedDataFrame\", data = cds_dims_data)
input_cds <- newCellDataSet(as.matrix(cds_counts_matrix), phenoData = sample_data, featureData = feature_data)
input_cds\@expressionFamily <- binomialff()
input_cds\@expressionFamily\@vfamily <- \"binomialff\"



set.seed(2017)
";


print R "# add cell data
	

pData(input_cds)\$cells <- NULL	
agg_cds <- aggregate_nearby_peaks(input_cds, distance = 10000)
agg_cds <- detectGenes(agg_cds)
agg_cds <- estimateSizeFactors(agg_cds)
agg_cds <- estimateDispersions(agg_cds)

#from here it is good
dimA<-t(as.matrix(cds_dims_data,rownames=F))
rownames(dimA)<-NULL
reducedDimA(agg_cds)<-dimA


#instead of this we can potentially add in the aggregate group kmers. Ok for now
agg_cds <- clusterCells(agg_cds, verbose = F,cores=10)

clustering_DA_sites <- differentialGeneTest(agg_cds,fullModelFormulaStr = '~Cluster')

#might not need to use this
# This takes a few minutes to run
#diff_timepoint <- differentialGeneTest(agg_cds,
#                  fullModelFormulaStr=\"~timepoint + num_genes_expressed\")

ordering_sites <- row.names(clustering_DA_sites)[order(clustering_DA_sites\$qval)][1:10000]

agg_cds <- setOrderingFilter(agg_cds, ordering_sites)

agg_cds <- reduceDimension(agg_cds, max_components = 2,
          residualModelFormulaStr=\"~num_genes_expressed\",
          reduction_method = \'DDRTree\')
agg_cds <- orderCells(agg_cds)





p<-plot_cell_trajectory(agg_cds, color_by = \"timepoint\")



";

if ($color_mapping !~ /none/i) {
	print R "
	p<-p+scale_colour_manual(values = c($color_mapping))";
}



print R "

ggsave(plot=p,filename=\"$opt{'O'}.annot_plot.png\",width=5,height=4,dpi=900)
ggsave(plot=p,filename=\"$opt{'O'}.annot_plot.pdf\",width=5,height=4);

p<-plot_cell_trajectory(agg_cds, color_by = \"State\")
ggsave(plot=p,filename=\"$opt{'O'}.state_plot.png\",width=5,height=4,dpi=900)
ggsave(plot=p,filename=\"$opt{'O'}.state_plot.pdf\",width=5,height=4);


cds<-agg_cds
#Overwrite monocle.CDS file with final analysis
saveRDS(cds,file=\"monocle.CDS.rds\")


#Save branch point coordinates
dp_mst <- minSpanningTree(cds)
reduced_dim_coords <- reducedDimK(cds)
ica_space_df <- data.frame(Matrix::t(reduced_dim_coords[1:2,]))
colnames(ica_space_df) <- c(\"prin_graph_dim_1\", \"prin_graph_dim_2\")
ica_space_df\$sample_name <- row.names(ica_space_df)
ica_space_df\$sample_state <- row.names(ica_space_df)
#edge_list <- as.data.frame(get.edgelist(dp_mst))
#colnames(edge_list) <- c(\"source\", \"target\")
#edge_df <- merge(ica_space_df, edge_list, by.x = \"sample_name\",by.y = \"source\", all = TRUE)
#edge_df <- plyr::rename(edge_df, c(prin_graph_dim_1 = \"source_prin_graph_dim_1\",prin_graph_dim_2 = \"source_prin_graph_dim_2\"))
#edge_df <- merge(edge_df, ica_space_df[, c(\"sample_name\",\"prin_graph_dim_1\", \"prin_graph_dim_2\")],by.x = \"target\", by.y = \"sample_name\", all = TRUE)
#edge_df <- plyr::rename(edge_df, c(prin_graph_dim_1 = \"target_prin_graph_dim_1\",prin_graph_dim_2 = \"target_prin_graph_dim_2\"))
#write.table(as.matrix(edge_df),file=\"$opt{'O'}/monocle_branchpoints.txt\",col.names=TRUE,row.names=FALSE,sep=\"\\t\",quote=FALSE)

p<-plot_cell_trajectory(agg_cds, color_by = \"timepoint\")
ggsave(plot=p,filename=\"$opt{'O'}.timepoint_plot.png\",width=5,height=4,dpi=900)
ggsave(plot=p,filename=\"$opt{'O'}.timepoint_plot.pdf\",width=5,height=4);




#add stuff here if you want to reroot agg_cds <- orderCells(agg_cds, root_state = \"D\")

p<-plot_cell_trajectory(agg_cds, color_by = \"Pseudotime\")
ggsave(plot=p,filename=\"$opt{'O'}.lambda_plot.png\",width=5,height=4,dpi=900)
ggsave(plot=p,filename=\"$opt{'O'}.lambda_plot.pdf\",width=5,height=4);

pData(input_cds)\$Pseudotime <- pData(agg_cds)[colnames(input_cds),]\$Pseudotime
pData(input_cds)\$State <- pData(agg_cds)[colnames(input_cds),]\$State


# writing out pseudotime etc so you can recreate everything
write.table(as.matrix(pData(agg_cds)),file=\"./$opt{'O'}_ddrt_aggragated_cells.txt\",col.names=TRUE,row.names=FALSE,sep=\"\\t\",quote=FALSE)
write.table(as.matrix(Biobase::exprs(agg_cds)),file=\"./$opt{'O'}_ddrt_aggragated_cells_norm_count.txt\",col.names=TRUE,row.names=TRUE,sep=\"\\t\",quote=FALSE)
write.table(as.matrix(fData(agg_cds)),file=\"./$opt{'O'}_ddrt_aggragated_features.txt\",col.names=TRUE,row.names=FALSE,sep=\"\\t\",quote=FALSE)
";



close R;

system("$Rscript $opt{'O'}.ddrt.r");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.ddrt.r");
}
}
1;
