package sci_commands::bam_tssenrich;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("bam_tssenrich");

sub bam_tssenrich {

    @ARGV = @_;

# Defaults:
    $bedtools = "bedtools";

    getopts("O:T:B:b:S:r:n:s:EX", \%opt);

$die2 = "
scitools bam_tssenrich [options] [duplicate removed and filtered bam file] [genome -- either 'hg38' or 'mm10' or 'dm6']

Options:
   -O   [STR]   Output prefix (default is bam file prefix)
   -T   [STR]   An alternative BED file for TSS signal (def = reference signal file for genome)
   -B   [STR]   An alternative BED file for background (def = reference background file for genome)
   -b   [STR]   Bedtools call (def = $bedtools)
   -E           Run using ENCODE ATAC-seq definitions to get a bulk signal
   -S   [STR]   An alternate BED file for unique TSS locations (def = reference tss file for genome)
   -r   [INT]   Range in bp for ENCODE ATAC-seq bulk signal (def = 1000)
   -n   [INT]   Bin in bp size for ENCODE ATAC-seq bulk signal (def = 5)
   -s   [STR]   Selection criteria for ENCODE ATAC-seq bulk signal. MID for middle constrained or MAX for max in range (def = MID)
   -X           Retain intermediate files, this includes individual value files for the TSS and background (def = remove)

";

    if (!defined $ARGV[0]) {die $die2};

    if (!defined $ARGV[1]) {die $die2};
    if ($ARGV[1] eq "hg38") {
        $tss_list = "/home/groups/oroaklab/refs/hg38/ensembl_tss/ensembl.hg38.uniqtss.chr.bed";
	$tss_signal = "/home/groups/oroaklab/refs/hg38/ensembl_tss/ensembl.hg38.tss.chr100bpWINDOW.bed";
        $bg_signal = "/home/groups/oroaklab/refs/hg38/ensembl_tss/ensembl.hg38.tss.chr.b1_b2.bed";
    } elsif ($ARGV[1] eq "mm10") {
        $tss_list = "/home/groups/oroaklab/refs/mm10/ensembl_tss/ensembl.mm10.uniqtss.chr.bed";
	$tss_signal = "/home/groups/oroaklab/refs/mm10/ensembl_tss/ensembl.mm10.tss.chr100bpWINDOW.bed";
        $bg_signal = "/home/groups/oroaklab/refs/mm10/ensembl_tss/ensembl.mm10.tss.chr.b1_b2.bed";
    } elsif ($ARGV[1] eq "dm6") {
        $tss_list = "/home/groups/oroaklab/refs/dm6/ensembl_tss/ensembl.hg38.uniqtss.chr.bed";
	$tss_signal = "/home/groups/oroaklab/refs/dm6/ensembl_tss/ensembl.dm6.tss.chr100bpWINDOW.bed";
        $bg_signal = "/home/groups/oroaklab/refs/dm6/ensembl_tss/ensembl.dm6.tss.chr.b1_b2.bed";
    } else {
        die $die2;
    }
    if (defined $opt{'S'}) {$tss_list = $opt{'S'}};
    if (defined $opt{'T'}) {$tss_signal = $opt{'T'}};
    if (defined $opt{'B'}) {$bg_signal = $opt{'B'}};

    if (defined $opt{'r'}) {$bulk_range = $opt{'r'};} else {$bulk_range = 1000;}
    if (defined $opt{'n'}) {$bulk_binsize = $opt{'n'};} else {$bulk_binsize = 5;}
    if (defined $opt{'s'}) {
	if ($opt{'s'} eq "MID" || $opt{'s'} eq "MAX") {
	    $bulk_select = $opt{'s'};
	} else {
	    die $die2;
	}
    } else {
	$bulk_select = "MID";
    }

    if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.bam$//};
    if (defined $opt{'b'}) {$bedtools = $opt{'b'}};

    if (defined $opt{'E'}) {
	system("cat $tss_list | awk -v range='$bulk_range' '\$2 > range+2' | awk -v range='$bulk_range' '{print \$1,\$2-range,\$3+range}' | tr ' ' '\t' > $opt{'O'}.temp_tss_range.bed");
	system("bedtools bamtobed -i $ARGV[0] | awk '{print \$1,\$2,\$2+1}' | bedops -e 1 - $opt{'O'}.temp_tss_range.bed | closest-features --closest --dist --delim \'\t\' - $tss_list > $opt{'O'}.temp.tss_range.closestfeatures");
	system("cat $opt{'O'}.temp.tss_range.closestfeatures | cut -f 7 | sort | uniq -c | sed -e 's/^[ ]*//g' | awk '{print \$2,\$1}' | tr ' ' '\t' | sort -n > $opt{'O'}.bulkTSSdist_$bulk_range.txt");

	open R, ">$opt{'O'}.bulkTSSenrich.r";
	print R "
library(ggplot2)
rangeR <- read.table(\"$opt{'O'}.bulkTSSdist_$bulk_range.txt\")
bins <- seq(from=min(as.numeric(rangeR\$V1)),to=max(as.numeric(rangeR\$V1)),by=as.numeric($bulk_binsize))

all <- rep(rangeR\$V1,rangeR\$V2)
groups <- table(cut(all,bins))

num_endbins <- 100/$bulk_binsize
average_endbins <- (sum(head(groups,n=num_endbins))+sum(tail(groups,n=num_endbins)))/(num_endbins*2)
num_middlebin <- round(length(bins)/2)
a <- num_middlebin-num_endbins
b <- num_middlebin+num_endbins

";

	if ($bulk_select eq "MID") {
	    print R "
TSS <- max(groups[a:b]/average_endbins)
TSSwhich <- which.max(groups[a:b]/average_endbins)
";
	} else {
	    print R "
TSS <- max(groups/average_endbins)
TSSwhich <- which.max(groups/average_endbins)
";
	}
	print R "
df <- data.frame(TSS)
rownames(df) <- \"bulk_TSS_enrichment\"
write.table(df,\"$opt{'O'}.bulkTSSenrich.log\",row.names=T,col.names=F,quote=F,sep=\"\t\",)

mainstr=paste(\"Binned Hist of Read Distances from TSS\\nTSS = \",round(TSS,5),\" with a binsize of $bulk_binsize at bin \",names(TSSwhich),sep=\"\")
binpos <- as.numeric(gsub(\"]\",\"\",gsub(\".*,\",\"\",names(TSSwhich))))
results <- hist(all,breaks=length(groups),plot=FALSE)
results_df <- data.frame(counts=as.numeric(results\$counts)/1000, mids=as.numeric(results\$mids))

png(\"$opt{'O'}.bulkTSSenrich.png\", width = 7, height = 4, units = \"in\", res = 300)
p <- ggplot(results_df, aes(x=mids,y=counts)) + geom_bar(stat=\"identity\",color=\"lightgrey\")
p + geom_vline(aes(xintercept=binpos), color=\"mediumorchid\", linetype=\"dashed\", size=1) + 
  geom_vline(aes(xintercept=0), color=\"black\", size=.1) + 
  theme_bw() + theme(panel.background = element_blank()) + 
  ggtitle(mainstr) + xlab(\"Distance to TSS (bp)\") + ylab(\"Counts in Thousands\") +
  theme(plot.title = element_text(face=\"bold\",size=14))
dev.off()
";

	close R;
	system("$Rscript $opt{'O'}.bulkTSSenrich.r");
	if (!defined $opt{'X'}) {
	    system("rm -f $opt{'O'}.temp.tss_range.closestfeatures $opt{'O'}.temp_tss_range.bed $opt{'O'}.bulkTSSenrich.r");
	}

    } else { 
        system("$bedtools intersect -bed -u -a $ARGV[0] -b $tss_signal | cut -f 4 | sed -e 's/:.*//g' | sort | uniq -c | sed -e 's/^[ ]*//g' | awk '{print \$2,\$1}' | tr ' ' '\t' > $opt{'O'}.tss_reads.value");
        system("$bedtools intersect -bed -u -a $ARGV[0] -b $bg_signal | cut -f 4 | sed -e 's/:.*//g' | sort | uniq -c | sed -e 's/^[ ]*//g' | awk '{print \$2,\$1}' | tr ' ' '\t' > $opt{'O'}.bg_reads.value");
        system("join -e \"0\" -a1 -a2 -o \"0,1.2,2.2\" $opt{'O'}.tss_reads.value $opt{'O'}.bg_reads.value | awk '{print \$1,\$2/(\$3+1)}' | tr ' ' '\t' > $opt{'O'}.TSSenrich.value");
        system("cat $opt{'O'}.TSSenrich.value | awk '{sum+=\$2}END{print \"Average TSS enrichment: \" sum/NR}' > $opt{'O'}.TSSenrich.log");
        if (!defined $opt{'X'}) {
            system("rm -f $opt{'O'}.tss_reads.value $opt{'O'}.bg_reads.value");
        }
    }

}

1;
