package sci_commands::atac_refcount;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("atac_refcount");

sub atac_refcount {

@ARGV = @_;

# Defaults:
$min_feature_size = 500;

getopts("O:b:s:X", \%opt);

$die2 = "
scitools atac-refcount [options] [Reference Genome or bed annotation] [duplicate removed and filtered bam file]

This script will take in a [duplicate removed and filtered bam file] and merge it with
known, annotated regions to generate a counts matrix. This is meant to serve as an alternative
to calling peaks de novo.

If [Reference Genome] is specified for argument 1, then will use promoter regions (2kb up from TSS, 500bp down) for that genome.
Will make a matrix of all cellID (columns) by all promoter regions (rows) with each value being read-count overlap.

If [Bed Annotation] is specified for argument 1, will make a matrix of all cellID (columns) 
by all unique annotation ids (rows) with each value being read-count overlap. 
MUST BE OF FORMAT:
<CHR><START><END><UNIQUE_ANNOTATION_ID>

[Reference Genome of bed annotation] shortcuts:
hg38  =   /home/groups/oroaklab/refs/hg38/refseq_tss_2kbup_1kbdown.filtered.txt
mm10  =  /home/groups/oroaklab/refs/mm10/mm10.RefGene.promoters_neg2000_pos500.bed.sorted.bed

Options:
   -O   [STR]   Output prefix (default is bam file prefix)
   -b   [STR]   Bedtools call (def = $bedtools)
   -s   [STR]   Samtools call (def = $samtools)
   -X           Retain intermediate files (def = remove)



";

if (!defined $ARGV[1]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.bam$//};
if (defined $opt{'s'}) {$samtools = $opt{'s'}};
if (defined $opt{'b'}) {$bedtools = $opt{'b'}};

if ($ARGV[0] eq "hg38") {$ARGV[0]="/home/groups/oroaklab/refs/hg38/refseq_tss_2kbup_1kbdown.filtered.txt"};
if ($ARGV[0] eq "mm10") {$ARGV[0]="/home/groups/oroaklab/refs/mm10/mm10.RefGene.promoters_neg2000_pos500.bed.sorted.bed"};

system("bedtools intersect -a $ARGV[0] -b $ARGV[1] -wa -wb | awk 'OFS=\"\\t\"{split(\$8,a,\":\"); print \$4,a[1]}' > $opt{'O'}.refcount.tmp");

open R, "$opt{'O'}.refcount.R";

print R
"
library(data.table)
  dat<-read.table(\"$opt{'O'}.refcount.tmp\")
  colnames(dat)<-c(\"ANNO\",\"CELLID\")
  dat_cast<-dcast(dat,ANNO~CELLID,fun.aggregate=length,value.var=\"CELLID\",fill=0)
  row.names(dat_cast)<-dat_cast\$ANNO
  dat_cast<-dat_cast[2:ncol(dat_cast)]
  write.table(dat_cast,\"$opt{'O'}.refcounts.matrix\",sep=\"\\t\",col.names=T,row.names=T,quote=F)
";
close R;

if (!defined $opt{'X'}) {
  system("rm -f $opt{'O'}.refcount.tmp $opt{'O'}.refcount.R");
}

}
1;