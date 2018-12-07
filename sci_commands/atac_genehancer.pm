package sci_commands::atac_genehancer;
use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("atac_genehancer");
sub atac_genehancer {
@ARGV = @_;
getopts("O:P:Xi", \%opt);
$die2 = "
scitools atac-genehancer [options] [peaks bed file]
   or    genehancer
This script is specific to files using the hg38 reference genome.
Will take a pead bed file, and merge with the GeneHancer data base.
Returns all peaks, but for those which overlap with known GeneHancer enhancers/promoters, will add extra fields.
New output format:
<chr><start><end><GenehancerSite><Sourece><Promoter/Enhancer><GeneHancer IDs><GeneHancer associated genes(comma-separated list)><Genehancer Scores>
Note: if multiple genehancer sites overlap with a peak, will make the tab-separated fields into a comma-separated ordered lists
Options:
   -O   [STR]   Output prefix (default is peaks prefix)
                (adds .GeneHancer.bed)
   -i 	[FLAG]	Will output the inverse of intersect (peaks not within Genehancer)
   -P   [STR]   Python call (def = $Pscript)
   -X   [FLAG]  Remove temp files
";
if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'P'}) {$opt{'P'}=$Pscript};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.bed//};

$bedtools_command="bedtools intersect -a $ARGV[0] -b /home/groups/oroaklab/refs/hg38/hg38_genehancer.txt -wb > $opt{'O'}.genehancer.temp";

if (defined $opt{'i'}){$bedtools_command="bedtools intersect -v -a $ARGV[0] -b /home/groups/oroaklab/refs/hg38/hg38_genehancer.txt -wb > $opt{'O'}.genehancer.inverse.bed";};
system($bedtools_command); 

if (!defined $opt{'i'}){
open py, ">$opt{'O'}.genehancer.py";
print py "

#bedtools intersect call:
#$bedtools_command


import csv
collapsed_dat={}
with open(\"$opt{'O'}.genehancer.temp\") as f:
	reader = csv.reader(f, delimiter=\"\\t\")
	d = list(reader)

for i in d:
	peak=i[0]+\"\\t\"+i[1]+\"\\t\"+i[2]
	genehancer_peak=i[4]+\"_\"+i[5]+\"_\"+i[6]
	genehancer_peak_type=i[7]
	genehancer_peak_conf=i[8]
	info=i[12].split(\";\")
	genehancer_id=info[0].split(\"=\")[1]
	genehancer_genes=[]
	genehancer_scores=[]
	for j in info:
		if j.startswith(\"connected_gene\"):
			genehancer_genes += [j.split(\"=\")[1]]
		if j.startswith(\"score\"):
			genehancer_scores += [j.split(\"=\")[1]]
	if peak not in collapsed_dat:
		collapsed_dat[peak]=[genehancer_peak,genehancer_peak_type,genehancer_peak_conf,genehancer_id,genehancer_genes,genehancer_scores]
	else:
		dat=collapsed_dat[peak]
		collapsed_dat[peak]=[dat[0].append(genehancer_peak),dat[1].append(genehancer_peak_type),dat[2].append(genehance_peak_conf),dat[3].append(genehancer_id),dat[4].append(genehancer_genes),dat[5].append(genehancer_scores)]
with open(\"$opt{'O'}.genehancer.bed\",\"w\") as fout:
	for peak in collapsed_dat:
		fout.write(peak+\"\t\")
		for i in collapsed_dat[peak]:
			if type(i) is list:
				fout.write(\",\".join(map(str, i))+\"\t\")
			else:
				fout.write(i+\"\\t\")
		fout.write(\"\\n\")
";

close py;
system("$opt{'P'} $opt{'O'}.genehancer.py");
if (!defined $opt{'X'}) {
        system("rm -f $opt{'O'}.genehancer.py");
        system("rm -f $opt{'O'}.genehancer.temp");
}
}
};
1;


