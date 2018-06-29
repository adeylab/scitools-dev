package sci_commands::fastq_simulate;

use sci_utils::general;
use sci_utils::modes;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("fastq_simulate");

sub fastq_simulate {

@ARGV = @_;
getopts("O:M:N:m:R:d:l:", \%opt);

# Defaults
$mode_name = "sci";
$num_reads = 100000;
$mut_rate = 1;
$modalityID = 0;
$read_length = 50;
@BASES = ("A", "C", "G", "T");
@BASES_E = @BASES; push @BASES_E, "N";

$die2 = "
scitools fastq-simulate [options] [annotation file]
   or    simulate-fastq

Simulates fastq files for split testing.

Read sequence will be randomers - this tool is currently only
designed for index split and mode testing.

Options:
   -O   [STR]   Output prefix (def = [annot prefix].simulated.[suffix])
   -M   [STR]   Mode (def = $mode_name)
   -N   [INT]   Number of reads to simulate (def = $num_reads)
   -m   [STR]   Run format specification file
         (def = $VAR{'sci_modes'})
   -R   [PCT]   Mutation rate (percentage, def = $mut_rate)
   -I   [STR]   Index files/directory - comma separated (req)
         (def = DIR=$VAR{'index_directory'})
   -d   [INT]   Modality ID (if multimodal mode, def = $modalityID)
   -l   [INT]   Length of 'read' sequence (def = $read_length)

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.annot$//;
if (defined $opt{'M'}) {$mode_name = $opt{'M'}};
if (defined $opt{'N'}) {$num_reads = $opt{'N'}};
if (defined $opt{'R'}) {$mut_rate = $opt{'R'}};
if (defined $opt{'l'}) {$read_length = $opt{'l'}};
if (!defined $opt{'I'}) {$opt{'I'} = "DIR=$VAR{'index_directory'}"};

# load mode
$mode = modes->new($mode_name,$VAR{'sci_modes'});
@MODALITIES = $mode->modalities();
if (!defined $MODALITIES[$modalityID]) {
	die "ERROR: Modality ID $modalityID is specified, but does not exist for mode $mode_name\n";
}
$modality = $MODALITIES[$modalityID];
@INDEXES = $mode->indexes($modality);

# read in indexes
read_indexdir($opt{'I'});
check_indexes();

# read annot file
read_annot($ARGV[0]);
# cellID array for random access
@CELLIDs = (); $cellCT = 0;
foreach $cellID (keys %CELLID_annot) {
	push @CELLIDs, $cellID;
	$cellCT++;
}

# open output files based on mode
open_outs();

# simulate reads based on index positions for specified modality
for ($readID = 0; $readID < $num_reads; $readID++) {
	# pick a random cellID
	$cellID = $CELLIDs[int(rand($cellCT))]; # random enough...
	if ($cellID =~ /-/) { # index names
		@IX_NAMES = split(/-/, $cellID);
	} else { # index seqs
		#########################
	}
	####################
}


#########################################################################################################

close_outs();
}

# SUBROUTINES #

sub check_indexes {
	$absent_index = 0;
	foreach $check_index (@INDEXES) {
		if (!defined $INDEX_TYPE_SEQ_seq{$check_index} && $check_index !~ /^umi$/i) {
			$absent_index++;
		}
		if ($INDEX_TYPE_length{$check_index} != $mode->index_length($modality,$check_index)) {
			die "ERROR: The index type ($check_index) specified has a different length in the mode configuration (".$mode->index_length($modality,$check_index).") than in the specified indexes ($INDEX_TYPE_length{$check_index}).\n";
		}
	}
	if ($absent_index>0) {die "ERROR: Index types specified in mode: $mode_name could not be found in the specified index file(s)/directory(s)!\n"};
}

sub open_outs {
	if ($mode->read_check(read1) eq "true") {open R1, "| gzip $opt{'O'}.simulated.R1.fq.gz"};
	if ($mode->read_check(read2) eq "true") {open R2, "| gzip $opt{'O'}.simulated.R2.fq.gz"};
	if ($mode->read_check(index1) eq "true") {open I1, "| gzip $opt{'O'}.simulated.I1.fq.gz"};
	if ($mode->read_check(index2) eq "true") {open I2, "| gzip $opt{'O'}.simulated.I2.fq.gz"};
}

sub close_outs {
	if ($mode->read_check(read1) eq "true") {close R1};
	if ($mode->read_check(read2) eq "true") {close R2};
	if ($mode->read_check(index1) eq "true") {close I1};
	if ($mode->read_check(index2) eq "true") {close I2};
}

1;