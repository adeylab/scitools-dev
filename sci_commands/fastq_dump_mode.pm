package sci_commands::fastq_dump_mode;

use sci_utils::general;
use sci_utils::modes;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("fastq_dump_mode");

sub fastq_dump_mode {

@ARGV = @_;
getopts("R:F:O:o:I:1:2:A:i:j:r:M:u:f:k:m:", \%opt);

# Defaults
$umi_length = 8;
$def_read_name_format = "barc";
$mode_name = "sci";

$die2 = "
scitools fastq-dump [options]
   or    dump-fastq

Takes sequencer fastq files (from bcl2fastq) and will format
them into fastq files with matched barcodes.

Options:
   -R   [STR]   Run name (preferred mode)
   -r   [STR]   Run ID (if additional sequencing for some
                libraries. Will add .[ID] to end of read
                names; optional)
   -M   [STR]   Mode. (def = $mode_name)
   -A   [STR]   Annotation file (will split fastqs)

Defaults:
   -F   [STR]   Fastq directory
         ($VAR{'fastq_input_directory'})
   -O   [STR]   Output fastq directory
         ($VAR{'SCI_fastq_directory'})
   -o   [STR]   Output prefix
         (def = run name)

Index Options:
   -m   [STR]   Run format specification file
         (def = $VAR{'sci_modes'})
   -I   [STR]   Index files/directory - comma separated
         Specify: DIR=[directory] to read all files present
           or   NAME=[file] to include files
         Def: DIR=$VAR{'index_directory'}
         Index file format, tab-delimited:
           col1 = <[PCR/Tn5/RNA/MET]_[Set_letter]_[i5/i7/r1]_[Row/Col_ID/well]>
           col2 = <Index_type>
           col3 = <Sequence>
   -f   [STR]   Read name format: barc OR name (def = $def_read_name_format)
   -k   [STR]   Index types to skip, comma separated
         Will set to a string of X's and every read will match.

To specify specific fastq files instead of defaults:
   (can be comma sep for multiple)
   -1   [STR]   Read 1 fastq
   -2   [STR]   Read 2 fastq
   -i   [STR]   Index 1 fastq
   -j   [STR]   Index 2 fastq

Will split with a hamming distance of 2

";

if (!defined $opt{'R'}) {die $die2};
if (!defined $opt{'F'}) {$opt{'F'} = $VAR{'fastq_input_directory'}};
if (!defined $opt{'O'}) {$opt{'O'} = $VAR{'SCI_fastq_directory'}};
if (!defined $opt{'o'}) {$opt{'o'} = $opt{'R'}};
if (defined $opt{'M'}) {$mode_name = $opt{'M'}};
if (defined $opt{'u'}) {$umi_length = $opt{'u'}};
if (!defined $opt{'I'}) {$opt{'I'} = "DIR=$VAR{'index_directory'}"};
if (!defined $opt{'f'}) {$opt{'f'} = $def_read_name_format};
if (defined $opt{'m'}) {$VAR{'sci_modes'} = $opt{'m'}};

# specify fastq files
if (!defined $opt{'1'}) {
	if (-e "$opt{'F'}/$opt{'R'}/Undetermined_S0_L001_R1_001.fastq.gz") {
		$read1 = "$opt{'F'}/$opt{'R'}/Undetermined_S0_L001_R1_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L002_R1_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L003_R1_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L004_R1_001.fastq.gz";
	} elsif (-e "$opt{'F'}/$opt{'R'}/Undetermined_S0_R1_001.fastq.gz") {
		$read1 = "$opt{'F'}/$opt{'R'}/Undetermined_S0_R1_001.fastq.gz";
	}
} else {
	$read1 = $opt{'1'};
}
if (!defined $opt{'2'}) {
	if (-e "$opt{'F'}/$opt{'R'}/Undetermined_S0_L001_R2_001.fastq.gz") {
		$read2 = "$opt{'F'}/$opt{'R'}/Undetermined_S0_L001_R2_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L002_R2_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L003_R2_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L004_R2_001.fastq.gz";
	} elsif (-e "$opt{'F'}/$opt{'R'}/Undetermined_S0_R2_001.fastq.gz") {
		$read2 = "$opt{'F'}/$opt{'R'}/Undetermined_S0_R2_001.fastq.gz";
	}
} else {
	$read2 = $opt{'2'};
}
if (!defined $opt{'i'}) {
	if (-e "$opt{'F'}/$opt{'R'}/Undetermined_S0_L001_I1_001.fastq.gz") {
		$index1 = "$opt{'F'}/$opt{'R'}/Undetermined_S0_L001_I1_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L002_I1_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L003_I1_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L004_I1_001.fastq.gz";
	} elsif (-e "$opt{'F'}/$opt{'R'}/Undetermined_S0_I1_001.fastq.gz") {
		$index1 = "$opt{'F'}/$opt{'R'}/Undetermined_S0_I1_001.fastq.gz";
	}
} else {
	$index1 = $opt{'i'};
}
if (!defined $opt{'j'}) {
	if (-e "$opt{'F'}/$opt{'R'}/Undetermined_S0_L001_I2_001.fastq.gz") {
		$index2 = "$opt{'F'}/$opt{'R'}/Undetermined_S0_L001_I2_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L002_I2_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L003_I2_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L004_I2_001.fastq.gz";
	} elsif (-e "$opt{'F'}/$opt{'R'}/Undetermined_S0_I2_001.fastq.gz") {
		$index2 = "$opt{'F'}/$opt{'R'}/Undetermined_S0_I2_001.fastq.gz";
	}
} else {
	$index2 = $opt{'j'};
}

if (!defined $read1) {die "ERROR: At a very minimum read1 must be specified or able to be detected.\n"};

# load in sci modes
read_mode($VAR{'sci_modes'});

if (!defined $MODE_GROUP_CLASS_PARTS{$ALIAS_mode{$mode_name}}) {
	die "ERROR: The specified mode or mode alias ($mode_name) cannot be found in the sci_modes config file ($VAR{'sci_modes'})\n";
}

# read in indexes
read_indexdir($opt{'I'});

# verify required index sets present and input fastq files
@OPENED_READS = ();
for ($mode_group = 0; $mode_group < @{$MODE_GROUP_CLASS_PARTS{$mode_name}}; $mode_group++) {
	@{$GROUP_index_set[$mode_group]} = ();
	foreach $item (keys %{$MODE_GROUP_CLASS_PARTS{$mode_name}[$mode_group]}) {
		if ($item ne "name" && $item ne "read_outputs") {
			# check index
			for ($item_position = 0; $item_position < @{$MODE_GROUP_CLASS_PARTS{$mode_name}[$mode_group]{$item}{'name'}}; $item_position++) {
				$index_type = $MODE_GROUP_CLASS_PARTS{$mode_name}[$mode_group]{$item}{'name'}[$item_position];
				if ($index_type ne "umi" && $index_type !~ /^read_/) {
					if (!defined $INDEX_TYPE_SEQ_seq{$index_type}) {
						die "ERROR: For mode: $mode_name, index type: $index_type is required, but was not found in the index file(s) provided ($opt{'I'})\n";
					} else {
						push @{$GROUP_index_set[$mode_group]}, $index_type;
					}
				}
			}
			if ($item eq "read1") {
				if (!defined $read1) {
					die "ERROR: read1 is required in the sci_modes confg file for mode: $mode_name and was not specified or cannot be detected.\n";
				} else {
					if (!defined $OPENED_status{'read1'}) {
						open R1, "$zcat $read1 |" || die "ERROR: Cannot open read 1 fastq files: $read1\n";
						open R1F, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.fail.read1.fq.gz";
						push @OPENED_READS, "read1";
						$OPENED_status{'read1'} = 1;
					}
				}
			} elsif ($item eq "read2") {
				if (!defined $read2) {
					die "ERROR: read2 is required in the sci_modes confg file for mode: $mode_name and was not specified or cannot be detected.\n";
				} else {
					if (!defined $OPENED_status{'read2'}) {
						open R2, "$zcat $read2 |" || die "ERROR: Cannot open read 2 fastq files: $read2\n";
						open R2F, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.fail.read2.fq.gz";
						push @OPENED_READS, "read2";
						$OPENED_status{'read2'} = 1;
					}
				}
			} elsif ($item eq "index1") {
				if (!defined $index1) {
					die "ERROR: index1 is required in the sci_modes confg file for mode: $mode_name and was not specified or cannot be detected.\n";
				} else {
					if (!defined $OPENED_status{'index1'}) {
						open I1, "$zcat $index1 |" || die "ERROR: Cannot open index 1 fastq files: $index1\n";
						open I1F, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.fail.index1.fq.gz";
						push @OPENED_READS, "index1";
						$OPENED_status{'index1'} = 1;
					}
				}
			} elsif ($item eq "index2") {
				if (!defined $index2) {
					die "ERROR: index2 is required in the sci_modes confg file for mode: $mode_name and was not specified or cannot be detected.\n";
				} else {
					if (!defined $OPENED_status{'index2'}) {
						open I2, "$zcat $index2 |" || die "ERROR: Cannot open index 2 fastq files: $index2\n";
						open I2F, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.fail.index2.fq.gz";
						push @OPENED_READS, "index2";
						$OPENED_status{'index2'} = 1;
					}
				}
			}
		}
	}
}

# make all-1-away hash
foreach $index_type (keys %INDEX_TYPE_SEQ_seq) {
	foreach $seq (keys %{$INDEX_TYPE_SEQ_seq{$index_type}}) {
		@TRUE = split(//, $seq);
		for ($i = 0; $i < @TRUE; $i++) {
			foreach $base (@BASES) {
				if ($TRUE[$i] ne $base) {
					@NEW = @TRUE; $NEW[$i] = $base;
					$new = join("", @NEW);
					if (!defined $INDEX_TYPE_SEQ_seq{$index_type}{$new}) {
						$INDEX_TYPE_SEQ_seq{$index_type}{$new} = $seq;
					}
				}
			}
		}
	}
}

# make all-2-away hash
foreach $index_type (keys %INDEX_TYPE_SEQ_seq) {
	foreach $id_seq (keys %{$INDEX_TYPE_SEQ_seq{$index_type}}) {
		$seq = $INDEX_TYPE_SEQ_seq{$index_type}{$id_seq};
		@TRUE = split(//, $seq);
		for ($i = 0; $i < @TRUE; $i++) {
			foreach $base (@BASES) {
				if ($TRUE[$i] ne $base) {
					@NEW = @TRUE; $NEW[$i] = $base;
					$new = join("", @NEW);
					if (!defined $INDEX_TYPE_SEQ_seq{$index_type}{$new}) {
						$INDEX_TYPE_SEQ_seq{$index_type}{$new} = $seq;
					}
				}
			}
		}
	}
} 

# Check in index exclusions
if (defined $opt{'k'}) {
	@INDEX_EXCLUSIONS = split(/,/, $opt{'k'});
	foreach $excluded_type (@INDEX_EXCLUSIONS) {
		if (!defined $INDEX_TYPE_SEQ_seq{$excluded_type}) {
			print STDERR "WARNING: index type: $excluded_type specified, but not found in loaded indexes.\n";
		} else {
			$xmer = "";
			for ($x_pos = 1; $x_pos < $INDEX_TYPE_length{$excluded_type}; $x_pos++) {
				$xmer .= "X";
			}
			$INDEX_EXCLUDE{$excluded_type} = $xmer;
			%{$INDEX_TYPE_SEQ_seq{$excluded_type}} = ();
			%{$INDEX_TYPE_SEQ_id{$excluded_type}} = ();
			$INDEX_TYPE_SEQ_seq{$excluded_type}{$xmer} = $xmer;
			$INDEX_TYPE_SEQ_id{$excluded_type}{$xmer} = "EXCLUDED";
			print STDERR "INFO: index type: $excluded_type will be excluded. Indexes in that position will auto-match to $xmer.\n";
		}
	}
}


$out1_handle = "$annot.1";
	open $out1_handle, "| $gzip > $opt{'O'}.$annot.1.fq.gz";
	$HANDLE1{$annot} = $out1_handle;

# Make output directory
system("mkdir $opt{'O'}/$opt{'o'}");

# Open passing output files
if (defined $opt{'A'}) {
	read_annot($opt{'A'});
	foreach $annot (keys %ANNOT_count) {
		for ($mode_group = 0; $mode_group < @{$MODE_GROUP_CLASS_PARTS{$mode_name}}; $mode_group++) {
			foreach $out_read_name (@{$MODE_GROUP_CLASS_PARTS{$mode_name}[$mode_group]{'read_outputs'}}) {
				$out_handle = $annot."-".$mode_group."-".$out_read_name;
				$mode_group_name = $MODE_GROUP_CLASS_PARTS{$mode_name}[$mode_group]{'name'};
				open $out_handle, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.$annot.$mode_group_name.$out_read_name.fq.gz";
			}
		}
	}
	system("cp $opt{'A'} $opt{'O'}/$opt{'o'}/$opt{'o'}.split.annot");
} else {
	foreach $out_read_name (@{$MODE_GROUP_CLASS_PARTS{$mode_name}[$mode_group]{'read_outputs'}}) {
		for ($mode_group = 0; $mode_group < @{$MODE_GROUP_CLASS_PARTS{$mode_name}}; $mode_group++) {
			$out_handle = $mode_group."-".$out_read_name;
			$mode_group_name = $MODE_GROUP_CLASS_PARTS{$mode_name}[$mode_group]{'name'};
			open $out_handle, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.$mode_group_name.$out_read_name.fq.gz";
		}
	}
}

# intialize counters
$totalCT = 0; $failCT = 0; $rnaCT = 0;
%READ_status = ();

# parse input fastq files
while ($r1tag = <R1>) {
	# set each read needing to be opened
	foreach $opened_read (@OPENED_READS) {
		$READ_status{$opened_read} = 0;
	}
	# parse read 1 as default
	chomp $r1tag; $r1seq = <R1>; $null = <R1>; $r1qual = <R1>; chomp $r1qual;
	$READ_status{'read1'} = 1;
	
	# go through each group
	for ($mode_group = 0; $mode_group < @{$MODE_GROUP_CLASS_PARTS{$mode_name}}; $mode_group++) {
		
		# set the status of each index to 0 (not found), then 1 (found), and 2 (passing)
		%INDEX_status = ();
		foreach $index_type (@{$GROUP_index_set[$mode_group]}) {
			$INDEX_status{$index_type} = 0;
		}
		
		# parse each read & open if necessary for non read1
		foreach $item (keys %{$MODE_GROUP_CLASS_PARTS{$mode_name}[$mode_group]}) {
			if ($item eq "read2") {
				if ($READ_status{'read2'} < 1) {
					$r2tag = <R2>; chomp $r2tag; $r2seq = <R2>; chomp $r2seq; $null = <R2>; $r2qual = <R2>; chomp $r2qual;
					$READ_status{'read2'} = 1;
				}
			}
			######################################################################################################
			if ($item ne "name" && $item ne "read_outputs") {
				for ($read_partID = 0; $read_partID < @{$MODE_GROUP_CLASS_PARTS{$mode_name}[$mode_group]{$item}{'name'}}; $read_partID++) {
					$part_name = $MODE_GROUP_CLASS_PARTS{$mode_name}[$mode_group]{$item}{'name'}[$read_partID];
					$offset = $MODE_GROUP_CLASS_PARTS{$mode_name}[$mode_group]{$item}{'offset'}[$read_partID];
					$length = $MODE_GROUP_CLASS_PARTS{$mode_name}[$mode_group]{$item}{'length'}[$read_partID];
					if ($part_name =~ /^read_/) {
						if ($length =~ /^\d+$/) {
							
						} else {
							
						}
					}
				}
			}
		}
		
	}
	######################################################################################################
	$r2tag = <R2>; chomp $r1tag; chomp $r2tag;
	$r1seq = <R1>; $r2seq = <R2>; chomp $r1seq; chomp $r2seq;
	$null = <R1>; $null = <R2>;
	$r1qual = <R1>; $r2qual = <R2>; chomp $r1qual; chomp $r2qual;
	
	$i1tag = <I1>; chomp $i1tag; $i2tag = <I2>; chomp $i2tag;
	$i1seq = <I1>; chomp $i1seq; $i2seq = <I2>; chomp $i2seq;
	$null = <I1>; $null = <I1>;
	$null = <I2>; $null = <I2>;
	
	$match = 0;
	
	if ($opt{'M'} =~ /^sci$/i) {
	
		$sci7 = substr($i1seq,0,$INDEX_TYPE_length{'sci_tn5_i7'}); if (defined $INDEX_EXCLUDE{'sci_tn5_i7'}) {$sci7 = $INDEX_EXCLUDE{'sci_tn5_i7'}};
		$pcr7 = substr($i1seq,$INDEX_TYPE_length{'sci_tn5_i7'},$INDEX_TYPE_length{'sci_pcr_i7'}); if (defined $INDEX_EXCLUDE{'sci_pcr_i7'}) {$pcr7 = $INDEX_EXCLUDE{'sci_pcr_i7'}};
		$sci5 = substr($i2seq,0,$INDEX_TYPE_length{'sci_tn5_i5'}); if (defined $INDEX_EXCLUDE{'sci_tn5_i5'}) {$sci5 = $INDEX_EXCLUDE{'sci_tn5_i5'}};
		$pcr5 = substr($i2seq,$INDEX_TYPE_length{'sci_tn5_i5'},$INDEX_TYPE_length{'sci_pcr_i5'}); if (defined $INDEX_EXCLUDE{'sci_pcr_i5'}) {$pcr5 = $INDEX_EXCLUDE{'sci_pcr_i5'}};
		
		if (defined $INDEX_TYPE_SEQ_seq{'sci_tn5_i7'}{$sci7} &&
		defined $INDEX_TYPE_SEQ_seq{'sci_pcr_i7'}{$pcr7} &&
		defined $INDEX_TYPE_SEQ_seq{'sci_tn5_i5'}{$sci5} &&
		defined $INDEX_TYPE_SEQ_seq{'sci_pcr_i5'}{$pcr5}) {
			$name = $INDEX_TYPE_SEQ_id{'sci_tn5_i7'}{$INDEX_TYPE_SEQ_seq{'sci_tn5_i7'}{$sci7}}."-".$INDEX_TYPE_SEQ_id{'sci_pcr_i7'}{$INDEX_TYPE_SEQ_seq{'sci_pcr_i7'}{$pcr7}}."-".$INDEX_TYPE_SEQ_id{'sci_tn5_i5'}{$INDEX_TYPE_SEQ_seq{'sci_tn5_i5'}{$sci5}}."-".$INDEX_TYPE_SEQ_id{'sci_pcr_i5'}{$INDEX_TYPE_SEQ_seq{'sci_pcr_i5'}{$pcr5}};
			$barc = $INDEX_TYPE_SEQ_seq{'sci_tn5_i7'}{$sci7}.$INDEX_TYPE_SEQ_seq{'sci_pcr_i7'}{$pcr7}.$INDEX_TYPE_SEQ_seq{'sci_tn5_i5'}{$sci5}.$INDEX_TYPE_SEQ_seq{'sci_pcr_i5'}{$pcr5};
			$match = 1;
		} else {
			$name = "sci_tn5_i7=".$sci7."-sci_pcr_i7=".$pcr7."-sci_tn5_i5=".$sci5."-sci_pcr_i5=".$pcr5;
			$barc = $sci7.$pcr7.$sci5.$pcr5;
		}
		
	} elsif ($opt{'M'} =~ /^rna$/i) {
		
		$umi = substr($r1seq,0,$umi_length);
		$fss = substr($r1seq,$umi_length,$INDEX_TYPE_length{'rna_fss_r1'}); if (defined $INDEX_EXCLUDE{'rna_fss_r1'}) {$fss = $INDEX_EXCLUDE{'rna_fss_r1'}};
		$pcr7 = substr($i1seq,0,$INDEX_TYPE_length{'std_pcr_i7'}); if (defined $INDEX_EXCLUDE{'std_pcr_i7'}) {$pcr7 = $INDEX_EXCLUDE{'std_pcr_i7'}};
		$pcr5 = substr($i2seq,0,$INDEX_TYPE_length{'std_pcr_i5'}); if (defined $INDEX_EXCLUDE{'std_pcr_i5'}) {$pcr5 = $INDEX_EXCLUDE{'std_pcr_i5'}};
		
		if (defined $INDEX_TYPE_SEQ_seq{'rna_fss_r1'}{$fss} &&
		defined $INDEX_TYPE_SEQ_seq{'std_pcr_i7'}{$pcr7} &&
		defined $INDEX_TYPE_SEQ_seq{'std_pcr_i5'}{$pcr5}) {
			$name = $INDEX_TYPE_SEQ_id{'rna_fss_r1'}{$INDEX_TYPE_SEQ_seq{'rna_fss_r1'}{$fss}}."-".$INDEX_TYPE_SEQ_id{'std_pcr_i7'}{$INDEX_TYPE_SEQ_seq{'std_pcr_i7'}{$pcr7}}."-".$INDEX_TYPE_SEQ_id{'std_pcr_i5'}{$INDEX_TYPE_SEQ_seq{'std_pcr_i5'}{$pcr5}};
			$barc = $INDEX_TYPE_SEQ_seq{'rna_fss_r1'}{$fss}.$INDEX_TYPE_SEQ_seq{'std_pcr_i7'}{$pcr7}.$INDEX_TYPE_SEQ_seq{'std_pcr_i5'}{$pcr5};
			$match = 2;
		} else {
			$name = "rna_fss_r1=".$fss."-std_pcr_i7=".$pcr7."-std_pcr_i5=".$pcr5."-umi=".$umi;
			$barc = $fss.$pcr7.$pcr5;
		}
		
	} elsif ($opt{'M'} =~ /^(met|nome)$/i) {
		
		$pcr7 = substr($i1seq,0,$INDEX_TYPE_length{'std_pcr_i7'}); if (defined $INDEX_EXCLUDE{'std_pcr_i7'}) {$pcr7 = $INDEX_EXCLUDE{'std_pcr_i7'}};
		$sci5 = substr($i2seq,0,$INDEX_TYPE_length{'met_tn5_i5'}); if (defined $INDEX_EXCLUDE{'met_tn5_i5'}) {$sci5 = $INDEX_EXCLUDE{'met_tn5_i5'}};
		$pcr5 = substr($i2seq,$INDEX_TYPE_length{'met_tn5_i5'},$INDEX_TYPE_length{'sci_pcr_i5'}); if (defined $INDEX_EXCLUDE{'sci_pcr_i5'}) {$pcr5 = $INDEX_EXCLUDE{'sci_pcr_i5'}};
		
		if (defined $INDEX_TYPE_SEQ_seq{'met_tn5_i5'}{$sci5} &&
		defined $INDEX_TYPE_SEQ_seq{'sci_pcr_i5'}{$pcr5} &&
		defined $INDEX_TYPE_SEQ_seq{'std_pcr_i7'}{$pcr7}) {
			$name = $INDEX_TYPE_SEQ_id{'std_pcr_i7'}{$INDEX_TYPE_SEQ_seq{'std_pcr_i7'}{$pcr7}}."-".$INDEX_TYPE_SEQ_id{'met_tn5_i5'}{$INDEX_TYPE_SEQ_seq{'met_tn5_i5'}{$sci5}}."-".$INDEX_TYPE_SEQ_id{'sci_pcr_i5'}{$INDEX_TYPE_SEQ_seq{'sci_pcr_i5'}{$pcr5}};
			$barc = $INDEX_TYPE_SEQ_seq{'std_pcr_i7'}{$pcr7}.$INDEX_TYPE_SEQ_seq{'met_tn5_i5'}{$sci5}.$INDEX_TYPE_SEQ_seq{'sci_pcr_i5'}{$pcr5};
			$match = 1;
		} else {
			$name = "std_pcr_i7=".$pcr7."-met_tn5_i5=".$sci5."-sci_pcr_i5=".$pcr5;
			$barc = $pcr7.$sci5.$pcr5;
		}
		
	} elsif ($opt{'M'} =~ /^i5o$/i) {
		
		$pcr7 = substr($i1seq,0,$INDEX_TYPE_length{'std_pcr_i7'}); if (defined $INDEX_EXCLUDE{'std_pcr_i7'}) {$pcr7 = $INDEX_EXCLUDE{'std_pcr_i7'}};
		$sci5 = substr($i2seq,0,$INDEX_TYPE_length{'i5o_tn5_i5'}); if (defined $INDEX_EXCLUDE{'i5o_tn5_i5'}) {$sci5 = $INDEX_EXCLUDE{'i5o_tn5_i5'}};
		$pcr5 = substr($i2seq,$INDEX_TYPE_length{'i5o_tn5_i5'},$INDEX_TYPE_length{'sci_pcr_i5'}); if (defined $INDEX_EXCLUDE{'sci_pcr_i5'}) {$pcr5 = $INDEX_EXCLUDE{'sci_pcr_i5'}};
		
		if (defined $INDEX_TYPE_SEQ_seq{'i5o_tn5_i5'}{$sci5} &&
		defined $INDEX_TYPE_SEQ_seq{'sci_pcr_i5'}{$pcr5} &&
		defined $INDEX_TYPE_SEQ_seq{'std_pcr_i7'}{$pcr7}) {
			$name = $INDEX_TYPE_SEQ_id{'std_pcr_i7'}{$INDEX_TYPE_SEQ_seq{'std_pcr_i7'}{$pcr7}}."-".$INDEX_TYPE_SEQ_id{'i5o_tn5_i5'}{$INDEX_TYPE_SEQ_seq{'i5o_tn5_i5'}{$sci5}}."-".$INDEX_TYPE_SEQ_id{'sci_pcr_i5'}{$INDEX_TYPE_SEQ_seq{'sci_pcr_i5'}{$pcr5}};
			$barc = $INDEX_TYPE_SEQ_seq{'std_pcr_i7'}{$pcr7}.$INDEX_TYPE_SEQ_seq{'i5o_tn5_i5'}{$sci5}.$INDEX_TYPE_SEQ_seq{'sci_pcr_i5'}{$pcr5};
			$match = 1;
		} else {
			$name = "std_pcr_i7=".$pcr7."-i5o_tn5_i5=".$sci5."-sci_pcr_i5=".$pcr5;
			$barc = $pcr7.$sci5.$pcr5;
		}
		
	} elsif ($opt{'M'} =~ /^(rna-atac|atac-rna)$/i) {
		
		$pcr7 = substr($i1seq,0,$INDEX_TYPE_length{'std_pcr_i7'}); if (defined $INDEX_EXCLUDE{'std_pcr_i7'}) {$pcr7 = $INDEX_EXCLUDE{'std_pcr_i7'}}; # both
		$umi = substr($r1seq,0,$umi_length); # rna
		$fss = substr($r1seq,$umi_length,$INDEX_TYPE_length{'rna_fss_r1'}); if (defined $INDEX_EXCLUDE{'rna_fss_r1'}) {$fss = $INDEX_EXCLUDE{'rna_fss_r1'}}; # rna
		$rna_pcr5 = substr($i2seq,0,$INDEX_TYPE_length{'std_pcr_i5'}); if (defined $INDEX_EXCLUDE{'std_pcr_i5'}) {$rna_pcr5 = $INDEX_EXCLUDE{'std_pcr_i5'}}; # rna
		$sci5 = substr($i2seq,0,$INDEX_TYPE_length{'i5o_tn5_i5'}); if (defined $INDEX_EXCLUDE{'i5o_tn5_i5'}) {$sci5 = $INDEX_EXCLUDE{'i5o_tn5_i5'}}; # atac
		$sci_pcr5 = substr($i2seq,$INDEX_TYPE_length{'i5o_tn5_i5'},$INDEX_TYPE_length{'sci_pcr_i5'}); if (defined $INDEX_EXCLUDE{'sci_pcr_i5'}) {$sci_pcr5 = $INDEX_EXCLUDE{'sci_pcr_i5'}}; # atac
		
		if (defined $INDEX_TYPE_SEQ_seq{'rna_fss_r1'}{$fss} &&
		defined $INDEX_TYPE_SEQ_seq{'std_pcr_i5'}{$rna_pcr5} &&
		defined $INDEX_TYPE_SEQ_seq{'std_pcr_i7'}{$pcr7}) { # rna
			$name = $INDEX_TYPE_SEQ_id{'rna_fss_r1'}{$INDEX_TYPE_SEQ_seq{'rna_fss_r1'}{$fss}}."-".$INDEX_TYPE_SEQ_id{'std_pcr_i7'}{$INDEX_TYPE_SEQ_seq{'std_pcr_i7'}{$pcr7}}."-".$INDEX_TYPE_SEQ_id{'std_pcr_i5'}{$INDEX_TYPE_SEQ_seq{'std_pcr_i5'}{$rna_pcr5}};
			$barc = $INDEX_TYPE_SEQ_seq{'rna_fss_r1'}{$fss}.$INDEX_TYPE_SEQ_seq{'std_pcr_i7'}{$pcr7}.$INDEX_TYPE_SEQ_seq{'std_pcr_i5'}{$rna_pcr5};
			$match = 2;
		} elsif (defined $INDEX_TYPE_SEQ_seq{'i5o_tn5_i5'}{$sci5} &&
		defined $INDEX_TYPE_SEQ_seq{'sci_pcr_i5'}{$sci_pcr5} &&
		defined $INDEX_TYPE_SEQ_seq{'std_pcr_i7'}{$pcr7}) { # atac
			$name = $INDEX_TYPE_SEQ_id{'std_pcr_i7'}{$INDEX_TYPE_SEQ_seq{'std_pcr_i7'}{$pcr7}}."-".$INDEX_TYPE_SEQ_id{'i5o_tn5_i5'}{$INDEX_TYPE_SEQ_seq{'i5o_tn5_i5'}{$sci5}}."-".$INDEX_TYPE_SEQ_id{'sci_pcr_i5'}{$INDEX_TYPE_SEQ_seq{'sci_pcr_i5'}{$sci_pcr5}};
			$barc = $INDEX_TYPE_SEQ_seq{'std_pcr_i7'}{$pcr7}.$INDEX_TYPE_SEQ_seq{'i5o_tn5_i5'}{$sci5}.$INDEX_TYPE_SEQ_seq{'sci_pcr_i5'}{$sci_pcr5};
			$match = 1;
		} else {
			$name = "rna_fss_r1=".$fss."-std_pcr_i7=".$pcr7."-std_pcr_i5=".$rna_pcr5."--std_pcr_i7=".$pcr7."-i5o_tn5_i5=".$sci5."-sci_pcr_i5=".$sci_pcr5;
			$barc = $fss.$pcr7.$rna_pcr5."--".$pcr7.$sci5.$sci_pcr5;
		}
	}
	
	if ($opt{'f'} =~ /barc/i) {$tag = $barc} else {$tag = $name};
	
	if ($match == 1) {
		
		$totalCT++;
		
		if (defined $opt{'r'}) {
			$r1out = "\@$tag:$totalCT.$opt{'r'}#0/1\n$r1seq\n\+\n$r1qual";
			$r2out = "\@$tag:$totalCT.$opt{'r'}#0/2\n$r2seq\n\+\n$r2qual";
		} else {
			$r1out = "\@$tag:$totalCT#0/1\n$r1seq\n\+\n$r1qual";
			$r2out = "\@$tag:$totalCT#0/2\n$r2seq\n\+\n$r2qual";
		}
		
		if (!defined $opt{'A'}) {
			print R1OUT "$r1out\n";
			print R2OUT "$r2out\n";
		} else {
			if (defined $CELLID_annot{$barc}) {
				$annot = $CELLID_annot{$barc};
				$ANNOT_count{$annot}++;
				$out1_handle = $HANDLE1{$annot}; $out2_handle = $HANDLE2{$annot};
				print $out1_handle "$r1out\n";
				print $out2_handle "$r2out\n";
			} elsif (defined $CELLID_annot{$name}) {
				$annot = $CELLID_annot{$name};
				$ANNOT_count{$annot}++;
				$out1_handle = $HANDLE1{$annot}; $out2_handle = $HANDLE2{$annot};
				print $out1_handle "$r1out\n";
				print $out2_handle "$r2out\n";
			} else {
				$non_annot_count++;
				print O1 "$r1out\n";
				print O2 "$r2out\n";
			}
		}
		
	} elsif ($match == 2) {
		
		$rnaCT++;
		
		if (defined $opt{'r'}) {
			$rna_out = "\@$tag:$totalCT.$opt{'r'}#$umi\n$r2seq\n\+\n$r2qual";
		} else {
			$rna_out = "\@$tag:$totalCT#$umi\n$r2seq\n\+\n$r2qual";
		}
		
		if (!defined $opt{'A'}) {
			print RNAOUT "$rna_out\n";
		} else {
			if (defined $CELLID_annot{$barc}) {
				$annot = $CELLID_annot{$barc};
				$ANNOT_count_rna{$annot}++;
				$out_rna_handle = $HANDLE_RNA{$annot};
				print $out_rna_handle "$r2out\n";
			} elsif (defined $CELLID_annot{$name}) {
				$annot = $CELLID_annot{$name};
				$ANNOT_count_rna{$annot}++;
				$out_rna_handle = $HANDLE_RNA{$annot};
				print $out_rna_handle "$r2out\n";
			} else {
				$non_annot_rna_count++;
				print OR "$r1out\n";
			}
		}
		
	} else {
	
		$failCT++;
		
		if (defined $opt{'r'}) {
			print R1FAIL "\@$tag:F_$failCT.$opt{'r'}#0/1\n$r1seq\n\+\n$r1qual\n";
			print R2FAIL "\@$tag:F_$failCT.$opt{'r'}#0/1\n$r2seq\n\+\n$r2qual\n";
		} else {
			print R1FAIL "\@$tag:F_$failCT#0/1\n$r1seq\n\+\n$r1qual\n";
			print R2FAIL "\@$tag:F_$failCT#0/1\n$r2seq\n\+\n$r2qual\n";
		}
		
	}
}

close R1; close R2; close I1; close I2;

if (defined $opt{'A'}) {
	foreach $annot (keys %ANNOT_count) {
		$out1_handle = $HANDLE1{$annot}; $out2_handle = $HANDLE2{$annot};
		close $out1_handle; close $out2_handle;
		print STDERR "Annot: $annot, count = $ANNOT_count{$annot}\n";
	}
	close O1; close O2;
}

if (defined $opt{'A'}) {
	read_annot($opt{'A'});
	foreach $annot (keys %ANNOT_count) {
		if ($opt{'M'} !~ /^rna$/i) {
			$out1_handle = $HANDLE1{$annot}; $out2_handle = $HANDLE2{$annot};
			close $out1_handle; close $out2_handle;
			print STDERR "Annot: $annot, count = $ANNOT_count{$annot}\n";
		}
		if ($opt{'M'} =~ /rna/i) {
			$out_rna_handle = $HANDLE_RNA{$annot};
			close $out_rna_handle;
			print STDERR "Annot: $annot, RNA count = $ANNOT_count_rna{$annot}\n";
		}
	}
	if ($opt{'M'} !~ /^rna$/i) {
		close O1; close O2;
	}
	if ($opt{'M'} =~ /rna/i) {
		close OR;
	}
} else {
	if ($opt{'M'} !~ /^rna$/i) {
		close R1OUT; close R2OUT;
	}
	if ($opt{'M'} =~ /rna/i) {
		close RNAOUT;
	}
}

close R1FAIL; close R2FAIL;

}
1;
