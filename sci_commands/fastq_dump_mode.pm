package sci_commands::fastq_dump_mode;

use sci_utils::general;
use sci_utils::modes;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("fastq_dump_mode");

sub fastq_dump_mode {

@ARGV = @_;
getopts("R:F:O:o:I:1:2:A:i:j:r:M:f:k:m:x", \%opt);

# Defaults
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
   -x           If -A is specified, reads matching barcodes but not an
                annotation will be sent to fail (def = to unmatched)

To specify specific fastq files instead of defaults:
   (can be comma sep for multiple)
   -1   [STR]   Read 1 fastq
   -2   [STR]   Read 2 fastq
   -i   [STR]   Index 1 fastq
   -j   [STR]   Index 2 fastq

Will split with a hamming distance of 2

";

if (!defined $opt{'R'} && !defined $opt{'1'}) {die $die2};
if (!defined $opt{'F'}) {$opt{'F'} = $VAR{'fastq_input_directory'}};
if (!defined $opt{'O'}) {$opt{'O'} = $VAR{'SCI_fastq_directory'}};
if (!defined $opt{'o'}) {$opt{'o'} = $opt{'R'}};
if (defined $opt{'M'}) {$mode_name = $opt{'M'}};
if (!defined $opt{'I'}) {$opt{'I'} = "DIR=$VAR{'index_directory'}"};
if (!defined $opt{'f'}) {$opt{'f'} = $def_read_name_format};
if (defined $opt{'m'}) {$VAR{'sci_modes'} = $opt{'m'}};

# detect & specify fastq files
detect_input_fastqs();
if (!defined $read1) {die "ERROR: At a very minimum read1 must be specified or able to be detected.\n"};

# load in sci mode
$mode = modes->new($mode_name,$VAR{'sci_modes'});
# get set of modalities
@MODALITIES = $mode->modalities();
# get ordered index list for the mode
for ($modalityID = 0; $modalityID < @MODALITIES; $modalityID++) {
	$modality = $MODALITIES[$modalityID];
	@{$MODALITY_INDEXES{$modality}} = $mode->indexes($modality);
}

# read in indexes
read_indexdir($opt{'I'});
check_indexes();

# build 2-hamming distance hash & set index exclusions
make_hamming_hash();

# open input fastqs that are present int he mode & were detected
open_input_fastqs();

# Make output directory
system("mkdir $opt{'O'}/$opt{'o'}");

# open output fastqs for each modality
open_outs();
open_fail_fastqs();

# start parsing fastq files
$read_number = 0;
while ($tag = <R1>) {
	# pull reads and build read_set object
	$read_set = pull_reads();
	
	$any_passing = 0;
	$out_barc = ""; $out_name = "";
	$umi = "null";
	
	for ($modalityID = 0; $modalityID < @MODALITIES; $modalityID++) {
		$modality = $MODALITIES[$modalityID]; $modality_pass = 1;
		for ($index_pos = 0; $index_pos < @{$MODALITY_INDEXES{$modality}}; $index_pos++) {
			$index_name = $MODALITY_INDEXES{$modality}[$index_pos];
			if ($index_name =~ /umi/) {
				# special case for umi's
				$umi = $mode->pull_index($modality,"umi",$read_set);
			} else {
				# regular index
				$index_seq = $mode->pull_index($modality,$index_name,$read_set);
				if (defined $INDEX_TYPE_SEQ_seq{$index_name}{$index_seq}) {
					$out_barc .= "$INDEX_TYPE_SEQ_seq{$index_name}{$index_seq}";
					$out_name .= "$INDEX_TYPE_SEQ_id{$index_name}{$index_seq}-";
				} else {
					$modality_pass = 0;
					$out_barc .= $index_seq;
					$out_name .= "$index_name=$index_seq-";
				}
			}
		}
		
		if (defined $opt{'A'}) {
			if (defined $CELLID_annot{$out_barc}) {
				$annot = $CELLID_annot{$out_barc}
			} else {
				$annot = "unmatched";
			}
			if (defined $opt{'x'} && $annot eq "unmatched") {
				$modality_pass = 0;
			}
		}
		
		if ($modality_pass > 0) {
			$any_passing++;
			if ($umi ne "null") {
				$out_barc .= ":UMI=$umi";
				$out_name .= ":UMI=$umi";
			}
			# print passing reads to this modality
			print_outs();
		}
	}
	
	if ($any_passing < 1) { # failed in all modalities - print to fail files
		print_failing();
	} elsif ($any_passing >= 2) {
		print STDERR "WARNING: Read $tag passed in more than one modality.\n";
	}
	
}

close_outs();

}

#### SUBROUTINES ####

sub check_indexes {
	$absent_index = 0;
	foreach $check_modality (@MODALITIES) {
		foreach $check_index (@{$MODALITY_INDEXES{$check_modality}}) {
			if (!defined $INDEX_TYPE_SEQ_seq{$check_index} && $check_index !~ /^umi$/i) {
				$absent_index++;
			}
			if ($INDEX_TYPE_length{$check_index} != $mode->index_length($check_modality,$check_index)) {
				die "ERROR: The index type ($check_index) specified has a different length in the mode configuration (".$mode->index_length($check_modality,$check_index).") than in the specified indexes ($INDEX_TYPE_length{$check_index}).\n";
			}
		}
	}
	if ($absent_index>0) {die "ERROR: Index types specified in mode: $mode_name could not be found in the specified index file(s)/directory(s)!\n"};
}

sub make_hamming_hash {
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
}

sub detect_input_fastqs {
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
}

sub open_input_fastqs {
	if ($mode->read_check(read1) eq "true") {open R1, "zcat $read1 |" || die "ERROR: Cannot open read1!\n"};
	if ($mode->read_check(read2) eq "true") {open R2, "zcat $read2 |" || die "ERROR: Cannot open read2!\n"};
	if ($mode->read_check(index1) eq "true") {open I1, "zcat $index1 |" || die "ERROR: Cannot open index1!\n"};
	if ($mode->read_check(index2) eq "true") {open I2, "zcat $index2 |" || die "ERROR: Cannot open index2!\n"};
}

sub open_outs {
	%OUT_HANDLES = ();
	for ($modalityID = 0; $modalityID < @MODALITIES; $modalityID++) {
		$modality = $MODALITIES[$modalityID];
		@{$MODALITY_OUTPUTS{$modality}} = $mode->outputs($modality);
		if (defined $opt{'A'}) {
			read_annot($opt{'A'});
			foreach $annot (keys %ANNOT_count) {
				foreach $out_type (@{$MODALITY_OUTPUTS{$modality}}) {
					$handle = "$annot.$modality.$out_type";
					$OUT_HANDLES{$handle} = 1;
					if ($mode->check_multimodal() =~ /true/) {
						open $handle, "| gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.$handle.fq.gz";
					} else {
						open $handle, "| gzip > $opt{'O'}/$opt{'o'}.$annot.$out_type.fq.gz";
					}
				}
			}
			if (!defined $opt{'x'}) {
				foreach $out_type (@{$MODALITY_OUTPUTS{$modality}}) {
					$handle = "unmatched.$modality.$out_type";
					$OUT_HANDLES{$handle} = 1;
					if ($mode->check_multimodal() =~ /true/) {
						open $handle, "| gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.$handle.fq.gz";
					} else {
						open $handle, "| gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.unmatched.$out_type.fq.gz";
					}
				}
			}
		} else {
			foreach $out_type (@{$MODALITY_OUTPUTS{$modality}}) {
				$handle = "$modality.$out_type";
				$OUT_HANDLES{$handle} = 1;
				if ($mode->check_multimodal() =~ /true/) {
					open $handle, "| gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.$handle.fq.gz";
				} else {
					open $handle, "| gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.$out_type.fq.gz";
				}
			}
		}
	}
}

sub print_outs {
	foreach $out_type (@{$MODALITY_OUTPUTS{$modality}}) {
		if (defined $opt{'A'}) {
			$handle = "$annot.$modality.$out_type";
		} else {
			$handle = "$modality.$out_type";
		}
		if (defined $OUT_HANDLES{$handle}) {
			if ($opt{'f'} =~ /barc/) {$tag = $out_barc} else {$tag = $out_name};
			$seq = $mode->pull_seq($modality,$out_type,$read_set);
			$qual = $mode->pull_qual($modality,$out_type,$read_set);
			print $handle "$tag\n$seq\n\+\n$qual\n";
		} else {
			print STDERR "WARNING: $handle does not exist for read $tag\n";
		}
	}
}

sub close_outs {
	foreach $handle (keys %OUT_HANDLES) {
		close $handle;
	}
	if ($mode->read_check(read1) eq "true") {close R1F};
	if ($mode->read_check(read2) eq "true") {close R2F};
	if ($mode->read_check(index1) eq "true") {close I1F};
	if ($mode->read_check(index2) eq "true") {close I2F};
}

sub open_fail_fastqs {
	if ($mode->read_check(read1) eq "true") {open R1F, "| gzip > $opt{'O'}/$opt{'o'}.fail.read1.fq.gz"};
	if ($mode->read_check(read2) eq "true") {open R2F, "| gzip > $opt{'O'}/$opt{'o'}.fail.read2.fq.gz"};
	if ($mode->read_check(index1) eq "true") {open I1F, "| gzip > $opt{'O'}/$opt{'o'}.fail.index1.fq.gz"};
	if ($mode->read_check(index2) eq "true") {open I2F, "| gzip > $opt{'O'}/$opt{'o'}.fail.index2.fq.gz"};
}

sub pull_reads {
	chomp $tag; $seq = <R1>; chomp $seq; $null = <R1>; $qual = <R1>; chomp $qual;
	$read_pull->{read1}->{tag} = $tag; $read_pull->{read1}->{seq} = $seq; $read_pull->{read1}->{qual} = $qual;
	if ($mode->read_check(read2) eq "true") {
		$tag = <R2>; chomp $tag; $seq = <R2>; chomp $seq; $null = <R2>; $qual = <R2>; chomp $qual;
		$read_pull->{read2}->{tag} = $tag; $read_pull->{read2}->{seq} = $seq; $read_pull->{read2}->{qual} = $qual;
	}
	if ($mode->read_check(index1) eq "true") {
		$tag = <I1>; chomp $tag; $seq = <I1>; chomp $seq; $null = <I1>; $qual = <I1>; chomp $qual;
		$read_pull->{index1}->{tag} = $tag; $read_pull->{index1}->{seq} = $seq; $read_pull->{index1}->{qual} = $qual;
	}
	if ($mode->read_check(index2) eq "true") {
		$tag = <I2>; chomp $tag; $seq = <I2>; chomp $seq; $null = <I2>; $qual = <I2>; chomp $qual;
		$read_pull->{index2}->{tag} = $tag; $read_pull->{index2}->{seq} = $seq; $read_pull->{index2}->{qual} = $qual;
	}
	return $read_pull;
}

sub print_failing {
	if ($mode->read_check(read1) eq "true") {print R1F $read_set->{read1}->{tag}."\n".$read_set->{read1}->{seq}."\n\+\n".$read_seq->{read1}->{qual}."\n"};
	if ($mode->read_check(read2) eq "true") {print R2F $read_set->{read2}->{tag}."\n".$read_set->{read2}->{seq}."\n\+\n".$read_seq->{read2}->{qual}."\n"};
	if ($mode->read_check(index1) eq "true") {print I1F $read_set->{index1}->{tag}."\n".$read_set->{index1}->{seq}."\n\+\n".$read_seq->{index1}->{qual}."\n"};
	if ($mode->read_check(index2) eq "true") {print I2F $read_set->{index2}->{tag}."\n".$read_set->{index2}->{seq}."\n\+\n".$read_seq->{index2}->{qual}."\n"};
}

1;
