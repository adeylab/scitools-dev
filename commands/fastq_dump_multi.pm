package commands::fastq_dump_multi;

use commands::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("fastq_dump_multi");

sub fastq_dump_multi {

@ARGV = @_;
getopts("R:F:O:o:I:1:2:A:i:j:r:M:u:f:k:", \%opt);

# Defaults
$umi_length = 8;
$def_read_name_format = "barc";

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
   -M   [STR]   Mode. Options: (def = sci)
                  sci (sci-ATAC-seq, sci-DNA-seq, sci-GCC)
                  met/nome (sci-MET, sci-NOMe)
                  rna-atac/atac-rna (sci-ATAC/RNA-seq)
                  rna (sci-RNA-seq)
                  i5o (sci- with i5 indexes only)
   -A   [STR]   Annotation file (will split fastqs)

Defaults:
   -F   [STR]   Fastq directory
         ($VAR{'fastq_input_directory'})
   -O   [STR]   Output fastq directory
         ($VAR{'SCI_fastq_directory'})
   -o   [STR]   Output prefix
         (def = run name)
   -I   [STR]   Index files/directory - comma separated
         Specify: DIR=[directory] to read all files present
           or   NAME=[file] to include files
         Def: DIR=$VAR{'index_directory'}
         Index file format, tab-delimited:
           col1 = <[PCR/Tn5/RNA/MET]_[Set_letter]_[i5/i7/r1]_[Row/Col_ID/well]>
           col2 = <Index_type>
           col3 = <Sequence>
   -u   [INT]   UMI length (def = $umi_length)
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
if (!defined $opt{'M'}) {$opt{'M'} = "sci"};
if (defined $opt{'u'}) {$umi_length = $opt{'u'}};
if (!defined $opt{'I'}) {$opt{'I'} = "DIR=$VAR{'index_directory'}"};
if (!defined $opt{'f'}) {$opt{'f'} = $def_read_name_format};

# read in indexes
read_indexdir($opt{'I'});
print STDERR "INFO: $indexes_loaded indexes read.\n";
foreach $type (keys %INDEX_TYPE_SEQ_seq) {
	print STDERR "\tIndex type: $type found.\n";
}

# verify mode and indexes
if ($opt{'M'} =~ /^sci$/i) {
	if (!defined $INDEX_TYPE_SEQ_seq{'sci_tn5_i5'} ||
		!defined $INDEX_TYPE_SEQ_seq{'sci_tn5_i7'} ||
		!defined $INDEX_TYPE_SEQ_seq{'sci_pcr_i5'} ||
		!defined $INDEX_TYPE_SEQ_seq{'sci_pcr_i7'}) {
			die "ERROR: Mode $opt{'M'} requires indexes of the following types: sci_tn5_i5, sci_tn5_i7, sci_pcr_i5, and sci_pcr_i7.\n";
	}
} elsif ($opt{'M'} =~ /^rna$/i) {
	if (!defined $INDEX_TYPE_SEQ_seq{'rna_fss_r1'} ||
		!defined $INDEX_TYPE_SEQ_seq{'std_pcr_i5'} ||
		!defined $INDEX_TYPE_SEQ_seq{'std_pcr_i7'}) {
			die "ERROR: Mode $opt{'M'} requires indexes of the following types: rna_fss_r1, std_pcr_i5, and std_pcr_i7.\n";
	}
} elsif ($opt{'M'} =~ /^(met|nome)$/i) {
	if (!defined $INDEX_TYPE_SEQ_seq{'met_tn5_i5'} ||
		!defined $INDEX_TYPE_SEQ_seq{'met_pcr_i5'} ||
		!defined $INDEX_TYPE_SEQ_seq{'std_pcr_i7'}) {
			die "ERROR: Mode $opt{'M'} requires indexes of the following types: met_tn5_i5, met_pcr_i5, and std_pcr_i7.\n";
	}
} elsif ($opt{'M'} =~ /^(rna-atac|atac-rna)$/i) {
	if (!defined $INDEX_TYPE_SEQ_seq{'i5o_tn5_i5'} ||
		!defined $INDEX_TYPE_SEQ_seq{'rna_fss_r1'} ||
		!defined $INDEX_TYPE_SEQ_seq{'sci_pcr_i5'} ||
		!defined $INDEX_TYPE_SEQ_seq{'std_pcr_i5'} ||
		!defined $INDEX_TYPE_SEQ_seq{'std_pcr_i7'}) {
			die "ERROR: Mode $opt{'M'} requires indexes of the following types: i5o_tn5_i5, rna_fss_r1, sci_pcr_i5, std_pcr_i5, and std_pcr_i7.\n";
	}
} elsif ($opt{'M'} =~ /i5o/i) {
	if (!defined $INDEX_TYPE_SEQ_seq{'i5o_tn5_i5'} ||
		!defined $INDEX_TYPE_SEQ_seq{'sci_pcr_i5'} ||
		!defined $INDEX_TYPE_SEQ_seq{'std_pcr_i7'}) {
			die "ERROR: Mode $opt{'M'} requires indexes of the following types: i5o_tn5_i5, sci_pcr_i5, and std_pcr_i7.\n";
	}
} else {
	die "ERROR: Cannot interpret mode: $opt{'M'}\n";
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

# Identify fastq files
if (!defined $opt{'1'}) {
	if (-e "$opt{'F'}/$opt{'R'}/Undetermined_S0_L001_R1_001.fastq.gz") {
		$r1 = "$opt{'F'}/$opt{'R'}/Undetermined_S0_L001_R1_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L002_R1_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L003_R1_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L004_R1_001.fastq.gz";
		$r2 = "$opt{'F'}/$opt{'R'}/Undetermined_S0_L001_R2_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L002_R2_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L003_R2_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L004_R2_001.fastq.gz";
		$i1 = "$opt{'F'}/$opt{'R'}/Undetermined_S0_L001_I1_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L002_I1_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L003_I1_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L004_I1_001.fastq.gz";
		$i2 = "$opt{'F'}/$opt{'R'}/Undetermined_S0_L001_I2_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L002_I2_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L003_I2_001.fastq.gz $opt{'F'}/$opt{'R'}/Undetermined_S0_L004_I2_001.fastq.gz";
	} elsif (-e "$opt{'F'}/$opt{'R'}/Undetermined_S0_R1_001.fastq.gz") {
		$r1 = "$opt{'F'}/$opt{'R'}/Undetermined_S0_R1_001.fastq.gz";
		$r2 = "$opt{'F'}/$opt{'R'}/Undetermined_S0_R2_001.fastq.gz";
		$i1 = "$opt{'F'}/$opt{'R'}/Undetermined_S0_I1_001.fastq.gz";
		$i2 = "$opt{'F'}/$opt{'R'}/Undetermined_S0_I2_001.fastq.gz";
	}
} else {
	$r1 = "$opt{'1'}";
	$r2 = "$opt{'2'}";
	$i1 = "$opt{'i'}";
	$i2 = "$opt{'j'}";
}

# Make output directory
system("mkdir $opt{'O'}/$opt{'o'}");

# Open output files
if (defined $opt{'A'}) {
	read_annot($opt{'A'});
	foreach $annot (keys %ANNOT_count) {
		if ($opt{'M'} !~ /^rna$/i) {
			$out1_handle = "$annot.1";
			open $out1_handle, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.$annot.1.fq.gz";
			$HANDLE1{$annot} = $out1_handle;
			$out2_handle = "$annot.2";
			open $out2_handle, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.$annot.2.fq.gz";
			$HANDLE2{$annot} = $out2_handle;
		}
		if ($opt{'M'} =~ /rna/i) {
			$out_rna_handle = "$annot.rna";
			open $out_rna_handle, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.$annot.rna.fq.gz";
			$HANDLE_RNA{$annot} = $out_rna_handle;
		}
	}
	if ($opt{'M'} !~ /^rna$/i) {
		open O1, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.unassigned.1.fq.gz";
		open O2, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.unassigned.2.fq.gz";
	}
	if ($opt{'M'} =~ /rna/i) {
		open OR, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.unassigned.rna.fq.gz";
	}
	system("cp $opt{'A'} $opt{'O'}/$opt{'o'}/$opt{'o'}.split.annot");
} else {
	if ($opt{'M'} !~ /^rna$/i) {
		open R1OUT, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.1.fq.gz";
		open R2OUT, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.2.fq.gz";
	}
	if ($opt{'M'} =~ /rna/i) {
		open RNAOUT, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.rna.fq.gz";
	}
}

open R1FAIL, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.fail.1.fq.gz";
open R2FAIL, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.fail.2.fq.gz";

$totalCT = 0; $failCT = 0; $rnaCT = 0;

open R1, "$zcat $r1 |" || die "ERROR: Cannot open read 1 fastq files: $r1\n";
open R2, "$zcat $r2 |" || die "ERROR: Cannot open read 2 fastq files: $r2\n";
open I1, "$zcat $i1 |" || die "ERROR: Cannot open index 1 fastq files: $i1\n";
open I2, "$zcat $i2 |" || die "ERROR: Cannot open index 2 fastq files: $i2\n";


while ($r1tag = <R1>) {

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
