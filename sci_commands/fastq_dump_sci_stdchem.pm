package sci_commands::fastq_dump_sci_stdchem;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = "fastq_dump_sci_stdchem";

sub fastq_dump_sci_stdchem {

@ARGV = @_;

getopts("R:F:O:o:1:2:A:i:j:r:N", \%opt);

# defaults
$hd_s = 2;
$hd_x = 1;

$die2 = "
scitools fastq-dump-sci-stdchem [options]
   or    dump-fastq-sci-stdchem

Takes sequencer fastq files (from bcl2fastq) and will format
them into fastq files with matched barcodes.

NOTE: This is an alternative fastq dump mode for standard chemistry sci libraries wherein the sci index is included in R2 rather than the multi index layout.
This works for the current version of sciWGS and sciMET/NOME

Options:
   -R   [STR]   Run name (preferred mode)
   -r   [STR]   Run ID (if additional sequencing for some
                libraries. Will add .[ID] to end of read
                names; optional)
   -A   [STR]   Annotation file (will split fastqs)
   -N           Run is not combinatorial indexing. (def = run is treated as combinatorial indexing)
   -D   [INT]   Hamming distance for sample & sci indexes
                (max 2, def = $hd_s)
   -d   [INT]   Hamming distance for sci index
                (max 2, def = $hd_x)

Defaults:
   -F   [STR]   Fastq directory
         ($VAR{'fastq_input_directory'})
   -O   [STR]   Output fastq directory
         ($VAR{'SCI_fastq_directory'})
   -o   [STR]   Output prefix
         (def = run name)


To specify specific fastq files instead of defaults:
   (can be comma sep for multiple)
   -1   [STR]   Read 1 fastq
   -2   [STR]   Read 2 fastq
   -i   [STR]   Index 1 fastq (opt)
   -j   [STR]   Index 2 fastq (opt)

";

if (!defined $opt{'R'}) {die $die2};
if (!defined $opt{'F'}) {$opt{'F'} = $VAR{'fastq_input_directory'}};
if (!defined $opt{'O'}) {$opt{'O'} = $VAR{'SCI_fastq_directory'}};
if (!defined $opt{'o'}) {$opt{'o'} = $opt{'R'}};
if (!defined $opt{'I'}) {$opt{'I'} = $VAR{'SCI_stdchem_index_file'}};


open IN, "$opt{'I'}";
while ($l = <IN>) {
	chomp $l;
	($id,$pos,$seq) = split(/\t/, $l);
	$POS_SEQ_seq{$pos}{$seq} = $seq;
	$POS_length{$pos} = length($seq);
} close IN;
if (defined $opt{'N'}) {$POS_SEQ_seq{'3'}{'null'} = ""};

# make all-1-away hash
foreach $pos (keys %POS_SEQ_seq) {
	if (($pos != 2 && $hd_s >= 1) || ($pos == 2 && $hd_x >= 1)) {
		foreach $seq (keys %{$POS_SEQ_seq{$pos}}) {
			@TRUE = split(//, $seq);
			for ($i = 0; $i < @TRUE; $i++) {
				if ($TRUE[$i] =~ /A/i) {
					@NEW = @TRUE; $NEW[$i] = "C"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
					@NEW = @TRUE; $NEW[$i] = "G"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
					@NEW = @TRUE; $NEW[$i] = "T"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
					@NEW = @TRUE; $NEW[$i] = "N"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
				} elsif ($TRUE[$i] =~ /C/i) {
					@NEW = @TRUE; $NEW[$i] = "A"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
					@NEW = @TRUE; $NEW[$i] = "G"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
					@NEW = @TRUE; $NEW[$i] = "T"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
					@NEW = @TRUE; $NEW[$i] = "N"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
				} elsif ($TRUE[$i] =~ /G/i) {
					@NEW = @TRUE; $NEW[$i] = "C"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
					@NEW = @TRUE; $NEW[$i] = "A"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
					@NEW = @TRUE; $NEW[$i] = "T"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
					@NEW = @TRUE; $NEW[$i] = "N"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
				} elsif ($TRUE[$i] =~ /T/i) {
					@NEW = @TRUE; $NEW[$i] = "C"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
					@NEW = @TRUE; $NEW[$i] = "G"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
					@NEW = @TRUE; $NEW[$i] = "A"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
					@NEW = @TRUE; $NEW[$i] = "N"; $new = join("", @NEW); $POS_SEQ_seq{$pos}{$new} = $seq;
				}
			}
		}
	}
}

# make all-2-away hash
foreach $pos (keys %POS_SEQ_seq) {
	if (($pos != 2 && $hd_s >= 2) || ($pos == 2 && $hd_x >= 2)) {
		foreach $id_seq (keys %{$POS_SEQ_seq{$pos}}) {
			$seq = $POS_SEQ_seq{$pos}{$id_seq};
			@TRUE = split(//, $seq);
			for ($i = 0; $i < @TRUE; $i++) {
				if ($TRUE[$i] =~ /A/i) {
					@NEW = @TRUE; $NEW[$i] = "C"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
					@NEW = @TRUE; $NEW[$i] = "G"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
					@NEW = @TRUE; $NEW[$i] = "T"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
					@NEW = @TRUE; $NEW[$i] = "N"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
				} elsif ($TRUE[$i] =~ /C/i) {
					@NEW = @TRUE; $NEW[$i] = "A"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
					@NEW = @TRUE; $NEW[$i] = "G"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
					@NEW = @TRUE; $NEW[$i] = "T"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
					@NEW = @TRUE; $NEW[$i] = "N"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
				} elsif ($TRUE[$i] =~ /G/i) {
					@NEW = @TRUE; $NEW[$i] = "C"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
					@NEW = @TRUE; $NEW[$i] = "A"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
					@NEW = @TRUE; $NEW[$i] = "T"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
					@NEW = @TRUE; $NEW[$i] = "N"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
				} elsif ($TRUE[$i] =~ /T/i) {
					@NEW = @TRUE; $NEW[$i] = "C"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
					@NEW = @TRUE; $NEW[$i] = "G"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
					@NEW = @TRUE; $NEW[$i] = "A"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
					@NEW = @TRUE; $NEW[$i] = "N"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq};
				}
			}
		}
	}
}


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

system("mkdir $opt{'O'}/$opt{'o'}");

if (defined $opt{'A'}) {
	read_annot($opt{'A'});
	foreach $annot (keys %ANNOT_count) {
		$out1_handle = "$annot.1";
		open $out1_handle, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.$annot.1.fq.gz";
		$HANDLE1{$annot} = $out1_handle;
		$out2_handle = "$annot.2";
		open $out2_handle, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.$annot.2.fq.gz";
		$HANDLE2{$annot} = $out2_handle;
	}
	open O1, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.unassigned.1.fq.gz";
	open O2, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.unassigned.2.fq.gz";
	system("cp $opt{'A'} $opt{'O'}/$opt{'o'}/$opt{'o'}.split.annot");
} else {
	open R1OUT, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.1.fq.gz";
	open R2OUT, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.2.fq.gz";
}

open R1FAIL, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.fail.1.fq.gz";
open R2FAIL, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.fail.2.fq.gz";

$totalCT = 0; $failCT = 0;

open R1, "$zcat $r1 |";
open R2, "$zcat $r2 |";
open I1, "$zcat $i1 |";
open I2, "$zcat $i2 |";

	
while ($r1tag = <R1>) {

	$r2tag = <R2>; chomp $r1tag; chomp $r2tag;
	$r1seq = <R1>; $r2seq = <R2>; chomp $r1seq; chomp $r2seq;
	$null = <R1>; $null = <R2>;
	$r1qual = <R1>; $r2qual = <R2>; chomp $r1qual; chomp $r2qual;
	
	$i1tag = <I1>; chomp $i1tag; $i2tag = <I2>; chomp $i2tag;
	$i1seq = <I1>; chomp $i1seq; $i2seq = <I2>; chomp $i2seq;
	$null = <I1>; $null = <I1>;
	$null = <I2>; $null = <I2>;
	
	$ix1 = substr($i1seq,0,$POS_length{'1'});
	$ix2 = substr($i2seq,0,$POS_length{'2'});
	if (!defined $opt{'N'}) {
		$ix3 = substr($r2seq,0,$POS_length{'3'});
		$r2dna = substr($r2seq,($POS_length{'3'}+20));
		$r2qtrm = substr($r2qual,($POS_length{'3'}+20));
		$r2qual = $r2qtrm;
	} else {
		$ix3 = "null";
		$r2dna = $r2seq;
		
	}
	
	if (defined $POS_SEQ_seq{'1'}{$ix1} &&
		defined $POS_SEQ_seq{'2'}{$ix2} &&
		defined $POS_SEQ_seq{'3'}{$ix3}) {
		
		$barc = $POS_SEQ_seq{'1'}{$ix1}.$POS_SEQ_seq{'2'}{$ix2}.$POS_SEQ_seq{'3'}{$ix3};
		
		$totalCT++;
		
		if (defined $opt{'r'}) {
			$r1out = "\@$barc:$totalCT.$opt{'r'}#0/1\n$r1seq\n\+\n$r1qual";
			$r2out = "\@$barc:$totalCT.$opt{'r'}#0/2\n$r2dna\n\+\n$r2qual";
		} else {
			$r1out = "\@$barc:$totalCT#0/1\n$r1seq\n\+\n$r1qual";
			$r2out = "\@$barc:$totalCT#0/2\n$r2dna\n\+\n$r2qual";
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
			} else {
				$non_annot_count++;
				print O1 "$r1out\n";
				print O2 "$r2out\n";
			}
		}
		
	} else {
	
		$barc = $ix1.$ix2.$ix3;
		
		if (defined $opt{'r'}) {
			print R1FAIL "\@$barc:F_$failCT.$opt{'r'}#0/1\n$r1seq\n\+\n$r1qual\n";
			print R2FAIL "\@$barc:F_$failCT.$opt{'r'}#0/1\n$r2dna\n\+\n$r2qual\n";
		} else {
			print R1FAIL "\@$barc:F_$failCT#0/1\n$r1seq\n\+\n$r1qual\n";
			print R2FAIL "\@$barc:F_$failCT#0/1\n$r2dna\n\+\n$r2qual\n";
		}
		
		$failCT++;
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

close R1FAIL; close R2FAIL;
}
1;
