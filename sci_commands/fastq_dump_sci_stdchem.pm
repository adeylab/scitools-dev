package sci_commands::fastq_dump_sci_stdchem;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = "fastq_dump_sci_stdchem";

sub fastq_dump_sci_stdchem {

%REVCOMP = ("A" => "T", "C" => "G", "G" => "C", "T" => "A", "N" => "N");
sub revcomp {
	@INSEQ = split(//, uc($_[0]));
	$revcomp = "";
	for ($pos = (@INSEQ-1); $pos >= 0; $pos--) {
		$revcomp .= $REVCOMP{$INSEQ[$pos]};
	}
	return $revcomp;
}

@ARGV = @_;

getopts("R:F:O:o:1:2:A:i:j:r:NVI:Cc:U:vXs", \%opt);

# defaults
$hd_s = 2;
$hd_x = 1;
$umi_len = 6;

$die2 = "
scitools fastq-dump-sci-stdchem [options]
   or    dump-fastq-sci-stdchem
   or    fastq-dump-std
   or    fastq-std

Takes sequencer fastq files (from bcl2fastq) and will format them into
fastq files with matched barcodes.

NOTE: This is an alternative fastq dump mode for standard chemistry sci
libraries wherein the sci index is included in R2 rather than the multi
index layout. This works for the current version of sciWGS and sciMET/NOME
and s3 s4 chemistries.

Options:
   -R   [STR]   Run name (preferred mode)
   -r   [STR]   Run ID (if additional sequencing for some
                libraries. Will add .[ID] to end of read
                names; optional)
   -A   [STR]   Annotation file (will split fastqs)                           # NOTE: RNA component for coassay will not split out properly yet - in development!
   -N           Run is not combinatorial indexing.
                (def = run is treated as combinatorial indexing)
   -D   [INT]   Hamming distance for sample & sci indexes
                (max 2, def = $hd_s)
   -d   [INT]   Hamming distance for sci index
                (max 2, def = $hd_x)
   -C           Coassay - also look for RNA reads (read 2 index)
                Also uses -D as hamming distance
   -s           Single-channel RNA (ie just do RNA split, forces -C)
   -U   [INT]   UMI length (preceeds RNA index on Read 2, def = $umi_len)
   -V 	        If flagged will use reverse compliment for index 2 (i5). 
                This is used for sequencing platforms NovaSeq 6000,
                MiSeq, HiSeq 2500, or HiSeq 2000 System. (Default: no)
   -v           If flagged will take reverse compliment for the
                coassay RNA index (def = no)
   -X           Debug

Defaults:
   -F   [STR]   Fastq directory
         ($VAR{'fastq_input_directory'})
   -O   [STR]   Output fastq directory
         ($VAR{'SCI_fastq_directory'})
   -o   [STR]   Output prefix
         (def = run name)

Specify Index Sets:
   -I   [STR]   Sci index file (def = $VAR{'SCI_stdchem_index_file'})
   -c   [STR]   RNA index file for coassay (def = $VAR{'RNA_Indexes'})

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
if (!defined $opt{'c'}) {$opt{'c'} = $VAR{'RNA_Indexes'}};
if (defined $opt{'s'} && !defined $opt{'C'}) {$opt{'C'} = 1};
if (defined $opt{'N'} && defined $opt{'C'}) {
	die "ERROR: Coassay demultiplexing (-C) is not compatible with non combinatorial indexing splitting (-N)\n";
}
if (!defined $opt{'U'}) {$opt{'U'} = $umi_len};

open IN, "$opt{'I'}";
while ($l = <IN>) {
	chomp $l;
	($id,$pos,$seq) = split(/\t/, $l);
	$POS_SEQ_seq{$pos}{$seq} = $seq;
	$POS_length{$pos} = length($seq);
} close IN;
if (defined $opt{'N'}) {$POS_SEQ_seq{'3'}{'null'} = ""};

if (defined $opt{'C'}) { # also look for RNA component
	open IN, "$opt{'c'}";
	while ($l = <IN>) {
		chomp $l;
		($id,$pos,$seq) = split(/\t/, $l);
		$pos = "R";
		$POS_SEQ_seq{$pos}{$seq} = $seq;
		$POS_length{$pos} = length($seq);
	} close IN;
}

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
					@NEW = @TRUE; $NEW[$i] = "C"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq} elsif ($POS_SEQ_seq{$pos}{$new} ne $seq) {$POS_SEQ_seq{$pos}{$new} = "NA"; $hamming_collisions++};
					@NEW = @TRUE; $NEW[$i] = "G"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq} elsif ($POS_SEQ_seq{$pos}{$new} ne $seq) {$POS_SEQ_seq{$pos}{$new} = "NA"; $hamming_collisions++};
					@NEW = @TRUE; $NEW[$i] = "T"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq} elsif ($POS_SEQ_seq{$pos}{$new} ne $seq) {$POS_SEQ_seq{$pos}{$new} = "NA"; $hamming_collisions++};
					@NEW = @TRUE; $NEW[$i] = "N"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq} elsif ($POS_SEQ_seq{$pos}{$new} ne $seq) {$POS_SEQ_seq{$pos}{$new} = "NA"; $hamming_collisions++};
				} elsif ($TRUE[$i] =~ /C/i) {
					@NEW = @TRUE; $NEW[$i] = "A"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq} elsif ($POS_SEQ_seq{$pos}{$new} ne $seq) {$POS_SEQ_seq{$pos}{$new} = "NA"; $hamming_collisions++};
					@NEW = @TRUE; $NEW[$i] = "G"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq} elsif ($POS_SEQ_seq{$pos}{$new} ne $seq) {$POS_SEQ_seq{$pos}{$new} = "NA"; $hamming_collisions++};
					@NEW = @TRUE; $NEW[$i] = "T"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq} elsif ($POS_SEQ_seq{$pos}{$new} ne $seq) {$POS_SEQ_seq{$pos}{$new} = "NA"; $hamming_collisions++};
					@NEW = @TRUE; $NEW[$i] = "N"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq} elsif ($POS_SEQ_seq{$pos}{$new} ne $seq) {$POS_SEQ_seq{$pos}{$new} = "NA"; $hamming_collisions++};
				} elsif ($TRUE[$i] =~ /G/i) {
					@NEW = @TRUE; $NEW[$i] = "C"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq} elsif ($POS_SEQ_seq{$pos}{$new} ne $seq) {$POS_SEQ_seq{$pos}{$new} = "NA"; $hamming_collisions++};
					@NEW = @TRUE; $NEW[$i] = "A"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq} elsif ($POS_SEQ_seq{$pos}{$new} ne $seq) {$POS_SEQ_seq{$pos}{$new} = "NA"; $hamming_collisions++};
					@NEW = @TRUE; $NEW[$i] = "T"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq} elsif ($POS_SEQ_seq{$pos}{$new} ne $seq) {$POS_SEQ_seq{$pos}{$new} = "NA"; $hamming_collisions++};
					@NEW = @TRUE; $NEW[$i] = "N"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq} elsif ($POS_SEQ_seq{$pos}{$new} ne $seq) {$POS_SEQ_seq{$pos}{$new} = "NA"; $hamming_collisions++};
				} elsif ($TRUE[$i] =~ /T/i) {
					@NEW = @TRUE; $NEW[$i] = "C"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq} elsif ($POS_SEQ_seq{$pos}{$new} ne $seq) {$POS_SEQ_seq{$pos}{$new} = "NA"; $hamming_collisions++};
					@NEW = @TRUE; $NEW[$i] = "G"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq} elsif ($POS_SEQ_seq{$pos}{$new} ne $seq) {$POS_SEQ_seq{$pos}{$new} = "NA"; $hamming_collisions++};
					@NEW = @TRUE; $NEW[$i] = "A"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq} elsif ($POS_SEQ_seq{$pos}{$new} ne $seq) {$POS_SEQ_seq{$pos}{$new} = "NA"; $hamming_collisions++};
					@NEW = @TRUE; $NEW[$i] = "N"; $new = join("", @NEW); if (!defined $POS_SEQ_seq{$pos}{$new}) {$POS_SEQ_seq{$pos}{$new} = $seq} elsif ($POS_SEQ_seq{$pos}{$new} ne $seq) {$POS_SEQ_seq{$pos}{$new} = "NA"; $hamming_collisions++};
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
		if (!defined $opt{'s'}) {
			$out1_handle = "$annot.1";
			open $out1_handle, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.$annot.1.fq.gz";
			$HANDLE1{$annot} = $out1_handle;
			$out2_handle = "$annot.2";
			open $out2_handle, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.$annot.2.fq.gz";
			$HANDLE2{$annot} = $out2_handle;
		}
		if (defined $opt{'C'}) { # open read 1 only for RNA
			$outR_handle = "$annot.R";
			open $outR_handle, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.$annot.RNA.fq.gz";
			$HANDLER{$annot} = $outR_handle;
		}
	}
	if (!defined $opt{'s'}) {
		open O1, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.unassigned.1.fq.gz";
		open O2, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.unassigned.2.fq.gz";
	}
	system("cp $opt{'A'} $opt{'O'}/$opt{'o'}/$opt{'o'}.split.annot");
} else {
	if (!defined $opt{'s'}) {
		open R1OUT, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.1.fq.gz";
		open R2OUT, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.2.fq.gz";
	}
	if (defined $opt{'C'}) { # open read 1 only for RNA
		open RNAOUT, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.RNA.fq.gz";
	}
}

open R1FAIL, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.fail.1.fq.gz";
open R2FAIL, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.fail.2.fq.gz";
open I1FAIL, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.fail.I1.fq.gz";
open I2FAIL, "| $gzip > $opt{'O'}/$opt{'o'}/$opt{'o'}.fail.I2.fq.gz";

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
	$null = <I1>; $null = <I2>;
	$i1qual = <I1>; $i2qual = <I2>; chomp $i1qual; chomp $i2qual;
	
	if (!defined $opt{'s'}) {
		if (defined $opt{'V'}) {
			$rev_i2seq = revcomp($i2seq);
			$i2seq = $rev_i2seq;
		}
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
	}
	if (defined $opt{'C'}) { # Pull RNA index and UMI
		$umi = substr($r2seq,0,$opt{'U'});
		$ixR = substr($r2seq,$opt{'U'},$POS_length{'R'});
		if (defined $opt{'v'}) {
			$rev_ixR = revcomp($ixR);
			$ixR = $rev_ixR;
		}
		if (defined $opt{'X'}) {
			print STDERR "DEBUG: -C flag, umi length = $opt{'U'}, ixR length = $POS_length{'R'}\nDEBUG:    Read 2 = $r2seq\nDEBUG:    UMI = $umi\nDEBUG:    ixR = $ixR\n";
		}
	}
	
	if (!defined $opt{'s'} &&
		defined $POS_SEQ_seq{'1'}{$ix1} &&
		defined $POS_SEQ_seq{'2'}{$ix2} &&
		defined $POS_SEQ_seq{'3'}{$ix3} &&
		$POS_SEQ_seq{'1'}{$ix1} ne "NA" &&
		$POS_SEQ_seq{'2'}{$ix2} ne "NA" &&
		$POS_SEQ_seq{'3'}{$ix3} ne "NA") {
		
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
		
	} elsif (defined $opt{'C'} &&             # RNA coassay check
		defined $POS_SEQ_seq{'1'}{$ix1} &&
		defined $POS_SEQ_seq{'2'}{$ix2} &&
		defined $POS_SEQ_seq{'R'}{$ixR} &&             # RNA coassay check
		$POS_SEQ_seq{'1'}{$ix1} ne "NA" &&
		$POS_SEQ_seq{'2'}{$ix2} ne "NA" &&
		$POS_SEQ_seq{'R'}{$ixR} ne "NA") {
		
		$barc = $POS_SEQ_seq{'1'}{$ix1}.$POS_SEQ_seq{'2'}{$ix2}.$POS_SEQ_seq{'R'}{$ixR};
		
		$totalCT++;
		
		if (defined $opt{'r'}) {
			$r1out = "\@$barc:UMI=$umi.$totalCT.$opt{'r'}#0/1\n$r1seq\n\+\n$r1qual";
		} else {
			$r1out = "\@$barc:UMI=$umi.$totalCT#0/1\n$r1seq\n\+\n$r1qual";
		}
		
		print RNAOUT "$r1out\n";
		
		# ADD IN ANNOT SUPPORT HERE!
		
	} else { # flag as failed
		$r1tag =~ s/\#.+$//;
		print R1FAIL "$r1tag#0/1\n$r1seq\n\+\n$r1qual\n";
		print R2FAIL "$r1tag#0/4\n$r2seq\n\+\n$r2qual\n";
		print I1FAIL "$r1tag#0/2\n$i1seq\n\+\n$i1qual\n";
		print I2FAIL "$r1tag#0/3\n$i2seq\n\+\n$i2qual\n";
		
		$failCT++;
	}

}

# ADD IN RNA FILE CLOSING

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