package sci_commands::bam_rmdup;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("bam_rmdup");

sub bam_rmdup {

@ARGV = @_;
getopts("s:O:xm:H:e:c:Cnrt:XNS", \%opt);

$memory = "2G";
$sort_threads = 1;
$chr_list = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22";

$die2 = "
scitools bam-rmdup [options] [sorted bam file] (additonal sorted bam file(s)) ...
   or    rmdup

Will produce a barcode-based duplicate removed bam and associated
complexity file. Will merge bams if multipe are specified.

Uses header from the first bam provided. If -O is not specified, the name
of the first bam provided will be used as the output prefix.

Will exclude /(M|Y|L|K|G|Un|un|random|alt|Random|Alt)/ chroms

Options:
   -O   [STR]   Output prefix (default is bam file prefix)
                 It is highly recommended to specify -O if
                 multiple bams are used as the input.
   -x           If multiple bams, do not include BAMID field
   -e   [PAT]   Patterns for chromosomes to exclude (comma sep)
                (def = chrM,chrY,chrUn,random,alt  can be set to 'none')
   -m   [MEM]   Samtools sort max memory K/M/G (def = $memory)
                 (only for multiple bams)
   -t   [INT]   Threads for sorting (def = $sort_threads)
   -C           By chromosome (for large files)
                 Bams must be idnexed (if not will make index)
                 Will ignore -e and just use chroms in list -c
   -c   [STR]   List chr to include:
                 def is chr1 to chr22 and chrX
                 provide as a comma sep list, eg: chr1,chr2,...
   -n           Perform name-sorted rmdup
                 Requires a name sorted bam for a single bam file
                 or multiple bam files in any sort.
   -N           If a single bam is provided and not name sorted,
                 then name sort first
   -r           If -n, retain name sort (def = chr/pos sort)
   -S           Discard reads with the same coordinates but
                 opposite strands (def = retain)
   -X           Retain intermediate merged nsrt bam
   -H   [BAM]   Use header from this bam instead.
   -s   [STR]   Samtools call (def = $samtools)

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.bam$//};
if (defined $opt{'s'}) {$samtools = $opt{'s'}};
if (defined $opt{'m'}) {$memory = $opt{'m'}};
if (!defined $ARGV[1]) {$opt{'x'} = 1};
if (defined $opt{'e'}) {@CHR_FILT = split(/,/, $opt{'e'})};
if (defined $opt{'c'}) {$chr_list = $opt{'c'}};
if (defined $opt{'t'}) {$sort_threads = $opt{'t'}};
if (defined $opt{'N'}) {$opt{'n'} = 1};

if (!defined $opt{'n'}) {
	if (!defined $ARGV[1]) {
		open OUT, "| $samtools view -bS - > $opt{'O'}.bbrd.q10.bam 2>/dev/null";
	} else {
		$out_prefix = "$opt{'O'}.bbrd.q10";
		open OUT, "| $samtools view -bSu - | $samtools sort -@ $sort_threads -m $memory -T $out_prefix.TMP - > $out_prefix.bam";
	}
	if (!defined $opt{'H'}) {$opt{'H'} = $ARGV[0]};
	open H, "$samtools view -H $opt{'H'} |";
	while ($l = <H>) {print OUT $l};
	close H;
	$reads_q10_to_other_chr = 0;
}

if (defined $opt{'C'}) {
	for ($bamID = 0; $bamID < @ARGV; $bamID++) {
		if (-e "$ARGV[$bamID].bai") {
			print STDERR "INFO: Bam $ARGV[$bamID] index found!\n";
		} else {
			print STDERR "INFO: Indexing $ARGV[$bamID] ...\n";
			system("$samtools index $ARGV[$bamID]");
		}
	}
}

if (defined $opt{'C'}) { # by chromosome
	@CHR_LIST = split(/,/, $chr_list);
	for ($chrID = 0; $chrID < @CHR_LIST; $chrID++) {
		%KEEP = (); %BARC_POS_ISIZE = (); %OBSERVED = ();
		for ($bamID = 0; $bamID < @ARGV; $bamID++) {
			open IN, "$samtools view -q 10 $ARGV[$bamID] $CHR_LIST[$chrID] |";
			while ($l = <IN>) {
				$q10_reads++;
				chomp $l;
				@P = split(/\t/, $l);
				if (!defined $opt{'x'}) {
					$P[0] =~ s/#.+$//;
					$P[0] .= ":BAMID=$bamID";
					$l = join("\t", @P);
				}
				$barc = $P[0]; $barc =~ s/:.+$//;
				if (defined $KEEP{$P[0]}) {
					print OUT "$l\n";
					$BARC_total{$barc}++;
					$BARC_kept{$barc}++;
					$total_kept++;
				} else {
					$BARC_total{$barc}++;
					if (defined $opt{'S'}) {
						$pos_isize = "$P[3]:$P[8]";
					} else {
						if ($P[1] & 16) {
							$pos_isize = "$P[3]:$P[8]:r";
						} else {
							$pos_isize = "$P[3]:$P[8]:f";
						}
					}
					if (!defined $BARC_POS_ISIZE{$barc}{$pos_isize} && !defined $OBSERVED{$P[0]}) {
						$BARC_POS_ISIZE{$barc}{$pos_isize} = 1;
						$KEEP{$P[0]} = 1;
						print OUT "$l\n";
						$BARC_kept{$barc}++;
						$total_kept++;
						$BAMID_included{$bamID}++;
					}
					$OBSERVED{$P[0]} = 1;
				}
			} close IN;
			%KEEP = ();
		}
	}
} elsif (defined $opt{'n'}) { # by namesort
	$opt{'O'} =~ s/\.nsrt//;
	if (defined $ARGV[1] || defined $opt{'N'}) { # multiple bams - qfilt, add bamID and nsort
		open OUT, "| $samtools sort -@ $sort_threads -m $memory -T $opt{'O'}.merged.nsrt.TMP -n - > $opt{'O'}.merged.nsrt.bam";
		open HEAD, "$samtools view -H $ARGV[0] |";
		while ($l = <HEAD>){print OUT "$l"};
		close HEAD;
		for ($bamID = 0; $bamID < @ARGV; $bamID++) {
			open IN, "$samtools view -q 10 $ARGV[$bamID] |";
			while ($l = <IN>) {
				chomp $l;
				@P = split(/\t/, $l);
				if (!defined $opt{'x'}) {
					$P[0] =~ s/#.+$//;
					$P[0] .= ":BAMID=$bamID";
					$l = join("\t", @P);
				}
				print OUT "$l\n";
			}
			close IN;
		} close OUT;
		open IN, "$samtools view -q 10 $opt{'O'}.merged.nsrt.bam |";
	} else { # single bam
		open IN, "$samtools view -q 10 $ARGV[0] |";
	}
	
	if (defined $opt{'r'}) {
		open OUT, "| $samtools view -bS - > $opt{'O'}.bbrd.q10.nsrt.bam";
	} else {
		open OUT, "| $samtools view -bSu - | $samtools sort -@ $sort_threads -m $memory -T $opt{'O'}.TMP - > $opt{'O'}.bbrd.q10.bam";
	}
	
	open HEAD, "$samtools view -H $ARGV[0] |";
	while ($l = <HEAD>){print OUT "$l"};
	close HEAD;
	
	$currentBarc = "null";
	while ($l = <IN>) {
		$q10_reads++;
		chomp $l;
		@P = split(/\t/, $l);
		$barc = $P[0]; $barc =~ s/:.+$//;
		if ($currentBarc ne $barc) {
			# set hashes
			%KEEP = (); %OBSERVED = (); %POS_ISIZE = ();
			$currentBarc = $barc;
		}
		# process read
		if (defined $KEEP{$P[0]}) {
			print OUT "$l\n";
			$BARC_total{$barc}++;
			$BARC_kept{$barc}++;
			$total_kept++;
		} elsif ($P[1] & 4) {} else {
			$filt_chrom = 0;
			if (!defined $opt{'e'} && $P[2] =~ /chrM|chrY|chrUn|random|alt/) {
				$filt_chrom+=10;
			} elsif ($opt{'e'} ne "none") {
				foreach $pattern (@CHR_FILT) {
					if ($P[2] =~ /$pattern/) {
						$filt_chrom+=10;
					}
				}
			}
			
			if ($filt_chrom < 1) {
				$BARC_total{$barc}++;
				if (defined $opt{'S'}) {
					$pos_isize = "$P[2]:$P[3]:$P[8]";
				} else {
					if ($P[1] & 16) {
						$pos_isize = "$P[2]:$P[3]:$P[8]:r";
					} else {
						$pos_isize = "$P[2]:$P[3]:$P[8]:f";
					}
				}
				if (!defined $POS_ISIZE{$pos_isize} && !defined $OBSERVED{$P[0]}) {
					$POS_ISIZE{$pos_isize} = 1;
					$KEEP{$P[0]} = 1;
					print OUT "$l\n";
					$BARC_kept{$barc}++;
					$total_kept++;
					$BAMID_included{$bamID}++;
				}
				$OBSERVED{$P[0]} = 1;
			} else {
				$reads_q10_to_other_chr++;
			}
		}
	} close IN;
	
	if (defined $ARGV[1] && !defined $opt{'X'}) {system("rm -f $opt{'O'}.merged.nsrt.bam")};
	
} else { # standard
	for ($bamID = 0; $bamID < @ARGV; $bamID++) {
		open IN, "$samtools view -q 10 $ARGV[$bamID] |";
		while ($l = <IN>) {
			$q10_reads++;
			chomp $l;
			@P = split(/\t/, $l);
			if (!defined $opt{'x'}) {
				$P[0] =~ s/#.+$//;
				$P[0] .= ":BAMID=$bamID";
				$l = join("\t", @P);
			}
			$barc = $P[0]; $barc =~ s/:.+$//;
			if (defined $KEEP{$P[0]}) {
				print OUT "$l\n";
				$BARC_total{$barc}++;
				$BARC_kept{$barc}++;
				$total_kept++;
			} elsif ($P[1] & 4) {} else {
				$filt_chrom = 0;
				if (!defined $opt{'e'} && $P[2] =~ /chrM|chrY|chrUn|random|alt/) {
					$filt_chrom+=10;
				} elsif ($opt{'e'} ne "none") {
					foreach $pattern (@CHR_FILT) {
						if ($P[2] =~ /$pattern/) {
							$filt_chrom+=10;
						}
					}
				}
				
				if ($filt_chrom < 1) {
					$BARC_total{$barc}++;
					if (defined $opt{'S'}) {
						$pos_isize = "$P[2]:$P[3]:$P[8]";
					} else {
						if ($P[1] & 16) {
							$pos_isize = "$P[2]:$P[3]:$P[8]:r";
						} else {
							$pos_isize = "$P[2]:$P[3]:$P[8]:f";
						}
					}
					if (!defined $BARC_POS_ISIZE{$barc}{$pos_isize} && !defined $OBSERVED{$P[0]}) {
						$BARC_POS_ISIZE{$barc}{$pos_isize} = 1;
						$KEEP{$P[0]} = 1;
						print OUT "$l\n";
						$BARC_kept{$barc}++;
						$total_kept++;
						$BAMID_included{$bamID}++;
					}
					$OBSERVED{$P[0]} = 1;
				} else {
					$reads_q10_to_other_chr++;
				}
			}
		} close IN;
		%KEEP = ();
	}
}

close OUT;

open OUT, ">$opt{'O'}.complexity.txt";
$rank = 1;
foreach $barc (sort {$BARC_kept{$b}<=>$BARC_kept{$a}} keys %BARC_kept) {
	$pct = sprintf("%.2f", ($BARC_kept{$barc}/$BARC_total{$barc})*100);
	print OUT "$rank\t$barc\t$BARC_total{$barc}\t$BARC_kept{$barc}\t$pct\n";
	$rank++;
} close OUT;

open OUT, ">$opt{'O'}.complexity.log";
$args = join(" ", @ARGV);
$ts = localtime(time);
$pct_off = sprintf("%.2f", ($reads_q10_to_other_chr/$q10_reads)*100);
$pct_kept = sprintf("%.2f", ($total_kept/$q10_reads)*100);
print OUT "$ts bam-rmdup completed. Args:\t$args
Q10 reads = $q10_reads
Aligning to other chroms: $reads_q10_to_other_chr ($pct_off %)
Total reads kept: $total_kept ($pct_kept %)
By Bam\n";
for ($bamID = 0; $bamID < @ARGV; $bamID++) {
	print OUT "\tBAMID=$bamID\t$ARGV[$bamID]\t$BAMID_included{$bamID} included reads.\n";
}
close OUT;

}
1;
