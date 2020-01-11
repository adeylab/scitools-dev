package sci_commands::bam_iswitch;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("bam_iswitch");

sub bam_iswitch {

@ARGV = @_;
$| = 1;

# defaults
$max_in_set = 100;
$increment = 10000000;
$filter = "T";

getopts("O:M:T:r:s:VN:", \%opt);

$die2 = "

scitools bam-iswitch [rmdup position-sorted bam file]
   or    iswitch

Detect & remove index switch reads in standard sci bam files

Options:
   -O   [STR]   Output prefix (def = input bam)
   -M   [INT]   Max allowed reads at identical position (def = $max_in_set)
   -T   [STR]   Switch filter(s), comma separated: (def = $filter)
                  T  = Matching both PCR, one or both Tn5 mismatches
                  P  = Matching both Tn5, one or both PCR index mismatches
                  7  = i5 side matching both, one or both i7 mismatch
                  5  = i7 side matching both, one or both i5 mismatch
                  A  = All of the above are filtered out (specify alone)
   -N   [INT]   Only test N reads and report (no bam output; def = filter all)
   -s   [STR]   Samtools call (def = $samtools)
   -V           Print progress to STDERR (def = only final)
   -r   [INT]   Reads to report progress at (forces -V, def = $increment)

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.bam$//};
if (defined $opt{'M'}) {$max_in_set = $opt{'M'}};
if (defined $opt{'r'}) {$increment = $opt{'r'}; $opt{'V'} = 1};
if (defined $opt{'s'}) {$samtools = $opt{'s'}};
if (defined $opt{'T'}) {$filter = $opt{'T'}};

# To Do: Add in filter for matching both PCR and either Tn5 with other Tn5 mismatch and vice versa

$full_filters = "T,P,7,5";
if ($filter =~ /A/i) {$filter = $full_filters};
@FILT_LIST = split(/,/, $filter);
foreach $filt (@FILT_LIST) {
	$FILT_CHECK{$filt} = 1;
}
@FULL_FILT_LIST = split(/,/, $full_filters);
foreach $filt (@FULL_FILT_LIST) {
	$FILT_COUNT{$filt} = 0;
}

$prev_pos = "";
@I5_P = (); @I5_T = (); @I7_P = (); @I7_T = ();
@BARC = (); @READS = (); @PASS = (); @R12 = ();
$reads_processed = 0;
$report = $increment;
$positions_exceeding_max = 0;
$reads_within_jackpot_positions = 0;
$fail_by_switching = 0;

if (!defined $opt{'N'}) {
	open OUT, "| $samtools view -bS - > $opt{'O'}.iSwitch_scrubbed.bam";

	open H, "$samtools view -H $ARGV[0] |";
	while ($l = <H>) {print OUT "$l"};
	close H;
}

open IN, "$samtools view $ARGV[0] |";

while ($l = <IN>) {
	if (defined $opt{'N'} && $total_input >= $opt{'N'}) {
		last;
	}
	chomp $l;
	@P = split(/\t/, $l);
	$barc = $P[0]; $barc =~ s/:.+$//;
	$BARC_input_reads{$barc}++;
	$total_input++;
	$pos = $P[2].":".$P[3];
	if ($pos eq $prev_pos) {
		if ($read_set_count<=$max_in_set) {
			$i7_t = substr($barc,0,8);
			$i7_p = substr($barc,8,10);
			$i5_t = substr($barc,18,8);
			$i5_p = substr($barc,26,10);
			if ($P[1] & 64) {push @R12, 1} else {push @R12, 2};
			push @I7_T, $i7_t; push @I7_P, $i7_p;
			push @I5_T, $i5_t; push @I5_P, $i5_p;
			push @BARC, $barc;
			push @READS, $l;
			push @PASS, 0;
		}
		$read_set_count++;
	} else {
		if ($read_set_count>1.5&&$read_set_count<=$max_in_set) {
			for ($i = 0; $i < @PASS; $i++) {
				for ($j = ($i+1); $j < @PASS; $j++) {
					if ($R12[$i] eq $R12[$j]) {
						# T filter:
						if (($I7_P[$i] eq $I7_P[$j] &&
							 $I5_P[$i] eq $I5_P[$j]) &&
							($I7_T[$i] ne $I7_T[$j] ||
							 $I5_T[$i] ne $I5_T[$j])) {
							$FILT_COUNT{'T'}++;
							#$BARC_FILT_COUNT{$BARC[$i]}{'T'}++;
							#$BARC_FILT_COUNT{$BARC[$j]}{'T'}++;
							if (defined $FILT_CHECK{'T'}) {
								$PASS[$i] = 1; $PASS[$j] = 1;
							}
						}
						# P filter:
						if (($I7_T[$i] eq $I7_T[$j] &&
							 $I5_T[$i] eq $I5_T[$j]) &&
							($I7_P[$i] ne $I7_P[$j] ||
							 $I5_P[$i] ne $I5_P[$j])) {
							$FILT_COUNT{'P'}++;
							#$BARC_FILT_COUNT{$BARC[$i]}{'P'}++;
							#$BARC_FILT_COUNT{$BARC[$j]}{'P'}++;
							if (defined $FILT_CHECK{'P'}) {
								$PASS[$i] = 1; $PASS[$j] = 1;
							}
						}
						# 7 filter:
						if (($I5_P[$i] eq $I5_P[$j] &&
							 $I5_T[$i] eq $I5_T[$j]) &&
							($I7_T[$i] ne $I7_T[$j] ||
							 $I7_P[$i] ne $I7_P[$j])) {
							$FILT_COUNT{'7'}++;
							#$BARC_FILT_COUNT{$BARC[$i]}{'7'}++;
							#$BARC_FILT_COUNT{$BARC[$j]}{'7'}++;
							if (defined $FILT_CHECK{'7'}) {
								$PASS[$i] = 1; $PASS[$j] = 1;
							}
						}
						# 5 filter:
						if (($I7_P[$i] eq $I7_P[$j] &&
							 $I7_T[$i] eq $I7_T[$j]) &&
							($I5_T[$i] ne $I5_T[$j] ||
							 $I5_P[$i] ne $I5_P[$j])) {
							$FILT_COUNT{'5'}++;
							#$BARC_FILT_COUNT{$BARC[$i]}{'5'}++;
							#$BARC_FILT_COUNT{$BARC[$j]}{'5'}++;
							if (defined $FILT_CHECK{'5'}) {
								$PASS[$i] = 1; $PASS[$j] = 1;
							}
						}
					}
				}
				if ($PASS[$i]==0) {
					if (!defined $opt{'N'}) {print OUT "$READS[$i]\n"};
					$BARC_passing_reads{$BARC[$i]}++;
					$total_passing++;
				} else {
					$fail_by_switching++;
				}
			}
		} elsif ($read_set_count < 2) {
			if (!defined $opt{'N'}) {print OUT "$l\n"};
			$BARC_passing_reads{$barc}++;
			$total_passing++;
		} elsif ($read_set_count>$max_in_set) {
			$positions_exceeding_max++;
			$reads_within_jackpot_positions+=$read_set_count;
		}
		@I5_P = (); @I5_T = (); @I7_P = (); @I7_T = ();
		@BARC = (); @READS = (); @PASS = (); @R12 = ();
		$read_set_count = 0;
	}
	$prev_pos = $pos;
	$reads_processed++;
	if ($reads_processed >= $report && defined $opt{'V'}) {
		$ts = localtime(time);
		$percent_passing = sprintf("%.2f", $total_passing/$total_input);
		$percent_jackpot = sprintf("%.2f", $reads_within_jackpot_positions/$total_input);
		$percent_switch = sprintf("%.2f", $fail_by_switching/$total_input);
		print STDERR "
$ts	$reads_processed reads processed
	Current position = $pos
	Matching PCR, Mismatching one or both Tn5 = $FILT_COUNT{'T'}
	Matching Tn5, Mismatching one or more PCR = $FILT_COUNT{'P'}
	Matching all i7, Mismatching one or both i5 = $FILT_COUNT{'5'}
	Matching all i5, Mismatching one or both i7 = $FILT_COUNT{'7'}
	Total reads filtered as index switching = $fail_by_switching ($percent_switch)
	Positions exceeding max ($max_in_set) = $positions_exceeding_max
	Reads in jackpot positions = $reads_within_jackpot_positions ($percent_jackpot)
	Total retained = $total_passing ($percent_passing)
";
		$report+=$increment;
	}
} close IN;
if (!defined $opt{'N'}) {close OUT};

open LOG, ">$opt{'O'}.iSwitch_summary.txt";
$ts = localtime(time);
$percent_passing = sprintf("%.2f", $total_passing/$total_input);
print LOG "
$ts	$reads_processed reads processed
	Matching PCR, Mismatching one or both Tn5 = $FILT_COUNT{'T'}
	Matching Tn5, Mismatching one or more PCR = $FILT_COUNT{'P'}
	Matching all i7, Mismatching one or both i5 = $FILT_COUNT{'5'}
	Matching all i5, Mismatching one or both i7 = $FILT_COUNT{'7'}
	Total reads filtered as index switching = $fail_by_switching ($percent_switch)
	Positions exceeding max ($max_in_set) = $positions_exceeding_max
	Reads in jackpot positions = $reads_within_jackpot_positions ($percent_jackpot)
	Total retained = $total_passing ($percent_passing)\n";
close LOG;

open OUT, ">$opt{'O'}.iSwitch_stats.txt";
foreach $barc (sort {$BARC_passing_reads{$b}<=>$BARC_passing_reads{$a}} keys %BARC_passing_reads) {
	$percent_passing = sprintf("%.2f", $BARC_passing_reads{$barc}/$BARC_input_reads{$barc});
	print OUT "$barc\t$BARC_input_reads{$barc}\t$BARC_passing_reads{$barc}\t$percent_passing\n";
} close OUT;

}
1;
