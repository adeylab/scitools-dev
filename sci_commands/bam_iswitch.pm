package sci_commands::bam_iswitch;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("bam_iswitch");

sub bam_iswitch {

@ARGV = @_;
$| = 1;

# defaults
$max_in_set = 20;
$increment = 10000000;

getopts("O:M:Tr:s:V", \%opt);

$die2 = "

scitools bam-iswitch [rmdup position-sorted bam file]
   or    iswitch

Detect & remove index switch reads in standard sci bam files

Options:
   -O   [STR]   Output prefix (def = input bam)
   -M   [INT]   Max allowed reads at identical position (def = $max_in_set)
   -T           Include Tn5 barcodes (def = just test PCR (outer) indexes)
   -s   [STR]   Samtools call (def = $samtools)
   -V           Print progress to STDERR (def = only final)
   -r   [INT]   Reads to report progress at (forces -V, def = $increment)

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.bam$//};
if (defined $opt{'M'}) {$max_in_set = $opt{'M'}};
if (defined $opt{'r'}) {$increment = $opt{'r'}; $opt{'V'} = 1};
if (defined $opt{'s'}) {$samtools = $opt{'s'}};

$prev_pos = "";
@I5_IX = (); @I7_IX = (); @BARC = (); @READS = (); @PASS = (); @R12 = ();
$reads_processed = 0;
$report = $increment;
$i5_length = -1*$i5_length;
$positions_exceeding_max = 0;
$reads_within_jackpot_positions = 0;

open OUT, "| $samtools view -bS - > $opt{'O'}.iSwitch_scrubbed.bam";

open H, "$samtools view -H $ARGV[0] |";
while ($l = <H>) {print OUT "$l"};
close H;

open IN, "$samtools view $ARGV[0] |";

while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$barc = $P[0]; $barc =~ s/:.+$//;
	$BARC_input_reads{$barc}++;
	$total_input++;
	$pos = $P[2].":".$P[3];
	if ($pos eq $prev_pos) {
		if ($read_set_count<=$max_in_set) {
			if (defined $opt{'T'}) {
				$i5_ix = substr($barc,18,18);
				$i7_ix = substr($barc,0,18);
			} else {
				$i5_ix = substr($barc,28,10);
				$i7_ix = substr($barc,8,10);
			}
			if ($P[1] & 64) {push @R12, 1} else {push @R12, 2};
			push @I5_IX, $i5_ix;
			push @I7_IX, $i7_ix;
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
						if ($I5_IX[$i] eq $I5_IX[$j] && $I7_IX[$i] ne $I7_IX[$j]) {
							$BARC_i5_switch{$BARC[$i]}++;
							$BARC_i5_switch{$BARC[$j]}++;
							$i5_switches++;
							$PASS[$i] = 1; $PASS[$j] = 1;
						} elsif ($I7_IX[$i] eq $I7_IX[$j] && $I5_IX[$i] ne $I5_IX[$j]) {
							$BARC_i7_switch{$BARC[$i]}++;
							$BARC_i7_switch{$BARC[$j]}++;
							$i7_switches++;
							$PASS[$i] = 1; $PASS[$j] = 1;
						}
					}
				}
				if ($PASS[$i]==0) {
					print OUT "$READS[$i]\n";
					$BARC_passing_reads{$BARC[$i]}++;
					$total_passing++;
				}
			}
		} elsif ($read_set_count < 2) {
			print OUT "$l\n";
			$BARC_passing_reads{$barc}++;
			$total_passing++;
		} elsif ($read_set_count>$max_in_set) {
			$positions_exceeding_max++;
			$reads_within_jackpot_positions+=$read_set_count;
		}
		@I5_IX = (); @I7_IX = (); @BARC = (); @READS = (); @PASS = (); @R12 = ();
		$read_set_count = 0;
	}
	$prev_pos = $pos;
	$reads_processed++;
	if ($reads_processed >= $report && defined $opt{'V'}) {
		$ts = localtime(time);
		$percent_passing = sprintf("%.2f", $total_passing/$total_input);
		print STDERR "
$ts	$reads_processed reads processed
	position = $pos
	$i5_switches i5 switches, $i7_switches i7 switches
	Total retained = $total_passing ($percent_passing)
	Positions exceeding max ($max_in_set) = $positions_exceeding_max
	Reads in jackpot positions = $reads_within_jackpot_positions
";
		$report+=$increment;
	}
} close IN; close OUT;

open LOG, ">$opt{'O'}.iSwitch_progress.log";
$ts = localtime(time);
$percent_passing = sprintf("%.2f", $total_passing/$total_input);
print LOG "
$ts	COMPLETE! $reads_processed total reads processed
	position = $pos
	$i5_switches i5 switches, $i7_switches i7 switches
	Total retained = $total_passing ($percent_passing)
	Positions exceeding max ($max_in_set) = $positions_exceeding_max
	Reads in jackpot positions = $reads_within_jackpot_positions\n";
close LOG;

open OUT, ">$opt{'O'}.iSwitch_stats.txt";
foreach $barc (sort {$BARC_passing_reads{$b}<=>$BARC_passing_reads{$a}} keys %BARC_passing_reads) {
	$percent_passing = sprintf("%.2f", $BARC_passing_reads{$barc}/$BARC_input_reads{$barc});
	print OUT "$barc\t$BARC_input_reads{$barc}\t$BARC_passing_reads{$barc}\t$percent_passing\n";
} close OUT;

}
1;
