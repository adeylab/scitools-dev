package sci_commands::annot_make;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("annot_make");

sub annot_make {

@ARGV = @_;
getopts("O:I:P:phD:dx", \%opt);

# DEFAULTS
@LETTERS = ("0", "A", "B", "C", "D", "E", "F", "G", "H");
%LETTER_NUM = ("A"=>"1", "B"=>"2", "C"=>"3", "D"=>"4", "E"=>"5", "F"=>"6", "G"=>"7", "H"=>"8");

$die2 = "
scitools annot-make [options] [annotation_description_1] [annotaiton_description_2] ...
   or    make-annot

Options:
   -O   [STR]   Output annotation file (default = STDOUT)
   -I   [STR]   Index file
         (default = $VAR{'SCI_index_file'})
         (Index names must be in form of: [Tier]_[set]_[i5/i7]_[A-H/1-12])
   -P   [STR]   Plate descriptor file (instead of written descriptors)
   -p           Print out sample plate file for modificaiton and exit
                (ExamplePlateDescriptor.csv)
   -D   [STR]   Dense plate descriptor
   -d           Print out example dense plate descriptor file
                (ExamplePlateDescriptor.dense.csv)
   -U   [STR]   Unique Dual Index Tn5 (UDI-Tn5) setup
   -u           Print example UDI-Tn5 sample sheet file
                (Example_UDI_Tn5_Sample_Sheet.csv)
   -h           More detailed description of plate / combo specification

";

$die3 = "

Annotation descriptors are provided as:

[Annotation_Name]+[Transposase or PCR descriptor]+[Transposase or PCR descriptor]+[etc...]
  must provide at least 1 transposase descriptor and at least 1 PCR descriptor for
  each annotation.

Transposase descriptor:
[TN5],[TN5 index set combination]=[row],[row]:[columns],[row],etc...
 or [NEX]

PCR Descriptor:
[PCR],[PCR index set combination]=[row]:[columns],[row],[row]:[columns],etc...

Before the \"=\" are two comma separated fields. The first is TN5/NEX or PCR to
define what stage of indexing is specified. The second is the index set IDs in
the order of i5 and then i7, (e.g. AA).

After the \"=\" are a series of comma separated fields, where each field
corresponds to a column of a plate, where columns are numbered 1-12. The column
field can be further specified with a subset of rows after a colon, where rows
are A-H.

Note: each index stage can be specified multiple times. The result is an all by
all of the transposase and PCR index sets that are specified.

Example:
My_Sample_1+NEX,AA=ALL+PCR,AC=1-8+PCR,AD=1:A-D,2-5,6:ACDFH,7-12 My_Sample_2+NEX,AA=ALL+NEX,BB=ALL+PCR,GE=1-8

   row specifications can be listed as a range OR letters with no spacing (e.g. ABCGH)
   commas in the column specification are ONLY for separating out columns NOT rows
      (i.e. 1,2,3:A,B would NOT be OK because of the comma between row letters)

";

if (defined $opt{'h'}) {die $die2.$die3};

if (defined $opt{'d'}) {
open OUT, ">ExamplePlateDescriptor.dense.csv";
print OUT "#INFO, This is the dense format for 'sci' plate descriptions.
#INFO, <- 'INFO' lines are ignored and can be used to comment.
#INFO, Each NEX section must start with '#NEX' then have a comma and then
#INFO, the two letter i5i7 combo for the index set. It then must be followed
#INFO, by a complete description of all wells, with the sample name entered
#INFO, into each well. ANy empty wells can be listed as 'empty' or 'null'
#INFO, ALL NEX sections must precede any #PCR sections!
#INFO, All NEX plates are assumed to be used for all PCR plates.
#INFO, A #PCR section has the header, then the 2-letter combo, then either
#INFO, 'all' for all wells of the plate, in which case no further lines are
#INFO, needed in the section, OR 'partial' and then 8 more rows with 12
#INFO, columns each, with a '0' for empty wells and '1' for used wells.
#INFO,
#NEX,AA
sample1,sample2,sample3,sample4,sample1,sample2,sample3,sample4
sample1,sample2,sample3,sample4,sample1,sample2,sample3,sample4
sample1,sample2,sample3,sample4,sample1,sample2,sample3,sample4
sample1,sample2,sample3,sample4,sample1,sample2,sample3,sample4
sample1,sample2,sample3,sample4,sample1,sample2,sample3,sample4
sample1,sample2,sample3,sample4,sample1,sample2,sample3,sample4
sample1,sample2,sample3,sample4,sample1,sample2,sample3,sample4
sample1,sample2,sample3,sample4,sample1,sample2,sample3,sample4
#PCR,AC,all
#PCR,AD,partial
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0";
close OUT;
exit;
}

if (defined $opt{'p'}) {
open OUT, ">ExamplePlateDescriptor.csv";
	print OUT "#NEX,MySampleID1,AA,Partial
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
#NEX,MySampleID1,BB,All
#PCR,MySampleID1,CE,Partial
1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,1,1,1,1,1,1
0,0,0,0,0,0,0,0,0,0,0,0
#NEX,MySampleID2,AA,Partial
0,0,0,0,0,0,1,1,1,1,1,1
0,0,0,0,0,0,1,1,1,1,1,1
0,0,0,0,0,0,1,1,1,1,1,1
0,0,0,0,0,0,1,1,1,1,1,1
0,0,0,0,0,0,1,1,1,1,1,1
0,0,0,0,0,0,1,1,1,1,1,1
0,0,0,0,0,0,1,1,1,1,1,1
0,0,0,0,0,0,1,1,1,1,1,1
#PCR,MySampleID2,EE,All
#PCR,MySampleID2,DF,Partial
1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,1,1,1,1,1,1
0,0,0,0,0,0,0,0,0,0,0,0
0,0,0,0,0,0,0,0,0,0,0,0
0,0,0,0,0,0,0,0,0,0,0,0
0,0,0,0,0,0,0,0,0,0,0,0\n";
close OUT;
exit;
}

if (defined $opt{'u'}) {
open OUT, ">Example_UDI_Tn5_Sample_Sheet.csv";
	print OUT "#INFO, this format assumes that i5 and i7 indexes are used
#INFO, going from A-H (i5) and 1-8 (i7). So each combination
#INFO, of index letters has 8 total samples, where the first listed
#INFO, is (i5 A, i7 1); then (i5 B, i7 2), etc... positions.
#INFO, PCR indexes are the plate combination, then a comma-separated
#INFO, list of the wells on the plate that are used.
#Nex,AA
sample1,sample2,sample3,sample4,sample5,sample6,sample7,sample8
#Nex,BB
sample9,sample10,sample11,empty,empty,empty,empty,empty
#PCR,CD
A1,A2\n";
close OUT;
exit;
}

if (!defined $ARGV[0] && !defined $opt{'P'} && !defined $opt{'D'} && !defined $opt{'U'}) {die $die2};

# Read in index file

if (!defined $opt{'I'}) {
    open IN, $VAR{'SCI_index_file'};
} else {
    open IN, $opt{'I'};
}

while ($l = <IN>) {
	chomp $l;
	($ID,$pos,$seq) = split(/\t/, $l);
	($tier,$set,$side,$wells) = split(/_/, $ID);
	if ($tier =~ /(Tn5|Nex)/i) {
		if ($side =~ /i5/) {
			$TN5SET_i5WELLS_seq{$set}{$wells} = $seq;
		} else {
			$TN5SET_i7WELLS_seq{$set}{$wells} = $seq;
		}
		$tier = "Tn5";
	} else {
		if ($side =~ /i5/) {
			$PCRSET_i5WELLS_seq{$set}{$wells} = $seq;
		} else {
			$PCRSET_i7WELLS_seq{$set}{$wells} = $seq;
		}
	}
} close IN;

if (defined $opt{'O'}) {open OUT, ">$opt{'O'}"};

# PROCESS ANNOT DESCRIPTOR FILE / INFO

if (defined $opt{'U'}) { # Tn5-UDI format
	open IN, "$opt{'U'}";
	while ($l = <IN>) {
		chomp $l;
		if ($l =~ /^#/) {
			@P = split(/,/, $l);		
			if ($P[0] =~ /(Tn5|Nex)/i) {
				($i5_set,$i7_set) = split(//, $P[1]);
				$l = <IN>; chomp $l;
				@SAMPLES = split(/,/, $l); unshift @SAMPLES, "0";
				for ($pos = 1; $pos <= 8; $pos++) {
					$annot = $SAMPLES[$pos];
					$rowLetter = $LETTERS[$pos];
					if ($annot ne "empty" && $annot ne "null" && $annot ne "0" && $annot ne "") {
						$pair = "$TN5SET_i5WELLS_seq{$i5_set}{$rowLetter},$TN5SET_i7WELLS_seq{$i7_set}{$pos}";
						$ANNOT_Tn5_pairs{$annot}{$pair} = 1;
					}
				}
			} elsif ($P[0] =~ /pcr/i) {
				($i5_set,$i7_set) = split(//, $P[1]);
				$l = <IN>; chomp $l;
				@PCR_SET = split(/,/, $l);
				for ($i = 0; $i < @PCR_SET; $i++) {
					$ix4 = $PCRSET_i5WELLS_seq{$i5_set}{substr($PCR_SET[$i],0,1)};
					$ix2 = $PCRSET_i7WELLS_seq{$i7_set}{substr($PCR_SET[$i],1)};
					foreach $annot (keys %ANNOT_Tn5_pairs) {
						foreach $pair (keys %{$ANNOT_Tn5_pairs{$annot}}) {
							($ix3,$ix1) = split(/,/, $pair);
							if (defined $opt{'O'}) {
								print OUT "$ix1$ix2$ix3$ix4\t$annot\n";
							} else {
								print "$ix1$ix2$ix3$ix4\t$annot\n";
							}
						}
					}
				}
			}
		}
	} close IN;
if (defined $opt{'O'}) {close OUT};
exit;
}

if (defined $opt{'D'}) { # dense plate descriptor format
	open IN, "$opt{'D'}";
	while ($l = <IN>) {
		chomp $l;
		if ($l =~ /^#/) {
			@P = split(/,/, $l);		
			if ($P[0] =~ /(Tn5|Nex)/i) {
				($i5_set,$i7_set) = split(//, $P[1]);
				for ($rowNum = 1; $rowNum <= 8; $rowNum++) {
					$row = <IN>; chomp $row; $rowLetter = $LETTERS[$rowNum];
					@ROW_COLS = split(/,/, $row); unshift @ROW_COLS, "0";
					for ($colNum = 1; $colNum <= 12; $colNum++) {
						$annot = $ROW_COLS[$colNum];
						if ($annot ne "empty" && $annot ne "null" && $annot ne "0" && $annot ne "") {
							$pair = "$TN5SET_i5WELLS_seq{$i5_set}{$rowLetter},$TN5SET_i7WELLS_seq{$i7_set}{$colNum}";
							$ANNOT_Tn5_pairs{$annot}{$pair} = 1;
							if (defined $opt{'x'}) {
								print STDERR "DEBUG: Row: $rowLetter, Col: $colNum, included as annot $annot\n";
							}
						} elsif (defined $opt{'x'}) {
							print STDERR "DEBUG: Row: $rowLetter, Col: $colNum, EXCLUDED ($annot)\n";
						}
					}
				}
			} elsif ($P[0] =~ /pcr/i) {
				($i5_set,$i7_set) = split(//, $P[1]);
				if ($P[2] =~ /^a/i) { # all PCR wells
					for ($rowNum = 1; $rowNum <= 8; $rowNum++) {
						$rowLetter = $LETTERS[$rowNum];
						for ($colNum = 1; $colNum <= 12; $colNum++) {
							$ix4 = $PCRSET_i5WELLS_seq{$i5_set}{$rowLetter};
							$ix2 = $PCRSET_i7WELLS_seq{$i7_set}{$colNum};
							foreach $annot (keys %ANNOT_Tn5_pairs) {
								foreach $pair (keys %{$ANNOT_Tn5_pairs{$annot}}) {
									($ix3,$ix1) = split(/,/, $pair);
									if (defined $opt{'O'}) {
										print OUT "$ix1$ix2$ix3$ix4\t$annot\n";
									} else {
										print "$ix1$ix2$ix3$ix4\t$annot\n";
									}
								}
							}
						}
					}
				} else { # partial PCR wells
					for ($rowNum = 1; $rowNum <= 8; $rowNum++) {
						$row = <IN>; chomp $row; $rowLetter = $LETTERS[$rowNum];
						@ROW_COLS = split(/,/, $row); unshift @ROW_COLS, "0";
						for ($colNum = 1; $colNum <= 12; $colNum++) {
							if ($ROW_COLS[$colNum]>0) {
								$ix4 = $PCRSET_i5WELLS_seq{$i5_set}{$rowLetter};
								$ix2 = $PCRSET_i7WELLS_seq{$i7_set}{$colNum};
								foreach $annot (keys %ANNOT_Tn5_pairs) {
									foreach $pair (keys %{$ANNOT_Tn5_pairs{$annot}}) {
										($ix3,$ix1) = split(/,/, $pair);
										if (defined $opt{'O'}) {
											print OUT "$ix1$ix2$ix3$ix4\t$annot\n";
										} else {
											print "$ix1$ix2$ix3$ix4\t$annot\n";
										}
									}
								}
							}
						}
					}
				}
			} # else it is description info & skip it
		} else {
			print STDERR "WARNING! A line was read that has no header associated with it! line: $l\n";
		}
	} close IN;
	if (defined $opt{'O'}) {close OUT};
	exit;
}

if (defined $opt{'P'}) { # plate descriptor, older version
	%NEX_ID_i5_i7_pair = ();
	%PCR_ID_i5_i7_pair = ();
	open IN, "$opt{'P'}";
	while ($l = <IN>) {
		chomp $l;
		if ($l =~ /^#/) {
			($class,$annot,$combo,$subset) = split(/,/, $l);;
			$class =~ s/^#//; $annot =~ s/ /_/g;
			if (!defined $ANNOT_flag{$annot}) {
				$ANNOT_flag{$annot} = $class;
			} else {
				$ANNOT_flag{$annot} .= ",$class";
			}
			($i5_set,$i7_set) = split(//, $combo);
			if ($subset =~ /all/i) {
				for ($rowNum = 1; $rowNum <= 8; $rowNum++) {
					$rowLetter = $LETTERS[$rowNum];
					for ($colNum = 1; $colNum <= 12; $colNum++) {
						if ($class =~ /(Tn5|Nex)/i) {
							$pair = "$TN5SET_i5WELLS_seq{$i5_set}{$rowLetter},$TN5SET_i7WELLS_seq{$i7_set}{$colNum}";
							$NEX_ID_i5_i7_pair{$annot}{$pair} = 1;
						} else {
							$pair = "$PCRSET_i5WELLS_seq{$i5_set}{$rowLetter},$PCRSET_i7WELLS_seq{$i7_set}{$colNum}";
							$PCR_ID_i5_i7_pair{$annot}{$pair} = 1;
						}
					}
				}
			} else {
				for ($rowNum = 1; $rowNum <= 8; $rowNum++) {
					$row = <IN>; chomp $row; $rowLetter = $LETTERS[$rowNum];
					@ROW_COLS = split(/,/, $row); unshift @ROW_COLS, "0";
					for ($colNum = 1; $colNum <= 12; $colNum++) {
						if ($ROW_COLS[$colNum]>0) {
							if ($class =~ /(Tn5|Nex)/i) {
								$pair = "$TN5SET_i5WELLS_seq{$i5_set}{$rowLetter},$TN5SET_i7WELLS_seq{$i7_set}{$colNum}";
								$NEX_ID_i5_i7_pair{$annot}{$pair} = 1;
							} else {
								$pair = "$PCRSET_i5WELLS_seq{$i5_set}{$rowLetter},$PCRSET_i7WELLS_seq{$i7_set}{$colNum}";
								$PCR_ID_i5_i7_pair{$annot}{$pair} = 1;
							}
						}
					}
				}
			}
		}
	} close IN;
	
	foreach $annot (keys %ANNOT_flag) {
		if ($ANNOT_flag{$annot} =~ /(Tn5|Nex)/i && $ANNOT_flag{$annot} =~ /pcr/i) {
			print STDERR "Printing $annot index combinations.\n";
			foreach $NEX_pair (keys %{$NEX_ID_i5_i7_pair{$annot}}) {
				($ix3,$ix1) = split(/,/, $NEX_pair);
				foreach $PCR_pair (keys %{$PCR_ID_i5_i7_pair{$annot}}) {
					($ix4,$ix2) = split(/,/, $PCR_pair);
					if (defined $opt{'O'}) {
						print OUT "$ix1$ix2$ix3$ix4\t$annot\n";
					} else {
						print "$ix1$ix2$ix3$ix4\t$annot\n";
					}
				}
			}
		} else {
			print STDERR "\nWARNING: Transposase (Nex/Tn5) AND PCR specifications must both be included for each annotation!\nBoth were not found for $annot! - SKIPPING!\n";
		}
	}
	if (defined $opt{'O'}) {close OUT};
	exit;
}

if (!defined $opt{'U'} && !defined $opt{'D'} && !defined $opt{'P'}) { # Longform annots
	foreach $annot_descriptor (@ARGV) {

		@DESCRIPTORS = split(/\+/, $annot_descriptor);
		$annot = shift(@DESCRIPTORS);
		%NEX_i5_i7_pair = ();
		%PCR_i5_i7_pair = ();
		
		print STDERR "\n##### PARSING ANNOT $annot #####\n";
		
		for ($i = 0; $i < @DESCRIPTORS; $i++) {
		
			print STDERR "\n##### Input Parse Descriptor: $DESCRIPTORS[$i] #####\n";
			($class,$columns) = split(/=/, $DESCRIPTORS[$i]);
			print STDERR "\tClass=$class, columns=$columns\n";
			($tier,$combo) = split(/,/, $class);
			print STDERR "\tTier=$tier, combo=$combo\n";
			($i5_set,$i7_set) = split(//, $combo);
			print STDERR "\ti5=$i5_set, i7=$i7_set\n";
			@COLUMN_DESCRIPTORS = split(/,/, $columns);
			for ($j = 0; $j < @COLUMN_DESCRIPTORS; $j++) {
				
				if ($COLUMN_DESCRIPTORS[$j] =~ /ALL/i) {$COLUMN_DESCRIPTORS[$j] = "1-12"};
				
				# determine if a subset of rows are specified (after :)
				if ($COLUMN_DESCRIPTORS[$j] =~ /:/) {
					($col,$rows) = split(/:/, $COLUMN_DESCRIPTORS[$j]);
				} else {
					$col = $COLUMN_DESCRIPTORS[$j];
					$rows = "A-H";
				}
				
				print STDERR "\tRows=$rows, which is:";
				# Add rows to the row list for the column(s)
				@ROW_LIST = split(//, $rows);
				@ROW_INCLUDE = ();
				for ($k = 0; $k < @ROW_LIST; $k++) {
					if ($ROW_LIST[($k+1)] eq "-") {
						$rowStart = $LETTER_NUM{$ROW_LIST[$k]};
						$k++; $k++;
						$rowEnd = $LETTER_NUM{$ROW_LIST[$k]};
						for ($l = $rowStart; $l <= $rowEnd; $l++) {
							push @ROW_INCLUDE, $LETTERS[$l];
							print STDERR " $LETTERS[$l]";
						}
					} else {
						push @ROW_INCLUDE, $ROW_LIST[$k];
						print STDERR " $ROW_LIST[$k]";
					}
				}
				
				# Go through the column(s)
				print STDERR "\n\tCols=$col, which is:";
				@COL_INCLUDE = ();
				if ($col =~ /-/) {
					($startCol,$endCol) = split(/-/, $col);
					for ($l = $startCol; $l <= $endCol; $l++) {
						push @COL_INCLUDE, $l;
						print STDERR " $l";
					}
				} else {
					push @COL_INCLUDE, $col;
					print STDERR " $col";
				}
				
				print STDERR "\n\t  --> Adding all pairs\n";
				# Now add all relevant barcodes to their respective groupings
				foreach $add_col (@COL_INCLUDE) {
					foreach $add_row (@ROW_INCLUDE) {
						if ($class =~ /(Tn5|Nex)/i) {
							$pair = "$TN5SET_i5WELLS_seq{$i5_set}{$add_row},$TN5SET_i7WELLS_seq{$i7_set}{$add_col}";
							$NEX_i5_i7_pair{$pair} = 1;
						} else {
							$pair = "$PCRSET_i5WELLS_seq{$i5_set}{$add_row},$PCRSET_i7WELLS_seq{$i7_set}{$add_col}";
							$PCR_i5_i7_pair{$pair} = 1;
						}
					}
				}
			}
		}

		foreach $NEX_pair (keys %NEX_i5_i7_pair) {
			($ix3,$ix1) = split(/,/, $NEX_pair);
			foreach $PCR_pair (keys %PCR_i5_i7_pair) {
				($ix4,$ix2) = split(/,/, $PCR_pair);
				if (defined $opt{'O'}) {
					print OUT "$ix1$ix2$ix3$ix4\t$annot\n";
				} else {
					print "$ix1$ix2$ix3$ix4\t$annot\n";
				}
			}
		}
	}
	if (defined $opt{'O'}) {close OUT};
	exit;
}

}
1;
