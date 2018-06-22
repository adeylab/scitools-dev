package sci_commands::annot_make_mode;

use sci_utils::general;
use sci_utils::modes;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("annot_make_mode");

sub annot_make_mode {

@ARGV = @_;
getopts("I:pM:m:", \%opt);

# DEFAULTS
@LETTERS = ("0", "A", "B", "C", "D", "E", "F", "G", "H");
%LETTER_NUM = ("A"=>"1", "B"=>"2", "C"=>"3", "D"=>"4", "E"=>"5", "F"=>"6", "G"=>"7", "H"=>"8");
# build the 96 well plate mapping for 1-96, 01-96, A1, and A01 well names
%PLATE_FORMAT = (); $wellID = 1;
for ($letter_pos = 1; $letter_pos < @LETTERS; $letter_pos++) {
	$letter = $LETTERS[$letter_pos];
	for ($num = 1; $num <= 12; $num++) {
		$well_name = $letter.$num;
		$PLATE_FORMAT{$well_name} = $wellID;
		if ($num<=9) {
			$well_name = $letter."0".$num;
			$PLATE_FORMAT{$well_name} = $wellID;
		}
		$wellID++;
	}
}
for ($num = 1; $num <= 96; $num++) {
	$PLATE_FORMAT{$num} = $num;
	if ($num<=9) {$PLATE_FORMAT{$num} = "0".$num};
}

$mode_name = "sci";
$def_read_name_format = "barc";

$die2 = "
scitools annot-make [options] [plate descriptor .csv file] [output prefix]
   or    make-annot

Options:
   -I   [STR]   Index files/directory - comma separated (req)
         (def = DIR=$VAR{'index_directory'})
   -M   [STR]   Mode (def = $mode_name)
   -m   [STR]   Run format specification file
         (def = $VAR{'sci_modes'})
   -f   [STR]   Read name format: barc OR name (def = $def_read_name_format)

   -p           Print out sample plate file for modificaiton and exit
                (ExamplePlateDescriptor.csv)

";

if (defined $opt{'p'}) {
open OUT, ">ExamplePlateDescriptor.csv";
	print OUT "#sci_tn5,MySampleID1,AA,Partial
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
#sci_tn5,MySampleID1,BB,All
#sci_pcr,MySampleID1,CE,Partial
1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,1,1,1,1,1,1
0,0,0,0,0,0,0,0,0,0,0,0
#sci_tn5,MySampleID2,AA,Partial
0,0,0,0,0,0,1,1,1,1,1,1
0,0,0,0,0,0,1,1,1,1,1,1
0,0,0,0,0,0,1,1,1,1,1,1
0,0,0,0,0,0,1,1,1,1,1,1
0,0,0,0,0,0,1,1,1,1,1,1
0,0,0,0,0,0,1,1,1,1,1,1
0,0,0,0,0,0,1,1,1,1,1,1
0,0,0,0,0,0,1,1,1,1,1,1
#sci_pcr,MySampleID2,EE,All
#sci_pcr,MySampleID2,DF,Partial
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

# parse options
if (!defined $ARGV[1]) {die $die2};
if (!defined $opt{'I'}) {$opt{'I'} = "DIR=$VAR{'index_directory'}"};
if (defined $opt{'M'}) {$mode_name = $opt{'M'}};
if (!defined $opt{'f'}) {$opt{'f'} = $def_read_name_format};
if (defined $opt{'m'}) {$VAR{'sci_modes'} = $opt{'m'}};

# load in sci mode
$mode = modes->new($mode_name,$VAR{'sci_modes'});
# get set of modalities
@MODALITIES = $mode->modalities();
# get ordered index list for the mode
for ($modalityID = 0; $modalityID < @MODALITIES; $modalityID++) {
	$modality = $MODALITIES[$modalityID];
	@{$MODALITY_INDEXES{$modality}} = $mode->indexes($modality);
}

# Read in indexes
read_indexdir($opt{'I'});
check_indexes();

# Read in plate descriptor file
$annotation_sets = load_plate_descriptions($ARGV[0]);








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
			if (!defined $INDEX_TYPE_class{$check_index}) {
				die "ERROR: The index type ($check_index) does not have a class specified! This is specified in the index files as a header.\n";
			}
			if (!defined $INDEX_TYPE_format{$check_index}) {
				die "ERROR: The index type ($check_index) does not have a format specified! This is specified in the index files as a header.\n";
			}
		}
	}
	if ($absent_index>0) {die "ERROR: Index types specified in mode: $mode_name could not be found in the specified index file(s)/directory(s)!\n"};
}

sub load_plate_descriptions { # eg: #NEX,MySampleID1,AA,Partial
	
	%ANNOT_SETS = ();
	open CSV, "$_[0]";
	while ($l = <CSV>) {
		chomp $l;
		if ($l =~ /^#/) {
		
			($class,$annot,$combo,$subset) = split(/,/, $l);
			$class =~ s/^#//; $annot =~ s/ /_/g;
			@CLASS_PARTS = split(/_/, $class);
			$class = $CLASS_PARTS[0]."_".$CLASS_PARTS[1];
			
			if (!defined $INDEX_CLASS_format{$class}) {
				die "ERROR: There is no index format specified for index class $class.\n\tThe plate descriptor files must match the index files for the class of the indexes.\n\t(for at least the first two '_' separated fields)\n";
			}
			
			# determine if it is a plate based or row/col based index set
			if ($INDEX_CLASS_format{$class} =~ /96|all|plate/) {
				if ($subset =~ /all/i) {
					for ($wellID = 1; $wellID <= 96; $wellID++) {
						$ANNOT_SETS{$annot}{$class}{#####################################################
					}
				} else {
					##################################################
				}
			} else {
				($i5_set,$i7_set) = split(//, $combo);
				
				
			}
			
			if ($subset =~ /all/i) {
				for ($rowNum = 1; $rowNum <= 8; $rowNum++) {
					$rowLetter = $LETTERS[$rowNum];
					for ($colNum = 1; $colNum <= 12; $colNum++) {
						
					}
				}
			} else {
				for ($rowNum = 1; $rowNum <= 8; $rowNum++) {
					$row = <CSV>; chomp $row; $rowLetter = $LETTERS[$rowNum];
					@ROW_COLS = split(/,/, $row); unshift @ROW_COLS, "0";
					for ($colNum = 1; $colNum <= 12; $colNum++) {
						if ($ROW_COLS[$colNum]>0) {
							
						}
					}
				}
			}
			
		}
	} close CSV;
	
	return \%ANNOT_SETS;
}


#### OLD


if (defined $opt{'P'}) {
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
	
}

}
1;