package sci_commands::annot_make_mode;

use sci_utils::general;
use sci_utils::modes;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("annot_make_mode");

sub annot_make_mode {

@ARGV = @_;
getopts("I:pM:m:v", \%opt);

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
%WELLID_letter = (); %WELLID_number = ();
$let = 1; $num = 1;
for ($wellID = 1; $wellID <= 96; $wellID++) {
	$WELLID_letter{$wellID} = $LETTERS[$let];
	$WELLID_number{$wellID} = $num;
	$num++;
	if ($num == 13) {
		$let++;
		$num = 1;
	}
}

$mode_name = "sci";
$def_read_name_format = "barc";

$die2 = "
scitools annot-make [options] [plate descriptor .csv file] [output prefix]
   or    make-annot

Note: If a multimodal format is specified, an index file will be created for
each modality.

Options:
   -I   [STR]   Index files/directory - comma separated (req)
         (def = DIR=$VAR{'index_directory'})
   -M   [STR]   Mode (def = $mode_name)
   -m   [STR]   Run format specification file
         (def = $VAR{'sci_modes'})
   -f   [STR]   Read name format: barc OR name (def = $def_read_name_format)
   -v           Verbose (for debugging)

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
if (defined $opt{'v'}) {print "INFO: Loading in mode $mode_name from config file $VAR{'sci_modes'}\n"};
$mode = modes->new($mode_name,$VAR{'sci_modes'});
# get set of modalities
@MODALITIES = $mode->modalities();
# get ordered index list for the mode
for ($modalityID = 0; $modalityID < @MODALITIES; $modalityID++) {
	$modality = $MODALITIES[$modalityID];
	@{$MODALITY_INDEXES{$modality}} = $mode->indexes($modality);
}

# Read in indexes
if (defined $opt{'v'}) {print "INFO: Reading in Index Directory: $opt{'I'}\n"};
read_indexdir($opt{'I'});
if (defined $opt{'v'}) {print "INFO: Checking Indexes that were Loaded\n"};
check_indexes();

# Read in plate descriptor file
if (defined $opt{'v'}) {print "INFO: Loading in plate file $ARGV[0]\n"};
load_plate_descriptions($ARGV[0]);

# Go through modalities then annots
foreach $modality (@MODALITIES) {
	if (defined $opt{'v'}) {print "INFO: Building matrix and printing annotations for modality: $modality\n"};
	if ($mode->check_multimodal() =~ /true/) {
		open OUT, ">$ARGV[1].$modality.annot";
	} else {
		open OUT, ">$ARGV[1].annot";
	}
	foreach $annot (keys %ANNOT_SETS) {
		if (defined $opt{'v'}) {print "INFO:\tBuilding for annot $annot\n"};
		@POS_INDEXES = (); @POS_IDS =();
		for ($index_pos = 0; $index_pos < @{$MODALITY_INDEXES{$modality}}; $index_pos++) {
			$index_type = $MODALITY_INDEXES{$modality}[$index_pos];
			$class = $INDEX_TYPE_class{$index_type};
			@{$POS_INDEXES[$index_pos]} = ();
			@{$POS_IDS[$index_pos]} = ();
			if (defined $opt{'v'}) {print "INFO:\t\tIndex position $index_pos, type=$index_type, class=$class\n"};
			if (defined $ANNOT_SETS{$annot}{$class}) {
				foreach $combo (keys %{$ANNOT_SETS{$annot}{$class}}) {
					if (defined $opt{'v'}) {print "INFO:\t\t\tCombo = $combo\n"};
					foreach $wellID (keys %{$ANNOT_SETS{$annot}{$class}{$combo}}) {
						if (defined $opt{'v'}) {print "INFO:\t\t\t\tWellID = $wellID\n"};
						if ($INDEX_CLASS_format{$class} =~ /96|all|plate/) {
							push @{$POS_INDEXES[$index_pos]}, $CLASS_COMBO_WELLID_seq{$class}{$combo}{$wellID};
							push @{$POS_IDS[$index_pos]}, $CLASS_COMBO_WELLID_id{$class}{$combo}{$wellID};
						} else {
							($i5_set,$i7_set) = split(//, $combo);
							if ($index_type =~ /i5/) {
								$well_name = $WELLID_letter{$wellID};
								push @{$POS_INDEXES[$index_pos]}, $CLASS_I5_COMBO_LETTER_seq{$class}{$combo}{$well_name};
								push @{$POS_IDS[$index_pos]}, $CLASS_I5_COMBO_LETTER_id{$class}{$combo}{$well_name};
							} elsif ($index_type =~ /i7/) {
								$well_name = $WELLID_number{$wellID};
								push @{$POS_INDEXES[$index_pos]}, $CLASS_I7_COMBO_NUMBER_seq{$class}{$combo}{$well_name};
								push @{$POS_IDS[$index_pos]}, $CLASS_I7_COMBO_NUMBER_id{$class}{$combo}{$well_name};
							}
						}
					}
				}
			}
		}
		if ($opt{'f'} =~ /barc/) {
			print_barcs($annot);
		} else {
			print_names($annot);
		}
	}
	close OUT;
}

#### SUBROUTINES ####

sub print_barcs {
	$annot_name = $_[0];
	%NEW = ();
	foreach $seq (@{$POS_INDEXES[0]}) {$NEW{$seq} = 1};
	if (defined $MODALITY_INDEXES{$modality}[1]) {
		for ($pos = 0; $pos < @{$MODALITY_INDEXES{$modality}}; $pos++) {
			%RUNNING = %NEW; %NEW = ();
			foreach $in_progress (keys %RUNNING) {
				foreach $seq (@{$POS_INDEXES[$pos]}) {
					$new = $in_progress.$seq;
					$NEW{$new} = 1;
				}
			}
		}
	}
	foreach $barc (keys %NEW) {
		print OUT "$barc\t$annot_name\n";
	}
}

sub print_names {
	$annot_name = $_[0];
	%NEW = ();
	foreach $seq (@{$POS_IDS[0]}) {$NEW{$seq} = 1};
	if (defined $MODALITY_INDEXES{$modality}[1]) {
		for ($pos = 0; $pos < @{$MODALITY_INDEXES{$modality}}; $pos++) {
			%RUNNING = %NEW; %NEW = ();
			foreach $in_progress (keys %RUNNING) {
				foreach $seq (@{$POS_IDS[$pos]}) {
					$new = $in_progress."-".$seq;
					$NEW{$new} = 1;
				}
			}
		}
	}
	foreach $barc (keys %NEW) {
		print OUT "$barc\t$annot_name\n";
	}
}

sub check_indexes {
	$absent_index = 0;
	foreach $check_modality (@MODALITIES) {
		if (defined $opt{'v'}) {print "INFO:\tChecking indexes for modality: $check_modality\n"};
		foreach $check_index (@{$MODALITY_INDEXES{$check_modality}}) {
			if (defined $opt{'v'}) {print "INFO:\t\tVerifying properties of index $check_index with mode $mode_name\n"};
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
			
			# build index matrix
			if (defined $opt{'v'}) {print "INFO:\t\tBuilding the index matrix for this index type's sequences:\n"};
			foreach $index_seq (keys %{$INDEX_TYPE_SEQ_id{$check_index}}) {
				if (defined $opt{'v'}) {print "INFO:\t\t\t$index_seq ... "};
				$class = $INDEX_TYPE_class{$check_index};
				$index_id = $INDEX_TYPE_SEQ_id{$check_index}{$index_seq};
				($id,$set,$side,$well_name) = split(/_/, $index_id);
				if ($INDEX_TYPE_format{$check_index} =~ /96|all|plate/) { # 96 well plate index set
					$wellID = $PLATE_FORMAT{$well_name};
					$combo = $set;
					$CLASS_COMBO_WELLID_seq{$class}{$combo}{$wellID} = $index_seq;
					$CLASS_COMBO_WELLID_id{$class}{$combo}{$wellID} = $index_id;
					if (defined $opt{'v'}) {print "96-well format, ID=$index_id, combo=$combo, wellID=$wellID\n"};
				} else { # rows by columns index set
					if ($side =~ /i5/) {
						$combo = $set;
						$CLASS_I5_COMBO_LETTER_seq{$class}{$combo}{$well_name} = $index_seq;
						$CLASS_I5_COMBO_LETTER_id{$class}{$combo}{$well_name} = $index_id;
						if (defined $opt{'v'}) {print "split format - i5 side, ID=$index_id, combo=$combo, well_name=$well_name\n"};
					} elsif ($side =~ /i7/) {
						$combo = $set;
						$CLASS_I7_COMBO_NUMBER_seq{$class}{$combo}{$well_name} = $index_seq;
						$CLASS_I7_COMBO_NUMBER_id{$class}{$combo}{$well_name} = $index_id;
						if (defined $opt{'v'}) {print "split format - i7 side, ID=$index_id, combo=$combo, well_name=$well_name\n"};
					} else { # undefined side
						die "ERROR: Cannot interpret the index position (i5 or i7) for the index $index_id.\nMake sure the format is properly specified in the index file and the index names are properly formatted:\n\t[ID]_[set]_[i5/i7]_[well]\n";
					}
				}
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
			$class =~ s/^#//;
			if ($class =~ /NEX/i) {$class = "sci_tn5"}; # for reverse compatability
			if ($class =~ /PCR/i) {$class = "sci_pcr"}; # for reverse compatability
			$annot =~ s/ /_/g;
			if (defined $opt{'v'}) {print "INFO:\tLoading $class, $annot, $combo, $subset\n"};
			@CLASS_PARTS = split(/_/, $class);
			$class = $CLASS_PARTS[0]."_".$CLASS_PARTS[1];
			
			if (!defined $INDEX_CLASS_format{$class}) {
				die "ERROR: There is no index format specified for index class $class.\n\tThe plate descriptor files must match the index files for the class of the indexes.\n\t(for at least the first two '_' separated fields)\n";
			}
			
			if ($subset =~ /all/i) {
				for ($wellID = 1; $wellID <= 96; $wellID++) {
					$ANNOT_SETS{$annot}{$class}{$combo}{$wellID} = 1;
					if (defined $opt{'v'}) {print "INFO:\t\tAdding wellID $wellID (annot=$annot,class=$class,combo=$combo)\n"};
				}
			} else {
				for ($rowNum = 1; $rowNum <= 8; $rowNum++) {
					$row = <CSV>; chomp $row;
					@ROW_COLS = split(/,/, $row); unshift @ROW_COLS, "0";
					for ($colNum = 1; $colNum <= 12; $colNum++) {
						if ($ROW_COLS[$colNum]!=0) {
							$well_name = $LETTERS[$rowNum].$colNum;
							$wellID = $PLATE_FORMAT{$well_name};
							$ANNOT_SETS{$annot}{$class}{$combo}{$wellID} = 1;
							if (defined $opt{'v'}) {print "INFO:\t\tAdding wellID $wellID (annot=$annot,class=$class,combo=$combo)\n"};
						}
					}
				}
			}
			
		}
	} close CSV;
	
}

}
1;