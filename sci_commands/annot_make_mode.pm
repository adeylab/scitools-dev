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

Note: If a multimodal format is specified, an index file will be created for
each modality.

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
load_plate_descriptions($ARGV[0]);

# Go through modalities then annots
foreach $modality (@MODALITIES) {
	if ($mode->check_multimodal() =~ /true/) {
		open OUT, ">$ARGV[1].$modality.annot";
	} else {
		open OUT, ">$ARGV[1].annot";
	}
	foreach $annot (keys %ANNOT_SETS) {
		@POS_INDEXES = ();
		for ($index_pos = 0; $index_pos < @{$MODALITY_INDEXES{$modality}}; $index_pos++) {
			$index_type = $MODALITY_INDEXES{$modality}[$index_pos];
			$class = $INDEX_TYPE_class{$index_type};
			if (defined $ANNOT_SETS{$annot}{$class}) {
				foreach $combo (keys %{$ANNOT_SETS{$annot}{$class}}) {
					
				}
			}
		}
	}
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
			if (!defined $INDEX_TYPE_class{$check_index}) {
				die "ERROR: The index type ($check_index) does not have a class specified! This is specified in the index files as a header.\n";
			}
			if (!defined $INDEX_TYPE_format{$check_index}) {
				die "ERROR: The index type ($check_index) does not have a format specified! This is specified in the index files as a header.\n";
			}
			
			# build index matrix
			foreach $index_seq (keys %{$INDEX_TYPE_SEQ_id{$check_index}}) {
				$index_id = $INDEX_TYPE_SEQ_id{$check_index}{$index_seq};
				###########################################################
			}
		}
	}
	if ($absent_index>0) {die "ERROR: Index types specified in mode: $mode_name could not be found in the specified index file(s)/directory(s)!\n"};
}

sub plate2seq {
	
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
			@CLASS_PARTS = split(/_/, $class);
			$class = $CLASS_PARTS[0]."_".$CLASS_PARTS[1];
			
			if (!defined $INDEX_CLASS_format{$class}) {
				die "ERROR: There is no index format specified for index class $class.\n\tThe plate descriptor files must match the index files for the class of the indexes.\n\t(for at least the first two '_' separated fields)\n";
			}
			
			
			if ($subset =~ /all/i) {
				for ($wellID = 1; $wellID <= 96; $wellID++) {
					$ANNOT_SETS{$annot}{$class}{$combo}{$wellID} = 1;
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
						}
					}
				}
			}
			
		}
	} close CSV;
	
}

}
1;