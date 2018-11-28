package sci_commands::annot_collapse;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("annot_collapse");

sub annot_collapse {

@ARGV = @_;

getopts("O:i", \%opt);

$die2 = "
scitools annot-collapse [options] [annot file] [group1=annot1,annot2,etc...] (group2=annotN) etc...
   
Will collapse clusters based on specified groupings.
Note - if the groups and annots are 1:1 it works as renaming

Options:
   -O   [STR]   Output prefix (def = annot file .collapsed.annot)
   -i           If an annot is not specified, include it as-is
                (def = exclude from output)


";

if (!defined $ARGV[1]) {die $die2};
if (!defined $opt{'O'}) {
	$opt{'O'} = $ARGV[0];
	$opt{'O'} =~ s/annot/collapsed/;
}

# parse group specifications
for ($i = 1; $i < @ARGV; $i++) {
	@ANNOT_LIST = split(/[=,]/, $ARGV[$1]);
	$groupID = shift(@ANNOT_LIST);
	foreach $annot (@ANNOT_LIST) {
		$ANNOT_group{$annot} = $groupID;
	}
}

if (-e "$opt{'O'}.annot") {die "ERROR: $opt{'O'}.annot already exists - will not overwrite!\n"};
open OUT, ">$opt{'O'}.annot";
open IN, $ARGV[0];
while ($l = <IN>) {
	($cellID,$annot) = split(/\t/, $l);
	if (defined $opt{'i'} && !defined $ANNOT_group{$annot}) {
		$ANNOT_group{$annot} = $annot;
	}
	if (defined $ANNOT_group{$annot}) {
		print OUT "$cellID\t$ANNOT_group{$annot}\n";
	}
}
close IN; close OUT;

}

1;