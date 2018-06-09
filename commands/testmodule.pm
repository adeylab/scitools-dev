package commands::testmodule;

use commands::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = "testmodule";

sub testmodule {
@ARGV = @_;
getopts("T:", \%opt);

$die2 = "
   this is the test module fail message. test option = -T [STR].
   ARGV0 is required as a string
";

if (!defined $ARGV[0]) {die $die2};

print "\nTEST MODULE ARGV0 = $ARGV[0]";
if (defined $opt{'T'}) {
	print ", opt T included as $opt{'T'}\n";
} else {
	print "\n";
}

# test using another scitools module here
#($COLOR_GRADIENT, $word) = load_gradient_defaults();
load_gradient_defaults();
load_defaults();
if (!defined $COLOR_GRADIENT{'PuOr'}) {print "DID NOT LOAD PuOr Gradient!\n"};
print "THIS IS PuOr gradient: $COLOR_GRADIENT{'PuOr'}\ncolor mapping is $color_mapping.\n";
####

read_annot($ARGV[1]);

print "ANNOT hello is $CELLID_annot{'hello'} and annot count is $annot_count\n";


exit;
}
1;
