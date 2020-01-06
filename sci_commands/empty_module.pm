package sci_commands::empty_module;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("empty_module");

sub empty_module {

@ARGV = @_;

# defaults

getopts("O:", \%opt);

$die2 = "
scitools [command] [arguments]

[Description]

Options:
   -O   [STR]   Output prefix (def = [def])

";

}
1;
