#!/usr/bin/perl

print "CONFIGURING SCITOOLS, PWD = $ENV{'PWD'}\n";

open SCITOOLS, "$ENV{'PWD'}/scitools-src" || die "ERROR: Cannot find 'scitools' in $ENV{'PWD'}\n";
print "   configuring scitools executable ...\n";
open MODDED, ">$ENV{'PWD'}/scitools";
while ($l = <SCITOOLS>) {
	chomp $l;
	if ($l =~ /#LIB#$/) {
		$l = "use lib \"$ENV{'PWD'}\"; #LIB#";
	} elsif ($l =~ /#CONFIG#/) {
		$l = "\$SCITOOLS_DEFAULTS = \"$ENV{'PWD'}/scitools.cfg\"; #CONFIG#";
	}
	print MODDED "$l\n";
} close SCITOOLS; close MODDED;
system("chmod +x $ENV{'PWD'}/scitools");

print "   configuring modules within $ENV{'PWD'}/commands ...\n";

opendir COMMANDS, "$ENV{'PWD'}/commands" || die "ERROR: Cannot open directory: $ENV{'PWD'}/commands!\n";
while (readdir COMMANDS) {
	if ($_ =~ /\.pm$/) {
		print "      $_ ...\n";
		open MODULE, "$ENV{'PWD'}/commands/$_";
		open MODDED, ">$ENV{'PWD'}/commands/$_.cfg";
		while ($l = <MODULE>) {
			chomp $l;
			if ($l =~ /#LIB#$/) {
				$l = "use lib \"$ENV{'PWD'}\"; #LIB#";
			}
			print MODDED "$l\n";
		} close MODULE; close MODDED;
		system("rm -f $ENV{'PWD'}/commands/$_ && mv $ENV{'PWD'}/commands/$_.cfg $ENV{'PWD'}/commands/$_");
	}
} closedir COMMANDS;

print "CONFIGURATION COMPLETE!\n";