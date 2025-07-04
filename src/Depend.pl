#!/usr/bin/perl

# Syntax: Depend.pl files  (*.f)

# designed to be called from a master Makefile
# searches each F90 file for lines of form:  
#   (whitespace) USE module (whitespace)    or
#   (whitespace) use module (whitespace)
# takes all found modules and converts it into string of form:
#   file.o: module1.o module2.o ...
# writes that string into Makefile.depend
# end result is Makefile.depend with all dependencies for all F90 files

if (@ARGV == 0) {
    print "Syntax: Depend.pl files\n";
    exit;
}

open(MAKE,">Makefile.depend");

foreach $file (@ARGV) {

    # slurp entire file into @lines

    open(FILE,$file) || die "Could not open file $file\n";
    @lines = <FILE>;
    close(FILE);

    # extract all module names of form "USE module" or "use module"

    @modules = ();
    foreach (@lines) {
	if (/^\s*USE (\w*)\s*/) {push (@modules, $1)}
	if (/^\s*use (\w*)\s*/) {push (@modules, $1)}
    }

    # delete duplicate module names in @module
    # from Perl Cookbook, p 102, ++ operator adds unseen to hash

    %seen = ();
    @unique = ();
    foreach $module (@modules) {
	push (@unique, $module) unless $seen{$module}++;
    }

    # write one-line dependency string into Makefile.depend

    $file =~ s/\.f//;
    print MAKE "$file.o:\t";
    foreach $module (@unique) {print MAKE "$module.o ";};
    print MAKE "\n";
}
