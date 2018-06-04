#!/usr/bin/perl
use strict;
use warnings;

my $fn = @ARGV == 1 ? $ARGV[0] : '<in>';
my $ln = 1;

while (<>) {
  if (/^[ ]*\/\/\/ \\verboutput ([a-z0-9_]+)((?:[ ]+[a-z0-9_-]+)*)$/) {
    # \verboutput some_cmd -> output of some_cmd wrapped in verbatim
    if (-x $1) {
      my $out = `2>&1 ./$1$2`;
      $out =~ s/^/\/\/\/     /gm;
      print "\/\/\/\n", $out;
    } else {
      print STDERR "$fn:$ln: error: cannot execute `$1$2'\n";
    }
  } elsif (/^[ ]*\/\/\/ \\listoutput ([a-z0-9_]+)((?:[ ]+[a-z0-9_-]+)*)$/) {
    # \listoutput some_cmd -> output of some_cmd wrapped into list
    if (-x $1) {
      my $out = `2>&1 ./$1$2`;
      $out =~ s/^/\/\/\/ /gm;
      $out =~ s/^\/\/\/   (\S.*)$/\/\/\/ * $1\n\/\/\//gm;
      $out =~ s/^\/\/\/ (\S.*):$/\/\/\/ ## $1/gm;
      print "\/\/\/\n", $out;
    } else {
      print STDERR "$fn:$ln: error: cannot execute `$1$2'\n";
    }
  } elsif(/^[ ]*\/\//) {
    # //// -> empty line
    s/^[ ]*\/\/\/\/$//;
    print;
  } else {
    # _ -> CUDA compatible remark
    s/(^|[ ])_[ ]/ \/*!\\remark Compatible with CUDA*\//;
    print;
  }
  ++$ln;
}
