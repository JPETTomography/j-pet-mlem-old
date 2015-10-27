#!/usr/bin/perl
use strict;
use warnings;

while (<>) {
  if (/^[ ]*\/\/\/ \\verboutput ([a-z0-9_]+)$/) {
    # \verboutput some_cmd -> output of some_cmd wrapped in verbatim
    if (-x $1) {
      my $out = `2>&1 ./$1`;
      $out =~ s/^/\/\/\/     /gm;
      print "\/\/\/\n", $out;
    } else {
      print STDERR "error: Cannot execute `$1'.\n";
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
}
