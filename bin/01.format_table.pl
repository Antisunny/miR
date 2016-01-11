#!/usr/bin/env perl
%TBL=();
while (<>) {
  s/\r?\n|\r//;
  @A = split/\t/;
  push $TBL{$A[1]},$A[0];
}
foreach  (keys %TBL) {
  print "".join'|',$TBL{$_}."\t$_\t0\t".length();
}
