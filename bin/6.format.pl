#!/usr/bin/perl

 open IN, $ARGV[0] or die "Open Error\n";
 %TBL = ();
 while (<IN>) {
    @fd = split/\t/ if /^ath/;
    $a = $fd{0} =~ s/-[0-9]p$//r;# =~ s/[a-z](\.[0-9])?$//r;
    $TBL{$a} = [$fd[1],$fd[2],$fd[3]] if ! $TBL{$a};
    if 
 }
