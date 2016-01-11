#!/usr/bin/perl

open IN,$ARGV[0];
%sd=();
while(<IN>){
	@fd = split/\t/;
	push @{$sd{$fd[1]}}, $fd[0];
}

print "".(join "/",@{$sd{$_}})."\t$_\t0\t".(length($_))."\n" for (keys %sd);
