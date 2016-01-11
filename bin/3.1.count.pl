#!/usr/bin/perl
%All=();
while(<>){
	s/\r?\n|\r//;
	$seq=<>;
	$seq =~ s/\r?\n|\r//;
	$All{length($seq)}++;
}
for(sort keys %All){
	print "$_\t$All{$_}\n";
}
