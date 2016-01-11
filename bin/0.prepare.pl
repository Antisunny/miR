#!/usr/bin/perl

use Getopt::Long;
use Pod::Usage;

GetOptions(
"mature|m=s" => \my $mature_in,
"hairpin|h=s" => \my $hairpin_in,
"outdir|o=s" => \my $out_dir,
"help" => \my $help
) or die pod2usage(-verbose=>1);
# pre-pro
die "\e[01;31moutpout dir is required\e[00m\n" unless $out_dir;
$mature_out = $out_dir."/athMature.fa";
$hairpin_out = "$out_dir/athHairpin.txt";
$hairpin_fa = "$out_dir/athHairpin.fa";

# get athMature.txt from mature.fa downloaded from miRbase, which contains mature miRNA seqs of several species.
# the new file contains only seq of $sepcies , FASTA format
open MIN,       $mature_in  or die "\e[01;31mError in reading $mature_in\n$!\e[00m\n";
open MOUT, '>', $mature_out or die "\e[01;31mError writing to $mature_out\n$!\e[00m\n";
$species = "Arabidopsis thaliana";
$pattern = qr/^>[A-Za-z0-9\- \.]+$species/;
while (<MIN>) {
    if(/$pattern/){
		$header = $_ =~ s/miR/MIR/r;
        print MOUT $header;
        my $seq  = <MIN>;
        print MOUT $seq;
    }
}
close MIN;
close MOUT;


# get athHairpin.txt from hairpin.fa 
# make hairpin.fa in FASTA format seq in one line, not multiline as in fasta format
open HIN,       $hairpin_in  or die "\e[01;31mError in reading $hairpin_in\n$!\e[00m\n";
open HOUT, '>', $hairpin_out or die "\e[01;31mError writing to $hairpin_out\n$!\e[00m\n";
open HOFA, ">", $hairpin_fa or die "\e[01;31mError writing to $hairpin_fa\n$!\e[00m\n";
$fas_seq = '';
$fas_fa = '';
$hit =1;
while (<HIN>){
   if (/$pattern/ ){
       print HOUT "$fas_seq\n" if $fas_seq;
	   print HOFA $fas_fa if $fas_fa;
       print HOUT $_;
	   print HOFA $_;
       $hit = 0;
       $fas_seq = '';
	   $fas_fa = '';
       next
   }
   if ($_ !~ /^>/ && $hit == 0){
	   $fas_fa .= $_;
       s/\r?\n|\r//;
       $fas_seq .= $_;
   }
}
close HIN;
close HOUT;

=head1
will get
athMature.fa # only of a.th.
athHairpin.fa # only of a.th.
athHairpin.txt # only of a.th. in one-line
from
mature.fa  # a fasta file containing all mature miRNAs of all species
hairpin.fa # a fasta format file containing all miRNAs= hairpin of all species
downloaded from
miRbase.org
=end
