#!/usr/bin/perl

use Getopt::Long;
my ($seqFile,$outdir,$infile,$outfile,$prefix) = ('','','','','');
GetOptions(
  'seq|s=s'     => \$seqFile,
  'outdir|o=s'  => \$outdir,
  'prefix|p=s'  => \$prefix,
  'help'        => \$help
) or die "help\n";

$infile= $ARGV[0];
# outdir
die 'output directory is required\n' unless $outdir;
use File::Path;
rmtree($outdir) if -d $outdir;
mkdir $outdir;
# prefix
unless($prefix){
    $prefix = [split/[\.-]/,[split/\//,$infile]->[-1]]->[0];
}
$outfile = "$outdir/${prefix}_removed_ids.txt";

open IN, $infile or die "Error $infile:$!\n";
open OUT, ">", $outfile or die "Error $outfile:$!\n";

%names=();
while(<IN>){
    $names{[split /\t/]->[1]}++
}

print OUT "all".(scalar keys %names)."\n";
print OUT "$_\n" foreach (keys %names);
close IN;
close OUT;

# extract fasta-formated data of non-trsnoRNA from `step1 result data`
$seqOut = "$outdir/${prefix}_wo_trsno.fa";
open SEQIN,       $seqFile or die "Error $seqFile: $!\n";
open SEQOUT, ">", $seqOut  or die "Error $seqOut: $!\n";

while(<SEQIN>){
    s/\r?\n|\r//;s/^>//;
    $raw_seq = <SEQIN>;
    print SEQOUT "$_\n$raw_seq" unless exists($names{$_});
}
close SEQIN;
close SEQOUT;
