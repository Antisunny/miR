#!/usr/bin/perl

use Getopt::Long;
use Pod::Usage qw/pod2usage/;

my ($seqFile,$outdir,$infile,$outfile,$prefix) = ('','','','','');
GetOptions(
  'seq|s=s'     => \$seqFile,
  'outdir|o=s'  => \$outdir,
  'prefix|p=s'  => \$prefix,
  'help'        => \$help
) or pod2usage(-verbose => 1);
pod2usage(-verbose => 1) if $help or @ARGV == 0;

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
    print SEQOUT ">$_\n$raw_seq" unless exists($names{$_});
}
close SEQIN;
close SEQOUT;

__END__
=head1 SYNOPSIS

  perl 1.extract_seq_no_trsno.pl [options] fasta_file

  Options:

  --seq or -s         specifies the fasta file gotten from
                      1.remove_adapter_w_lenqua.pl

  --prefix or -p      specifies the unique name for this inout file. This will be
                      used as the prefix for all output files

  --outdir or -o      specifies the output dir, which will includes all the output
                      files of this step

  --help or -h        will print this help info page.

  fasta_file          the removed file generated from 2.discard-reads-of-trsno.#!/usr/bin/env perl

=cut
