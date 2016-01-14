#!/usr/bin/perl
use Pod::Usage qw/pod2usage/;
use Getopt::Long;
my ($seqFile,$outDir,$inFile,$outFile,$prefix) = ('','','','','');
GetOptions(
  'outdir|o=s' => \$outDir,
  'prefix|p=s' => \$prefix,
  'ref|r=s'    => \$refHairpin,
  'help|h'       => \$help
) or pod2usage(-verbose=>1);
pod2usage(-verbose => 1) if $help or @ARGV == 0;

$inFile = $ARGV[0];
die "prefix is required\n"     unless $prefix;
die "output dir is required\n" unless $outDir;

use File::Path;
rmtree($outDir) if -d $outDir;
mkdir $outDir;

$fmtFile = "$outDir/${prefix}_no_trsno_n.fa";
open IN,     $inFile  or die "Error $inFile : $!\n";
open FMT,">",$fmtFile or die "Error $fmtFile : $!\n";
while(<IN>){
    $seq = <IN>;
    $inf=[split/-/]->[0];
    $seq = $seq =~ s/T/U/gr =~ s/\r?\n|\r//r;
    print FMT ">sample-$inf-$seq--0\n$seq\n";
}
close(IN);
close(FMT);

print STDERR "Starting map to hairpin\n";
$mappedOut  = "$outDir/${prefix}_mapped_to_hairpin.sam";
$cmd_bowtie = "bowtie2 -f -p 20 -x $refHairpin -U $fmtFile -S $mappedOut";
system $cmd_bowtie;
$cmd_sort   = "samtools sort -T ${prefix}_temp -O 'sam' -o $outDir/${prefix}.hairpin.std.sam -@ 10 $mappedOut";
system $cmd_sort;

print STDERR "Starting map to mature\n";
$refMature = $refHairpin =~ s/Hairpin/Mature/r;
$mappedOut  = "$outDir/${prefix}_mapped_to_mature.sam";
$cmd_bowtie = "bowtie2 -f -p 20 -x $refMature -U $fmtFile -S $mappedOut";
system $cmd_bowtie;
$cmd_sort   = "samtools sort -T ${prefix}_temp -O 'sam' -o $outDir/${prefix}.mature.std.sam -@ 10 $mappedOut";
system $cmd_sort;

__END__
=head1 SYNOPSIS

  perl 4.remapping.pl [options]

  Options:

  --outdir or -o      specifies the output dir, which will includes all the output
                      files of this step

  --prefix or -p      specifies the unique name for this inout file. This will be
                      used as the prefix for all output files

  --ref or -r         the reference sequences of miRNA hairpin and mature from 0.prepare.pl
                      both should be indexed using bowtie2 first.

  --help or -h        will print this help info page.

=cut
