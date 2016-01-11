#!/usr/bin/perl -w

#########################################
# map reads of potential miRNA to hairpin
# and trim 1bp per time 10 times for unmapped
# and merge all the mapped together
# and sort and igv index
use Getopt::Long;
my ($seqFile,$outDir,$inFile,$outFile,$prefix) = ('','','','','');
GetOptions(
  'outdir|o=s' => \$outDir,
  'prefix|p=s' => \$prefix,
  'ref|r=s'      => \$refHairpin,

) or die "help\n";

$inFile = $ARGV[0];
die "prefix is required\n"     unless $prefix;
die "output dir is required\n" unless $outDir;
# make sure $outDir is empty
use File::Path;
rmtree($outDir) if -d $outDir;
mkdir $outDir;

$fmtFile = "$outDir/${prefix}_no_trsno_n.fa";
open IN,     $inFile  or die "Error $inFile : $!\n";
open FMT,">",$fmtFile or die "Error $fmtFile : $!\n";
while(<IN>){
    $seq = <IN>;
    $inf=[split/-/]->[0];
    $seq =~ s/T/U/g;
    print FMT ">sample-$inf-$seq--0\n$seq";
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
