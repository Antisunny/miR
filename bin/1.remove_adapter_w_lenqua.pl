#!/usr/bin/perl

use feature ":5.16";
use Getopt::Long qw/GetOptions/;
use Pod::Usage qw/pod2usage/;
no warnings "all";

GetOptions(
    'adapter|a=s' => \my $adapter,
    'min-length|mil=i' => \my $min_length,
    'max-length|mal=i' => \my $max_length,
    'min-quality|maq=i' => \my $max_quality,
    'outdir|o=s'  => \my $outdir,
    'prefix|p=s'  => \my $prefix,
    'help|h'        => \my $help
 ) or pod2usage(-verbose => 2);

pod2usage( -verbose =>1) if $help or @ARGV == 0;
die "adapter is required\n" unless $adapter;
my $infile = $ARGV[0];
if ($ARGV[0] =~ /\.(gz|bz2)$/){
  die "compressed fastq is not supported. uncompress it."
}
unless($prefix){
	$prefix = [split /\//, $infile]->[-1];
	$prefix =~ m/([A-Za-z0-9_\.-]+)\.(fastq|fq)$/;
	$prefix = $1;
}

$outfile = $prefix."_armd.fa";
open IN,$infile or die "Error in reading $infile\n";
open OUT,">>","$outdir/$outfile" or die "$!\nError in writing results to $outdir/$outfile\n";
open IN,  '<',   $infile                     or die "Error $infile\n";
open OUT, '>>', "$outdir/$outfile"           or die "Error write results to $outfile\n";
open LOG, '>>', "$outdir/${prefix}_run.log"  or die "$!\n";
open DBT ,'>>', "$outdir/${prefix}_dist.txt" or die "$!\n";
my $cc = 0;
my $all_reads_num=0;
my $subseq_num=0;
my $clean_reads_num=0;
my %PhredScore=();
my @range;
push @range, $min_length ? $min_length : 17;
push @range, $max_length ? $max_length : 30;
my $min_quality = $min_quality ? $min_quality : 20;
my %length_distribution=();
my %length_distribution_A=();
my %length_distribution_G=();
my %length_distribution_C=();
my %length_distribution_T=();

while(my $raw1=<IN>) {
	chomp(my $raw=<IN>);
	chomp(my $raw3=<IN>);
	chomp(my $raw4=<IN>);
	chomp($raw1);
	if($raw=~/$adapter/g) {
		my $match_position=pos($raw);
		my $mirna='';
		my $quality='';
		my $score='';
		if($match_position>length($adapter)) {
			$mirna = substr $raw, 0, $match_position-length($adapter);
			$quality = substr $raw4, 0, $match_position-length($adapter);
			$score = &MeanPhredScore($quality);
			my $readslen = length($mirna);
			if( $readslen >= $range[0] && $readslen <= $range[1]){
				$subseq_num++;
				$PhredScore{$score}++;
				if($score >= $min_quality){
					$length_distribution{$readslen}++;
					next if $mirna =~ /N/i;
					my $startswith = substr $mirna,0,1;
					given ($startswith) {
						when (/[Aa]/) { $length_distribution_A{$readslen}++; }
						when (/[Gg]/) { $length_distribution_G{$readslen}++; }
						when (/[Cc]/) { $length_distribution_C{$readslen}++; }
						when (/[Tt]/) { $length_distribution_T{$readslen}++; }
						default {next}
					}
					print OUT ">${prefix}-$clean_reads_num-score:$score\n$mirna\n" if $mirna !~ /N/i;
					$clean_reads_num++;
				}
			}
		}
	}
	$all_reads_num++;
	if($all_reads_num%10000 == 0){
		$cc += 1;
		print LOG "*";
		if ($cc%80 == 0) {
			print LOG "\t$all_reads_num\n";
			print "$all_reads_num\n";
		}
	}
}
close(IN);
close(OUT);

# LOG summary of the reads
print LOG "\n".("*" x 40)."\n$name\ntotal: $all_reads_num\n$min_length-$max_length: $subseq_num\n and score>=$min_quality: $clean_reads_num\n";
print LOG "mean phred score distribution\n";
my @ss = keys %PhredScore;
for (sort {$a cmp $b} @ss) {
	print LOG "$_\t$PhredScore{$_}\n";
}
close(LOG);
my $calculated_max_length = &max(keys %length_distribution);
my $calculated_min_length = &min(keys %length_distribution);
$min_length = $range[0] > $calculated_min_length ? $range[0] : $calculated_min_length;
$max_length = $range[1] > $calculated_max_length ? $calculated_max_length : $range[1];

# DBT save all the count into a tablular
my $title = "$prefix\t";
for ($min_length..$max_length){
  $title .= $_."\t";
}
$title =~ s/\t$/\n/;
print DBT $title;
sub save_recoreds{
	my ($mark,%a) = @_;
	my $record = "";
	foreach ($min_length..$max_length){
		$record .= $a{$_}."\t";
	}
  $record =~ s/\t$/\n/;
	print DBT $record;
}
&save_recoreds('all',%length_distribution);
&save_recoreds('A',  %length_distribution_A);
&save_recoreds('T',  %length_distribution_T);
&save_recoreds('C',  %length_distribution_C);
&save_recoreds('G',  %length_distribution_G);
close DBT;

sub MeanPhredScore{
	my @inf = split//, shift;
	my $score=0;
	for(my $i=0;$i<scalar(@inf);$i++){
		$score += ord($inf[$i]);
	}
	int($score/length(@inf)+0.5);
}
sub max{
	my $max = shift;
	for (@_) {$max = $_ if ($_ > $max);}
	$max
}
sub min{
    my $min = shift;
    for (@_) {$min = $_ if ($_ < $min)}
    $min
}
__END__
=head1 SYNOPSIS

  perl 1.remove_adapter_w_lenqua.pl [options] fastq-file

  Options:
    --adapter or -a specifies the adapter sequence. If it is more than 15bg long,
                    only the first 15 bp is kept and used for indexing.

    --min-length    specifies the minimal length of the kept seq.[17]

    --max-length    specifies the maximal length of the kept seq.[30]

    --min-quality   specifies the mininal quality of the seq to be kept.[20]

    --prefix or -p  specifies the unique name for this inout file. This will be
                    used as the prefix for all output files, ex. $prefix."_armd.fa"
                    is the kept seq in fasta format, $prefix."_run.log" is the log
                    file, $prefix."_dist.fa" is the statistics of the reads length
                    after trimming.
                    if not specified, the name of the input file will be used w/o
                    extension.

    --outdir or -o  specifies the output dir, which will includes all the output
                    files of this step

    --help or -h    will print this help info page.

    fastq-file      specifies the input fastq file to be filtered. The only format
                    supported is fastq, not compressed (gziped or bzip2ed is not
                    supported).

=back
=head1 DESCRIPTION

  B<This program> will read the fastq-formatted file with options to get trimmed
  reads that contains adapter sequence and is larger than the $min_length and less than
  the $max_length in length and the mean quality of which is more than the $min_quality.

=cut
