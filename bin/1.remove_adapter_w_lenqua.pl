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
  # get basename
	$prefix = [split /\//, $infile]->[-1];
  # remove extension and the left as $prefix
	$prefix =~ m/([A-Za-z0-9_\.-]+)\.(fastq|fq)$/;
	$prefix = $1;
}

# add file handler
$outfile = $prefix."_armd.fa";
open IN,$infile or die "Error in reading $infile\n";
open OUT,">>","$outdir/$outfile" or die "$!\nError in writing results to $outdir/$outfile\n";
open IN,  '<',   $infile                     or die "Error $infile\n";
open OUT, '>>', "$outdir/$outfile"           or die "Error write results to $outfile\n";
open LOG, '>>', "$outdir/${prefix}_run.log"  or die "$!\n";
open DBT ,'>>', "$outdir/${prefix}_dist.txt" or die "$!\n";

# global variables
my $cc = 0;
my $all_reads_num=0;    # all reads of the fastq file
my $subseq_num=0;       # reads between $min_length and $max_length
my $clean_reads_num=0;  # reads in length range and quality > 20
my %dist_qlt=();      # key as the mean quality of a read, value as the count of this quality
my @range;              # the reads longer than $min_length or shorter than $max_length kept
push @range, $min_length // 17;
push @range, $max_length // 30;
my $adapter_length = length($adapter);
my %kept_dist=();
my %dist_A=();
my %dist_G=();
my %dist_C=();
my %dist_T=();
my $too_long=0;
my $too_short=0;
my $no_adapter=0;

while(<IN>) {
  # remove the new line character
	my $sqs=<IN> =~ s/\r?\n|\r//r;
	my $pls=<IN>;
	my $qlt=<IN> =~ s/\r?\n|\r//r;
	if($sqs=~/$adapter/g) {
		my $match_position=$+[0];
		my $subseq='';
		my $quality='';
		my $score='';
    # make sure the seq left is longer than 0bp
		if($match_position > $adapter_length) {
			$subseq = substr $sqs, 0, $match_position-$adapter_length;
			my $readslen = length($subseq);
			if( $readslen >= $range[0] && $readslen <= $range[1]){
        # get mean quality of the reads
        $score = &mean_qlt(substr $qlt, 0, $match_position - $adapter_length);
				# $subseq_num++;
				$dist_qlt{$score}++;
				if($score >= $min_quality){
					$kept_dist{$readslen}++;
					next if $subseq =~ /N/i;
					my $startswith = substr $subseq,0,1;
					given ($startswith) {
						when (/[Aa]/) { $dist_A{$readslen}++; }
						when (/[Gg]/) { $dist_G{$readslen}++; }
						when (/[Cc]/) { $dist_C{$readslen}++; }
						when (/[Tt]/) { $dist_T{$readslen}++; }
						default {next}
					}
					print OUT ">${prefix}-$clean_reads_num-score:$score\n$subseq\n";
					$clean_reads_num++;
				}
			}elsif($readslen < $range[0]){
        $too_long++
      }elsif($readslen > $range[1]){
        $too_short++
      }
		}
	}else{
    $no_adapter++;
    say STDERR "no adapter in the read";
  }
	$all_reads_num++;
}
close(IN);
close(OUT);

# LOG summary of the reads
#$too_long = $all_reads_num - $subseq_num - $too_short;
$clean_ratio     = sprintf("%.2f%%", $clean_reads_num/$all_reads_num);
$too_long_ratio  = sprintf("%.2f%%", $too_long/$all_reads_num);
$too_short_ratio = sprintf("%.2f%%", $too_short/$all_reads_num);
$no_adapter_ratio= sprintf("%.2f%%", $no_adapter_ratio/$all_reads_num);
print LOG "sample name\t$name\ntotal reads\t$all_reads_num\nkept reads([$min_length,$max_length] and > $min_quality)\t$clean_reads_num\ntoo long\t$too_long ($too_long_ratio)\ntoo short\t$too_short($too_short_ratio)\n";
print LOG "mean phred score distribution\n";
print LOG &tablulate(\%dist_qlt);
my $clt_max_length = &max(keys %kept_dist);
my $clt_min_length = &min(keys %kept_dist);
$min_length = $range[0] > $clt_min_length ? $range[0] : $clt_min_length;
$max_length = $range[1] > $clt_max_length ? $clt_max_length : $range[1];

# DBT save all the count into a tablular
print DBT &tablulate(\%kept_dist,1);
print DBT &tablulate(\%dist_A,0);
print DBT &tablulate(\%dist_T,0);
print DBT &tablulate(\%dist_C,0);
print DBT &tablulate(\%dist_G,0);
# my $title = "$prefix\t";
# $title .= "$_\t" for ($min_length..$max_length);
# $title =~ s/\t$/\n/;
# print DBT $title;
#
# sub save_recoreds{
# 	my ($mark,%a) = @_;
# 	my $record = "";
# 	foreach ($min_length..$max_length){
# 		$record .= $a{$_}."\t";
# 	}
#   $record =~ s/\t$/\n/;
# 	print DBT $record;
# }
# #&save_recoreds('all',%kept_dist);
# &save_recoreds('A',  %dist_A);
# &save_recoreds('T',  %dist_T);
# &save_recoreds('C',  %dist_C);
# &save_recoreds('G',  %dist_G);
close DBT;

sub mean_qlt{
	my @inf = split//, shift;
	my $score=0;
	for(my $i=0;$i<scalar(@inf);$i++){
		$score += ord($inf[$i])-33;
	}
	int($score/scalar(@inf)+0.5);
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
sub tablulate(){
  # hash is_header_print
  $in = shift;
  $is_header = shift // 1;
  die "a hash type required" if HASH == ref %$in;
  @head_fields = sort {$a cmp $b} keys %$in;
  $header = join "\t", @head_fields;
  $value .= "$_\t" for (@head_fields);
  $value =~ s/\t$//;
  return $is_header ? "$header\n$value\n" : "$value\n";
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
