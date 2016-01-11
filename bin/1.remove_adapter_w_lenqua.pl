#!/usr/bin/perl

use feature ":5.16";
use Getopt::Long;
use Pod::Usage;
no warnings "all";

my ($adapter,$outdir,$prefix,$help) = ("","pre-check",'',0);
GetOptions(
    'adapter|a=s' => \$adapter,
    'outdir|o=s'  => \$outdir,
    'prefix|p=s'  => \$prefix,
    'help'        => \$help
 ) or pod2usage(1);

die "adapter is required\n" unless $adapter;
my $infile = $ARGV[0];
unless($prefix){
	$prefix = [split /\//, $infile]->[-1];
	$prefix =~ m/([A-Za-z0-9_\.-]+)\.(fastq|fq)$/;
	$prefix = $1;
}

$outfile = $prefix."_armd.fa";

open IN,$infile or die "Error in reading $infile\n";
open OUT,">>","$outdir/$outfile" or die "$!\nError in writing results to $outdir/$outfile\n";

my $cc = 0;
open IN,  '<',   $infile                     or die "Error $infile\n";
open OUT, '>>', "$outdir/$outfile"           or die "Error write results to $outfile\n";
open LOG, '>>', "$outdir/${prefix}_run.log"  or die "$!\n";
open DBT ,'>>', "$outdir/${prefix}_dist.txt" or die "$!\n";

my $all_reads_num=0; # all that reads
my $subseq_num=0; # 12-30
my $clean_reads_num=0; # 12-30 and >20
my %PhredScore=();
my @range = (12,30);
my $minQual = 20;
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
				if($score >= $minQual){
					$length_distribution{$readslen}++;
					# check the leading base distribution
					next if $mirna =~ /N/i;
					my $startswith = substr $mirna,0,1;
					given ($startswith) {
						when (/[Aa]/) { $length_distribution_A{$readslen}++; }
						when (/[Gg]/) { $length_distribution_G{$readslen}++; }
						when (/[Cc]/) { $length_distribution_C{$readslen}++; }
						when (/[Tt]/) { $length_distribution_T{$readslen}++; }
						default { print STDERR "Impossible\n" }
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
print LOG "\n".("*" x 40)."\n$name\ntotal: $all_reads_num\n12-30: $subseq_num\n and score>=20: $clean_reads_num\n";
print LOG "mean phred score distribution\n";
my @ss = keys %PhredScore;
for (sort {$a cmp $b} @ss) {
	print LOG "$_\t$PhredScore{$_}\n";
}
close(LOG);

my $min_length = $range[0];
my $max_length = $range[1];

# DBT save all the count into a tablular 
my $title = "$prefix\t";
for ($min_length..$max_length){
	if ($_ == $max_length) {
        $title .= $_."\n";
    }else{
        $title .= $_."\t";
    }
}
print DBT $title;
sub save_recoreds{
	my ($mark,%a) = @_;
	my $record = "";
	foreach ($min_length..$max_length){
		if ($_ != $max_length){$record .= $a{$_}."\t";}
		else{$record .= $a{$_}."\n";}
	}
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
