#!/usr/bin/perl

use Getopt::Long;
use Pod::Usage qw/pod2usage/;

GetOptions(
	'sam|m=s' => \my $samFile,
	'table|t=s' => \my $tblFile,
	'outdir|o=s' => \my $outDir,
	'prefix|p=s' => \my $prefix,
	'help|h'  => \my $help
) or pod2usage(-verbose => 1);

pod2usage(-verbose => 1) if $help;
$cate  = $samFile =~ s|[^/]+/||gr =~ s/$prefix\.//r =~ s/\.std\.sam//r;
$out_file = "$outDir/$prefix-$cate.count";
die "\e[01;31moutpout dir is required\e[00m\n" unless $outDir;
die "\e[01;31mprefix is required\e[00m\n" unless $prefix;
open IN,  $tblFile       or die "\e[01;31mError reading $tblFile\n$!\e[00m\n";
open SAM, $samFile       or die "\e[01;31mError reading $samFile\n$!\e[00m\n";
open OUT, '>',$out_file  or die "\e[01;31mError writing to $out_file\n$!\e[00m\n";
print OUT "name\ttotal\ttrimmed\tintack\n";

# format SAM file to samller
@SAM =();
while (<SAM>) {
    next if /^@/;
    my ($a,$b,$c) = (split/\t/)[0,2,3];
    next if $b eq '*';
    push @SAM, [$a,$b,$c];
}
close SAM;

while (<IN>){
	s/\r?\n|\r//;
  &dist((split/\t/)[0,2,3],\@SAM);
}

sub dist{
    my ($miRNA,$start,$length,$lib) = @_;
    my $total=0;
    my $trimmedOnly=0;
    my $intack=0;
	print "$miRNA ...\n";
    for (@$lib){
        if ($$_[1] =~ m#$miRNA#i){
            if (abs($$_[2]-$start) <= 2) {
                @heads = split /-+/, $$_[0];
                if(length($heads[2]) < $length){
                    $trimmedOnly++;
                    print STDERR "*";
                }
                if(length($heads[2]) == $length){
                    $intack++;
                    print STDERR "#";
                }
            }
            $total++;
        }
    }
	print STDERR "\n";
    return -1 if ! $total;
    print OUT "$miRNA\t$total\t$trimmedOnly\t$intack\n";
}

__END__
=head1 SYNOPSIS

	perl 5.count.pl [options]

	Options:

	--sam or -s      specifies the 4.remapping SAM file of hairpin or mature(not both).

	--table or -t    specifies the 00.info_table tbl file hairpin or mature accordingly.

	--outdir or -o   specifies the output dir, which will includes all the output
                      files of this step

	--prefix or -p   specifies the unique name for this input file. This will be
								   used as the prefix for all output files

	--help or -h     will print this help info page.

=cut
