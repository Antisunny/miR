#!/usr/bin/perl

use Getopt::Long;
use Pod::Usage;

GetOptions(
	'sam|m=s' => \my $samFile,
	'table|t=s' => \my $tblFile,
	'outdir|o=s' => \my $outDir,
	'prefix|p=s' => \my $prefix,
	'htlp|h'  => \my $help
) or pod2usage(-verbose => 1);

pod2usage(-verbose => 1) if $help;
$cate  = $samFile =~ s#[^/]+/##gr =~ s/$prefix\.//r =~ s/\.std\.sam//r;
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
    #&miRNA_coverage((split /\t/)[0,2]);
    &hairpin_dist((split/\t/)[0,2,3],\@SAM);
}

sub miRNA_coverage{
	my($miRNA,$start) = @_;
	my %total = ();
	print STDERR "$miRNA\n";
    # $miRNA like `ath-miR447a.2-3p`
	print "$miRNA\n";
	my $hairpin = join '-',(split(/-|\./,$miRNA))[0,1] =~s/mi/MI/r;
	seek(SAM,0,SEEK_SET);
	while(<SAM>){
        next if /^@/;
		if(/$hairpin/){
			my @flds = split /\t/;
			my $seq =(split /-+/, $flds[0])[0];
			if(abs($flds[3]-$start) <= 2 ){
                # allow 2bp shift 
				for(my $i = $flds[3]; $i <= $flds[3]+length($seq)-1; $i++){ ##total coverage
					if ($i >= $start){
                        $total{$i}++;
                    }
				}
			}
		}
	}			
	close SAM;
    
    if (! %total){
        print "$miRNA empty %total\n";
        return 1;
    }
    print OUT "position\t";
	print OUT "$_\t"  for (sort {$a<=>$b} (keys %total));
	print OUT "\n$miRNA\t";
	print OUT "$total{$_}\t" for (sort {$a<=>$b} (keys %total));
	print OUT "\n";
}
sub hairpin_dist{
    my ($miRNA,$start,$length,$lib) = @_;
    my $total=0;
    my $trimmedOnly=0;
    my $intack=0;
	print "$miRNA ...\n";
    #format ath-MIR447a
	# seek(SAM,0,SEEK_SET);
    # while(<SAM>){
    #     s/\r?\n|\r//;
    #     next if /^@/;
    #     if (/$miRNA/i){
    #         @flds = split/\t/;
    #         if(abs($flds[3]-$start) <= 2){
    #             @heads = split /-+/, $flds[0];
    #             # format Sample-SRR505135-TGACAGAAGAGAG--0 
    #             if(length($heads[2]) < $length){
    #                 $trimmedOnly++;
    #                 print STDERR "*";
    #             }
    #             if(length($heads[2]) == $length){
    #                 $intack++;
    #                 print STDERR "#";
    #             }
    #         }
    #         $total++;
    #     }
    # }
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

