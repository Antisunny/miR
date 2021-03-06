#!/usr/bin/perl

use feature ":5.16";
use Getopt::Long;
use Pod::Usage qw/pod2usage/;
no warnings "all";

GetOptions(
    "discard_file|d=s" =>\my $discard_file,
    "prefix|p=s" => \my $prefix,
    "outdir|o=s" => \my $outDir,
    "help|h"     => \my $help
) or pod2usage(-verbose=>1);

pod2usage(-verbose=>1) if $help or @ARGV == 0;
$sam_file  = $ARGV[0];
die "discard file is required\n" unless $discard_file;
die "prefix|outdir is required\n" unless $prefix && $outDir;
die "sam File is required\n" unless $sam_file;
$removed = "$outDir/${prefix}_removed.dat";
$kept    = "$outDir/${prefix}_kept.dat";

open DIF, $discard_file or die "\033[01;31mError in reading $discard_file\n$!\033[00m\n";
@chr1 = ();
@chr2 = ();
@chr3 = ();
@chr4 = ();
@chr5 = ();
@chrc = ();
@chrm = ();
while(<DIF>){
    next if /^@/;
    @field = split/\t/;
    given($field[0]){
        when('Chr1'){
            push @chr1, [$field[3]-30,$field[4]+30];
        }
        when('Chr2'){
            push @chr2, [$field[3]-30,$field[4]+30];
        }
        when('Chr3'){
            push @chr3, [$field[3]-30,$field[4]+30];
        }
        when('Chr4'){
            push @chr4, [$field[3]-30,$field[4]+30];
        }
        when('Chr5'){
            push @chr5, [$field[3]-30,$field[4]+30];
        }
        when('ChrC'){
            push @chrc, [$field[3]-30,$field[4]+30];
        }
        when('ChrM'){
            push @chrm, [$field[3]-30,$field[4]+30];
        }
    }
}

open SAM, $sam_file      or "\033[01;31mError in reading $sam_file\n$!\033[00m\n";
open KEPT, ">", $kept    or "\033[01;31mError in writing to $kept\n$!\033[00m\n";
open REMV, ">", $removed or "\033[01;31mError in writing to $removed\n$!\033[00m\n";
print REMV "Chromesome\tindex\tstart_pos\tlength\trange_start\trange_end\n";
while(<SAM>){
	next if /^@/;
    @field = split/\t/;
    $removed=0;
    given($field[2]){
        when('Chr1'){
            for $enpair (@chr1){
                if ($field[3] >= $$enpair[0] && $field[3]+length($field[9]) <= $$enpair[1]){
                    $removed +=1;
                    print REMV "Chr1\t$field[0]\t$field[3]\t".length($field[9])."\t$$enpair[0]\t$$enpair[1]\n";
					last
                }
            }
        }
        when('Chr2'){
            for $enpair (@chr2){
                if ($field[3] >= $$enpair[0] && $field[3]+length($field[9]) <= $$enpair[1]){
                    $removed +=1;
                    print REMV "Chr2\t$field[0]\t$field[3]\t".length($field[9])."\t$$enpair[0]\t$$enpair[1]\n";
					last
                }
            }
        }
        when('Chr3'){
            for $enpair (@chr3){
                if ($field[3] >= $$enpair[0] && $field[3]+length($field[9]) <= $$enpair[1]){
                    $removed +=1;
                    print REMV "Chr3\t$field[0]\t$field[3]\t".length($field[9])."\t$$enpair[0]\t$$enpair[1]\n";
					last
                }
            }
        }
        when('Chr4'){
            for $enpair (@chr4){
                if ($field[3] >= $$enpair[0] && $field[3]+length($field[9]) <= $$enpair[1]){
                    $removed +=1;
                    print REMV "Chr4\t$field[0]\t$field[3]\t".length($field[9])."\t$$enpair[0]\t$$enpair[1]\n";
					last
                }
            }
        }
        when('Chr5'){
            for $enpair (@chr5){
                if ($field[3] >= $$enpair[0] && $field[3]+length($field[9]) <= $$enpair[1]){
                    $removed +=1;
                    print REMV "Chr5\t$field[0]\t$field[3]\t".length($field[9])."\t$$enpair[0]\t$$enpair[1]\n";
					last
                }
            }
        }
        when('ChrC'){
            for $enpair (@chrc){
                if ($field[3] >= $$enpair[0] && $field[3]+length($field[9]) <= $$enpair[1]){
                    $removed +=1;
                    print REMV "ChrC\t$field[0]\t$field[3]\t".length($field[9])."\t$$enpair[0]\t$$enpair[1]\n";
					last
                }
            }
        }
        when('ChrM'){
            for $enpair (@chrm){
                if ($field[3] >= $$enpair[0] && $field[3]+length($field[9]) <= $$enpair[1]){
                    $removed +=1;
                    print REMV "ChrM\t$field[0]\t$field[3]\t".length($field[9])."\t$$enpair[0]\t$$enpair[1]\n";
					last
                }
            }
        }
    }
    unless($removed > 0){
        print KEPT $_;
	}
}

=head1 SYNOPSIS

  perl 2.discard-reads-of-trsno.pl [options] sam_file

  Options:

  --discard_file or -d  specifies the file in GFF3 format containing all the
                        sites to be discarded.
                        By default, the file has been give, named tair-t-r-sn-sno.txt
                        ,which is extraced from TAIR10 Annotation file for A.th.

  --prefix or -p        specifies the unique name for this inout file. This will be
                        used as the prefix for all output files, ex.
                        ${prefix}_removed.dat contains info of reads matching sites
                        specified in $discard_file, ${prefix}_kept.dat contains the
                        reads kept for later use.

  --outdir or -o        specifies the output dir, which will includes all the output
                        files of this step

  --help or -h          will print this help info page.

  sam_file              SAM file generated from mapping to TAIR1O genome using bowtie2
                        . First index the Genome using bowtie2-build, and map
                        before this step.
=cut
