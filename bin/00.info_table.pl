#!/usr/bin/perl
use Getopt::Long;
use Pod::Usage;

GetOptions(
"mature|m=s" => \my $mature,
"hairpin|h=s" => \my $hairpin_oneline,
"matureGFF|g=s" => \my $mature_gff,
"outdir|o=s" => \my $out_dir,
"help" => \my $help
) or pod2usage(-verbose => 1);

#pre
pod2usage( -verbose =>1) if $help;
die "\e[01;31moutpout dir is required\e[00m\n" unless $out_dir;
$out_file1 = "$out_dir/all-hairpin-start-length.tbl";
$out_file2 = "$out_dir/all-mature-start-length.tbl";
open MUG, $mature_gff       or die "\e[01;31mError reading $mature_gff\n$!\e[00m\n";
open MUF, $mature           or die "\e[01;31mError reading $mature\n$!\e[00m\n";
open HPF, $hairpin_oneline  or die "\e[01;31mError reading $hairpin_oneline\n$!\e[00m\n";
open OUT, '>',$out_file1    or die "\e[01;31mError writing to $out_file1\n$!\e[00m\n";
open OUT2,'>',$out_file2    or die "\e[01;31mError writing to $out_file2\n$!\e[00m\n";

#process
while(<MUG>){
    next if /^#/;
    s/\r?\n|\r//;
    $start_pos =0;
    @flds = split /\s+/;
    next if $flds[2] ne "miRNA";
    %attrs = split /;|=/, $flds[8];
    $flag_fas = 1;
    $seq_fas  = '';
    seek(MUF,0,SEEK_SET);
    # after the loop
    while( $raw = <MUF>){
        $raw =~ s/\r?\n|\r//;
        if ($raw =~ /^>/){
            $raw = substr $raw, 1;
            @beads = split /[ ]/, $raw;
            if (lc($beads[0]) eq lc($attrs{Name})){
                $seq_fas = <MUF>;
                $seq_fas =~ s/\r?\n|\r//;
                next
           }
        }
    }
    # $seq_fas some empty
    # get start position
    seek(HPF,0,SEEK_SET);
    while($rawh = <HPF>){
        chomp($rawh);
        $name = $attrs{Name} =~ s/mi/MI/r;
        if ($rawh =~ /^>/ && $rawh =~ /$name/){
            chomp($ref = <HPF>);
            $start_pos = index $ref, $seq_fas; # 1-based

            next if $start_pos == -1;
            if ($start_pos != -1){
              if ($start_pos-2 >= 0){
                $mature_shift2 = substr $ref,$start_pos-2,length($seq_fas)+4;
              }

            }
        }
    }
    next if $start_pos == -1;
    $start_pos -=1; # covert to 0-based
    $name = $attrs{Name} =~ s/mi/MI/r;
    print OUT  "$name\t$seq_fas\t$start_pos\t".($flds[4]-$flds[3]+1)."\n";
    print OUT2 "$name\t$mature_shift2\t0\t".(length($mature_shift2))."\n";
}

=head1 SYNOPSIS
 will get
  all-hairpin-start-length.tbl
  all-mature-start-length.tbl
 from
  step-0 result files
 and
  gff3 file
   from miRbase should be given
=cut
