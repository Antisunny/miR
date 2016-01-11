#!/bin/bash
miR="$HOME/Bin/miR"
# prepare
cd $miR/bin
perl 0.prepare.pl -m ../mature.fa -h ../hairpin.fa -o ..
perl 00.info_table.pl -m ../athMature.fa -h ../athHairpin.txt -g ../ath.gff3 -o ..

domiR(){
	mark=$1
	mkdir -p $mark/1
	perl $miR/bin/1.remove_adapter_w_lenqua.pl -a AAGATCGGAAGAGCA -p $mark -o ${mark}/1 ${mark}.fq
	echo "["`date +%m/%d/%Y-%H:%M:%S`"] Done with step 1; Saved to $mark/1"
	bowtie2 -f -p 20 -x /mnt/data2/BT2/tair10_bt2 -U $mark/1/${mark}_armd.fa -S $mark/1/${mark}.sam >& $mark/1/bowtie2.log
	echo "["`date +%m/%d/%Y-%H:%M:%S`"] Done with mapping; Saved to $mark/1/${mark}.sam"
	samtools sort -T ${mark}_temp -O 'sam' -o $mark/1/${mark}.std.sam -@ 10 $mark/1/${mark}.sam
	mkdir -p $mark/2
	perl $miR/bin/2.discard-reads-of-trsno.pl -d $miR/tair-t-r-sn-sno.txt -p ${mark} -o ${mark}/2 ${mark}/1/${mark}.std.sam
	echo "["`date +%m/%d/%Y-%H:%M:%S`"] Done with step 2; Saved to $mark/2"
	mkdir -p $mark/3
	perl $miR/bin/3.extract_seq_no_trsno.pl -s $mark/1/${mark}_armd.fa -o $mark/3 -p ${mark} $mark/2/${mark}_removed.dat
	echo -e "["`date +%m/%d/%Y-%H:%M:%S`"] Done with step 3\nStart remapping"
	mkdir -p $mark/4
	perl $miR/bin/4.remapping.pl -o $mark/4 -p $mark -r $miR/athHairpin $mark/3/${mark}_wo_trsno.fa >& $mark/4/bowtie2.log
	echo "["`date +%m/%d/%Y-%H:%M:%S`"] Done with step 4 remapping to hairpin and mature"
	echo "["`date +%m/%d/%Y-%H:%M:%S`"] Starting hairpin count";
	mkdir -p $mark/5
	perl $miR/bin/5.count.pl -o $mark/5 -m ${mark}/4/${mark}.hairpin.std.sam -t $miR/all-hairpin-start-length.tbl -p $mark
	echo "["`date +%m/%d/%Y-%H:%M:%S`"] Start mature count"
	perl $miR/bin/5.count.pl -o $mark/5 -m ${mark}/4/${mark}.mature.std.sam -t $miR/all-mature-start-length.tbl3 -p $mark
	echo "["`date +%m/%d/%Y-%H:%M:%S`"] Done with step 5 counting";
}
domiRs(){
	for pf in $@;do
		echo "Processing $pf"
		domiR $pf
	done
}
domiRs ilp1_1 ilp1_2 Col-0_1 Col-0_2
