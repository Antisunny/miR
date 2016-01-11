#!/usr/bin/py3
####################################################
# >tair-t-r-sn-sno.txt contains all the 4 interfering
# RNA that are supposed to be removed from the mapped
# reader sooner or later
# >samname is the reads that has mapped to the genome TAIR10
# via bowtie2 and bowtie2-build
# 1. all reads of 4 RNA(t,r,sn,sno) saved to ${name}removed-reads.txt
# in the form of "chr chr name_of_reads `reads-start` start_pos `length` length_of_seq `in` start_pos end_pos"
####################################################
import sys,os
from pathlib import Path
chr1 = list()
chr2 = list()
chr3 = list()
chr4 = list()
chr5 = list()
chrc = list()
chrm = list()

from optparse import OptionParser
parser =OptionParser(usage="%prog [options] SAM_FILE",version="%prog version 1.1", description='discard reads in the range specified by -d option')
parser.add_option('-d','--discard',metavar='FILE',dest='discard_file',help="specifies the ranges to be discarded in GFF3 format")
parser.add_option('-o','--outdir',metavar="PATH",dest='outdir',help="output directory")
parser.add_option('-p','--prefix',metavar="Str",dest='prefix',help='prefix for all output files')
(option,args) = parser.parse_args()
samfile = args[0]
outdir = option.outdir
if not option.prefix:
    if samfile.find('/'):
        prefix = samfile.split('/')[-1].split('.')[0]
    else:
        prefix = samfile.split('.')[0]
else:
    prefix = option.prefix
removed_reads_output_path = '{}/{}_removed_reads.txt'.format(outdir,prefix)
kept_reads_output_path = '{}/{}_kept.fq'.format(outdir,prefix)
if not Path(outdir).is_dir():
    os.mkdir(outdir)

for record in open(option.discard_file):
    if record.startswith('@'):
        continue
    fields = record.split("\t")

    if   fields[0] == 'Chr1':
        se = (int(fields[3])-30,int(fields[4])+30)
        chr1.append(se)
    elif fields[0] == 'Chr2':
        se = (int(fields[3])-30,int(fields[4])+30)
        chr2.append(se)
    elif fields[0] == 'Chr3':
        se = (int(fields[3])-30,int(fields[4])+30)
        chr3.append(se)
    elif fields[0] == 'Chr4':
        se = (int(fields[3])-30,int(fields[4])+30)
        chr4.append(se)
    elif fields[0] == 'Chr5':
        se = (int(fields[3])-30,int(fields[4])+30)
        chr5.append(se)
    elif fields[0] == 'ChrC':
        se = (int(fields[3])-30,int(fields[4])+30)
        chrc.append(se)
    elif fields[0] == 'ChrM':
        se = (int(fields[3])-30,int(fields[4])+30)
        chrm.append(se)

for line in open(samfile):
    field = line.split("\t")
    removed=0
    if field[2] == 'Chr1':
        for enpair in chr1:
            if int(field[3]) >= int(enpair[0]) and int(field[3])+len(field[9]) <= int(enpair[1]):
                removed+=1
                print('Chr1',field[0],field[3],str(len(field[9])),enpair[0],enpair[1],sep='\t',file=open(removed_reads_output_path,'a'))
            else:
                print(line,file=open(kept_reads_output_path,'a'))
        continue
    if field[2] == 'Chr2':
        for enpair in chr2:
            if int(field[3]) >= int(enpair[0]) and int(field[3])+len(field[9]) <= int(enpair[1]):
                removed+=1
                print('Chr2',field[0],field[3],str(len(field[9])),enpair[0],enpair[1],sep='\t',file=open(removed_reads_output_path,'a'))
            else:
                print(line,file=open(kept_reads_output_path,'a'))
        continue
    if field[2] == 'Chr3':
        for enpair in chr1:
            if int(field[3]) >= int(enpair[0]) and int(field[3])+len(field[9]) <= int(enpair[1]):
                removed+=1
                print('Chr3',field[0],field[3],str(len(field[9])),enpair[0],enpair[1],sep='\t',file=open(removed_reads_output_path,'a'))
            else:
                print(line,file=open(kept_reads_output_path,'a'))
        continue
    if field[2] == 'Chr4':
        for enpair in chr1:
            if int(field[3]) >= int(enpair[0]) and int(field[3])+len(field[9]) <= int(enpair[1]):
                removed+=1
                print('Chr4',field[0],field[3],str(len(field[9])),enpair[0],enpair[1],sep='\t',file=open(removed_reads_output_path,'a'))
            else:
                print(line,file=open(kept_reads_output_path,'a'))
        continue
    if field[2] == 'Chr5':
        for enpair in chr1:
            if int(field[3]) >= int(enpair[0]) and int(field[3])+len(field[9]) <= int(enpair[1]):
                removed+=1
                print('Chr5',field[0],field[3],str(len(field[9])),enpair[0],enpair[1],sep='\t',file=open(removed_reads_output_path,'a'))
            else:
                print(line,file=open(kept_reads_output_path,'a'))
        continue
    if field[2] == 'ChrC':
        for enpair in chr1:
            if int(field[3]) >= int(enpair[0]) and int(field[3])+len(field[9]) <= int(enpair[1]):
                removed+=1
                print('ChrC',field[0],field[3],str(len(field[9])),enpair[0],enpair[1],sep='\t',file=open(removed_reads_output_path,'a'))
            else:
                print(line,file=open(kept_reads_output_path,'a'))
        continue
    if field[2] == 'ChrM':
        for enpair in chr1:
            if int(field[3]) >= int(enpair[0]) and int(field[3])+len(field[9]) <= int(enpair[1]):
                removed+=1
                print('ChrM',field[0],field[3],str(len(field[9])),enpair[0],enpair[1],sep='\t',file=open(removed_reads_output_path,'a'))
            else:
                print(line,file=open(kept_reads_output_path,'a'))
        continue
    if removed > 0:
        print(field[2],sep='\t',file=open(removed_reads_output_path,'a'))

# LOG
# those to be removed are allowed a 30bp shift either upstream or upstream
