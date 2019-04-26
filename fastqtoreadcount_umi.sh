#! /bin/bash

# Nagendra Palani - University of Minnesota Genomics Center - 04/2018 - MN, USA
# Revision - 04/2019

#Pipeline for TnSeq data preprocessing, includes UMI deduplication.

#This script is designed for use on the University of Minnesota HPC facility. You would (most probably) need to customize
# for your environment. All tools are available through Conda install. When running at UMN MSI, Run with MINIMUM RAM of 32G.

#Tool Dependencies - bbtools, Hisat2, GNU Parallel, Bioawk.
#Tools used are multi-threaded or use GNU Parallel for multi-core execution.

#File dependencies - ReferenceOrganism.fasta, TAgenemap.txt (this file needs to be alphabetically sorted).

# NOTE: MAKE SURE THAT R1 AND R2 FILES ARE GUNZIPPED. THIS WORKFLOW WILL WORK ONLY WITH UNCOMPRESSED FASTQ FILES.
# DO NOT USE FASTQ.GZ. LEADS TO ISSUES WITH SOFTWARE USED FOR UMI ANALYSIS.
# run this in the directory containing the gzipped files before launching this workflow.
# > parallel gunzip {} ::: *.fastq.gz


#base folder containing sequencing data
basedr=~/Downloads/SeqFolder



#sub-directories to be created for outputs

mkdir $basedr/00_umi_filter
mkdir $basedr/0_umi_trim
mkdir $basedr/1_umi_clump
mkdir $basedr/2_umi_header
mkdir $basedr/3_umi_filter_r1

mkdir $basedr/RefIndex
mkdir $basedr/0_filter
mkdir $basedr/1_cutadapt
mkdir $basedr/2_sizecut
mkdir $basedr/3_TAonly
mkdir $basedr/4_alignfiles
mkdir $basedr/5_sam016
mkdir $basedr/6_idxbamfiles
mkdir $basedr/7_readfreqs
mkdir $basedr/8_mappedinserts

#Place fasta & TAposition files from matlab output in reffiles folder within base directory
reffasta=ReferenceOrganism.fasta
alignidxname="$(basename "${reffasta%.fasta*}")_idx"

TAgenemap_sort=TAgenemap.txt


module load hisat2
#--------------------------------------------------------------------------

#Build Hisat2 aligner index

hisat2-build -f $basedr/reffiles/$reffasta $basedr/RefIndex/$alignidxname

#--------------------------------------------------------------------------
# hisat2-build has conflicts running within tnseqenv


module load python2
source activate tnseq_env
module load parallel
module load samtools
module load hisat2


#hisat2 threads - set threadnum equal to number of processor cores available.
threadnum=12

#samtools threads (1 less than available cores)
samthread=11

# Run with MINIMUM RAM of 32G
#--------------------------------------------------------------------------

#UMI based deduplication

for r2file in $basedr/*_R2_001.fastq;
do
bbduk.sh -Xmx30g in="$r2file" outm=$basedr/00_umi_filter/"$(basename "${r2file%_R2*}" )_filtered.fastq" literal=CTACAAGAGCGGTGAG k=16 hdist=3 rcomp=f maskmiddle=f
done


parallel "reformat.sh -Xmx3g in={} out=$basedr/0_umi_trim/{/.}_trim.fastq forcetrimright=50" ::: $basedr/00_umi_filter/*_filtered.fastq

for trimfile in $basedr/0_umi_trim/*_trim.fastq;
do
clumpify.sh -Xmx30g in="$trimfile" out=$basedr/1_umi_clump/"$(basename "${trimfile%_trim*}" )_clump.fastq" dedupe=t passes=3 subs=3
done

bawk_cond='{print $name;}'
parallel "bioawk -c fastx '$bawk_cond' {} > $basedr/2_umi_header/{/.}_header.txt" ::: $basedr/1_umi_clump/*_clump.fastq

for headerfile in $basedr/*_R1_001.fastq;
do
filterbyname.sh -Xmx30g in="$headerfile" out=$basedr/3_umi_filter_r1/"$(basename "${headerfile%_R1_001*}" )_umir1.fastq" names=$basedr/2_umi_header/"$(basename "${headerfile%_R1_001*}" )_filtered_clump_header.txt" include=t
done

#--------------------------------------------------------------------------

#filter sequences for valid transposons - bbduk (bbtools package) - specify filtering seq in literal

#sleeping beauty literal = TAAACTTCCGACTTCAACTG
#mariner transposon literal = GGTCTTATCATCCAACCTGT

for seqfile in $basedr/3_umi_filter_r1/*_umir1.fastq;
do
bbduk.sh -Xmx30g in="$seqfile" outm=$basedr/0_filter/"$(basename "${seqfile%_umir1*}" )_filtered.fastq" literal=GGTCTTATCATCCAACCTGT k=20 hdist=3 rcomp=f maskmiddle=f

done
#--------------------------------------------------------------------------

#trim 5' adapters - bbduk (bbtools package) - specify adapter seq in literal

for seqfilefilt in $basedr/0_filter/*.fastq;
do
bbduk.sh -Xmx30g in="$seqfilefilt" out=$basedr/1_cutadapt/"$(basename "${seqfilefilt%_filter*}" )_cutadapt.fastq" literal=GGTCTTATCATCCAACCTGT ktrim=l k=20 mink=11 hdist=3 rcomp=f

done

#--------------------------------------------------------------------------

#cut sequences to length X bp -bioawk - change startposition & length in substr($seq,start,length)

bawk_sizecut='{ print "@"$name" "$comment; print substr($seq,1,70); print "+"; print substr($qual,1,70);}'

parallel "bioawk -c fastx '$bawk_sizecut' {} > $basedr/2_sizecut/{/.}_sizecut.fastq" ::: $basedr/1_cutadapt/*.fastq

#--------------------------------------------------------------------------

#keep only sequences that start with TA - bioawk

bawk_keepta='{ if ($seq ~ /^TA/) { print "@"$name" "$comment; print $seq; print "+"; print $qual;}}'

parallel "bioawk -c fastx '$bawk_keepta' {} > $basedr/3_TAonly/{/.}_TAonly.fastq" ::: $basedr/2_sizecut/*.fastq

#--------------------------------------------------------------------------

#run aligner to generate SAM files - build index in RefIndex directory and fill-in index name here.

for alreffile in $basedr/3_TAonly/*.fastq;
do
hisat2 -q -x $basedr/RefIndex/$alignidxname -U "$alreffile" -S $basedr/4_alignfiles/"$(basename "${alreffile%_cutadapt*}" )_aligned.sam" --rdg 10000,10000 --rfg 10000,10000 --no-unal --no-spliced-alignment --no-softclip -p $threadnum
done
#replace this step with Bowtie2 or BWA-mem if that's your preference.

# filter for samflags 0 (forward read) or 16 (reverse read)
bawk_samfilter='$flag==0 || $flag==16'
parallel "bioawk -c sam -H '$bawk_samfilter' {} > $basedr/5_sam016/{/.}_filt.sam" ::: $basedr/4_alignfiles/*.sam

# Generate indexed bam

for samtobam in $basedr/5_sam016/*_filt.sam;
do
samtools sort -o $basedr/6_idxbamfiles/"$(basename "${samtobam%_aligned*}" ).bam" -@ $samthread $samtobam
done

parallel "samtools index {}" ::: $basedr/6_idxbamfiles/*.bam

#--------------------------------------------------------------------------
# Find read frequencies for each insertion

bawk_samposition='{ if($flag==0) print $pos+1; if($flag==16)  print $pos+length($seq)-1; }'

parallel "bioawk -c sam '$bawk_samposition' {} | cut -f 1 | LANG=en_EN sort -g -k1 | uniq -c | sed 's/^ *//' | LANG=en_EN sort -b -k2 > $basedr/7_readfreqs/{/.}_readfreq.txt" ::: $basedr/4_alignfiles/*.sam

# Map reads to TA positions
 
parallel "join -1 2 -2 1 {} $basedr/reffiles/$TAgenemap_sort | sort -n -r -k2 | tr ' ' '\t' > $basedr/8_mappedinserts/{/.}_mapped.txt" ::: $basedr/7_readfreqs/*_readfreq.txt





