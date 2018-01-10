#! /bin/bash

# Nagendra Palani - University of Minnesota Genomics Center - 01/2018 - MN, USA
# contact : nagendra [AT] umn //dot// edu

#Read through https://github.com/jbadomics/tnseq for verbose explanation of analysis. Hat tip to jbadomics for Bioawk.

#Pipeline for TnSeq data preprocessing - setup for UMN MSI Mesabi HPC. 
#Uses Read_1 files from paired end sequencing to look for transposon junctions & extracts adjacent genomic sequence for mapping.

#Input - demultiplexed Fastq (Mutiple Tnseq fastq files from a single sequencing run)
#Output - mapped TA positions with raw readcounts

#Tool Dependencies - bbtools, Hisat2, GNU Parallel, Bioawk - Script is fully parallelized.

#File dependencies - ReferenceOrganismGenome.fasta, list of TA positions in Reference Genome
#Create a folder 'reffiles' within base directory & place dependency files within reffiles.
 

#modules to be loaded
module load parallel
module load hisat2

#base folder containing sequencing data - copy fastq files from UMN datarelease folder to a local profile and change path of base directory appropriately.
basedr=~/Downloads/tnseqtemp


reffasta=ReferenceOrganismGenome.fasta
alignidxname="$(basename "${reffasta%.fasta*}")_idx"

TAgenemap=RefOrganism_TA_position.txt
#Make sure to sort TA list alphabetically (shell command : sort -b )

#hisat2 threads - set threadnum equal to number of processor cores available / requested in PBS job.
threadnum=20

#sub-directories to be created for outputs

mkdir $basedr/RefIndex
mkdir $basedr/0_filter
mkdir $basedr/1_cutadapt
mkdir $basedr/2_sizecut
mkdir $basedr/3_TAonly
mkdir $basedr/4_alignfiles
mkdir $basedr/5_readfreqs
mkdir $basedr/6_mappedinserts

#--------------------------------------------------------------------------
# Build required files

#Build Hisat2 aligner index

hisat2-build -f $basedr/reffiles/$reffasta $basedr/RefIndex/$alignidxname

#--------------------------------------------------------------------------

#filter sequences for valid transposons - bbduk (bbtools package) - specify filtering seq in literal (transposon inverted repeat sequence from 5'end of read)

for seqfile in $basedr/*.fastq;
do
bbduk.sh -Xmx4g in="$seqfile" outm=$basedr/0_filter/"$(basename "${seqfile%_R1_001*}" )_filtered.fastq" literal=GGACTTATCAGCCAACCTGT k=20 hdist=3 rcomp=f maskmiddle=f

done
#--------------------------------------------------------------------------

#trim 5' adapters - bbduk (bbtools package) - specify adapter seq in literal (transposon inverted repeat sequence from 5'end of read) - can be merged with previous step for same literal - separated here for code flexibility.

for seqfilefilt in $basedr/0_filter/*.fastq*;
do
bbduk.sh -Xmx4g in="$seqfilefilt" out=$basedr/1_cutadapt/"$(basename "${seqfilefilt%_filter*}" )_cutadapt.fastq" literal=GGACTTATCAGCCAACCTGT ktrim=l k=20 mink=11 hdist=3 rcomp=f

done

#--------------------------------------------------------------------------

#cut genomic junction sequences to length X bp -bioawk - change startposition & length in substr($seq,start,length)

bawk_sizecut='{ print "@"$name" "$comment; print substr($seq,1,25); print "+"; print substr($qual,1,25);}'

parallel "bioawk -c fastx '$bawk_sizecut' {} > $basedr/2_sizecut/{/.}_sizecut.fastq" ::: $basedr/1_cutadapt/*.fastq

#--------------------------------------------------------------------------

#keep only sequences that start with TA - bioawk

bawk_keepta='{ if ($seq ~ /^TA/) { print "@"$name" "$comment; print $seq; print "+"; print $qual;}}'

parallel "bioawk -c fastx '$bawk_keepta' {} > $basedr/3_TAonly/{/.}_TAonly.fastq" ::: $basedr/2_sizecut/*.fastq

#--------------------------------------------------------------------------

#run aligner to generate SAM files - build index in RefIndex directory and fill-in index name here.

for alreffile in $basedr/3_TAonly/*.fastq;
do
hisat2 -q -x $basedr/RefIndex/$alignidxname -U "$alreffile" -S $basedr/4_alignfiles/"$(basename "${alreffile%_cutadapt*}" )_aligned.sam" --no-unal -p $threadnum

done

#--------------------------------------------------------------------------
# Find read frequencies for each insertion

bawk_samposition='{ if($flag==0) print $pos+1; if($flag==16)  print $pos+length($seq)-1; }'

parallel "bioawk -c sam '$bawk_samposition' {} | cut -f 1 | sort -g -k1 | uniq -c | sed 's/^ *//' | sort -b -k2 > $basedr/5_readfreqs/{/.}_readfreq.txt" ::: $basedr/4_alignfiles/*.sam
#--------------------------------------------------------------------------
# Map reads to TA positions

parallel "join -1 2 -2 1 {} $basedr/reffiles/$TAgenemap | sort -n -r -k2 | tr ' ' '\t' > $basedr/6_mappedinserts/{/.}_mapped.txt" ::: $basedr/5_readfreqs/*_readfreq.txt




