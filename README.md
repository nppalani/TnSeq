# TnSeq
Scripts to process Mariner transposon based Tnseq sequencing files to readcounts per insertion position. Works for sequencing libraries prepared by TraDIS or INSeq. Can be adapted for Tn-Seq (MmeI based) with very minor modification.

The shell script ['fastqtoreadcount.sh'](https://github.com/nppalani/TnSeq/blob/master/fastqtoreadcount.sh) takes Tnseq sequencing data (fastq) files as input and outputs raw readcounts per insertion position. See [sample output](https://github.com/nppalani/TnSeq/blob/master/Sample_Output_aligned_readfreq_mapped.txt)

[Instructions](https://github.com/nppalani/TnSeq/blob/master/TA_position_list.md) to create a file containing all TA positions in a genome 
