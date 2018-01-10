### Generate a list of all TA positions for a given organism.


#### The following MATLAB commands will read a genome sequence in fasta format, find positions of all TA sites, and write the output to a txt file. You can do this in your desired programming language.

```matlab
reffasta=fastaread('ReferenceOrganismGenome.fasta'); % download from NCBI.
TApositions = transpose(strfind(upper(reffasta.Sequence),'TA'))+1 ;
writetable(TApositions,'Ref_org_TAposition.txt','WriteVariableNames',0,'delimiter','\t');
```

#### Sort the file alphabetically. 
```shell
sort -b -k1 Ref_org_TAposition.txt > Ref_org_TAposition_sorted.txt
```
#### Use the sorted file with the fastqtoreads.sh script. 
