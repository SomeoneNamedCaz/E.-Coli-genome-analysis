# NOT Working yet
# $1 is conda ragtag env name, $2 is conda gsalign env name, $3 is megaCats environment nam, $4 is refGenome.fasta, $5 is refGenome.gbff, $6 is all of the genomes to analyze (fastas), $7 metadata in a tsv.
# must be used in folder that has all of the phython code
conda activate $1 && mkdir ragtagOutputs && for fileName in $6 && do ragtag.py scaffold $4 $fileName -t 16 -o ./ragtagOutputs/$fileName/ && done && python getLongestContig.py . && cd ragtagOutputs && cd longestScaffoldFiles && conda activate $2 && mkdir gsAlignOutputs && for fileName in $6 && do gsAlign -r  $4 -q $fileName -o ./gsAlignOutputs/$fileName -t 16 && done || echo "requires conda ragtag env name, conda gsalign env name, megaCats environment name (just need correct r and perl versions), refGenome.fasta, refGenome.gbff, all of the genomes to analyze (fastas), metadata in a tsv format (needs to addd scaffold_ to the front and .vcf to the back of the genome names)."
cd ..; cd ..
python gsalignSNPsToMultipleAlign.py "./ragtagOutputs/longestScaffoldFiles/gsAlignOutputs/*.vcf" "./SnpAlignment/allSnps.afa"
conda activate $3 # megaCats
perl metadata_parser.pl $7 "./SnpAlignment/allSnps.afa"
python parsingMetaCatsREsults 1* ./SnpAlignment/allSnps.afa ./SnpAlignment/allSnpsIndexes.afa $4 $7