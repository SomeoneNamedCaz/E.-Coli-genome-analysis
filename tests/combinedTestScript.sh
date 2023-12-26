#!/bin/sh



cd /Users/cazcullimore/Documents/ericksonLabCode/tests/tenAssembliesFromEachCategory;
conda activate gsAlign; mkdir gsAlignOutputs; for fileName in *.fasta; do gsAlign -r   /Users/cazcullimore/Documents/ericksonLabCode/refGenomes/k-12.fasta -q $fileName -o scaffold_$fileName -t 16; done;
mv *.vcf ./gsAlignOutputs/;
rm *.maf;
time python gsalignSNPsToMultipleAlign.py DATA_DIR + "gsAlignOutputs/*.vcf" ./SnpAlign/carefulAndFixedUsingOppositeStrand.afa;

time python megaCatsPythonVersion.py ./SnpAlign/carefulAndFixedUsingOppositeStrand.afa /Users/cazcullimore/Documents/ericksonLabCode/metaDataForMetaCatsWithExtraMastitis.tsv ./megaCatsTest.txt True;
time python /Users/cazcullimore/Documents/ericksonLabCode/parsingMegaCatsResults.py
time python convertSNPAlignmentToNormalAlignment.py ./SnpAlign/carefulAndFixedUsingOppositeStrand.afa ./SnpAlign/carefulAndFixedUsingOppositeStrandIndexes.txt /Users/cazcullimore/Documents/ericksonLabCode/refGenomes/k-12.gbff carefulNormalAlign.afa;

time python testingMultiAlignusingAnnotatedFiles.py