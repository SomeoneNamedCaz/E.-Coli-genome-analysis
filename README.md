# E.-Coli-genome-analysis

code to analyze E. coli genome differences across species by their animal host or by pathogenicity.

This code can be downloaded and each python file is run as its own program using the terminal or the ide. Most programs require commandline arguments to specify the input files. However, most of the files are divided into two projects: one looking at gene presence or absence in the pan genome, and another that looks at snps. There are also some file management programs and programs that check that things are working properly.

# SortNatureAssemblies.py
This takes the assemblies from https://doi.org/10.1038/s41467-021-20988-w and groups them into pathogenic and commensal strains. This doesn't support commandline arguments.

# downloadingXMLFileOfNucIDs.py
This program runs in python 3.8 (maybe other versions too, but it doesn't work on some versions of python). It uses an XML file that has the IDs of the Nucleotide entries in NCBI, and downloads the associated files. You can get the XML file through using the eutils for NCBI and using esearch and requesting all of the entries and then downloading that page using your browser.

# downloadingFromSequenceSetBrowser.py
It also downloads the genomes using bioproject names. This program used a csv from ncbi to download the mastitis genomes (this file is in the repo as WGSProjectNamesWithExtraStuffRemoved.tsv which comes from allMastitisBioProjects.csv)

# unannotatedGbFilesToFasta.py
This takes files downloaded from ncbi using the downloading python programs and combines them into a fasta file for each genome. This is done because the files directly downloaded only contained a single contig of the genome.

# MakingMetadataSheetForMetaCats.py
This program makes the metadata that is required for megacats (https://github.com/bpickett/megaCATS). It doesn't support commandline arguments.

# getLongestContig.py
This program takes a scaffold, gets the top contig and outputs this as a file. It was used to get the longest contig on a ragtag scaffold. This allows us to use a regular single alignment program (we used gsAlign). This works with gsAlignToMultipleAlign.py to create a multiple alignment without using Mauve. This step loses some of the DNA, but we only lost ~0.5Mb of the E. coli genomes.

# gsAlignToMultipleAlign.py
This takes .vcf files from gsAlign and outputs a fake "snp" genome that contains all of the snps that are present in ten or more genomes. It also outputs the original indexes of the snps (so it is possible to map the snp genome back to the original one) and the original indexes of frame shifts. This is to bypass the multiple alignment process of mauve that would take too long to run on the almost 1,000 genomes that we would need them to run.

The snp genome can be used as the input for megaCats.

# parsingMegaCatsResults.py
This takes the statistics from megaCats and interprets them. It outputs a file of the significant snps that looks like the input megaCats stats. It also outputs two files for each metadata category: one looks at the genes that have the most snps in them; the other looks at each snp individually, ranks them by significance, and determines the gene that it is in. Both files look at the type of mutations that the snps provide.

This program doesn't look at snps that aren't within a gene.

# fixingPathogenicityFile.py
When I ran megaCats the pathogenicity correlated snps had a third group. Since I didn't know the cause I worked around it. However, I have now discovered the cause (no newline at the end of the file), so it probably isn't needed anymore.

# findingOddGenomeOut.py
This program discovers which genome is causing the problem with the third pathogenicity group.

# functions.py
This is just a python file with the functions that are shared across the repository.

# concatenateContigs.py
This turns an assembly into an assembly with all of the contigs added together into a single contig. This was an attempt to make Mauve faster.

# pValueHistogram.py
This makes a histogram from the pvalues from a megaCats output file.

# turningGsAlignIntoMultipleAlign.py
This was an attempt that didn't work out because it wasn't able to recover enough of the genome. It used the aligned sequence files from gsAlign instead of the vcf snp files.

# parallelgsalignSNPsToMultipleAlign.py
this was an attempt at multithreading, but it didn't result in much speed up (maybe needed to copy the set over). The single threaded version is pretty fast right now in comparison (~2hrs).


