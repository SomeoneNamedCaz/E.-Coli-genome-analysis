# E.-Coli-genome-analysis

code to analyze E. coli genome differences by their animal host or by pathogenicity.

This code can be downloaded and each python file is run as its own program using the terminal or the ide. Most programs require commandline arguments to specify the input files. However, most of the files are divided into two projects: one looking at gene presence or absence in the pan genome, and another that looks at snps. There are also some file management programs and programs that check that things are working properly.

# SortNatureAssemblies.py
This takes the assemblies from https://doi.org/10.1038/s41467-021-20988-w and groups them into pathogenic and commensal strains. This doesn't support commandline arguments.

# downloadingXMLFileOfNucIDs.py
This program runs in python 3.8 (maybe other versions too, but it doesn't work on some versions of python). It uses an XML file that has the IDs of the Nucleotide entries in NCBI, and downloads the associated files. You can get the XML file through using the eutils for NCBI and using esearch and requesting all of the entries and then downloading that page using your browser.

# downloadingFromSequenceSetBrowser.py
It also downloads the genomes using bioproject names. This program used a csv from ncbi to download the mastitis genomes (this file is in the repo as WGSProjectNamesWithExtraStuffRemoved.tsv which comes from allMastitisBioProjects.csv)

# unannotatedGbFilesToFasta.py
This takes files downloaded from ncbi using the downloading python programs and combines them into a fasta file for each genome. This is done because the files directly downloaded only contained a single contig of the genome.
  - combiningFastaFiles.py is the version of this program that combines fasta contig files instead of gb files.

# MakingMetadataSheetForMetaCats.py
This program makes the metadata that is required for megacats (https://github.com/bpickett/megaCATS). It doesn't support commandline arguments.

# getLongestContig.py
This program takes a scaffold, gets the top contig and outputs this as a file. It was used to get the longest contig on a ragtag scaffold. This allows us to use a regular single alignment program (we used gsAlign). This works with gsAlignToMultipleAlign.py to create a multiple alignment without using Mauve. This step loses some of the DNA, but we only lost ~0.5Mb of the E. coli genomes.

# gsAlignToMultipleAlign.py
This takes .vcf files from gsAlign and outputs a fake "snp" genome that contains all of the snps that are present in ten or more genomes. It also outputs the original indexes of the snps (so it is possible to map the snp genome back to the original one) and the original indexes of frame shifts. This is to bypass the multiple alignment process of mauve that would take too long to run on >1,000 genomes that we would need them to run.

The snp genome can be used as the input for megaCats, and runs in ~10 minutes.

This program works by putting all of the snps in a list with some information. This list is sorted, and looped through. For each snp the snp genome associated with that snp is appended the snp. Once all of the snps at a certain position have been entered, the snps are aligned. This alignment process works by adding the first nucleotide of the reference nucleotides corresponding to the snp (this is only more than one if there is a deletion at the position. Then dashes are added equal to the length of the longest insert or the length of the logest insert minus the length of the insert for that genome. Then refernce nucleotides are addded to make each genome the same length. If a snp at the next position will be within a snp at the current position, the alignment stops before the next snp's position so it can be added.

# parsingMegaCatsResults.py
This takes the statistics from megaCats and interprets them: It requires the combined megaCats statistics, the file with the snp genomes (from gsAlignToMultipleAlign.py), the indexes also from gsAlignToMultipleAlign.py, a suffix to be added to the output file names, the reference annotated genome (.gb or .gbff). Optional commandline arguments are a flag to remove the sparce entries (according to megaCats), which is on by default, and the output directory, which is the current directory by default. It outputs a file of the significant snps that looks like the input megaCats stats. It also outputs two files for each metadata category: one looks at the genes that have the most snps in them; the other looks at each snp individually, ranks them by significance, and determines the gene that it is in. Both files look at the type of mutations that the snps provide.

This program includes snps that are outside of genes in the sorted snp files. It doesn't inlcude them in the file that looks at the genes with the most snps

# functions.py
This is just a python file with the functions that are shared across the repository.

# concatenateContigs.py
This turns an assembly into an assembly with all of the contigs added together into a single contig. This was an attempt to make Mauve faster.

# pValueHistogram.py
This makes a histogram from the pvalues from a megaCats output file.

# AnnotatedGenomeAnalysis:
This is a file that looks at the frequencies of the snps that were identified here https://doi.org/10.1038/s41467-021-20988-w for our genomes. 
