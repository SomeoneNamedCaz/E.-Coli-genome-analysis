# SnpAlign

SnpAlign is a tool for running genome-wide association studies (GWAS) from assembly files.
It uses ragtag for scaffolding and gsAlign to get the variants.
It used to require MegaCats to do the chi-square test but now does that internally.


## Installation
    usage: 
        snpalign [options] --ref <reference genbank file> --metadata <tsv with metadata categories> --fastas <fasta genome assembly files
    optional:
            --name: specifies a name to add to each output file
            --outdir: specify the output directory
    

## scripts to download from NCBI
[DownloadingFromNCBI](src%2Fsnpalign%2FDownloadingFromNCBI) contains the python scripts to download the assemblies used for the analysis.
A zip folder with these files can also be found at: https://drive.google.com/file/d/1zYs9sBaFS4PnkYI4bdSdWTbiVavh3cdP/view?usp=sharing
- downloadingFromSequenceSetBrowser.py was used for downloading the mastitis strains
- downloadingFromXMLFileOfNucIDs.py was used for downloading the commensal strains

## Descriptions of the Most Important Files

### MakingMetadataSheetForMetaCats.py
This program makes the metadata that is required for the chisquare analysis (https://github.com/bpickett/megaCATS). It doesn't support commandline arguments.

### alignVcfSnps.py
This takes .vcf files from gsAlign and outputs a fake "snp" genome that contains all snps that are present in ten or more genomes. It also outputs the original indexes of the snps (so it is possible to map the snp genome back to the original one) and the original indexes of frame shifts. This is to bypass the multiple alignment process of mauve that would take too long to run on >1,000 genomes that we would need them to run.

The snp genome can be used as the input for megaCats, and runs in ~10 minutes.

This program works by putting all the snps in a list. This list is sorted, and looped through. For each snp the snp genome associated with that snp is appended the snp. Once all the snps at a certain position have been entered, the snps are aligned. This alignment process works by adding the first nucleotide of the reference nucleotides corresponding to the snp (this is only more than one if there is a deletion at the position. Then dashes are added equal to the length of the longest insert or the length of the longest insert minus the length of the insert for that genome. Then reference nucleotides are added to make each genome the same length. If a snp at the next position will be within a snp at the current position, the alignment stops before the next snp's position, so it can be added.

### parsingMegaCatsResults.py
This takes the statistics from megaCats and interprets them: It requires the combined megaCats statistics, the file with the snp genomes (from gsAlignToMultipleAlign.py), the indexes also from gsAlignToMultipleAlign.py, a suffix to be added to the output file names, the reference annotated genome (.gb or .gbff). Optional commandline arguments are a flag to remove the sparce entries (according to megaCats), which is on by default, and the output directory, which is the current directory by default. It outputs a file of the significant snps that looks like the input megaCats stats. It also outputs two files for each metadata category: one looks at the genes that have the most snps in them; the other looks at each snp individually, ranks them by significance, and determines the gene that it is in. Both files look at the type of mutations that the snps provide.

This program includes snps that are outside of genes in the sorted snp files. It doesn't include them in the file that looks at the genes with the most snps

### functions.py
This is just a python file with the functions that are shared across the repository.
