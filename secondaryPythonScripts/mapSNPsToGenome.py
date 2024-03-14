# maybe something like this
#       |gene A start       geneAStop|
# ACTAGATAGAGCATAGAGTAACTAGATAGAGCATAGAGTA
# ACTAGATAGAGCATAGAGTACCTAGATAGAGCATAGAGTA
#                     |SNP misSense T->K
from functions import *
from reconstructNormalAlignment import *


# if there's a dash in the genome that doesn't have any snps then we know there was an insert so if we account
# for all the inserts we can figure out where the genes and snps should be

def mapsSnps(refGenomePath, snpIndexes, alignedSnps, significantSnps):
	refGenomeContigs = getContigs(refGenomePath)
	if len(refGenomeContigs) > 1:
		raise Exception("ref genome has more than one contig, this is not a problem as long as the first contig is the chromosome")
		
	refSeq = refGenomeContigs[0]
	refGenes = getGenesOnContigsByPosition(refGenomePath, refGenomeContigs)
	print("pre align")
	alignedGenomes = reconstructNormalAlignmentHelper(refSeq, snpIndexes, {"genomeWithoutAnySnps":alignedSnps["genomeWithoutAnySnps"]})
	print("post align")
	
	# significantSnps
	
	insertIndexes = []
	numShifts = 0
	indexThatMatchesOriginalSeq = 0 # this is index not counting the inserts
	nextGeneIndex = 0
	annotationLine = ""#[] # change to list of chars if slow
	snpLine = ""
	index = 0
	indexOfNextSnpLocation = 0
	for nuc in alignedGenomes["genomeWithoutAnySnps"]:
		if indexThatMatchesOriginalSeq == refGenes[nextGeneIndex].startPos:
			stuffToAdd = " " * (index - len(annotationLine) - 2 - len(refGenes[nextGeneIndex].name)) + refGenes[nextGeneIndex].name + "->|<-"
			if refGenes[nextGeneIndex].isForward:
				annotationLine += stuffToAdd + " start"
			else:
				annotationLine += stuffToAdd + " end (gene is on complement strand)"
		elif indexThatMatchesOriginalSeq == refGenes[nextGeneIndex].stopPos:
			stuffToAdd = " " * (index - len(annotationLine) - 2 - len(refGenes[nextGeneIndex].name)) + refGenes[nextGeneIndex].name + "->|<-"
			if refGenes[nextGeneIndex].isForward:
				annotationLine += stuffToAdd + " end"
			else:
				annotationLine += stuffToAdd + " start (gene is on complement strand)"
			if nextGeneIndex != len(refGenes) - 1:
				nextGeneIndex += 1
		if indexThatMatchesOriginalSeq == significantSnps[indexOfNextSnpLocation][10]:
			snpLine += " " * (index - len(snpLine) - 2) + "->|<-" + " ".join(significantSnps[indexOfNextSnpLocation])
			indexOfNextSnpLocation += 1
		if nuc != "-":
			indexThatMatchesOriginalSeq += 1
		index += 1
	
	return annotationLine, alignedGenomes, snpLine
	
	# for index in range(len()):

if __name__ == "__main__":
	try:
		snpGenomePath = sys.argv[1]
		snpIndexesPath = sys.argv[2]
		print("snp genome path", snpGenomePath)
		print("index path", snpIndexesPath)
		
		refGenomePath = sys.argv[3]  # .gbff
		
		outputPath = sys.argv[4]
	except IndexError:
		print("\nplease provide the path to the snp genome, the corresponding indexes, the reference sequence (.gb)")
		print("and the output path output path")
		exit(1)
	reconstructNormalAlignment(snpGenomePath, snpIndexesPath, refGenomePath, outputPath)






