# maybe something like this
#       |gene A start       geneAStop|
# ACTAGATAGAGCATAGAGTAACTAGATAGAGCATAGAGTA
# ACTAGATAGAGCATAGAGTACCTAGATAGAGCATAGAGTA
#                     |SNP misSense T->K
from functions import *
from reconstructNormalAlignment import *
from parsingMegaCatsResults import getMutationType


# if there's a dash in the genome that doesn't have any snps then we know there was an insert so if we account
# for all the inserts we can figure out where the genes and snps should be

def format(originalString, newString, endLengthOfString):
	availableSpace = endLengthOfString - len(originalString) - 2
	if availableSpace < len(newString): # add new string to the right
		return originalString + " " * (endLengthOfString - len(originalString) - 2) + "->|<-" + newString
	else: # add it to the right
		return originalString + " " * (endLengthOfString - len(originalString) - 2 - len(newString)) + newString + "->|<-"
def moveToNextGene(nextGeneIndex,refGenes, onlyForward=True ):
	if nextGeneIndex != len(refGenes) - 1:
		nextGeneIndex += 1
	if onlyForward:
		while nextGeneIndex != len(refGenes) - 1 and not refGenes[nextGeneIndex].isForward:
			nextGeneIndex += 1
	else:
		while nextGeneIndex != len(refGenes) - 1 and refGenes[nextGeneIndex].isForward:
			nextGeneIndex += 1
	return nextGeneIndex

def mapsSnps(refGenomePath, snpIndexes, alignedSnps, significantSnps):
	refGenomeContigs = getContigs(refGenomePath)
	if len(refGenomeContigs) > 1:
		raise Exception("ref genome has more than one contig, this is not a problem as long as the first contig is the chromosome")
		
	refSeq = refGenomeContigs[0]
	refGenes = getGenesOnContigsByPosition(refGenomePath, refGenomeContigs)
	print("len(refGenes)",len(refGenes))
	print("pre align")
	alignedGenomes = reconstructNormalAlignmentHelper(refSeq, snpIndexes, {"genomeWithoutAnySnps":alignedSnps["genomeWithoutAnySnps"]})
	print("len align", len(alignedGenomes["genomeWithoutAnySnps"]))
	print("post align")
	
	# significantSnps
	
	insertIndexes = []
	numShifts = 0
	indexThatMatchesOriginalSeq = 0 # this is index not counting the inserts
	nextForwardGeneIndex = 0
	nextReverseGeneIndex = 0
	forwardGeneLine = ""#[] # change to list of chars if slow
	reverseGeneLine = ""
	snpLine = ""
	index = 0
	estimatedSNPPosition = 0
	# indexOfNextSnpLocation = 0
	for nuc in alignedGenomes["genomeWithoutAnySnps"]:
		# print(indexThatMatchesOriginalSeq)
		if indexThatMatchesOriginalSeq == refGenes[nextForwardGeneIndex].startPos :#or indexThatMatchesOriginalSeq == refGenes[nextGeneIndex + 1].startPos:
			forwardGeneLine = format(forwardGeneLine, refGenes[nextForwardGeneIndex].name + " start", index)
		elif indexThatMatchesOriginalSeq == refGenes[nextForwardGeneIndex].stopPos:
			forwardGeneLine = format(forwardGeneLine, refGenes[nextForwardGeneIndex].name + " end", index)
			nextForwardGeneIndex = moveToNextGene(nextForwardGeneIndex,refGenes,onlyForward=True)
			# print(nextForwardGeneIndex)
		if indexThatMatchesOriginalSeq == refGenes[nextReverseGeneIndex].startPos:
			reverseGeneLine = format(reverseGeneLine, refGenes[nextReverseGeneIndex].name + " end (gene is on complement strand)", index)
		elif indexThatMatchesOriginalSeq == refGenes[nextReverseGeneIndex].stopPos:
			reverseGeneLine = format(reverseGeneLine,refGenes[nextReverseGeneIndex].name + " start (gene is on complement strand)", index)
			nextReverseGeneIndex = moveToNextGene(nextReverseGeneIndex, refGenes, onlyForward=False)
		# try:
		for indexOfNextSnpLocation in range(max(estimatedSNPPosition - 5,0), min(estimatedSNPPosition + 5,len(significantSnps))):
			if indexThatMatchesOriginalSeq == int(significantSnps[indexOfNextSnpLocation][11]) + int(significantSnps[indexOfNextSnpLocation][2]):
				# print("added a snp")
				indexOfPathogenGroup = 8
				indexOfCommensalGroup = 5
				if significantSnps[indexOfNextSnpLocation][1] == "pathogen":
					indexOfPathogenGroup = 5
					indexOfCommensalGroup = 8
				# getMutationType()
				shortenedSnp = [significantSnps[indexOfNextSnpLocation][0],
				                                            significantSnps[indexOfNextSnpLocation][2],
				                                            "pathogen nuc:", significantSnps[indexOfNextSnpLocation][indexOfPathogenGroup],
				                                            "commensal nuc:", significantSnps[indexOfNextSnpLocation][indexOfCommensalGroup],
				                                            significantSnps[indexOfNextSnpLocation][-2]]
				snpLine = format(snpLine, " ".join(shortenedSnp), index)
				# indexOfNextSnpLocation += 1
				estimatedSNPPosition += 1
		# except IndexError:
		# 	pass
		if nuc != "-":
			indexThatMatchesOriginalSeq += 1
		index += 1
	
	return forwardGeneLine,reverseGeneLine, alignedGenomes["genomeWithoutAnySnps"], snpLine
	
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






