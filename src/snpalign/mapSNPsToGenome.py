# maybe something like this
#       |gene A start       geneAStop|
# ACTAGATAGAGCATAGAGTAACTAGATAGAGCATAGAGTA
# ACTAGATAGAGCATAGAGTACCTAGATAGAGCATAGAGTA
#                     |SNP misSense T->K
try:
	from .reconstructNormalAlignment import *
except ImportError:
	from reconstructNormalAlignment import *

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

def mapsSnps(refGenomePath, snpIndexes, alignedSnps, significantSnps, numThreads=16):
	refGenomeContigs = getContigs(refGenomePath)
	if len(refGenomeContigs) > 1:
		raise Exception("ref genome has more than one contig, this is not a problem as long as the first contig is the chromosome")
		
	refSeq = refGenomeContigs[0]
	refGenes = getGenesOnContigsByPosition(refGenomePath, refGenomeContigs)
	
	alignedGenomes = reconstructNormalAlignmentHelper(refSeq, snpIndexes, alignedSnps,numThreads=numThreads)#{"genomeWithoutAnySnps":alignedSnps["genomeWithoutAnySnps"], "scaffold_k-12.fasta.vcf":alignedSnps["scaffold_k-12.fasta.vcf"]})
	
	
	indexThatMatchesOriginalSeq = 0 # this is index not counting the inserts
	nextForwardGeneIndex = moveToNextGene(-1,refGenes,onlyForward=True)
	nextReverseGeneIndex = moveToNextGene(-1, refGenes, onlyForward=False)
	forwardGeneLine = ""#[] # change to list of chars if slow
	reverseGeneLine = ""
	
	numSnpLines = 8
	snpLines = []
	for _ in range(numSnpLines):
		snpLines.append("")
	index = 0
	estimatedSNPPosition = 0
	
	for nuc in alignedGenomes[nameOfRefSnpGenome]:
		
		if indexThatMatchesOriginalSeq == refGenes[nextForwardGeneIndex].startPos :#or indexThatMatchesOriginalSeq == refGenes[nextGeneIndex + 1].startPos:
			forwardGeneLine = format(forwardGeneLine, refGenes[nextForwardGeneIndex].name + " start", index)
		elif indexThatMatchesOriginalSeq == refGenes[nextForwardGeneIndex].stopPos:
			forwardGeneLine = format(forwardGeneLine, refGenes[nextForwardGeneIndex].name + " end", index)
			nextForwardGeneIndex = moveToNextGene(nextForwardGeneIndex,refGenes,onlyForward=True)
		if indexThatMatchesOriginalSeq == refGenes[nextReverseGeneIndex].startPos:
			reverseGeneLine = format(reverseGeneLine, refGenes[nextReverseGeneIndex].name + " end (gene is on complement strand)", index)
		elif indexThatMatchesOriginalSeq == refGenes[nextReverseGeneIndex].stopPos:
			reverseGeneLine = format(reverseGeneLine,refGenes[nextReverseGeneIndex].name + " start (gene is on complement strand)", index)
			nextReverseGeneIndex = moveToNextGene(nextReverseGeneIndex, refGenes, onlyForward=False)
		# try:
		for indexOfNextSnpLocation in range(max(estimatedSNPPosition - 5,0), min(estimatedSNPPosition + 5,len(significantSnps))):
			if indexThatMatchesOriginalSeq == int(significantSnps[indexOfNextSnpLocation][11]) + int(significantSnps[indexOfNextSnpLocation][2]):

				indexOfPathogenGroup = 8
				indexOfCommensalGroup = 5
				if significantSnps[indexOfNextSnpLocation][1] == "pathogen":
					indexOfPathogenGroup = 5
					indexOfCommensalGroup = 8
				shortenedSnp = [significantSnps[indexOfNextSnpLocation][0],
				                                            significantSnps[indexOfNextSnpLocation][2],
				                                            "pathogenic nuc is", significantSnps[indexOfNextSnpLocation][indexOfPathogenGroup] + ",",
				                                            "commensal nuc is", significantSnps[indexOfNextSnpLocation][indexOfCommensalGroup] + ",",
				                                            significantSnps[indexOfNextSnpLocation][-2]]
				snpLines[estimatedSNPPosition % len(snpLines)] = format(snpLines[estimatedSNPPosition % len(snpLines)], " ".join(shortenedSnp), index)

				estimatedSNPPosition += 1

		if nuc != "-":
			indexThatMatchesOriginalSeq += 1
		index += 1
	
	geneLines = [forwardGeneLine,reverseGeneLine]
	# snpLines = [firstSnpLine, secondSnpLine]
	return geneLines, alignedGenomes, snpLines
	
def writeMappedSnps( namePrefix,refGenomePath,snpAlignPath, snpIndexPath,sigSnpsPath,metadataPath, outDir="./", debug=False, findM12andK12=False, numThreads=16):
	
	if debug:
		print(sigSnpsPath)
	alignedSnps = readInFastaAsDict(snpAlignPath)
	indexes = readIndexes(snpIndexPath)
	nameOfSigsnpsFile = "snpsSortedBySignificanceWithGenesContainingThem" + namePrefix
	for sigSnpsFile in glob(sigSnpsPath + nameOfSigsnpsFile + "*.tsv"):
		significantSnps = []
		with open(sigSnpsFile) as file:
			isFirstLine = True
			for line in file:
				if isFirstLine:
					isFirstLine = False
					continue
				lineData = line.split("\t")[:12]
				
				significantSnps.append(lineData)
		
		significantSnps.sort(key=lambda x: (int(x[11]) + int(x[2])))
		if debug:
			print(significantSnps[1])
		geneLines, alignedGenomes, snpLines = mapsSnps(refGenomePath=refGenomePath, snpIndexes=indexes, alignedSnps=alignedSnps,
		                   significantSnps=significantSnps, numThreads=16)
	
		metadata = readMetaDataAsDict(metadataPath)
		
		blacklist = {nameOfRefSnpGenome}
		if findM12andK12:
			blacklist.add("scaffold_k-12.fasta.vcf")
		allgenomesWithoutK12M12 = {re.sub("\..+", "",key): value for key, value
		                           in alignedGenomes.items() if key not in blacklist}
		longestGenomeNameLen = max([len(name + metadata[name][1]) for name in allgenomesWithoutK12M12.keys()])
		shift = longestGenomeNameLen * " "
		
		
		with open(outDir + namePrefix + re.sub(nameOfSigsnpsFile + "|\..+","",sigSnpsFile) + "mappedSnps.txt", "w") as file:
			nCharsToWrite = 150
			for x in range(0, len(alignedGenomes[nameOfRefSnpGenome]) - nCharsToWrite, nCharsToWrite):
				for geneLine in geneLines:
					file.write(shift + geneLine[x:x + nCharsToWrite] + "\n")
				for snpLine in snpLines:
					file.write(shift + snpLine[x:x + nCharsToWrite] + "\n")
				if findM12andK12:
					file.write(("M12:" + " " * (longestGenomeNameLen - len("M12:"))) + alignedGenomes[nameOfRefSnpGenome][x:x + nCharsToWrite] + "\n")
					file.write(("K12:" + " " * (longestGenomeNameLen - len("K12:"))) + alignedGenomes["scaffold_k-12.fasta.vcf"][x:x + nCharsToWrite] + "\n")
				else:
					file.write(
						(nameOfRefSnpGenome + ":" + " " * (longestGenomeNameLen - len(nameOfRefSnpGenome) + 1)) + alignedGenomes[nameOfRefSnpGenome][
						                                                        x:x + nCharsToWrite] + "\n")
				for genomeName, seq in allgenomesWithoutK12M12.items():
					ogGenomeName = deepcopy(genomeName)
					genomeName = re.sub("\..+", "",genomeName)
					genomeShift = genomeName + metadata[ogGenomeName][1] + " " * (longestGenomeNameLen - len(genomeName + metadata[ogGenomeName][1]))
					file.write(genomeShift + seq[x:x + nCharsToWrite] + "\n")
				file.write("#" * (nCharsToWrite + len(shift)) + "\n")
			
		if debug:
			print("wrote out data to:", outDir + "mappedSnps.txt")

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






