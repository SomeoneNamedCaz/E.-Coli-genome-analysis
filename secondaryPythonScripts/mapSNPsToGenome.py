# maybe something like this
#       |gene A start       geneAStop|
# ACTAGATAGAGCATAGAGTAACTAGATAGAGCATAGAGTA
# ACTAGATAGAGCATAGAGTACCTAGATAGAGCATAGAGTA
#                     |SNP misSense T->K
import re

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

def mapsSnps(refGenomePath, snpIndexes, alignedSnps, significantSnps, metadataPath=DATA_DIR + "metaDataForMetaCatsWithExtraMastitis.tsv"):
	refGenomeContigs = getContigs(refGenomePath)
	if len(refGenomeContigs) > 1:
		raise Exception("ref genome has more than one contig, this is not a problem as long as the first contig is the chromosome")
		
	refSeq = refGenomeContigs[0]
	refGenes = getGenesOnContigsByPosition(refGenomePath, refGenomeContigs)
	print("len(refGenes)",len(refGenes))
	print("pre align")
	alignedGenomes = reconstructNormalAlignmentHelper(refSeq, snpIndexes, alignedSnps)#{"genomeWithoutAnySnps":alignedSnps["genomeWithoutAnySnps"], "scaffold_k-12.fasta.vcf":alignedSnps["scaffold_k-12.fasta.vcf"]})
	print("len align", len(alignedGenomes["genomeWithoutAnySnps"]))
	print("post align")
	print(alignedGenomes)
	cpdBStart = -1
	cdpBEnd = -1
	# {name: seq[] for name, seq in alignedGenomes.items()}
	# exit(0)
	# significantSnps
	
	insertIndexes = []
	numShifts = 0
	indexThatMatchesOriginalSeq = 0 # this is index not counting the inserts
	nextForwardGeneIndex = moveToNextGene(-1,refGenes,onlyForward=True)
	nextReverseGeneIndex = moveToNextGene(-1, refGenes, onlyForward=False)
	forwardGeneLine = ""#[] # change to list of chars if slow
	reverseGeneLine = ""
	# firstSnpLine = ""
	# secondSnpLine = ""
	numSnpLines = 8
	snpLines = []
	for _ in range(numSnpLines):
		snpLines.append("")
	index = 0
	estimatedSNPPosition = 0
	# indexOfNextSnpLocation = 0
	for nuc in alignedGenomes["genomeWithoutAnySnps"]:
		# print(indexThatMatchesOriginalSeq)
		
		if indexThatMatchesOriginalSeq == refGenes[nextForwardGeneIndex].startPos :#or indexThatMatchesOriginalSeq == refGenes[nextGeneIndex + 1].startPos:
			if refGenes[nextForwardGeneIndex].name == "cpdB":
				cpdBStart = index
			forwardGeneLine = format(forwardGeneLine, refGenes[nextForwardGeneIndex].name + " start", index)
		elif indexThatMatchesOriginalSeq == refGenes[nextForwardGeneIndex].stopPos:
			if refGenes[nextForwardGeneIndex].name == "cpdB":
				cpdBEnd = index
				metadata = readMetaDataAsDict(metadataPath)
				blacklist = {"genomeWithoutAnySnps", "scaffold_k-12.fasta.vcf"}
				allgenomesWithoutK12M12 = {key: value for key, value
				                           in alignedGenomes.items() if key not in blacklist}
				# with open("cpdBsFromSnpAlign", "w") as outFile:
				# 	for name, seq in allgenomesWithoutK12M12.items():
				# 		outFile.write(">"+name.split()[-1] + "_" +  metadata[name][1] + "\n")
				# 		outFile.write(seq[cpdBStart - 500:cpdBEnd+500] + "\n")
				# exit(0)
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
				                                            "pathogenic nuc is", significantSnps[indexOfNextSnpLocation][indexOfPathogenGroup] + ",",
				                                            "commensal nuc is", significantSnps[indexOfNextSnpLocation][indexOfCommensalGroup] + ",",
				                                            significantSnps[indexOfNextSnpLocation][-2]]
				snpLines[estimatedSNPPosition % len(snpLines)] = format(snpLines[estimatedSNPPosition % len(snpLines)], " ".join(shortenedSnp), index)
				# if estimatedSNPPosition % 2 == 0:
				# 	firstSnpLine = format(firstSnpLine, " ".join(shortenedSnp), index)
				# else:
				# 	secondSnpLine = format(secondSnpLine, " ".join(shortenedSnp), index)
				# indexOfNextSnpLocation += 1
				estimatedSNPPosition += 1
		# except IndexError:
		# 	pass
		if nuc != "-":
			indexThatMatchesOriginalSeq += 1
		index += 1
	
	metadata = readMetaDataAsDict(metadataPath)

	shift = ""#len("commensal") * " "
	blacklist = {"genomeWithoutAnySnps","scaffold_k-12.fasta.vcf"}
	allgenomesWithoutK12M12 = [value for key, value in alignedGenomes.items() if key not in blacklist]
	
	geneLines = [forwardGeneLine,reverseGeneLine]
	# snpLines = [firstSnpLine, secondSnpLine]
	return geneLines, alignedGenomes, snpLines
	
def writeMappedSnps():
	namePrefix = "M12RefGenome"  # "M12RefGenome"
	if namePrefix == "K12RefGenome":
		refGenomePath = DATA_DIR + "refGenomes/k-12.gbff"
	elif namePrefix == "M12RefGenome":
		refGenomePath = DATA_DIR + "refGenomes/M12.gbk"
	snpAlignPath = DATA_DIR + "RedoingEverything/SnpAlign/" + namePrefix + ".afa"
	snpIndexPath = DATA_DIR + "RedoingEverything/SnpAlign/" + namePrefix + "Indexes.txt"
	sigSnpsPath = DATA_DIR + "RedoingEverything/snpsSortedBySignificanceWithGenesContainingThem" + namePrefix + "Pathogenicity.tsv"
	if namePrefix == "K12RefGenome":
		outFilePath = "k12RefGenomeMappedSNPs.txt"
	elif namePrefix == "M12RefGenome":
		outFilePath = "M12RefGenomeMappedSNPsAllGenomes.txt"
	
	metadataPath = DATA_DIR + "metaDataForMetaCatsWithExtraMastitis.tsv"
	
	print(sigSnpsPath)
	alignedSnps = readInFastaAsDict(snpAlignPath)
	indexes = readIndexes(snpIndexPath)
	
	significantSnps = []
	with open(sigSnpsPath) as file:
		isFirstLine = True
		for line in file:
			if isFirstLine:
				isFirstLine = False
				continue
			lineData = line.split("\t")[:12]
			pos = lineData[11]
			
			significantSnps.append(lineData)
	
	significantSnps.sort(key=lambda x: (int(x[11]) + int(x[2])))
	print(significantSnps[1])
	geneLines, alignedGenomes, snpLines = mapsSnps(refGenomePath=refGenomePath, snpIndexes=indexes, alignedSnps=alignedSnps,
	                   significantSnps=significantSnps)
	
	metadata = readMetaDataAsDict(metadataPath)
	
	patternToShortenGenomes = "scaffold_|.fasta|.vcf|.result"
	blacklist = {"genomeWithoutAnySnps", "scaffold_k-12.fasta.vcf"}
	allgenomesWithoutK12M12 = {key: value for key, value
	                           in alignedGenomes.items() if key not in blacklist}
	longestGenomeNameLen = max([len(re.sub(patternToShortenGenomes, "",name) + metadata[name][1]) for name in allgenomesWithoutK12M12.keys()])
	shift = longestGenomeNameLen * " "
	
	
	with open(outFilePath, "w") as file:
		nCharsToWrite = 300
		for x in range(0, len(alignedGenomes["genomeWithoutAnySnps"]) - nCharsToWrite, nCharsToWrite):
			for geneLine in geneLines:
				file.write(shift + geneLine[x:x + nCharsToWrite] + "\n")
			for snpLine in snpLines:
				file.write(shift + snpLine[x:x + nCharsToWrite] + "\n")
			file.write(("M12:" + " " * (longestGenomeNameLen - len("M12:"))) + alignedGenomes["genomeWithoutAnySnps"][x:x + nCharsToWrite] + "\n")
			file.write(("K12:" + " " * (longestGenomeNameLen - len("K12:"))) + alignedGenomes["scaffold_k-12.fasta.vcf"][x:x + nCharsToWrite] + "\n")
			for genomeName, seq in allgenomesWithoutK12M12.items():
				ogGenomeName = deepcopy(genomeName)
				genomeName = re.sub(patternToShortenGenomes, "",genomeName)
				genomeShift = genomeName + metadata[ogGenomeName][1] + " " * (longestGenomeNameLen - len(genomeName + metadata[ogGenomeName][1]))
				file.write(genomeShift + seq[x:x + nCharsToWrite] + "\n")
			file.write("#" * (nCharsToWrite + len(shift)) + "\n")
			
	# file.write(dataline + "\n")
	# for dataline in outData:
	# file.write(dataline[x:] + "\n")
	print("wrote out data to:", outFilePath)

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






