from secondaryPythonScripts.functions import *

def makeCovFile(metaDataFilePath, outFilePath, snpGenomeFilePath=None):
	if snpGenomeFilePath is not None:
		snpGenomeData = readInFastaAsDict(snpGenomeFilePath)
	with open(metaDataFilePath) as metaDataFile, open(outFilePath, "w") as outFile:
		outFile.write("#FID\tIID\tPAT\tMAT\tSEX\tPHENOTYPE\n")
		firstLine = True
		for line in metaDataFile:
			line = line.strip()
			if firstLine or line == "":
				firstLine = False
				continue
			cols = line.split("\t")
			name = cols[0]
			nameWithoutExtension = name.split(".")[0]
			lineData = [nameWithoutExtension,nameWithoutExtension, cols[1]]
			if snpGenomeFilePath is not None:
				if name in snpGenomeData.keys():
					outFile.write("\t".join(lineData) + "\n")
			else:
				outFile.write("\t".join(lineData) + "\n")
		
def makeMapFile(indexPath, outFilePath):
	with open(indexPath) as indexFile, open(outFilePath, "w") as outFile:
		i = 0
		for line in indexFile:
			line = line.strip()
			if line == "":
				continue
			lineData = ["0","snp" + str(i), "0", str(i)]#str(int(float(line) * 1000))]
			outFile.write("\t".join(lineData) + "\n")
			i += 1
def makePedFile(snpGenomeFilePath, metaDataFilePath, outFilePath):
	snpGenomeData = readInFastaAsDict(snpGenomeFilePath)
	seqToMetadata = readMetaDataAsDict(metaDataFilePath)
	with open(outFilePath, "w") as outFile:
		# FamilyID IndividualID FatherID motherID sex phenotype SNP1 SNP2
		for name, seq in snpGenomeData.items():
			doubleSeq = []
			for nuc in seq:
				doubleSeq.append(nuc)
				doubleSeq.append(nuc)
			nameWithoutExtension = name.split(".")[0]
			lineData = [nameWithoutExtension,nameWithoutExtension,"0",'0','0',"0"] + doubleSeq
			outFile.write("\t".join(lineData) + "\n")
			