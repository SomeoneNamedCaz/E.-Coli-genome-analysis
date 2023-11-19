from functions import *

def makeMapFile(metaDataFilePath, outFilePath):
	with open(metaDataFilePath) as metaDataFile, open(outFilePath, "w") as outFile:
		firstLine = True
		for line in metaDataFile:
			line = line.strip()
			if firstLine or line == "":
				firstLine = False
				continue
			cols = line.split("\t")
			lineData = cols[0] + ''.join(cols[1:])
			outFile.write("\t".join(lineData) + "\n")
		
		
def makePedFile(snpGenomeFilePath, outFilePath):
	snpGenomeData = readInFastaAsDict(snpGenomeFilePath)
	with open(outFilePath, "w") as outFile:
		# FamilyID IndividualID FatherID motherID sex phenotype SNP1 SNP2
		for name, seq in snpGenomeData.items():
			lineData = [name,name,"0",'0','0','0'] + list(seq)
			outFile.write("\t".join(lineData) + "\n")
			