import sys
sys.path.insert(1, '/Users/cazcullimore/dev/ericksonLabCode/')
sys.path.insert(1, '/Users/cazcullimore/dev/ericksonLabCode/secondaryPythonScripts')
import unittest
from secondaryPythonScripts.mapSNPsToGenome import *
from secondaryPythonScripts.functions import *


class MyTestCase(unittest.TestCase):
	def test_something(self):
		
		namePrefix = "K12RefGenome" # "M12RefGenome"
		if namePrefix == "K12RefGenome":
			refGenomePath = DATA_DIR + "refGenomes/k-12.gbff"
		elif namePrefix == "M12RefGenome":
			refGenomePath = DATA_DIR + "refGenomes/M12.gbk"
		snpAlignPath = DATA_DIR + "RedoingEverything/SnpAlign/" + namePrefix +".afa"
		snpIndexPath = DATA_DIR + "RedoingEverything/SnpAlign/" + namePrefix +"Indexes.txt"
		sigSnpsPath = DATA_DIR + "RedoingEverything/snpsSortedBySignificanceWithGenesContainingThem" + namePrefix + "Pathogenicity.tsv"
		if namePrefix == "K12RefGenome":
			outFilePath = "k12RefGenomeMappedSNPs.txt"
		elif namePrefix == "M12RefGenome":
			outFilePath = "M12RefGenomeMappedSNPs.txt"
		
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
				
		
		significantSnps.sort(key=lambda x:(int(x[11])+ int(x[2])))
		print(significantSnps[1])
		outData = mapsSnps(refGenomePath=refGenomePath, snpIndexes=indexes, alignedSnps=alignedSnps, significantSnps=significantSnps)
		print(outData)
		with open(outFilePath, "w") as file:
			nCharsToWrite = 100
			for x in range(0,len(outData[2]) - nCharsToWrite,nCharsToWrite):
				for dataline in outData:
					file.write(dataline[x:x+nCharsToWrite] + "\n")
					# file.write(dataline + "\n")
		print("wrote out data to:", outFilePath)

if __name__ == '__main__':
	unittest.main()
