import sys
sys.path.insert(1, '/Users/cazcullimore/dev/ericksonLabCode/')
sys.path.insert(1, '/Users/cazcullimore/dev/ericksonLabCode/secondaryPythonScripts')
import unittest
from secondaryPythonScripts.mapSNPsToGenome import *
from secondaryPythonScripts.functions import *


class MyTestCase(unittest.TestCase):
	def test_something(self):
		refGenomePath = DATA_DIR +  "refGenomes/M12.gbk"
		namePrefix = "M12RefGenome" #"K12RefGenome"
		snpAlignPath = DATA_DIR + "RedoingEverything/SnpAlign/" + namePrefix +".afa"
		snpIndexPath = DATA_DIR + "RedoingEverything/SnpAlign/" + namePrefix +"Indexes.txt"
		sigSnpsPath = DATA_DIR + "RedoingEverything/snpsSortedBySignificanceWithGenesContainingThemM12RefGenomePathogenicity.tsv"
		alignedSnps = readInFastaAsDict(snpAlignPath)
		indexes = readIndexes(snpIndexPath)
		
		significantSnps = []
		with open(sigSnpsPath) as file:
			for line in file:
				significantSnps.append(line.split("\t")[:12])
		print(significantSnps[1])
		annotationLine, alignedGenomes, snpLine = mapsSnps(refGenomePath=refGenomePath, snpIndexes=indexes, alignedSnps=alignedSnps, significantSnps=significantSnps)
		with open("out.txt", "w") as file:
			nCharsToWrite = 100
			for x in range(0,len(annotationLine) - 100,100):
				file.write(annotationLine[x:x+100] + "\n")
				file.write(alignedGenomes["genomeWithoutAnySnps"][x:x+100] + "\n")
				file.write(snpLine[x:x+100] + "\n")

if __name__ == '__main__':
	unittest.main()
