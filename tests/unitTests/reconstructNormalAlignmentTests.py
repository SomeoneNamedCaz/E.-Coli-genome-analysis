import sys
sys.path.insert(1, '/Users/cazcullimore/dev/ericksonLabCode/')
from scr.snpalign.reconstructNormalAlignment import *
import unittest
sys.path.insert(1, '/scr/snpalign/secondaryPythonScripts')
from scr.snpalign.functions import *

class MyTestCase(unittest.TestCase):
	def testAddRefNucs(self):
		# testFastaPath = TEST_DATA_DIR + "fakeFastas/test2.fasta"
		# refFastaPath = TEST_DATA_DIR + "fakeFastas/testRef.fasta"
		self.assertEqual(addRefNucs("ACTG", [1,2],"123456789"),"1AC456789")
		self.assertEqual(addRefNucs("ACTG", [1, 2,6,9], "0123456789"), "0AC345T78G")
		self.assertEqual(addRefNucs("ACTG", [1.0, 1.001,6,9],"0123456789"), "0AC2345T78G")
		self.assertEqual(addRefNucs("A---", [1.0,2.0,3.0,4.0], "0123456789"), "0A---56789")
		self.assertEqual(addRefNucs("AC--", [0.0,0.001,2,3], "0123456789"), "AC1--456789")

	def testReconstructNormalAlignment(self):
		self.maxDiff = None
		testFastaPath = TEST_DATA_DIR + "fakeFastas/test2.fasta"
		refFastaPath = TEST_DATA_DIR + "fakeFastas/testRef.fasta"
		snpAlignText = "A----T"
		snpGenomePath = "tempSnpAlign.afa"
		snpIndexPath = "tempSnpAlignIndexes.txt"
		normalAlignPath = "tempNormalAlignPath.afa"
		with open(snpGenomePath, "w") as tempFile:
			tempFile.write(snpAlignText + "\n")
		with open(snpIndexPath, "w") as tempFile:
			tempFile.write("309.0\n310.0\n311.0\n312.0\n313.0\n699.0\n")
		reconstructNormalAlignment(snpGenomePath=snpGenomePath, snpIndexesPath=snpIndexPath, refGenomePath=refFastaPath, outputPath=normalAlignPath)
		
		originalSeq = readInFastaAsList(testFastaPath)[1]
		reconstructedSeq = re.sub("-","",readInFastaAsList(normalAlignPath)[1])
		self.assertEqual(reconstructedSeq,originalSeq)
		
		os.remove(snpGenomePath)
		os.remove(snpIndexPath)
		os.remove(normalAlignPath)
		

if __name__ == '__main__':
	unittest.main()
