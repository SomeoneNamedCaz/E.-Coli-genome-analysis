import unittest
import sys
sys.path.insert(1, '/Users/cazcullimore/Documents/ericksonLabCode/secondaryPythonScripts')
from tests.functions import *

class MyTestCase(unittest.TestCase):

    def testReadInFasta(self):
        fastaPath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/unitTests/k-12.fasta"
        fileData = readInFastaAsList(fastaPath)
        print(len(fileData))
        with open(fastaPath) as file:
            entryData = ""
            for line in file:
                if line[0] == ">":
                    self.assertEqual(fileData[0], line.strip())
                    fileData.pop(0)
                    if entryData != "":
                        self.assertEqual(fileData[0], entryData)
                        fileData.pop(0)
                        entryData = ""
                else:
                    entryData += line.strip()
        if entryData != "":
            self.assertEqual(fileData[0], entryData)
            fileData.pop(0)
        self.assertEqual(fileData, [])# should be empty


    def testGetAssemblyContigs(self):
        gbPath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/unitTests/k-12.gbff"
        fastaPath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/unitTests/k-12.fasta"
        index = 0
        gbContigs = getContigs(gbPath)
        fastaContigs = readInFastaAsList(fastaPath)[1:][::2]
        self.assertEqual(len(gbContigs), len(fastaContigs))
        for gbContig in gbContigs:
            foundMatch = False
            for fastaContig in fastaContigs:
                if gbContig == fastaContig:
                    foundMatch = True
            if not foundMatch:
                raise Exception()

    def testGetFullGenomeContig(self):
        gbPath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/unitTests/k-12.gbff"
        fastaPath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/unitTests/k-12.fasta"
        index = 0
        self.assertEqual(getContigs(gbPath), readInFastaAsList(fastaPath)[1:][::2])

if __name__ == '__main__':
    unittest.main()
