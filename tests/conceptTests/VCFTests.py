import unittest
from secondaryPythonScripts.functions import *

class MyTestCase(unittest.TestCase):


    def testForWeirdThingsInVCFFiles(self):
        fastaPath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/unitTests/k-12.fasta"
        gbPath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/unitTests/k-12.gbff"
        vcfPaths = "/Users/cazcullimore/Documents/ericksonLabCode/k-12RefGenomeAnalysisFiles/AllAssemblies/APEC_assemblies/ragtagOutputs/longestScaffoldFiles/carefulGsAlign/*.vcf"

        indexOfRefColumn = 3
        indexOfPositionColumn = 1
        genomeDNA = readInFastaAsList(fastaPath)[1] # it's on the second line
        for filePath in glob(vcfPaths):
            with open(filePath) as file:
                for line in file:
                    line = line.strip()
                    if line == "" or line[0] == "#":
                        continue
                    cols = line.split("\t")
                    # print(indexOfPositionColumn)
                    snpPositionIndex = int(cols[indexOfPositionColumn]) - 1
                    refNucs = cols[indexOfRefColumn]
                    if genomeDNA[snpPositionIndex:len(refNucs) + snpPositionIndex] != refNucs:
                        print(filePath,"is lying")
                        print("vcf ref nucs", refNucs,"at", snpPositionIndex, "aren't", genomeDNA[snpPositionIndex:len(refNucs) + snpPositionIndex])



if __name__ == '__main__':
    unittest.main()