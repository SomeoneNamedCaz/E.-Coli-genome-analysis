import unittest
from scr.snpalign.alignVcfSnps import *

class MyTestCase(unittest.TestCase):


    def testForWeirdThingsInVCFFiles(self):
        fastaPath = DATA_DIR + "tests/unitTests/k-12.fasta"
        gbPath = DATA_DIR + "tests/unitTests/k-12.gbff"
        vcfPaths = DATA_DIR + "k-12RefGenomeAnalysisFiles/AllAssemblies/APEC_assemblies/ragtagOutputs/longestScaffoldFiles/carefulGsAlign/*.vcf"

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
    def testByParseVCFs(self):
        refFastaPath = DATA_DIR + "tests/unitTests/k-12.fasta"
        refGenomeSeq = readInFastaAsList(refFastaPath)[1]
        vcfs = glob(DATA_DIR + "k-12RefGenomeAnalysisFiles/AllAssemblies/APEC_assemblies/ragtagOutputs/longestScaffoldFiles/carefulGsAlign/*.vcf")
        snps = readInSnps(vcfs, refGenomeSeq, ignoreRefSeq=True)
        snps.sort(key=lambda a: a[1])
        lastPos = -1
        hasNs = False
        hasNoNs = False
        for snp in snps:
            snpType = snp[4]
            snpPos = snp[1]
            refNucs = snp[2]
            altNucs = snp[3]
            if snpPos == lastPos:
                if refNucs.count("N") != 0 or altNucs.count("N") != 0:
                    hasNs = True
                else:
                    hasNoNs = True
                if hasNoNs and hasNs:
                    print("has Ns and Doesn't have Ns")
            lastPos = snpPos
        readInSnps(vcfs, refGenomeSeq, ignoreRefSeq=True)

if __name__ == '__main__':
    unittest.main()