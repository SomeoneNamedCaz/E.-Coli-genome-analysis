import unittest
import sys
from parsingMegaCatsResults import getMutationType

from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
# from Bio
from Bio.Entrez import *
import re
sys.path.insert(1, '/Users/cazcullimore/dev/ericksonLabCode/secondaryPythonScripts')
from secondaryPythonScripts.functions import *

class MyTestCase(unittest.TestCase):

    def testReadIn1ContigFasta(self):
        fastaPath = TEST_DATA_DIR + "unitTests/k-12.fasta"
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

    def testReadInAssemblyFasta(self):
        fastaPath = TEST_DATA_DIR + "tenAssembliesFromEachCategory/1465_SS_220.fasta"
        fileData = readInFastaAsList(fastaPath)
        print(len(fileData))
        with open(fastaPath) as file:
            entryData = ""
            for line in file:
                if line[0] == ">":
                    if entryData != "":
                        self.assertEqual(fileData.pop(0), entryData)
                        entryData = ""
                    self.assertEqual(fileData.pop(0), line.strip())
                else:
                    entryData += line.strip()
        if entryData != "":
            self.assertEqual(fileData[0], entryData)
            fileData.pop(0)
        self.assertEqual(fileData, [])# should be empty


    def testGetAssemblyContigs(self):
        gbPath = TEST_DATA_DIR + "unitTests/k-12.gbff"
        fastaPath = TEST_DATA_DIR + "unitTests/k-12.fasta"
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
        gbPath = DATA_DIR + "unitTests/k-12.gbff"
        fastaPath = DATA_DIR + "unitTests/k-12.fasta"
        index = 0
        self.assertEqual(getContigs(gbPath), readInFastaAsList(fastaPath)[1:][::2])

    def testGetGenes(self): # test
        # this looks at if the genes are similar to what they should look like
        testGbPath = TEST_DATA_DIR + "scaffold_1465_SS_220.fasta.gbk"
        
        genes = getGenesOnContigs(testGbPath, getContigs(testGbPath))
        geneAASeq = ""
        inFeatures = False
        inCDS = False
        inTranslation = False
        with open(testGbPath) as testGb:
            geneName = "no name"
            for line in testGb:
                line = line.strip()
                cols = line.split()
                if len(cols) == 0:
                    continue
                if cols[0] == "FEATURES":
                    inFeatures = True
                elif line == 'ORIGIN':
                    inCDS = False
                elif inFeatures:  # if in annotations
                    if cols[0] == "CDS":
                        inCDS = True
                        inTranslation = False
                    elif cols[0] == "gene":
                        inCDS = False
                        inTranslation = False
                        try:
                            print("\n"+genes[geneName].sequence)
                            print(geneName, geneAASeq)
                            print(translate(genes[geneName].sequence))
                            self.assertEqual(geneAASeq[1:], translate(genes[geneName].sequence)[1:])
                            # exit(0)
                            geneName = "no name"
                        except KeyError:
                            pass
                    elif inCDS:
                        if line[:7] == '/gene="':
                            geneName = re.sub('/gene="', "", line)[:-1]
                        if "/translation=\"" in line:
                            inTranslation = True
                            geneAASeq = re.sub("(/translation=\")|\"", "", line)
                        elif inTranslation:
                            geneAASeq += re.sub('"', "",line)
                            if '"' in line:
                                inTranslation = False
        

    # def testTwoFilesHaveSameInfo(self):
    #     filesContainTheSameInformation()
    
    def testGetSnpInfo(self):
        print(getSnpInfo("pathogen(24_A,160_T)|commensal(2_A,86_T)"))
        self.assertEqual(("T","A",0), getSnpInfo("pathogen(24_A,160_T)|commensal(2_A,86_T)"))
        self.assertEqual(("T", "A", 0), getSnpInfo("pathogen(24_A,160_T)|commensal(2_A,86_T)"))
    
    def testGetMutationTypeSilent(self):
        testGene = Gene(0,10,"ACTGACTGACT","testGene", "testProduct", isForward=True)
        testSnp = SNP(2,"T","C",{},0,"group1", "testGene")
        indexesOfFrameshifts = set()
        self.assertEqual(SNP.mutationType.silent,getMutationType(testGene,testSnp,indexesOfFrameshifts))
    def testGetMutationTypeMisSense(self):
        testGene = Gene(0, 10, "ACTGACTGACT", "testGene", "testProduct", isForward=True)
        testSnp = SNP(1, "C", "A", {}, 0, "group1", "testGene")
        indexesOfFrameshifts = set()
        self.assertEqual(SNP.mutationType.misSense, getMutationType(testGene, testSnp, indexesOfFrameshifts))
        
    def testGetMutationTypeSilentBackwards(self):
        testGene = Gene(0, 10, "ACTGACTGACT", "testGene", "testProduct", isForward=False)
        testSnp = SNP(0, "A", "G", {}, 0, "group1", "testGene")
        indexesOfFrameshifts = set()
        self.assertEqual(SNP.mutationType.silent, getMutationType(testGene, testSnp, indexesOfFrameshifts))
        
    def testGetMutationTypeMisSenseBackwards(self):
        testGene = Gene(0, 10, "ACTGACTGACT", "testGene", "testProduct", isForward=False)
        testSnp = SNP(1, "G", "C", {}, 0, "group1", "testGene")
        indexesOfFrameshifts = set()
        self.assertEqual(SNP.mutationType.misSense, getMutationType(testGene, testSnp, indexesOfFrameshifts))
        
    def testGetMutationTypeFrameShiftBackwards(self):
        testGene = Gene(0, 10, "ACTGACTGACT", "testGene", "testProduct", isForward=False)
        testSnp = SNP(1, "G", "-", {}, 0, "group1", "testGene")
        indexesOfFrameshifts = {1}
        self.assertEqual(SNP.mutationType.frameShift, getMutationType(testGene, testSnp, indexesOfFrameshifts))

if __name__ == '__main__':
    unittest.main()
