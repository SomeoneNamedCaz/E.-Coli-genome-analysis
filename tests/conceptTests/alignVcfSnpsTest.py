import unittest
import sys
import os
sys.path.insert(1, '/Users/cazcullimore/Documents/ericksonLabCode/')
try: # so pycharm recognizes it
    from secondaryPythonScripts.functions import *
except ImportError:
    pass
from convertSNPAlignmentToNormalAlignment import *
from alignVcfSnps import alignVcfSnps

class MyTestCase(unittest.TestCase):
    def testSingleFile(self):
        filePath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/fake vcfs/test3.vcf"
        outPath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/snpAlign/singleFileSnps.afa"
        expectedIndexes = [308.0,309.0, 310.0,311.0, 312.0,313.0,318.0, 318.001]
        for i in range(319,319+15):
            expectedIndexes.append(float(i))
        for i in range(600,600+17):
            expectedIndexes.append(600.0 + i/1000)
        indexFilePath = outPath[:-4] + "Indexes.txt"
        alignVcfSnps(filePath, outFilePath=outPath, numSnpsRequired=1, ignoreRefSeq=True)
        self.assertEqual(readInFastaAsList(outPath)[1], "GATGGTACT--------------ATCCCCCCCCCCCCCCC")
        
        # test indexes
        indexes = readInIndexes(indexFilePath)
        self.assertListEqual(indexes, expectedIndexes)


    def testInsertOverlap(self):
        filePath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/fake vcfs/insertOverlap*.vcf"
        outPath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/snpAlign/insertOverlap.afa"
        indexFilePath = outPath[:-4] + "Indexes.txt"
        alignedSnps = alignVcfSnps(filePath, outFilePath=outPath, numSnpsRequired=1, ignoreRefSeq=True)
        self.assertEqual("ACCGTCGT",''.join(alignedSnps["insertOverlap2.vcf"]))
        self.assertEqual("ACCA----", ''.join(alignedSnps["insertOverlap1.vcf"]))

    def testInsertWithinDeletion(self):
        filePath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/fake vcfs/insertWithinDelete*.vcf"
        outPath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/snpAlign/insertWithinDelete.afa"
        indexFilePath = outPath[:-4] + "Indexes.txt"
        alignedSnps = alignVcfSnps(filePath, outFilePath=outPath, numSnpsRequired=1, ignoreRefSeq=True)
        # test alignment
        self.assertEqual("TAAAATTTTTAAAATTTTGT", ''.join(alignedSnps["insertWithinDelete2.vcf"]))
        self.assertEqual("A---------A---------", ''.join(alignedSnps["insertWithinDelete1.vcf"]))
        
        # test indexes
        indexes = readInIndexes(indexFilePath)
        expectedIndexes = [309.0,310.0,310.001,310.002,310.003,311.0,311.001,311.002,311.003,312.0,314.0,314.001,314.002,314.003,315.0,315.001,315.002,315.003,316.0,317.0,]
        self.assertListEqual(indexes, expectedIndexes)

    def testLotsOfFiles(self):
        filePath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/fake vcfs/test*.vcf"
        # filePath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/fakeFastas/gsAlignOutputs/test*.vcf"
        outPath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/snpAlign/lotsOfFIles.afa"
        indexFilePath = outPath[:-4] + "Indexes.txt"
        alignedSnps = alignVcfSnps(filePath, outFilePath=outPath, numSnpsRequired=1, ignoreRefSeq=True)
        # test alignment
        #                      308 stuff         |317                   | 318   |319           |340  |439           440| |442  |443
        self.assertEqual("A-----------------A----------------------A-------TGAAAAAAAACGAAATGGTTGA-----------------TGG-----TTGA----------------TTGTGTGTGTGTGTGTGTGTGTGTGT", ''.join(alignedSnps["test2.vcf"]))
        self.assertEqual("G---------ATGGTAAAA----------------------AC------T--------------TGGTTGA-----------------TGG-----TTGATCCCCCCCCCCCCCCCGTGTGTGTGTGTGTGTGTGTGTGTGT", ''.join(alignedSnps["test3.vcf"]))
        self.assertEqual("ACC-------A----AAAA----------------------ATGGGGGTCGAAAAAAAACGAAATGGTTGA----------------------------T----------------G------------------------T", ''.join(alignedSnps["test4.vcf"]))
        self.assertEqual("T---------A----AAAATTTTTTTTTTTTTTTTTTTTTTA-------TGAAAAAAAACGAAAT-----A-----------------TGGTTGGTTTGA----------------GTGTGTGTGTGTGTGTGTGTGTGTGC", ''.join(alignedSnps["test5.vcf"]))
        self.assertEqual("AGGGGGGGGGA----AAAA----------------------AC------TTAAAAAAAACGAAATGGTTGACCGGAGTTGGGAGCCTCTGG-----TTGA----------------GTGTGTGTGTGTGTGTGTGTGTGTGT", ''.join(alignedSnps["test6.vcf"]))
        self.assertEqual("A---------A----AAAA----------------------A-------TGAAAAAAAACGAAATGGTTGA-----------------TGG-----TTGA----------------GTGTGTGTGTGTGTGTGTGTGTGTGT", ''.join(alignedSnps["test7.vcf"]))
        self.assertEqual("A---------A----AAAA----------------------A-------TGAAAAAAAACGAAATGGTTGA-----------------TGG-----TTGA----------------GTGTGTGTGTGTGTGTGTGTGTGTGT", ''.join(alignedSnps["test8.vcf"]))
        
        # test indexes
        indexes = readInIndexes(indexFilePath)
        expectedIndexes = (InsertIndexes(308, 10)
                           + InsertIndexes(309, 5)
                           + DeleteIndexes(310, 3)
                           + InsertIndexes(317, 23)
                           + InsertIndexes(318, 8)
                           + DeleteIndexes(319, 15)
                           + DeleteIndexes(340, 6)
                           + InsertIndexes(439, 18)  #
                           + DeleteIndexes(440, 2)
                           + InsertIndexes(442, 6)
                           + DeleteIndexes(443, 3)
                           + InsertIndexes(599, 17)
                           + DeleteIndexes(699, 25)
                           + DeleteIndexes(3449, 1))
        print("expected")
        print(expectedIndexes)
        print("actual")
        print(indexes)
        self.assertListEqual(indexes, expectedIndexes)
        
    def testReferenceGenomeMatching(self):
        gsAlignVCFs = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/fake vcfs/*.vcf"
        # gsAlignVCFs = "/Users/cazcullimore/Documents/ericksonLabCode/AllAssembliesInOneFolder/ragtagOutputs/longestScaffoldFiles/gsAlignOutputs/*.vcf"#"/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/vcfsFromScaffolds/*.vcf"
        alignedSnpPath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/snpAlign/allSnps.afa"
        reRunAlignment = True

        refSeq = readInFastaAsList("/Users/cazcullimore/Documents/ericksonLabCode/refGenomes/k-12.fasta")[1]
        if not os.path.exists(alignedSnpPath) or reRunAlignment:
            alignVcfSnps(gsAlignVCFs,alignedSnpPath, numSnpsRequired=1)
        indexFilePath = alignedSnpPath[:-4] + "Indexes.txt"
        alignedSnps = readInFastaAsList(alignedSnpPath)[1:][::2]
        indexes = []
        insertIndexes = []
        with open(indexFilePath) as indexFile:
            fileLineIndex = -1
            for line in indexFile:
                fileLineIndex += 1
                line = line.strip()
                fLine = float(line)
                iLine = int(fLine)
                indexes.append(iLine)
                if fLine != iLine:
                    insertIndexes.append(fLine)
                else:
                    foundRefNuc = False
                    for snpGenome in alignedSnps:
                        if refSeq[iLine] == snpGenome[fileLineIndex]:

                            if not foundRefNuc:
                                print("foundRef",refSeq[iLine - 2:iLine + 2] + "  " + snpGenome[fileLineIndex])
                                print("did found ref nuc", iLine)
                                print("__")
                            foundRefNuc = True
                    if not foundRefNuc:
                        print(refSeq[iLine-2:iLine+2] +"  " +snpGenome[fileLineIndex])
                        print("didn't find ref nuc", iLine)
                        print("__")
        
    def testReconstructFakeGenomes(self):
        
        pathToRefGenomeFasta = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/fakeFastas/testRef.fasta"
        # pathToRefGenomeGb = "/Users/cazcullimore/Documents/ericksonLabCode/refGenomes/k-12.gbff"
        namePrefix = "reconstructFakeGenomes"  # the name to give to start all the files created
        genomeFastaFile = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/fakeFastas/*.fasta"
        vcfFilePath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/fakeFastas/gsAlignOutputs/*.vcf"
        
        normalAlignFastaPath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/annotatedNormalAlignFiles/" + namePrefix + ".fasta"
        
        snpAlignPath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/snpAlign/" + namePrefix + ".afa"
        snpIndexPath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/snpAlign/" + namePrefix + "Indexes.txt"
        
        # print("annotatedNormalAlignPath", annotatedNormalAlignPath)
        print("vcfFilePath", vcfFilePath)
        print("snpAlignPath", snpAlignPath)
        
        reAlignSnps = True
        redoNormalAligment = True
        typesOfSnpsToSkipDuringAlignment = []  # ["INSERT", "DELETE"]
        if reAlignSnps:
            alignVcfSnps(vcfFilePath, outFilePath=snpAlignPath, thingsToSkip=typesOfSnpsToSkipDuringAlignment,
                         numSnpsRequired=1,refSeqPath=pathToRefGenomeFasta, ignoreRefSeq=False)
        
        if redoNormalAligment:
            makeNormalAlignment(snpGenomePath=snpAlignPath, snpIndexesPath=snpIndexPath,
                                refGenomePath=pathToRefGenomeFasta, outputPath=normalAlignFastaPath)
        
        

        """maybe look for the snps that are in genes and the compare them to the annotated files but first just test that the ref seq matches some of the aligned seqs"""

def readInIndexes(fileName):
    indexes = []
    with open(fileName) as file:
        for line in file:
            line = line.strip()
            if line == "":
                continue
            indexes.append(float(line))
    return indexes
def InsertIndexes(startIndex, numIndexes):
    indexes = []
    for i in range(0, numIndexes):
        indexes.append(startIndex + i / 1000)
    return indexes
def addConsequtiveInserts(startIndex, numIndexes): ## for deletes
    indexes = []
    for i in range(startIndex, startIndex + numIndexes):
        indexes.append(float(i))
    return indexes
def DeleteIndexes(startIndex, numIndexes):
    return addConsequtiveInserts(startIndex, numIndexes)
if __name__ == '__main__':
    unittest.main()
