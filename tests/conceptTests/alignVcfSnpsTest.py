import re
import unittest
import sys
import os
sys.path.insert(1, '/Users/cazcullimore/dev/ericksonLabCode/')
try: # so pycharm recognizes it
    from secondaryPythonScripts.functions import *
except ImportError:
    pass
from reconstructNormalAlignment import *
from alignVcfSnps import alignVcfSnps

class MyTestCase(unittest.TestCase):
    def testSingleFile(self):
        filePath = TEST_DATA_DIR + "fake vcfs/test3.vcf"
        outPath = TEST_DATA_DIR + "snpAlign/singleFileSnps.afa"
        expectedIndexes = ([308.0] + InsertIndexes(309,5) +
                           InsertIndexes(318,2) +
                           DeleteIndexes(319, 15) +
                           InsertIndexes(599, 17))
        indexFilePath = outPath[:-4] + "Indexes.txt"
        alignVcfSnps(filePath, outFilePath=outPath, numSnpsRequired=1, ignoreRefSeq=True,
                     includeGenomeWithoutAnySnps=True)
        self.assertEqual(readInFastaAsDict(outPath)[filePath.split("/")[-1]], "GATGGTACT--------------ATCCCCCCCCCCCCCCC")
        
        # test indexes
        indexes = readInIndexes(indexFilePath)
        self.assertListEqual(indexes, expectedIndexes)


    def testInsertOverlap(self):
        filePath = TEST_DATA_DIR + "fake vcfs/insertOverlap*.vcf"
        outPath = TEST_DATA_DIR + "snpAlign/insertOverlap.afa"
        indexFilePath = outPath[:-4] + "Indexes.txt"
        alignedSnps = alignVcfSnps(filePath, outFilePath=outPath, numSnpsRequired=1, ignoreRefSeq=True)
        self.assertEqual("ACCGTCGT",''.join(alignedSnps["insertOverlap2.vcf"]))
        self.assertEqual("ACCA----", ''.join(alignedSnps["insertOverlap1.vcf"]))

    def testInsertWithinDeletion(self):
        filePath = TEST_DATA_DIR + "fake vcfs/insertWithinDelete*.vcf"
        outPath = TEST_DATA_DIR + "snpAlign/insertWithinDelete.afa"
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
        filePath = TEST_DATA_DIR + "fake vcfs/test*.vcf"
        # filePath = TEST_DATA_DIR + "fakeScaffolds/gsAlignOutputs/test*.vcf"
        outPath = TEST_DATA_DIR + "snpAlign/lotsOfFIles.afa"
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
        gsAlignVCFs = TEST_DATA_DIR + "fake vcfs/test*.vcf"
        # gsAlignVCFs = DATA_DIR + "AllAssembliesInOneFolder/ragtagOutputs/longestScaffoldFiles/gsAlignOutputs/*.vcf"#TEST_DATA_DIR + "vcfsFromScaffolds/*.vcf"
        alignedSnpPath = TEST_DATA_DIR + "snpAlign/allSnps.afa"
        reRunAlignment = True

        refSeq = readInFastaAsList(TEST_DATA_DIR + "fakeScaffolds/testRef.fasta")[1]
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
                        print(refSeq[iLine-2:iLine+2],[genome[fileLineIndex]for genome in alignedSnps])
                        print("didn't find ref nuc", iLine)
                        print("__")
                        raise Exception("didn't find ref nuc", iLine,refSeq[iLine-2:iLine+2],[genome[fileLineIndex]for genome in alignedSnps])
        
    def testReconstructFakeGenomes(self):
        
        pathToRefGenomeFasta = TEST_DATA_DIR + "fakeScaffolds/testRef.fasta"
        # pathToRefGenomeGb = DATA_DIR + "refGenomes/k-12.gbff"
        namePrefix = "reconstructFakeGenomes"  # the name to give to start all the files created
        genomeFastaFiles = TEST_DATA_DIR + "fakeScaffolds/*.fasta"
        vcfFilePath = TEST_DATA_DIR + "fakeScaffolds/gsAlignOutputs/*.vcf"
        
        normalAlignFastaPath = TEST_DATA_DIR + "annotatedNormalAlignFiles/" + namePrefix + ".fasta"
        
        snpAlignPath = TEST_DATA_DIR + "snpAlign/" + namePrefix + ".afa"
        snpIndexPath = TEST_DATA_DIR + "snpAlign/" + namePrefix + "Indexes.txt"
        
        # print("annotatedNormalAlignPath", annotatedNormalAlignPath)
        print("vcfFilePath", vcfFilePath)
        print("snpAlignPath", snpAlignPath)
        
        reAlignSnps = True
        redoNormalAligment = True
        typesOfSnpsToSkipDuringAlignment = []  # ["INSERT", "DELETE"]
        if reAlignSnps:
            alignVcfSnps(vcfFilePath, outFilePath=snpAlignPath, thingsToSkip=typesOfSnpsToSkipDuringAlignment,
                         numSnpsRequired=1, refSeqPath=pathToRefGenomeFasta, ignoreRefSeq=False)
        
        if redoNormalAligment:
            reconstructNormalAlignment(snpGenomePath=snpAlignPath, snpIndexesPath=snpIndexPath,
                                       refGenomePath=pathToRefGenomeFasta, outputPath=normalAlignFastaPath)
        
        normalAlignNameToSeq = readInFastaAsDict(normalAlignFastaPath)
        self.maxDiff = None
        for file in glob(genomeFastaFiles):
            self.assertEqual(readInFastaAsList(file)[1], re.sub("-","", normalAlignNameToSeq[file.split("/")[-1]+".vcf"]))

    def testReconstructingSingleAlignement(self):
        pathToRefGenomeFasta = DATA_DIR + "refGenomes/k-12.fasta"
        # pathToRefGenomeGb = DATA_DIR + "refGenomes/k-12.gbff"
        namePrefix = "reconstructScaffold"  # the name to give to start all the files created
        genomeFastaFile = TEST_DATA_DIR + "scaffolds/longestScaffoldFiles/scaffold_1465_SS_220.fasta"
        vcfFilePath = TEST_DATA_DIR + "vcfsFromScaffolds/scaffold_1465_SS_220.fasta.vcf"
        
        normalAlignFastaPath = TEST_DATA_DIR + "annotatedNormalAlignFiles/" + namePrefix + ".fasta"
        snpAlignPath = TEST_DATA_DIR + "snpAlign/" + namePrefix + ".afa"
        snpIndexPath = TEST_DATA_DIR + "snpAlign/" + namePrefix + "Indexes.txt"
        
        
        # print("annotatedNormalAlignPath", annotatedNormalAlignPath)
        print("vcfFilePath", vcfFilePath)
        print("snpAlignPath", snpAlignPath)
        
        reAlignSnps = True
        redoNormalAligment = True
        typesOfSnpsToSkipDuringAlignment = []  # ["INSERT", "DELETE"]
        if reAlignSnps:
            alignVcfSnps(vcfFilePath, outFilePath=snpAlignPath, thingsToSkip=typesOfSnpsToSkipDuringAlignment,
                         numSnpsRequired=1, refSeqPath=pathToRefGenomeFasta, ignoreRefSeq=False)
        
        if redoNormalAligment:
            reconstructNormalAlignment(snpGenomePath=snpAlignPath, snpIndexesPath=snpIndexPath,
                                       refGenomePath=pathToRefGenomeFasta, outputPath=normalAlignFastaPath)
        
        normalAlignNameToSeq = readInFastaAsDict(normalAlignFastaPath)
        genomeName = genomeFastaFile.split("/")[-1] + ".vcf"
        originalGenome = readInFastaAsList(genomeFastaFile)[1]
        reconstructedQueryGenome = re.sub("-", "", normalAlignNameToSeq[genomeName])
        reconstructedRefGenome =  re.sub("-", "", normalAlignNameToSeq["genomeWithoutAnySnps"])
        self.assertEqual(reconstructedRefGenome,readInFastaAsList(pathToRefGenomeFasta)[1])
        # never going to pass because if it doesn't align everything it won't reconstruct everything
        self.assertEqual(originalGenome,
                         reconstructedQueryGenome)
        """maybe look for the snps that are in genes and the compare them to the annotated files but first just test that the ref seq matches some of the aligned seqs"""


    def testDeletesReconstructGenome(self):
        ref = "ACCGAGGAGTAGAGAGC"
        indexes = [2, 4,5,6]
        query = list(ref)
        indexes.reverse()
        snps = []
        for index in indexes:
            snps.append(query[index])
            query[index] = "-"
        snps.reverse()
        query = "".join(query)
        alignedSnps = {"query": "-"*len(snps), "genomeWithNoSnps": "".join(snps)}
        indexes.reverse()
        reconstructedAlignments = reconstructNormalAlignmentHelper(refSeq=ref, snpIndexes=indexes, alignedSnps=alignedSnps)
        self.assertEqual(reconstructedAlignments["query"], query)
        self.assertEqual(reconstructedAlignments["genomeWithNoSnps"], ref)
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
