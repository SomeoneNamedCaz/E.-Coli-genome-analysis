import sys
sys.path.insert(1, '/Users/cazcullimore/Documents/ericksonLabCode/')
from alignVcfSnps import *
from megaCatsPythonVersion import *
from parsingMegaCatsResults import *
from reconstructNormalAlignment import *
import unittest
sys.path.insert(1, '/Users/cazcullimore/Documents/ericksonLabCode/secondaryPythonScripts')
try:
    from secondaryPythonScripts.functions import *
except ImportError:
    from functions import *

# this doesn't work because it the gbk file isn't aligned
class MyTestCase(unittest.TestCase):
    
    def testSingleAlignment(self):
        
        pathToRefGenomeFasta = "/Users/cazcullimore/Documents/ericksonLabCode/refGenomes/k-12.fasta"
        pathToRefGenomeGb = "/Users/cazcullimore/Documents/ericksonLabCode/refGenomes/k-12.gbff"
        namePrefix = "singleAlignment"  # the name to give to start all the files created
        # vcfFilePath = glob("/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/vcfsFromScaffolds/*.vcf")[0]
        genomeFastaFile = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/vcfsFromScaffolds/scaffold_1465_SS_220.fasta.vcf"
        vcfFilePath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/vcfsFromScaffolds/scaffold_1465_SS_220.fasta.vcf"
        scaffoldGb = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/annotatedNormalAlignFiles/scaffold_1465_SS_220.fasta.gbk"
        
        annotatedNormalAlignPath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/annotatedNormalAlignFiles/" + namePrefix+ ".gbk"
        normalAlignFastaPath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/annotatedNormalAlignFiles/" + namePrefix+ ".fasta"
        
        snpAlignPath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/snpAlign/" + namePrefix + ".afa"
        snpIndexPath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/snpAlign/" + namePrefix + "Indexes.txt"
    
        print("annotatedNormalAlignPath",annotatedNormalAlignPath)
        print("vcfFilePath", vcfFilePath)
        print("snpAlignPath",snpAlignPath)
    
        reAlignSnps = False
        redoNormalAligment = False
        reAnnotate = False
        typesOfSnpsToSkipDuringAlignment = ["INSERT", "DELETE"]
        if reAlignSnps:
            alignVcfSnps(vcfFilePath, outFilePath=snpAlignPath, thingsToSkip=typesOfSnpsToSkipDuringAlignment, numSnpsRequired=1)
        
        if redoNormalAligment:
            reconstructNormalAlignment(snpGenomePath=snpAlignPath, snpIndexesPath=snpIndexPath,
                                       refGenomePath=pathToRefGenomeGb, outputPath=normalAlignFastaPath)
        if reAnnotate:
            # need fasta of
            os.system("cd " + "/".join(normalAlignFastaPath.split("/")[:-1]) + "; for fileName in *.fasta; do conda run -n prokkaEnv prokka $fileName --force --centre X --compliant; done; cp */*.gbk " + annotatedNormalAlignPath)
        normalAlignGenes = getGenesOnContigsByPosition(annotatedNormalAlignPath, getContigs(annotatedNormalAlignPath))
        
        assemblyGenes = getGenesOnContigs(scaffoldGb, getContigs(scaffoldGb))
        alignedSeq = readInFastaAsList(snpAlignPath)[1]
        lineOfIndexfile = -1
        possibleScafIndexes = list(range(-4, 4))
        possibleMyAlignIndexes = list(range(-4, 4))
        with open(snpIndexPath) as indexFile:
            for index in indexFile:
                lineOfIndexfile += 1
                index = index.strip()
                snpPos = int(float(index))
                print("index", index)
                for gene in normalAlignGenes:
                    if gene.startPos <= snpPos and gene.stopPos >= snpPos:
                        print(gene.name + "\nnorm\t",end="")
                        printNearbyNucs(gene.sequence, snpPos - gene.startPos, 10)
                       
                        try:
                            print("scaff\t", end="")
                            printNearbyNucs(assemblyGenes[gene.name][1], snpPos - gene.startPos, 10)
                        except KeyError:
                            pass
                        
                        print("myalign\t", end="")
                        printNearbyNucs(alignedSeq, lineOfIndexfile, 10)
        
def printNearbyNucs(seq, pos, numFlankingNucs):
    print(seq[pos - numFlankingNucs + 1:pos],
          seq[pos],
          seq[pos + 1: pos + numFlankingNucs])
if __name__ == '__main__':
	unittest.main()
