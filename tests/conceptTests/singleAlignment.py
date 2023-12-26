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
        pass
        pathToRefGenomeFasta = "/Users/cazcullimore/Documents/ericksonLabCode/refGenomes/k-12.fasta"
        pathToRefGenomeGb = "/Users/cazcullimore/Documents/ericksonLabCode/refGenomes/k-12.gbff"
        namePrefix = "shortSingleAlignment"  # the name to give to start all the files created
        # vcfFilePath = glob("/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/vcfsFromScaffolds/*.vcf")[0]
        vcfFilePath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/filesToReplicateSingleAignError/shortenedScaffold_1465_SS_220.vcf"
        scaffoldDir = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/filesToReplicateSingleAignError/"
        scaffoldName = "shortenedScaffold_1465_SS_220.fasta.gbk"
        scaffoldGb = scaffoldDir + scaffoldName
        annotateScaffoldFile = True
        if annotateScaffoldFile:
            os.system(
                "cd " + scaffoldDir + "; for fileName in *.fasta; do conda run -n prokkaEnv prokka $fileName --force --centre X --compliant; cp */*.gbk " + scaffoldGb + "; rm -r PROKKA*/;done;")
            os.system("cd " + scaffoldDir + "; conda run -n gsAlign gsAlign -r /Users/cazcullimore/Documents/ericksonLabCode/refGenomes/k-12.fasta -q " + re.sub("\.gbk", "",scaffoldName) + " -unique -sen -one; mv output.vcf " + vcfFilePath)

        self.reconstructAnnotateAndCompareGenomes(vcfFilePath,scaffoldGbs=scaffoldGb, namePrefix=namePrefix,reAlignSnps=True,redoNormalAligment= True,reAnnotate= True)
        
    def testTwoGenomeAlignmnet(self):
        pathToRefGenomeFasta = "/Users/cazcullimore/Documents/ericksonLabCode/refGenomes/k-12.fasta"
        pathToRefGenomeGb = "/Users/cazcullimore/Documents/ericksonLabCode/refGenomes/k-12.gbff"
        namePrefix = "singleAlignment"  # the name to give to start all the files created
        # vcfFilePath = glob("/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/vcfsFromScaffolds/*.vcf")[0]
        genomeFastaFile = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/vcfsFromScaffolds/scaffold_1465_SS_220.fasta.vcf"
        vcfFilePath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/vcfsFromScaffolds/scaffold_1465_SS_220.fasta.vcf"
        scaffoldGb = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/scaffold_1465_SS_220.fasta.gbk"
        # reconstructAnnotateAndCompareGenomes(vcfFilePath, scaffoldGbs=scaffoldGb, namePrefix=namePrefix)
        # reconstructAnnotateAndCompareGenomes()
    def testMultiAlignment(self):
        pathToRefGenomeFasta = "/Users/cazcullimore/Documents/ericksonLabCode/refGenomes/k-12.fasta"
        pathToRefGenomeGb = "/Users/cazcullimore/Documents/ericksonLabCode/refGenomes/k-12.gbff"
        namePrefix = "multiAlignment"  # the name to give to start all the files created
        vcfFilePath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/*.vcf"
        scaffoldGb = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/scaffold_*.gbk"
        scaffoldDir = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/"
        annotateScaffoldFiles = False
        if annotateScaffoldFiles:
            os.system(
                "cd " + scaffoldDir + "; for fileName in *.fasta; do conda run -n prokkaEnv prokka $fileName --force --centre X --compliant; cp PROKKA*/*.gbk ./$fileName.gbk; rm -r PROKKA*/;done;")
            # os.system(
            #     "cd " + scaffoldDir + "; for fileName in *.fasta; do conda run -n gsAlign gsAlign -r /Users/cazcullimore/Documents/ericksonLabCode/refGenomes/k-12.fasta -q " +
            #         "${fileName%\".gbk\"} -unique -sen -one; mv output.vcf " + vcfFilePath)
        
        
        self.reconstructAnnotateAndCompareGenomes(vcfFilePath, scaffoldGbs=scaffoldGb, namePrefix=namePrefix, normalAlignPath="/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/annotatedNormalAlignFiles/multiAlignNormalFiles/",
                                                  reAlignSnps=False,redoNormalAligment= False,reAnnotate= False)

    def reconstructAnnotateAndCompareGenomes(self,vcfs, scaffoldGbs, namePrefix, normalAlignPath="/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/annotatedNormalAlignFiles/",
                                             reAlignSnps=True,redoNormalAligment= True,reAnnotate= True):
        
        annotatedRefGenomePath = "/Users/cazcullimore/Documents/ericksonLabCode/refGenomes/k-12.gbff"
        namePrefix = namePrefix  # the name to give to start all the files created
        # vcfs = glob("/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/vcfsFromScaffolds/*.vcf")[0]
        # scaffoldGb = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/annotatedNormalAlignFiles/scaffold_1465_SS_220.fasta.gbk"
        
        combinedNormalAlignFastaPath = normalAlignPath + namePrefix + ".fasta"
        
        snpAlignPath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/snpAlign/" + namePrefix + ".afa"
        snpIndexPath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/snpAlign/" + namePrefix + "Indexes.txt"
        
        print("scaffold gbs", scaffoldGbs)
        print("annotatedNormalAlignPath", combinedNormalAlignFastaPath)
        print("vcfs", vcfs)
        print("snpAlignPath", snpAlignPath)
        
        
        typesOfSnpsToSkipDuringAlignment = []#["INSERT", "DELETE"]
        if reAlignSnps:
            alignVcfSnps(vcfs, outFilePath=snpAlignPath, thingsToSkip=typesOfSnpsToSkipDuringAlignment,
                         numSnpsRequired=1, ignoreRefSeq=False)
        
        if redoNormalAligment:
            reconstructNormalAlignment(snpGenomePath=snpAlignPath, snpIndexesPath=snpIndexPath,
                                       refGenomePath=annotatedRefGenomePath, outputPath=combinedNormalAlignFastaPath)
            normalAlignData = readInFastaAsDict(combinedNormalAlignFastaPath)
            for name, seq in normalAlignData.items():
                name = name.strip()
                seq = seq.strip()
                with open(normalAlignPath + name + ".fasta", "w") as normalAlignFile:
                    normalAlignFile.write(">" + name + "\n")
                    normalAlignFile.write(seq)#re.sub("-", "", seq)+"\n'")
        
            os.system("rm " + combinedNormalAlignFastaPath)
        
        combinedNormalAlignFastaDir = "/".join(combinedNormalAlignFastaPath.split("/")[:-1])
        if reAnnotate:
            # os.system("rm " + combinedNormalAlignFastaDir + "/*.gbk; rm -r" + combinedNormalAlignFastaDir +"/PROKKA*/")
            os.system("cd " + combinedNormalAlignFastaDir + "; for fileName in *.fasta; do conda run -n prokkaEnv prokka $fileName --force --centre X --compliant && cp PROKKA*/*.gbk $fileName.gbk && rm -r PROKKA*/;done;")
        
        normalAlignGenes = {}
        for file in glob(normalAlignPath + "*.gbk"):
            normalAlignGenes[re.sub("(.fasta)+|(.vcf)+", "", file.split("/")[-1])] = []
        pool = ThreadPoolExecutor(8)
        futures = []
        for file in glob(normalAlignPath + "*.gbk"):
            genomeName = re.sub("(.fasta)+|(.vcf)+","",file.split("/")[-1])
            def getNormalAlignGenes(genomeName, file):
                normalAlignGenes[genomeName] = getGenesOnContigs(file, getContigs(file))
            futures.append(pool.submit(getNormalAlignGenes, genomeName, file))
        
        scaffoldGenesFromAllFIles = {}
        for file in glob(scaffoldGbs):
            genomeName = re.sub("(.fasta)+|(.vcf)+", "", file.split("/")[-1])
            def x(genomeName, file):
                scaffoldGenesFromAllFIles[genomeName] = getGenesOnContigs(file, getContigs(file))
            
            futures.append(pool.submit(x, genomeName, file))
        refGenes = []
        def x():
            refGenes.append(getGenesOnContigsByPosition(annotatedRefGenomePath, getContigs(annotatedRefGenomePath)))
        
        futures.append(pool.submit(x))
        
        snpAlign = readInFastaAsDict(snpAlignPath)
        
        snps = readInSnps(glob(vcfs), refGenomeSeq=getContigs(annotatedRefGenomePath)[0], ignoreRefSeq=False)
        snps.sort(key=lambda x: x[1])
        # block until done
        for future in futures:
            future.result()
        refGenes = refGenes[0]
        with open(snpIndexPath) as indexFile:
            self.CompareGenomes(indexFile, scaffoldGenesFromAllFIles, refGenes, normalAlignGenes, snps, snpAlign)
    
    def CompareGenomes(self, indexes, scaffoldGenesFromAllFIles, refGenes, normalAlignGenes, snps, snpAlign):
        for genomeNamePath, scaffoldGenes in scaffoldGenesFromAllFIles.items(): # genomeWithoutAnysnpsGenes, 1465_SS scaffoldGenes
            lineOfIndexfile = -1
            for index in indexes:
                lineOfIndexfile += 1
                index = index.strip()
                snpPos = int(float(index))
                genomeName = genomeNamePath.split("/")[-1]
                for i in range(len(refGenes)):
                    try:
                        refGene = refGenes[i]
                        scaffoldGene = scaffoldGenes[refGene.name]
                        
                        normalAlignGene = normalAlignGenes[genomeName][refGene.name]
                        if refGene.name.count("_") != 0 or refGene.name == "unnamed":
                            continue
                        if refGene.startPos <= snpPos and refGene.stopPos > snpPos:
                            print("index", index, "genome", genomeName)
                            
                            snpPosInGene = snpPos - refGene.startPos
                            refSeq = refGene.sequence
                            normalSeq = normalAlignGene.sequence
                            scaffoldSeq = scaffoldGene.sequence
                            if not scaffoldGene.isForward:
                                scaffoldSeq = reverseComplement(scaffoldSeq)
                            if not refGene.isForward:
                                refSeq = reverseComplement(refGene.sequence)
                                diffStart = len(scaffoldSeq) - len(refSeq)
                            else:
                                diffStart = len(refSeq) -  len(scaffoldSeq)
                            if not normalAlignGene.isForward:
                                normalSeq = reverseComplement(normalSeq)
                            
                            print(refGene.name + "\nnorm\n" + str(normalAlignGene.startPos) + "\t", end="")
                            printNearbyNucs(normalSeq, snpPosInGene, 10)
                            print("scaff\n" + str(scaffoldGene.startPos) + "\t", end="")
                            printNearbyNucs(scaffoldSeq,snpPosInGene  + diffStart, 10)
                            print("reverse scaff\n" + str(scaffoldGene.startPos) + "\t", end="")
                            printNearbyNucs(reverseComplement(scaffoldSeq),snpPosInGene, 10)
                            print("ref\n" + str(refGene.startPos) + "\t", end="")
                            printNearbyNucs(refSeq, snpPosInGene, 10)
                            # print("snpAlign\t\t", end="")
                            print("refGene.isForward",refGene.isForward)
                            print("normal",normalAlignGene.isForward)
                            print("scaffold",scaffoldGene.isForward)
                            
                            if normalSeq[snpPosInGene] != scaffoldSeq[snpPosInGene + diffStart]:
                                
                                print('something weird happened')
                                numDeletesBeforeSnp = 0
                                numDeletesAfterSnp = 0
                                for snp in snps:
                                   if snpPos > snp[1] and refGene.startPos < snp[1]: # does it work if on complement strand??
                                        if snp[4] != "SUBSTITUTE":
                                            numDeletesBeforeSnp += len(snp[2]) - len(snp[3])
                                   if snpPos < snp[1] and refGene.stopPos > snp[1]:
                                       if snp[4] != "SUBSTITUTE":
                                           numDeletesAfterSnp += len(snp[2]) - len(snp[3])
                                print("normalSeq v scaffold", geneSimilarity(normalAlignGene.sequence, scaffoldGene.sequence))
                                print("scaffold v ref", geneSimilarity(scaffoldGene.sequence, refGene.sequence))
                                print("normalSeq v ref", geneSimilarity(normalAlignGene.sequence, refGene.sequence))
                                print("numDeletesBeforeSnp", numDeletesBeforeSnp)
                                print("numDeletesAfterSnp",numDeletesAfterSnp)
                            # self.assertEqual(scaffoldGene.sequence[snpPosInGene], snpAlign[genomeNamePath.split("/")[-1]][lineOfIndexfile]) # if this test fails it means that there is a problem with gsAlign or with my indexes
                            print("--------------")
                    except IndexError:
                        print("INDEX ERROR")
                    except KeyError:
                        pass
    def testCompareGenes(self):
        pass
    def testSnpAlignMatchesScaffold(self): #doesn't work because doesn't start from the same place and any indels mess it up
        pathToRefGenomeFasta = "/Users/cazcullimore/Documents/ericksonLabCode/refGenomes/k-12.fasta"
        pathToRefGenomeGb = "/Users/cazcullimore/Documents/ericksonLabCode/refGenomes/k-12.gbff"
        namePrefix = "snpAlignMatchesScaffold"  # the name to give to start all the files created
        vcfFilePath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/output.vcf"
        scaffoldFastaPath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/scaffold_1465_SS_220.fasta"
        refSeq = readInFastaAsList(pathToRefGenomeFasta)[1]
        snpAlign, indexes, frameShifts = alignVCFSnpsHelper([vcfFilePath],refSeq=refSeq,thingsToSkip=[],ignoreRefSeq=False,numSnpsRequired=1)
        scaffoldGbPath = ""
        scaffoldSeq = readInFastaAsList(scaffoldFastaPath)[1]
        numIndels = 0
        for alignedNuc, index in zip(snpAlign["output.vcf"], indexes):
            if float(int(index)) != float(index):
                numIndels += 1
                print("indel", index, alignedNuc)
            elif alignedNuc == "-":
                numIndels -= 1
                print("indel", index, alignedNuc)
                continue
            print("index", index)
            printNearbyNucs(scaffoldSeq,int(index) - 16751 + 25 + numIndels, 5)
            print("\t",alignedNuc)
            # self.assertNotEquals(alignedNuc, scaffoldSeq[int(index) - 78847 + 18188])
        print("hi")
def printNearbyNucs(seq, pos, numFlankingNucs):
    try:
        print(seq[pos - numFlankingNucs + 1:pos],
              seq[pos],
              seq[pos + 1: pos + numFlankingNucs])
    except IndexError:
        print("out of range")
    
if __name__ == '__main__':
    unittest.main()