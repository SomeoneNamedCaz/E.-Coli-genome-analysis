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
    def __init__(self, arg):
        super().__init__(arg)
        self.annotatedRefGenomePath = "/Users/cazcullimore/Documents/ericksonLabCode/refGenomes/k-12.gbff"
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
        for genomeNamePath, scaffoldGenes in scaffoldGenesFromAllFIles.items():
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
                                diffStart = len(refSeq) - len(scaffoldSeq)
                            if not normalAlignGene.isForward:
                                normalSeq = reverseComplement(normalSeq)
                            # if normalSeq[snpPosInGene + diffStart] != scaffoldSeq[snpPosInGene + diffStart]:
                            #
                            #     numDeletesBeforeSnp = 0
                            #     numDeletesAfterSnp = 0
                            #     for snp in snps:
                            #         if snpPos > snp[1] and refGene.startPos < snp[
                            #             1]:  # does it work if on complement strand??
                            #             if snp[4] != "SUBSTITUTE":
                            #                 numDeletesBeforeSnp += len(snp[2]) - len(snp[3])
                            #         if snpPos < snp[1] and refGene.stopPos > snp[1]:
                            #             if snp[4] != "SUBSTITUTE":
                            #                 numDeletesAfterSnp += len(snp[2]) - len(snp[3])
                            #     # diffStart = numDeletesBeforeSnp
                            print(refGene.name + "\nnorm\n" + str(normalAlignGene.startPos) + "\t", end="")
                            printNearbyNucs(normalSeq, snpPosInGene - diffStart, 10)
                            print("scaff\n" + str(scaffoldGene.startPos) + "\t", end="")
                            printNearbyNucs(scaffoldSeq, snpPosInGene - diffStart, 10)
                            print("shifted scaf\n" + str(scaffoldGene.startPos) + "\t", end="")
                            printNearbyNucs(scaffoldSeq, snpPosInGene + diffStart, 10)
                            print("ref\n" + str(refGene.startPos) + "\t", end="")
                            printNearbyNucs(refSeq, snpPosInGene, 10)
                            # print("snpAlign\t\t", end="")
                            print("refGene.isForward",refGene.isForward)
                            print("normal",normalAlignGene.isForward)
                            print("scaffold",scaffoldGene.isForward)
                            
                            if normalSeq[snpPosInGene + diffStart] != scaffoldSeq[snpPosInGene + diffStart]:
                                
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
                                if refGene.isForward:
                                    # thought process: if it didn't match the first time it means that it has prepended
                                    # nucleotides and hopefully not appended nucleotides
                                    ref = "ACTGFGAGAGATAGAGAGARRGR"
                                    scaf =            "AGAGAGARRGR"
                                    
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
    def testCompareGenomesNhaR(self):
        indexes = ["18767", " 18779", " 18812", " 18830", " 18836", " 18839", " 18863", " 18867", " 18875", " 18878", " 18881", " 18893", " 18905", " 18911", " 18971", " 18974", " 18989", " 18995", " 19007", " 19013", " 19019", " 19037", " 19046", " 19049", " 19052", " 19055", " 19056", " 19058", " 19067", " 19070", " 19075", " 19081", " 19082", " 19085", " 19091", " 19097", " 19106", " 19109", " 19115", " 19121", " 19122", " 19124", " 19130", " 19133", " 19160", " 19175", " 19178", " 19181", " 19184", " 19196", " 19202", " 19206", " 19210", " 19226", " 19229", " 19241", " 19244", " 19247", " 19252", " 19262", " 19268", " 19271", " 19274", " 19280", " 19286", " 19289", " 19292", " 19293", " 19295", " 19296", " 19298", " 19304", " 19310", " 19313", " 19322", " 19328", " 19332", " 19334", " 19335", " 19337", " 19358", " 19361", " 19367", " 19373", " 19374", " 19376", " 19379", " 19394", " 19397", " 19409", " 19415", " 19420", " 19424", " 19427", " 19430", " 19433", " 19442", " 19445", " 19446", " 19454", " 19457", " 19458", " 19463", " 19466", " 19470", " 19471", " 19476", " 19478", " 19479", " 19481", " 19484", " 19487", " 19490", " 19493", " 19496", " 19499", " 19502", " 19508", " 19526", " 19529", " 19538", " 19541", " 19542", "19612"]
        scaffoldGenesFromAllFiles = {"scaffold_1466_SS_221.fasta.gbk":{"nhaR": Gene(5872, 6772, 'ATGTCTCATATCAATTACAACCACTTGTATTACTTCTGGCATGTCTACAAAGAAGGTTCTGTGGTTGGCGCAGCGGAGGCGCTTTATTTAACACCACAAACCATTACCGGGCAGATCCGGGCGCTGGAAGAGCGCCTGCAAGGGAAACTATTTAAGCGTAAAGGACGTGGTCTGGAACCCAGCGAACTGGGGGAACTGGTCTATCGCTATGCCGATAAAATGTTCACCTTAAGCCAGGAAATGCTGGATATCGTCAACTATCGCAAAGAGTCCAACTTATTGTTTGATGTTGGTGTGGCAGATGCACTTTCCAAACGTCTGGTCAGCAGTGTTCTGGATGCCGCAGTTGTGGAAGACGAGCAGATCCATCTACGCTGTTTCGAATCGACGCACGAGATGCTTTTAGAGCAGTTGAGTCAGCATAAACTGGATATGATCATCTCTGACTGTCCGATCGATTCCACTCAGCAGGAAGGGCTGTTTTCCATGAAAATTGGCGAATGTGGTGTCAGTTTCTGGTGCACTAACCCACTACCAGAAAAGCCGTTTCCTGCCTGTCTTGAAGAGCGTCGTTTACTTATTCCGGGGCGTCGCTCAATGTTGGGGCGTAAACTATTAAACTGGTTTAACTCCCAGGGCTTGAACGTCGAAATTTTGGGTGAGTTTGATGATGCTGCGTTGATGAAAGCCTTTGGGGCGACGCATAACGCTATTTTCGTTGCACCTTCGCTTTACGCTAATGATTTCTATAACGATGACTCGGTTGTGGAGATAGGCCGTGTTGAGAACGTGATGGAAGAGTACCACGCGATTTTTGCCGAAAGGATGATTCAGCACCCTGCAGTACAGCGTATCTGCAATACAGACTATTCTGCGCTGTTTACTCCAGCTTCAAAATAA',
                                 "nhaR", "product",True)}}
        refGenes = [Gene(18714,19620, 'ATGAGCATGTCTCATATCAATTACAACCACTTGTATTACTTCTGGCATGTCTATAAAGAAGGTTCCGTGGTTGGCGCAGCGGAGGCGCTTTATTTAACTCCACAAACCATTACCGGACAGATTCGAGCGCTGGAAGAGCGCCTGCAAGGCAAATTATTTAAACGCAAGGGACGTGGTCTCGAACCCAGCGAGCTGGGAGAACTGGTCTATCGCTATGCCGATAAAATGTTCACCTTAAGCCAGGAAATGCTGGATATTGTGAACTATCGCAAAGAATCCAATTTATTGTTTGACGTTGGCGTGGCTGATGCACTTTCCAAACGCCTGGTCAGTAGCGTACTTAACGCCGCAGTGGTAGAAGGCGAGCCCATTCATCTTCGCTGCTTCGAATCCACCCACGAAATGCTGCTGGAGCAATTAAGTCAGCATAAACTGGATATGATCATTTCTGACTGTCCGATAGACTCTACGCAGCAGGAAGGCCTGTTCTCCGTGAGAATTGGCGAATGTGGCGTGAGTTTCTGGTGTACAAATCCACCACCAGAAAAACCGTTCCCGGCTTGTCTGGAAGAACGGCGACTTTTGATTCCTGGGCGACGTTCAATGTTAGGGCGCAAATTGCTTAACTGGTTTAACTCCCAGGGATTAAACGTAGAAATCCTCGGCGAGTTTGATGATGCCGCTTTGATGAAAGCTTTTGGTGCGATGCACAATGCAATCTTCGTTGCCCCAACGCTTTATGCATATGACTTTTATGCCGATAAAACTGTCGTAGAAATTGGTCGCGTCGAGAATGTGATGGAAGAGTACCATGCTATTTTTGCTGAGCGGATGATTCAGCACCCGGCGGTACAGCGAATCTGCAATACGGATTATTCTGCGCTTTTTAGTCCAGCGGTGCGTTAA',
                                 "nhaR", "product",True)]
        normalAlignGenes = {"scaffold_1466_SS_221.fasta.gbk":{"nhaR": Gene(18715,19615, 'ATGTCTCATATCAATTACAACCACTTGTATTACTTCTGGCATGTCTACAAAGAAGGTTCTGTGGTTGGCGCAGCGGAGGCGCTTTATTTAACACCACAAACCATTACCGGGCAGATCCGGGCGCTGGAAGAGCGCCTGCAAGGGAAACTATTTAAGCGTAAAGGACGTGGTCTGGAACCCAGCGAACTGGGGGAACTGGTCTATCGCTATGCCGATAAAATGTTCACCTTAAGCCAGGAAATGCTGGATATCGTCAACTATCGCAAAGAGTCCAACTTATTGTTTGATGTTGGTGTGGCAGATGCACTTTCCAAACGTCTGGTCAGCAGTGTTCTGGATGCCGCAGTTGTGGAAGACGAGCAGATCCATCTACGCTGTTTCGAATCGACGCACGAGATGCTTTTAGAGCAGTTGAGTCAGCATAAACTGGATATGATCATCTCTGACTGTCCGATCGATTCCACTCAGCAGGAAGGGCTGTTTTCCATGAAAATTGGCGAATGTGGTGTCAGTTTCTGGTGCACTAACCCACTACCAGAAAAGCCGTTTCCTGCCTGTCTTGAAGAGCGTCGTTTACTTATTCCGGGGCGTCGCTCAATGTTGGGGCGTAAACTATTAAACTGGTTTAACTCCCAGGGCTTGAACGTCGAAATTTTGGGTGAGTTTGATGATGCTGCGTTGATGAAAGCCTTTGGGGCGACGCATAACGCTATTTTCGTTGCACCTTCGCTTTACGCTAATGATTTCTATAACGATGACTCGGTTGTGGAGATAGGCCGTGTTGAGAACGTGATGGAAGAGTACCACGCGATTTTTGCCGAAAGGATGATTCAGCACCCGGCGGTACAGCGAATCTGCAATACGGATTATTCTGCGCTTTTTAGTCCAGCGGTGCGTTAA',
                                 "nhaR", "product",True)}}
        vcfFilePath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/*.vcf"
        snps = readInSnps(glob(vcfFilePath), refGenomeSeq=getContigs(self.annotatedRefGenomePath)[0], ignoreRefSeq=False)
        self.CompareGenomes(indexes, scaffoldGenesFromAllFiles, refGenes, normalAlignGenes, snps, {})
    def testCompareGenomesApaH(self):
        snpIndexPath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/snpAlign/MultiAlignmentIndexes.txt"
        scaffoldGenesFromAllFiles = {"scaffold_1466_SS_221.fasta.gbk":{"apaH": Gene(36990, 37839, 'ATGGCGACATACCTTATTGGCGACGTTCATGGTTGTTACGATGAACTGATCGCATTGCTGCATAAAGTAGAATTTACCCCTGGGAAAGATACCCTCTGGCTGACGGGCGATCTGGTCGCGCGCGGCCCGGGTTCGCTGGATGTTCTGCGCTATGTGAAATCCTTAGGCGACAGCGTACGTCTGGTGCTGGGTAATCACGATCTGCATCTGCTGGCGGTATTTGCCGGGATCAGCCGCAATAAACCGAAAGATCGCCTGACACCGCTGCTGGAAGCGCCGGATGCCGACGAGCTGCTTAACTGGCTGCGTCGCCAGCCTTTGCTGCAAATCGACGAAGAGAAAAAGTTGGTGATGGCCCACGCCGGGATCACGCCGCAGTGGGATCTGCAGACCGCCAAAGAGTGCGCGCGCGATGTAGAAGCGGTGCTGTCGAGTGACTCCTATCCCTTCTTTCTTGATGCCATGTACGGCGATATGCCAAATAACTGGTCACCGGAATTGCGGGGGCTGGGAAGACTGCGTTTTATCACCAACGCCTTTACCCGTATGCGTTTTTGCTTCCCGAACGGTCAACTGGATATGTACAGCAAAGAATCGCCGGAAGAGGCCCCTGCCCCACTGAAACCGTGGTTTGCGATTCCTGGCCCCGTCGCTGAAGAGTACAACATCGCCTTTGGTCACTGGGCATCGCTGGAAGGCAAAGGTACGCCGGAAGGTATTTACGCGCTGGATACCGGCTGCTGCTGGGGTGGTACATTAACCTGCCTGCGCTGGGAAGATAAACAGTATTTTGTCCAGCCGTCGAACCGGCATAAGGATTTGAGTGAGGGAGAGGCGGTAGCGTCTTAA',
                                 "apaH", "product",False)}}
        refGenes = [Gene(50379,51222, 'ATGGCGACATACCTTATTGGCGACGTTCATGGTTGTTACGATGAACTGATCGCATTGCTGCATAAAGTAGAATTTACCCCTGGGAAAGATACCCTCTGGCTGACGGGCGATCTGGTCGCGCGCGGGCCGGGTTCGCTGGATGTTCTGCGCTATGTGAAATCCTTAGGCGACAGCGTACGTCTGGTGCTGGGCAATCACGATCTGCATCTGCTGGCGGTATTTGCCGGGATCAGCCGCAATAAACCGAAAGATCGCCTGACACCGCTGCTGGAAGCGCCGGATGCCGACGAGCTGCTTAACTGGCTGCGGCGCCAGCCTCTGCTGCAAATCGACGAAGAGAAAAAGCTGGTGATGGCCCACGCAGGGATCACGCCGCAGTGGGATCTGCAGACCGCCAAAGAGTGCGCACGCGATGTAGAAGCGGTGCTATCGAGTGACTCCTATCCCTTCTTTCTTGATGCCATGTACGGCGATATGCCAAATAACTGGTCACCGGAATTGCGGGGGCTGGGAAGACTGCGTTTTATCACCAACGCTTTTACCCGTATGCGTTTTTGCTTCCCGAACGGTCAACTGGATATGTACAGCAAAGAATCGCCGGAAGAGGCCCCTGCCCCACTGAAACCGTGGTTTGCGATTCCTGGCCCTGTCGCTGAAGAATACAGCATCGCCTTTGGTCACTGGGCATCGCTGGAGGGCAAAGGTACGCCGGAAGGTATATACGCGCTGGATACCGGCTGCTGCTGGGGTGGTACATTAACCTGCCTGCGCTGGGAAGATAAACAGTATTTTGTCCAGCCGTCGAACCGGCATAAGGATTTGGGCGAAGCGGCGGCGTCTTAA',
                                 "apaH", "product",False)]
        normalAlignGenes = {"scaffold_1466_SS_221.fasta.gbk":{"apaH": Gene(50377,51220, 'ATGGCGACATACCTTATTGGCGACGTTCATGGTTGTTACGATGAACTGATCGCATTGCTGCATAAAGTAGAATTTACCCCTGGGAAAGATACCCTCTGGCTGACGGGCGATCTGGTCGCGCGCGGCCCGGGTTCGCTGGATGTTCTGCGCTATGTGAAATCCTTAGGCGACAGCGTACGTCTGGTGCTGGGTAATCACGATCTGCATCTGCTGGCGGTATTTGCCGGGATCAGCCGCAATAAACCGAAAGATCGCCTGACACCGCTGCTGGAAGCGCCGGATGCCGACGAGCTGCTTAACTGGCTGCGTCGCCAGCCTTTGCTGCAAATCGACGAAGAGAAAAAGTTGGTGATGGCCCACGCCGGGATCACGCCGCAGTGGGATCTGCAGACCGCCAAAGAGTGCGCGCGCGATGTAGAAGCGGTGCTGTCGAGTGACTCCTATCCCTTCTTTCTTGATGCCATGTACGGCGATATGCCAAATAACTGGTCACCGGAATTGCGGGGGCTGGGAAGACTGCGTTTTATCACCAACGCCTTTACCCGTATGCGTTTTTGCTTCCCGAACGGTCAACTGGATATGTACAGCAAAGAATCGCCGGAAGAGGCCCCTGCCCCACTGAAACCGTGGTTTGCGATTCCTGGCCCCGTCGCTGAAGAGTACAACATCGCCTTTGGTCACTGGGCATCGCTGGAAGGCAAAGGTACGCCGGAAGGTATTTACGCGCTGGATACCGGCTGCTGCTGGGGTGGTACATTAACCTGCCTGCGCTGGGAAGATAAACAGTATTTTGTCCAGCCGTCGAACCGGCATAAGGATTTGGGCGAAGCGGCGGCGTCTTAA',
                                 "apaH", "product",False)}}
        vcfFilePath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/*.vcf"
        snps = readInSnps(glob(vcfFilePath), refGenomeSeq=getContigs(self.annotatedRefGenomePath)[0], ignoreRefSeq=False)
        with open(snpIndexPath) as indexes:
            self.CompareGenomes(indexes, scaffoldGenesFromAllFiles, refGenes, normalAlignGenes, snps, {})
    def testSnpAlignMatchesScaffold(self): #doesn't work because doesn't start from the same place and any indels mess it up
        pathToRefGenomeFasta = "/Users/cazcullimore/Documents/ericksonLabCode/refGenomes/k-12.fasta"
        pathToRefGenomeGb = "/Users/cazcullimore/Documents/ericksonLabCode/refGenomes/k-12.gbff"
        namePrefix = "snpAlignMatchesScaffold"  # the name to give to start all the files created
        vcfFilePath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/output.vcf"
        scaffoldFastaPath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/scaffold_1465_SS_220.fasta"
        refSeq = readInFastaAsList(pathToRefGenomeFasta)[1]
        snpAlign, indexes, frameShifts = alignVCFSnpsHelper([vcfFilePath],refSeq=refSeq,thingsToSkip=[],ignoreRefSeq=False,numSnpsRequired=1)
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