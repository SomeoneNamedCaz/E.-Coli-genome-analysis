import sys
sys.path.insert(1, '/Users/cazcullimore/dev/ericksonLabCode/')
from alignVcfSnps import *
from megaCatsPythonVersion import *
from parsingMegaCatsResults import *
from reconstructNormalAlignment import *
import unittest
sys.path.insert(1, '/Users/cazcullimore/dev/ericksonLabCode/secondaryPythonScripts')
try:
    from secondaryPythonScripts.functions import *
except ImportError:
    from functions import *
# hscC is weird
# this doesn't work because it the gbk phylogroupSnpFile isn't aligned
class MyTestCase(unittest.TestCase):
    def __init__(self, arg):
        super().__init__(arg)
        self.annotatedRefGenomePath = DATA_DIR + "refGenomes/k-12.gbff"
    def testSingleAlignment(self):
        pass
        pathToRefGenomeFasta = DATA_DIR + "refGenomes/k-12.fasta"
        pathToRefGenomeGb = DATA_DIR + "refGenomes/k-12.gbff"
        namePrefix = "shortSingleAlignment"  # the name to give to start all the files created
        # vcfFilePath = glob(TEST_DATA_DIR + "vcfsFromScaffolds/*.vcf")[0]
        vcfFilePath = TEST_DATA_DIR + "filesToReplicateSingleAignError/shortenedScaffold_1465_SS_220.vcf"
        scaffoldDir = TEST_DATA_DIR + "filesToReplicateSingleAignError/"
        scaffoldName = "shortenedScaffold_1465_SS_220.fasta.gbk"
        scaffoldGb = scaffoldDir + scaffoldName
        annotateScaffoldFile = False
        if annotateScaffoldFile:
            os.system(
                "cd " + scaffoldDir + "; for fileName in *.fasta; do conda run -n prokkaEnv prokka $fileName --force --centre X --compliant; cp */*.gbk " + scaffoldGb + "; rm -r PROKKA*/;done;")
            os.system("cd " + scaffoldDir + "; conda run -n gsAlign gsAlign -r /Users/cazcullimore/Documents/ericksonLabCode/refGenomes/k-12.fasta -q " + re.sub("\.gbk", "",scaffoldName) + " -unique -sen -one; mv output.vcf " + vcfFilePath)

        self.reconstructAnnotateAndCompareGenomes(vcfFilePath,scaffoldGbs=scaffoldGb, namePrefix=namePrefix,reAlignSnps=True,redoNormalAligment= True,reAnnotate= False)
        
    def testTwoGenomeAlignment(self):
        pathToRefGenomeFasta = DATA_DIR + "refGenomes/k-12.fasta"
        pathToRefGenomeGb = DATA_DIR + "refGenomes/k-12.gbff"
        namePrefix = "singleAlignment"  # the name to give to start all the files created
        # vcfFilePath = glob(TEST_DATA_DIR + "vcfsFromScaffolds/*.vcf")[0]
        genomeFastaFile = TEST_DATA_DIR + "vcfsFromScaffolds/scaffold_1465_SS_220.fasta.vcf"
        vcfFilePath = TEST_DATA_DIR + "vcfsFromScaffolds/scaffold_1465_SS_220.fasta.vcf"
        scaffoldGb = TEST_DATA_DIR + "scaffold_1465_SS_220.fasta.gbk"
        # reconstructAnnotateAndCompareGenomes(vcfFilePath, scaffoldGbs=scaffoldGb, namePrefix=namePrefix)
        # reconstructAnnotateAndCompareGenomes()
    def testMultiAlignment(self):
        pathToRefGenomeFasta = DATA_DIR + "refGenomes/k-12.fasta"
        pathToRefGenomeGb = DATA_DIR + "refGenomes/k-12.gbff"
        namePrefix = "multiAlignment"  # the name to give to start all the files created
        vcfDir = DATA_DIR + "k-12RefGenomeAnalysisFiles/AllAssemblies/allBovineScaffolds/gsAlignOutput/"
        vcfFilePath = vcfDir + "*.vcf"
        scaffoldDir = DATA_DIR + "k-12RefGenomeAnalysisFiles/AllAssemblies/allBovineScaffolds/"
        scaffoldGb = scaffoldDir + "*.gbk"
        annotateScaffoldFiles = False
        if annotateScaffoldFiles:
            os.system(
                "cd " + scaffoldDir + "; for fileName in *.fasta; do conda run -n prokkaEnv prokka $fileName --force --centre X --compliant; cp PROKKA*/*.gbk ./$fileName.gbk; rm -r PROKKA*/; done;")
            os.system(
                "cd " + scaffoldDir + "; for fileName in *.fasta; do conda run -n gsAlign gsAlign -r /Users/cazcullimore/dev/data/refGenomes/k-12.fasta -q " +
                    "${fileName%\".gbk\"} -unique -sen -one; mv output.vcf " + vcfDir + "$fileName.vcf; done")
        
        
        self.reconstructAnnotateAndCompareGenomes(vcfFilePath, scaffoldGbs=scaffoldGb, namePrefix=namePrefix, normalAlignPath=TEST_DATA_DIR + "annotatedNormalAlignFiles/multiAlignNormalFiles/",
                                                  reAlignSnps=False,redoNormalAligment=False,reAnnotate= False)

    def reconstructAnnotateAndCompareGenomes(self,vcfs, scaffoldGbs, namePrefix, normalAlignPath=TEST_DATA_DIR + "annotatedNormalAlignFiles/",
                                             reAlignSnps=True,redoNormalAligment= True,reAnnotate= True):
        
        annotatedRefGenomePath = "/Users/cazcullimore/dev/ericksonLabCode/tests/unitTests/k-12.fasta.gbk"
        namePrefix = namePrefix  # the name to give to start all the files created
        # vcfs = glob(TEST_DATA_DIR + "vcfsFromScaffolds/*.vcf")[0]
        # scaffoldGb = TEST_DATA_DIR + "annotatedNormalAlignFiles/scaffold_1465_SS_220.fasta.gbk"
        
        combinedNormalAlignFastaPath = normalAlignPath + namePrefix + ".fasta"
        
        snpAlignPath = TEST_DATA_DIR + "snpAlign/" + namePrefix + ".afa"
        snpIndexPath = TEST_DATA_DIR + "snpAlign/" + namePrefix + "Indexes.txt"
        
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
        normalFutures = []
        for file in glob(normalAlignPath + "*.gbk"):
            genomeName = re.sub("(.fasta)+|(.vcf)+","",file.split("/")[-1])
            def getNormalAlignGenes(genomeName, file):
                return getGenesOnContigs(file, getContigs(file))
            normalFutures.append((genomeName,pool.submit(getNormalAlignGenes, genomeName, file)))
        
        scaffoldGenesFromAllFIles = {}
        for file in glob(scaffoldGbs):
            genomeName = re.sub("(.fasta)+|(.vcf)+", "", file.split("/")[-1])
            scaffoldGenesFromAllFIles[genomeName] = []
        scaffoldFutures = []
        for file in glob(scaffoldGbs):
            genomeName = re.sub("(.fasta)+|(.vcf)+", "", file.split("/")[-1])
            def x(genomeName, file):
                return getGenesOnContigs(file, getContigs(file))
            
            scaffoldFutures.append((genomeName,pool.submit(x, genomeName, file)))
        refGenes = getGenesOnContigsByPosition(annotatedRefGenomePath, getContigs(annotatedRefGenomePath))
        
        
        snpAlign = readInFastaAsDict(snpAlignPath)
        
        snps = readInSnps(glob(vcfs), refGenomeSeq=getContigs(annotatedRefGenomePath)[0], ignoreRefSeq=False)
        snps.sort(key=lambda x: x[1])
        # block until done
        for future in scaffoldFutures:
            scaffoldGenesFromAllFIles[future[0]] = future[1].result()
        
        for future in normalFutures:
            normalAlignGenes[future[0]] = future[1].result()
            
            
        with open(snpIndexPath) as indexFile:
            self.CompareGenomes(indexFile, scaffoldGenesFromAllFIles, refGenes, normalAlignGenes, snps, snpAlign)
    
    def CompareGenomes(self, indexes, scaffoldGenesFromAllFiles, refGenes, normalAlignGenes, snps, snpAlign):
        lineOfIndexfile = -1
        genomeNameToMetadata = readMetaDataAsDict(DATA_DIR + "metaDataForMetaCatsWithExtraMastitis.tsv")
        for index in indexes:
            index = index.strip()
            snpPos = int(float(index))
            if snpPos < 4431321 + 569:#4433164 + 740:#4512411 + 630: #1570645: #
                continue
            lineOfIndexfile += 1
            geneDistributions = {"pathogen": {"A": 0, "T": 0, "C": 0, "G": 0, "noGene":0},
                                 "commensal": {"A": 0, "T": 0, "C": 0, "G": 0, "noGene":0}, "cow": {"A": 0, "T": 0, "C": 0, "G": 0},
                                 "chicken": {"A": 0, "T": 0, "C": 0, "G": 0}}
            numGenomesWithoutGene = 0
            numGenomesWithBadAlignment = 0
            for genomeNamePath, scaffoldGenes in scaffoldGenesFromAllFiles.items():
                
                genomeName = genomeNamePath.split("/")[-1]
                print("genomeName", genomeName)
                for i in range(len(refGenes)):
                    try:
                        refGene = refGenes[i]
                        if refGene.name.count("_") != 0 or refGene.name == "unnamed":
                            continue
                        if refGene.startPos <= snpPos and refGene.stopPos > snpPos:
                            print("index", index, "genome", genomeName)
                            try:
                                scaffoldGene = scaffoldGenes[refGene.name]
                                
                                normalAlignGene = normalAlignGenes[genomeName][refGene.name]
                            except KeyError:
                                print("scaffold or normal align doesn't have gene")
                                geneDistributions[
                                    genomeNameToMetadata[re.sub("(.fasta)|(.gbk)", "", genomeName) + ".fasta.vcf"][-1]][
                                    "noGene"] += 1
                                numGenomesWithoutGene += 1
                                continue
                            
                            snpPosInGene = snpPos - refGene.startPos
                            refSeq = refGene.sequence
                            normalSeq = normalAlignGene.sequence
                            scaffoldSeq = scaffoldGene.sequence
                            if not scaffoldGene.isForward:
                                scaffoldSeq = reverseComplement(scaffoldSeq)
                            if not refGene.isForward:
                                refSeq = reverseComplement(refGene.sequence)
                                diffScaffoldStart = len(refSeq) - len(scaffoldSeq)
                                diffNormalStart = len(refSeq) - len(normalSeq)
                            else:
                                diffScaffoldStart = len(refSeq) - len(scaffoldSeq)
                                diffNormalStart = len(refSeq) - len(normalSeq)
                            if not normalAlignGene.isForward:
                                normalSeq = reverseComplement(normalSeq)
                            print(refGene.name + " index within gene " + str(snpPosInGene) + "\nnorm\n" + str(normalAlignGene.startPos) + "\t", end="")
                            printNearbyNucs(normalSeq, snpPosInGene - diffNormalStart, 10)
                            print("scaff\n" + str(scaffoldGene.startPos) + "\t", end="")
                            printNearbyNucs(scaffoldSeq, snpPosInGene - diffScaffoldStart, 10)
                            # print("shifted scaf\n" + str(scaffoldGene.startPos) + "\t", end="")
                            # printNearbyNucs(scaffoldSeq, snpPosInGene + diffScaffoldStart, 10)
                            print("ref\n" + str(refGene.startPos) + "\t", end="")
                            printNearbyNucs(refSeq, snpPosInGene, 10)
                            # print("snpAlign\t\t", end="")
                            print("refGene.isForward",refGene.isForward)
                            print("normal",normalAlignGene.isForward)
                            print("scaffold",scaffoldGene.isForward)
                            print("diffScaffoldStart", diffScaffoldStart)
                            print("diffNormalStart", diffNormalStart)
                            # print("genomeNameToMetadata[genomeName]][scaffoldSeq[snpPosInGene - diffScaffoldStart]",
                            #       genomeNameToMetadata[re.sub(".gbk","",genomeName) + ".fasta.vcf"])
                            print("print success")
                            if normalSeq[snpPosInGene - diffNormalStart] == scaffoldSeq[snpPosInGene - diffScaffoldStart]:
                                print("diffStart")
                                geneDistributions[genomeNameToMetadata[re.sub("(.fasta)|(.gbk)","",genomeName) + ".fasta.vcf"][-1]][scaffoldSeq[snpPosInGene - diffScaffoldStart]] += 1
                            else:
                                print("not diffStart")
                                if normalSeq[snpPosInGene] == scaffoldSeq[snpPosInGene]:
                                    print("otherOneWorks")
                                    geneDistributions[genomeNameToMetadata[re.sub("(.fasta)|(.gbk)","",genomeName) + ".fasta.vcf"][-1]][scaffoldSeq[snpPosInGene]] += 1
                                else:
                                    numGenomesWithBadAlignment += 1
                                    print("broken alignment not added to counts")
                                print('something weird happened')
                            print("trying without shift")
                            print(refGene.name + " index within gene " + str(snpPosInGene) + "\nnorm\n" + str(
                                normalAlignGene.startPos) + "\t", end="")
                            printNearbyNucs(normalSeq, snpPosInGene, 10)
                            print("scaff\n" + str(scaffoldGene.startPos) + "\t", end="")
                            printNearbyNucs(scaffoldSeq, snpPosInGene, 10)
                            # print("shifted scaf\n" + str(scaffoldGene.startPos) + "\t", end="")
                            # printNearbyNucs(scaffoldSeq, snpPosInGene + diffScaffoldStart, 10)
                            print("ref\n" + str(refGene.startPos) + "\t", end="")
                            printNearbyNucs(refSeq, snpPosInGene, 10)
                            
                            numDeletesBeforeSnp = 0
                            numDeletesAfterSnp = 0
                            for snp in snps:
                                if snpPos > snp[1] and refGene.startPos < snp[1]: # does it work if on complement strand??
                                    if snp[4] != "SUBSTITUTE":
                                        numDeletesBeforeSnp += len(snp[2]) - len(snp[3])
                                if snpPos < snp[1] and refGene.stopPos > snp[1]:
                                    if snp[4] != "SUBSTITUTE":
                                       numDeletesAfterSnp += len(snp[2]) - len(snp[3])
                            print("shift num deletes before snp - ")
                            print(refGene.name + " index within gene " + str(snpPosInGene) + "\nnorm\n" + str(
                                normalAlignGene.startPos) + "\t", end="")
                            printNearbyNucs(normalSeq, snpPosInGene - numDeletesBeforeSnp, 10)
                            print("scaff\n" + str(scaffoldGene.startPos) + "\t", end="")
                            printNearbyNucs(scaffoldSeq, snpPosInGene - numDeletesBeforeSnp, 10)
                            print("ref\n" + str(refGene.startPos) + "\t", end="")
                            printNearbyNucs(refSeq, snpPosInGene, 10)
                            
                            print("shift num deletes before snp + ")
                            print(refGene.name + " index within gene " + str(snpPosInGene) + "\nnorm\n" + str(
                                normalAlignGene.startPos) + "\t", end="")
                            printNearbyNucs(normalSeq, snpPosInGene + numDeletesBeforeSnp, 10)
                            print("scaff\n" + str(scaffoldGene.startPos) + "\t", end="")
                            printNearbyNucs(scaffoldSeq, snpPosInGene + numDeletesBeforeSnp, 10)
                            print("ref\n" + str(refGene.startPos) + "\t", end="")
                            printNearbyNucs(refSeq, snpPosInGene, 10)
                            
                            print("shift num deletes after snp - ")
                            print(refGene.name + " index within gene " + str(snpPosInGene) + "\nnorm\n" + str(
                                normalAlignGene.startPos) + "\t", end="")
                            printNearbyNucs(normalSeq, snpPosInGene - numDeletesAfterSnp, 10)
                            print("scaff\n" + str(scaffoldGene.startPos) + "\t", end="")
                            printNearbyNucs(scaffoldSeq, snpPosInGene - numDeletesAfterSnp, 10)
                            print("ref\n" + str(refGene.startPos) + "\t", end="")
                            printNearbyNucs(refSeq, snpPosInGene, 10)
                            
                            print("shift num deletes after snp + ")
                            print(refGene.name + " index within gene " + str(snpPosInGene) + "\nnorm\n" + str(
                                normalAlignGene.startPos) + "\t", end="")
                            printNearbyNucs(normalSeq, snpPosInGene + numDeletesAfterSnp, 10)
                            print("scaff\n" + str(scaffoldGene.startPos) + "\t", end="")
                            printNearbyNucs(scaffoldSeq, snpPosInGene + numDeletesAfterSnp, 10)
                            print("ref\n" + str(refGene.startPos) + "\t", end="")
                            printNearbyNucs(refSeq, snpPosInGene, 10)
                            
                            print("normalSeq v scaffold", geneSimilarity(normalAlignGene.sequence, scaffoldGene.sequence))
                            print("scaffold v ref", geneSimilarity(scaffoldGene.sequence, refGene.sequence))
                            print("normalSeq v ref", geneSimilarity(normalAlignGene.sequence, refGene.sequence))
                            print("numDeletesBeforeSnp", numDeletesBeforeSnp)
                            print("numDeletesAfterSnp",numDeletesAfterSnp)
                            # self.assertEqual(scaffoldGene.sequence[snpPosInGene], snpAlign[genomeNamePath.split("/")[-1]][lineOfIndexfile]) # if this test fails it means that there is a problem with gsAlign or with my indexes
                            print("--------------")
                    except IndexError:
                        print("INDEX ERROR")
                    except KeyError as err:
                        print("KEY ERROR",err)
            print(geneDistributions)
            print("numGenomesWithoutGene",numGenomesWithoutGene)
            print("numGenomesWithBadAlignment", numGenomesWithBadAlignment)
    def testCompareGenomesNhaR(self):
        indexes = ["18767", " 18779", " 18812", " 18830", " 18836", " 18839", " 18863", " 18867", " 18875", " 18878", " 18881", " 18893", " 18905", " 18911", " 18971", " 18974", " 18989", " 18995", " 19007", " 19013", " 19019", " 19037", " 19046", " 19049", " 19052", " 19055", " 19056", " 19058", " 19067", " 19070", " 19075", " 19081", " 19082", " 19085", " 19091", " 19097", " 19106", " 19109", " 19115", " 19121", " 19122", " 19124", " 19130", " 19133", " 19160", " 19175", " 19178", " 19181", " 19184", " 19196", " 19202", " 19206", " 19210", " 19226", " 19229", " 19241", " 19244", " 19247", " 19252", " 19262", " 19268", " 19271", " 19274", " 19280", " 19286", " 19289", " 19292", " 19293", " 19295", " 19296", " 19298", " 19304", " 19310", " 19313", " 19322", " 19328", " 19332", " 19334", " 19335", " 19337", " 19358", " 19361", " 19367", " 19373", " 19374", " 19376", " 19379", " 19394", " 19397", " 19409", " 19415", " 19420", " 19424", " 19427", " 19430", " 19433", " 19442", " 19445", " 19446", " 19454", " 19457", " 19458", " 19463", " 19466", " 19470", " 19471", " 19476", " 19478", " 19479", " 19481", " 19484", " 19487", " 19490", " 19493", " 19496", " 19499", " 19502", " 19508", " 19526", " 19529", " 19538", " 19541", " 19542", "19612"]
        scaffoldGenesFromAllFiles = {"scaffold_1466_SS_221.fasta.gbk":{"nhaR": Gene(5872, 6772, 'ATGTCTCATATCAATTACAACCACTTGTATTACTTCTGGCATGTCTACAAAGAAGGTTCTGTGGTTGGCGCAGCGGAGGCGCTTTATTTAACACCACAAACCATTACCGGGCAGATCCGGGCGCTGGAAGAGCGCCTGCAAGGGAAACTATTTAAGCGTAAAGGACGTGGTCTGGAACCCAGCGAACTGGGGGAACTGGTCTATCGCTATGCCGATAAAATGTTCACCTTAAGCCAGGAAATGCTGGATATCGTCAACTATCGCAAAGAGTCCAACTTATTGTTTGATGTTGGTGTGGCAGATGCACTTTCCAAACGTCTGGTCAGCAGTGTTCTGGATGCCGCAGTTGTGGAAGACGAGCAGATCCATCTACGCTGTTTCGAATCGACGCACGAGATGCTTTTAGAGCAGTTGAGTCAGCATAAACTGGATATGATCATCTCTGACTGTCCGATCGATTCCACTCAGCAGGAAGGGCTGTTTTCCATGAAAATTGGCGAATGTGGTGTCAGTTTCTGGTGCACTAACCCACTACCAGAAAAGCCGTTTCCTGCCTGTCTTGAAGAGCGTCGTTTACTTATTCCGGGGCGTCGCTCAATGTTGGGGCGTAAACTATTAAACTGGTTTAACTCCCAGGGCTTGAACGTCGAAATTTTGGGTGAGTTTGATGATGCTGCGTTGATGAAAGCCTTTGGGGCGACGCATAACGCTATTTTCGTTGCACCTTCGCTTTACGCTAATGATTTCTATAACGATGACTCGGTTGTGGAGATAGGCCGTGTTGAGAACGTGATGGAAGAGTACCACGCGATTTTTGCCGAAAGGATGATTCAGCACCCTGCAGTACAGCGTATCTGCAATACAGACTATTCTGCGCTGTTTACTCCAGCTTCAAAATAA',
                                 "nhaR", "product",True)}}
        refGenes = [Gene(18714,19620, 'ATGAGCATGTCTCATATCAATTACAACCACTTGTATTACTTCTGGCATGTCTATAAAGAAGGTTCCGTGGTTGGCGCAGCGGAGGCGCTTTATTTAACTCCACAAACCATTACCGGACAGATTCGAGCGCTGGAAGAGCGCCTGCAAGGCAAATTATTTAAACGCAAGGGACGTGGTCTCGAACCCAGCGAGCTGGGAGAACTGGTCTATCGCTATGCCGATAAAATGTTCACCTTAAGCCAGGAAATGCTGGATATTGTGAACTATCGCAAAGAATCCAATTTATTGTTTGACGTTGGCGTGGCTGATGCACTTTCCAAACGCCTGGTCAGTAGCGTACTTAACGCCGCAGTGGTAGAAGGCGAGCCCATTCATCTTCGCTGCTTCGAATCCACCCACGAAATGCTGCTGGAGCAATTAAGTCAGCATAAACTGGATATGATCATTTCTGACTGTCCGATAGACTCTACGCAGCAGGAAGGCCTGTTCTCCGTGAGAATTGGCGAATGTGGCGTGAGTTTCTGGTGTACAAATCCACCACCAGAAAAACCGTTCCCGGCTTGTCTGGAAGAACGGCGACTTTTGATTCCTGGGCGACGTTCAATGTTAGGGCGCAAATTGCTTAACTGGTTTAACTCCCAGGGATTAAACGTAGAAATCCTCGGCGAGTTTGATGATGCCGCTTTGATGAAAGCTTTTGGTGCGATGCACAATGCAATCTTCGTTGCCCCAACGCTTTATGCATATGACTTTTATGCCGATAAAACTGTCGTAGAAATTGGTCGCGTCGAGAATGTGATGGAAGAGTACCATGCTATTTTTGCTGAGCGGATGATTCAGCACCCGGCGGTACAGCGAATCTGCAATACGGATTATTCTGCGCTTTTTAGTCCAGCGGTGCGTTAA',
                                 "nhaR", "product",True)]
        normalAlignGenes = {"scaffold_1466_SS_221.fasta.gbk":{"nhaR": Gene(18715,19615, 'ATGTCTCATATCAATTACAACCACTTGTATTACTTCTGGCATGTCTACAAAGAAGGTTCTGTGGTTGGCGCAGCGGAGGCGCTTTATTTAACACCACAAACCATTACCGGGCAGATCCGGGCGCTGGAAGAGCGCCTGCAAGGGAAACTATTTAAGCGTAAAGGACGTGGTCTGGAACCCAGCGAACTGGGGGAACTGGTCTATCGCTATGCCGATAAAATGTTCACCTTAAGCCAGGAAATGCTGGATATCGTCAACTATCGCAAAGAGTCCAACTTATTGTTTGATGTTGGTGTGGCAGATGCACTTTCCAAACGTCTGGTCAGCAGTGTTCTGGATGCCGCAGTTGTGGAAGACGAGCAGATCCATCTACGCTGTTTCGAATCGACGCACGAGATGCTTTTAGAGCAGTTGAGTCAGCATAAACTGGATATGATCATCTCTGACTGTCCGATCGATTCCACTCAGCAGGAAGGGCTGTTTTCCATGAAAATTGGCGAATGTGGTGTCAGTTTCTGGTGCACTAACCCACTACCAGAAAAGCCGTTTCCTGCCTGTCTTGAAGAGCGTCGTTTACTTATTCCGGGGCGTCGCTCAATGTTGGGGCGTAAACTATTAAACTGGTTTAACTCCCAGGGCTTGAACGTCGAAATTTTGGGTGAGTTTGATGATGCTGCGTTGATGAAAGCCTTTGGGGCGACGCATAACGCTATTTTCGTTGCACCTTCGCTTTACGCTAATGATTTCTATAACGATGACTCGGTTGTGGAGATAGGCCGTGTTGAGAACGTGATGGAAGAGTACCACGCGATTTTTGCCGAAAGGATGATTCAGCACCCGGCGGTACAGCGAATCTGCAATACGGATTATTCTGCGCTTTTTAGTCCAGCGGTGCGTTAA',
                                 "nhaR", "product",True)}}
        vcfFilePath = "/Users/cazcullimore/dev/ericksonLabCode/tests/testFiles/vcfsFromScaffolds/scaffold_1466_SS_221.fasta.vcf"
        snps = readInSnps(glob(vcfFilePath), refGenomeSeq=getContigs(self.annotatedRefGenomePath)[0], ignoreRefSeq=False)
        self.CompareGenomes(indexes, scaffoldGenesFromAllFiles, refGenes, normalAlignGenes, snps, {})
    def testCompareGenomesApaH(self):
        snpIndexPath = TEST_DATA_DIR + "snpAlign/MultiAlignmentIndexes.txt"
        scaffoldGenesFromAllFiles = {"scaffold_1466_SS_221.fasta.gbk":{"apaH": Gene(36990, 37839, 'ATGGCGACATACCTTATTGGCGACGTTCATGGTTGTTACGATGAACTGATCGCATTGCTGCATAAAGTAGAATTTACCCCTGGGAAAGATACCCTCTGGCTGACGGGCGATCTGGTCGCGCGCGGCCCGGGTTCGCTGGATGTTCTGCGCTATGTGAAATCCTTAGGCGACAGCGTACGTCTGGTGCTGGGTAATCACGATCTGCATCTGCTGGCGGTATTTGCCGGGATCAGCCGCAATAAACCGAAAGATCGCCTGACACCGCTGCTGGAAGCGCCGGATGCCGACGAGCTGCTTAACTGGCTGCGTCGCCAGCCTTTGCTGCAAATCGACGAAGAGAAAAAGTTGGTGATGGCCCACGCCGGGATCACGCCGCAGTGGGATCTGCAGACCGCCAAAGAGTGCGCGCGCGATGTAGAAGCGGTGCTGTCGAGTGACTCCTATCCCTTCTTTCTTGATGCCATGTACGGCGATATGCCAAATAACTGGTCACCGGAATTGCGGGGGCTGGGAAGACTGCGTTTTATCACCAACGCCTTTACCCGTATGCGTTTTTGCTTCCCGAACGGTCAACTGGATATGTACAGCAAAGAATCGCCGGAAGAGGCCCCTGCCCCACTGAAACCGTGGTTTGCGATTCCTGGCCCCGTCGCTGAAGAGTACAACATCGCCTTTGGTCACTGGGCATCGCTGGAAGGCAAAGGTACGCCGGAAGGTATTTACGCGCTGGATACCGGCTGCTGCTGGGGTGGTACATTAACCTGCCTGCGCTGGGAAGATAAACAGTATTTTGTCCAGCCGTCGAACCGGCATAAGGATTTGAGTGAGGGAGAGGCGGTAGCGTCTTAA',
                                 "apaH", "product",False)}}
        refGenes = [Gene(50379,51222, 'ATGGCGACATACCTTATTGGCGACGTTCATGGTTGTTACGATGAACTGATCGCATTGCTGCATAAAGTAGAATTTACCCCTGGGAAAGATACCCTCTGGCTGACGGGCGATCTGGTCGCGCGCGGGCCGGGTTCGCTGGATGTTCTGCGCTATGTGAAATCCTTAGGCGACAGCGTACGTCTGGTGCTGGGCAATCACGATCTGCATCTGCTGGCGGTATTTGCCGGGATCAGCCGCAATAAACCGAAAGATCGCCTGACACCGCTGCTGGAAGCGCCGGATGCCGACGAGCTGCTTAACTGGCTGCGGCGCCAGCCTCTGCTGCAAATCGACGAAGAGAAAAAGCTGGTGATGGCCCACGCAGGGATCACGCCGCAGTGGGATCTGCAGACCGCCAAAGAGTGCGCACGCGATGTAGAAGCGGTGCTATCGAGTGACTCCTATCCCTTCTTTCTTGATGCCATGTACGGCGATATGCCAAATAACTGGTCACCGGAATTGCGGGGGCTGGGAAGACTGCGTTTTATCACCAACGCTTTTACCCGTATGCGTTTTTGCTTCCCGAACGGTCAACTGGATATGTACAGCAAAGAATCGCCGGAAGAGGCCCCTGCCCCACTGAAACCGTGGTTTGCGATTCCTGGCCCTGTCGCTGAAGAATACAGCATCGCCTTTGGTCACTGGGCATCGCTGGAGGGCAAAGGTACGCCGGAAGGTATATACGCGCTGGATACCGGCTGCTGCTGGGGTGGTACATTAACCTGCCTGCGCTGGGAAGATAAACAGTATTTTGTCCAGCCGTCGAACCGGCATAAGGATTTGGGCGAAGCGGCGGCGTCTTAA',
                                 "apaH", "product",False)]
        normalAlignGenes = {"scaffold_1466_SS_221.fasta.gbk":{"apaH": Gene(50377,51220, 'ATGGCGACATACCTTATTGGCGACGTTCATGGTTGTTACGATGAACTGATCGCATTGCTGCATAAAGTAGAATTTACCCCTGGGAAAGATACCCTCTGGCTGACGGGCGATCTGGTCGCGCGCGGCCCGGGTTCGCTGGATGTTCTGCGCTATGTGAAATCCTTAGGCGACAGCGTACGTCTGGTGCTGGGTAATCACGATCTGCATCTGCTGGCGGTATTTGCCGGGATCAGCCGCAATAAACCGAAAGATCGCCTGACACCGCTGCTGGAAGCGCCGGATGCCGACGAGCTGCTTAACTGGCTGCGTCGCCAGCCTTTGCTGCAAATCGACGAAGAGAAAAAGTTGGTGATGGCCCACGCCGGGATCACGCCGCAGTGGGATCTGCAGACCGCCAAAGAGTGCGCGCGCGATGTAGAAGCGGTGCTGTCGAGTGACTCCTATCCCTTCTTTCTTGATGCCATGTACGGCGATATGCCAAATAACTGGTCACCGGAATTGCGGGGGCTGGGAAGACTGCGTTTTATCACCAACGCCTTTACCCGTATGCGTTTTTGCTTCCCGAACGGTCAACTGGATATGTACAGCAAAGAATCGCCGGAAGAGGCCCCTGCCCCACTGAAACCGTGGTTTGCGATTCCTGGCCCCGTCGCTGAAGAGTACAACATCGCCTTTGGTCACTGGGCATCGCTGGAAGGCAAAGGTACGCCGGAAGGTATTTACGCGCTGGATACCGGCTGCTGCTGGGGTGGTACATTAACCTGCCTGCGCTGGGAAGATAAACAGTATTTTGTCCAGCCGTCGAACCGGCATAAGGATTTGGGCGAAGCGGCGGCGTCTTAA',
                                 "apaH", "product",False)}}
        vcfFilePath = TEST_DATA_DIR + "*.vcf"
        snps = readInSnps(glob(vcfFilePath), refGenomeSeq=getContigs(self.annotatedRefGenomePath)[0], ignoreRefSeq=False)
        with open(snpIndexPath) as indexes:
            self.CompareGenomes(indexes, scaffoldGenesFromAllFiles, refGenes, normalAlignGenes, snps, {})
    def testSnpAlignMatchesScaffold(self): #doesn't work because doesn't start from the same place and any indels mess it up
        pathToRefGenomeFasta = DATA_DIR + "refGenomes/k-12.fasta"
        pathToRefGenomeGb = DATA_DIR + "refGenomes/k-12.gbff"
        namePrefix = "snpAlignMatchesScaffold"  # the name to give to start all the files created
        vcfFilePath = TEST_DATA_DIR + "output.vcf"
        scaffoldFastaPath = TEST_DATA_DIR + "scaffold_1465_SS_220.fasta"
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