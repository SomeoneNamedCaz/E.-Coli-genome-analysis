from secondaryPythonScripts.functions import *
from alignVcfSnps import *
from secondaryPythonScripts.makePlinkFiles import *
from megaCatsPythonVersion import *


# pathToAssemblies = "./AllAssembliesInOneFolder/*.fasta"
pathToAssemblies = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/tenAssembliesFromEachCategory/*.fasta"
pathToRefGenome = "/Users/cazcullimore/Documents/ericksonLabCode/refGenomes/k-12.fasta"

namePrefix = "tenFromEachCategory" # the name to give to start all the files created
# pathToGsAlignVCFs = "/".join(pathToAssemblies.split("/")[:-1]) + "/ragtagOutputs/longestScaffoldFiles/gsAlignOutputs/*.vcf"
# pathToGsAlignVCFs = "./AllAssembliesInOneFolder/carefulAssemblyGsAlign/*.vcf"
# snpAlignPath = "/Users/cazcullimore/Documents/ericksonLabCode/RedoingEverything/SnpAlign/all.afa"
pathToGsAlignVCFs = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/vcfsFromScaffolds/*.vcf"
snpAlignPath = "/Users/cazcullimore/Documents/ericksonLabCode/RedoingEverything/SnpAlign/" + namePrefix +".afa"
snpIndexPath = "/Users/cazcullimore/Documents/ericksonLabCode/RedoingEverything/SnpAlign/" + namePrefix +"Indexes.txt"

plinkOutputPath = "/Users/cazcullimore/Documents/ericksonLabCode/RedoingEverything/plinkFiles/"
metadataFilePath = "/Users/cazcullimore/Documents/ericksonLabCode/metaDataForMetaCatsWithExtraMastitis.tsv"
reRunScaffoldAndGsAlign = False
reAlignSnps = False
reDoPlinkFiles = True
runPlink = False
runMegaCats = True

if reRunScaffoldAndGsAlign:
	os.system("cd " + pathToAssemblies + "; conda activate ragtagEnv; mkdir ragtagOutputs;" +
	" for fileName in *.fasta; do ragtag.py scaffold " + pathToRefGenome +
	" $fileName -t 16 -o ./ragtagOutputs/$fileName/; done;)" +
	"python /Users/cazcullimore/Documents/ericksonLabCode/getLongestContig.py .; cd ragtagOutputs; cd longestScaffoldFiles; conda activate gsAlign; mkdir gsAlignOutputs; for fileName in *.fasta; do gsAlign -r   /Users/cazcullimore/Documents/ericksonLabCode/refGenomes/k-12.fasta -q $fileName -o ./gsAlignOutputs/$fileName -t 16; done;" +
	"cd /Users/cazcullimore/Documents/ericksonLabCode/")

if reAlignSnps:
	alignVcfSnps(pathToGsAlignVCFs,outFilePath=snpAlignPath, thingsToSkip=["INSERT", "DELETE"])


pedFilePath = plinkOutputPath + namePrefix + ".ped"
covFilePath = plinkOutputPath + namePrefix + ".cov"
mapFilePath = plinkOutputPath + namePrefix + ".map"

if not os.path.exists(plinkOutputPath):
	os.mkdir(plinkOutputPath)
if reDoPlinkFiles:
	makePedFile(snpGenomeFilePath=snpAlignPath, metaDataFilePath=metadataFilePath, outFilePath=pedFilePath)
	makeCovFile(metaDataFilePath=metadataFilePath, outFilePath=covFilePath, snpGenomeFilePath=snpAlignPath)
	makeMapFile(indexPath=snpIndexPath, outFilePath=mapFilePath)
if runPlink:
	os.system("cd RedoingEverything/plinkStats/; ./plink --chr-set -1 --allow-extra-chr --assoc --ped " + pedFilePath + " --allow-no-sex "#--all-pheno --pheno " + covFilePath
	          + " --map " + mapFilePath)

if runMegaCats:
	pass