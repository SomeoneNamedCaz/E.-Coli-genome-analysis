from secondaryPythonScripts.functions import *
from alignVcfSnps import *
from secondaryPythonScripts.makePlinkFiles import *
from megaCatsPythonVersion import *
from parsingMegaCatsResults import *
from reconstructNormalAlignment import *

pathToAssemblies = "./AllAssembliesInOneFolder/*.fasta"
# pathToAssemblies = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/tenAssembliesFromEachCategory/*.fasta"
pathToRefGenomeFasta = "/Users/cazcullimore/Documents/ericksonLabCode/refGenomes/k-12.fasta"
pathToRefGenomeGb = "/Users/cazcullimore/Documents/ericksonLabCode/refGenomes/k-12.gbff"


namePrefix = "fullScript" # the name to give to start all the files created
pathsToGsAlignVCFs = "/".join(pathToAssemblies.split("/")[:-1]) + "/ragtagOutputs/longestScaffoldFiles/gsAlignOutputs/*.vcf"
# pathsToGsAlignVCFs = "./AllAssembliesInOneFolder/carefulAssemblyGsAlign/*.vcf"
# snpAlignPath = "/Users/cazcullimore/Documents/ericksonLabCode/RedoingEverything/SnpAlign/all.afa"
# pathsToGsAlignVCFs = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/vcfsFromScaffolds/*.vcf"
snpAlignPath = "/Users/cazcullimore/Documents/ericksonLabCode/RedoingEverything/SnpAlign/" + namePrefix +".afa"
snpIndexPath = "/Users/cazcullimore/Documents/ericksonLabCode/RedoingEverything/SnpAlign/" + namePrefix +"Indexes.txt"

plinkOutputPath = "/Users/cazcullimore/Documents/ericksonLabCode/RedoingEverything/plinkFiles/"
metadataFilePath = "/Users/cazcullimore/Documents/ericksonLabCode/metaDataForMetaCatsWithExtraMastitis.tsv"
megaCatStatsFilePath = "/Users/cazcullimore/Documents/ericksonLabCode/RedoingEverything/metadataStats.tsv"

normalAlignPath = "/Users/cazcullimore/Documents/ericksonLabCode/RedoingEverything/" + namePrefix + "normalAlign.afa"

reRunScaffoldAndGsAlign = False
reAlignSnps = True
reDoPlinkFiles = True
runPlink = False
runMegaCats = True
reDoMegaCatsStats = True
runNormalAlignment = True

typesOfSnpsToSkipDuringAlignment = []#["INSERT", "DELETE"]

if reRunScaffoldAndGsAlign:
	os.system("cd " + pathToAssemblies + "; conda activate ragtagEnv; mkdir ragtagOutputs;" +
	" for fileName in *.fasta; do ragtag.py scaffold " + pathToRefGenomeFasta +
	" $fileName -t 16 -o ./ragtagOutputs/$fileName/; done;)" +
	          "python /Users/cazcullimore/Documents/ericksonLabCode/getLongestContig.py .; cd ragtagOutputs; cd longestScaffoldFiles; conda activate gsAlign; mkdir gsAlignOutputs; for fileName in *.fasta; do gsAlign -r   /Users/cazcullimore/Documents/ericksonLabCode/refGenomes/k-12.fasta -q $fileName -o ./gsAlignOutputs/$fileName -t 16; done;" +
	          "cd /Users/cazcullimore/Documents/ericksonLabCode/")

if reAlignSnps:
	alignVcfSnps(pathsToGsAlignVCFs, outFilePath=snpAlignPath, thingsToSkip=typesOfSnpsToSkipDuringAlignment, ignoreRefSeq=False)

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
	if reDoMegaCatsStats:
		doMegaCatsStats(alignedFilePath=snpAlignPath, metadataFilePath=metadataFilePath, outFilePath=megaCatStatsFilePath,
		                matchMegacatsStyle=True)
	parseMegaCatsFile(megaCatsFile=megaCatStatsFilePath, snpGenomePath=snpAlignPath, snpIndexesPath=snpIndexPath,
	                  suffix=namePrefix, metaDataFilePath=metadataFilePath, annotatedRefGenomePath=pathToRefGenomeGb,
	                  removeSparce=True, outputDirectory="./RedoingEverything/")
	
if runNormalAlignment:
	reconstructNormalAlignment(snpGenomePath=snpAlignPath, snpIndexesPath=snpIndexPath, refGenomePath=pathToRefGenomeGb, outputPath=normalAlignPath)