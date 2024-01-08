from secondaryPythonScripts.functions import *
from alignVcfSnps import *
from secondaryPythonScripts.makePlinkFiles import *
from megaCatsPythonVersion import *
from parsingMegaCatsResults import *
from reconstructNormalAlignment import *
from secondaryPythonScripts.makeMLFiles import *
from secondaryPythonScripts.getSignificantClusters import *
from secondaryPythonScripts.divideGenomesByPhylogroup import *

pathToAssemblies = DATA_DIR + "AllAssembliesInOneFolder/*.fasta"
# pathToAssemblies = TEST_DATA_DIR + "tenAssembliesFromEachCategory/*.fasta"
pathToRefGenomeFasta = DATA_DIR + "refGenomes/k-12.fasta"
pathToRefGenomeGb = DATA_DIR + "refGenomes/k-12.gbff"


namePrefix = "justMastitis" # the name to give to start all the files created
pathsToGsAlignVCFs = DATA_DIR + "k-12RefGenomeAnalysisFiles/justCowVCFs/*.vcf"

snpAlignPath = DATA_DIR + "RedoingEverything/SnpAlign/" + namePrefix +".afa"
snpIndexPath = DATA_DIR + "RedoingEverything/SnpAlign/" + namePrefix +"Indexes.txt"

genomeNameToPhylogroupPath = "Phylogroups.txt"
phylogroupSnpAlignPrefix = DATA_DIR + "RedoingEverything/SnpAlign/" + namePrefix

plinkOutputPath = DATA_DIR + "RedoingEverything/plinkFiles/"
metadataFilePath = DATA_DIR + "metaDataForMetaCatsWithExtraMastitis.tsv"
megaCatStatsFilePath = DATA_DIR + "RedoingEverything/" + namePrefix + "megaCatsStats.tsv"

normalAlignPath = DATA_DIR + "RedoingEverything/" + namePrefix + "normalAlign.afa"

reRunScaffoldAndGsAlign = False
reAlignSnps = False
reRunPhylogroups = False
divideByPhylogroup = False
reDoPlinkFiles = False
runPlink = False
runMegaCats = True
reDoMegaCatsStats = False
runNormalAlignment = False
runMLAnalysis = True

typesOfSnpsToSkipDuringAlignment = []#["INSERT", "DELETE"]

if reRunScaffoldAndGsAlign:
	os.system("cd " + pathToAssemblies + "; mkdir ragtagOutputs;" +
	" for fileName in *.fasta; do conda run -n ragtagEnv ragtag.py scaffold " + pathToRefGenomeFasta +
	" $fileName -t 16 -o ./ragtagOutputs/$fileName/; done;)" +
	          "python /Users/cazcullimore/Documents/ericksonLabCode/getLongestContig.py .; cd ragtagOutputs; cd longestScaffoldFiles; mkdir gsAlignOutputs; for fileName in *.fasta; do conda run -n gsAlign gsAlign -r " + pathToRefGenomeFasta + " -q $fileName -o ./gsAlignOutputs/$fileName -t 16; done;")

if reAlignSnps:
	alignVcfSnps(pathsToGsAlignVCFs, outFilePath=snpAlignPath, thingsToSkip=typesOfSnpsToSkipDuringAlignment, ignoreRefSeq=False)

pedFilePath = plinkOutputPath + namePrefix + ".ped"
covFilePath = plinkOutputPath + namePrefix + ".cov"
mapFilePath = plinkOutputPath + namePrefix + ".map"

if reRunPhylogroups:
	os.system("rm " + genomeNameToPhylogroupPath + "; for fileName in " + pathToAssemblies + "; do conda run -n ezclermontEnv ezclermont $fileName >> " + genomeNameToPhylogroupPath + "; done;")

if divideByPhylogroup:
	divideGenomesByPhylogroup(snpAlignPath, phylogroupsPath=genomeNameToPhylogroupPath, metadataPath=metadataFilePath, outDir=phylogroupSnpAlignPrefix)


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
		calculateMegaCatsStats(alignedFilePath=snpAlignPath, metadataFilePath=metadataFilePath, outFilePath=megaCatStatsFilePath,
		                       matchMegacatsStyle=True)
		
		for phylogroupSnpFile in glob(phylogroupSnpAlignPrefix + "Phylogroup*.afa"):
			fileName = phylogroupSnpFile.split("/")[-1]
			megaCatStatsPhylogroupPath = "/".join(megaCatStatsFilePath.split("/")[:-1]) + re.sub("\.afa","",fileName) + ".tsv"
			print("megaCatStatsPhylogroupPath",megaCatStatsPhylogroupPath)
			calculateMegaCatsStats(alignedFilePath=phylogroupSnpFile, metadataFilePath=metadataFilePath,
			                       outFilePath=megaCatStatsPhylogroupPath,
			                       matchMegacatsStyle=True)
			
	parseMegaCatsFile(megaCatsFile=megaCatStatsFilePath, snpGenomePath=snpAlignPath, snpIndexesPath=snpIndexPath,
	                  suffix=namePrefix, metaDataFilePath=metadataFilePath, annotatedRefGenomePath=pathToRefGenomeGb,
	                  removeSparce=True, outputDirectory=DATA_DIR + "RedoingEverything")
	
	for phylogroupSnpFile in glob(phylogroupSnpAlignPrefix + "Phylogroup*.afa"):
		fileName = re.sub("\.afa", "",phylogroupSnpFile.split("/")[-1])
		megaCatStatsPhylogroupPath = "/".join(megaCatStatsFilePath.split("/")[:-1]) + fileName + ".tsv"
		try:
			parseMegaCatsFile(megaCatsFile=megaCatStatsPhylogroupPath, snpGenomePath=phylogroupSnpFile, snpIndexesPath=snpIndexPath,
			                  suffix=fileName, metaDataFilePath=metadataFilePath,
			                  annotatedRefGenomePath=pathToRefGenomeGb,
			                  removeSparce=True, outputDirectory=DATA_DIR + "RedoingEverything")
		except Exception:
			print("no data found exception handled")
	
if runNormalAlignment:
	reconstructNormalAlignment(snpGenomePath=snpAlignPath, snpIndexesPath=snpIndexPath, refGenomePath=pathToRefGenomeGb, outputPath=normalAlignPath)
	
if runMLAnalysis:
	snpAlign = readInFastaAsDict(snpAlignPath)
	metadata = readMetaDataAsDict(metadataFilePath)
	for key in metadata.keys():
		metadata[key] = metadata[key][1:]
	mlData = makeMLData(snpAlign, metadata)
	trainRFModel(mlData)