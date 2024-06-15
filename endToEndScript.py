from secondaryPythonScripts.functions import *
from alignVcfSnps import *
from secondaryPythonScripts.makePlinkFiles import *
from megaCatsPythonVersion import *
from parsingMegaCatsResults import *
from reconstructNormalAlignment import *
from secondaryPythonScripts.makeMLFiles import *
from secondaryPythonScripts.getSignificantClusters import *
from secondaryPythonScripts.divideGenomesByPhylogroup import *
import sys

# assemblyDir = "/".join(pathToAssemblies.split("/")[:-1])
# # pathToAssemblies = TEST_DATA_DIR + "tenAssembliesFromEachCategory/*.fasta"
# pathToRefGenomeFasta = DATA_DIR + "refGenomes/M12 closed genome (1).fasta"
# pathToRefGenomeGb = DATA_DIR + "refGenomes/M12.gbk"
#
#
# namePrefix = "M12RefGenome" # the name to give to start all the files created
# pathsToGsAlignVCFs = assemblyDir + "/ragtagOutputs/longestScaffoldFiles/gsAlignOutputs/*.vcf"
# pathToAssemblies = DATA_DIR + "m-12RefGenomeAnalysisFiles/allBovineScaffolds/*.fasta"

if "--name" in sys.argv:
	try:
		nameIdx = sys.argv.index("--name") + 1
		namePrefix = sys.argv[nameIdx]
		sys.argv.pop(nameIdx)
		sys.argv.pop(nameIdx - 1)
	except IndexError:
		print("name flag given without name")
		exit(1)
else:
	namePrefix = ""
	
if "--ref" in sys.argv:
	idx = sys.argv.index("--ref") + 1
	pathToRefGenomeGb = sys.argv[idx]
	
	sys.argv.pop(idx)
	sys.argv.pop(idx - 1)
else:
	print("no annotated genbank file given for the reference sequence. usage: --ref <path to referenceSequence>")
	exit(1)
# namePrefix = "M12RefGenome"# the name to give to start all the files created

# make ref genome fasta
pathToRefGenomeFasta = re.sub("\..+$","",pathToRefGenomeGb) + ".fasta"
with open(pathToRefGenomeFasta, "w") as fastaFile:
	i = 0
	for contig in getContigs(pathToRefGenomeGb):
		fastaFile.write(">refContig"+ str(i))
		fastaFile.write(contig)

pathToAssemblies = sys.argv[1]
pathOfAnnotatedScaffolds = "/".join(pathToAssemblies.split("/")[:-1]) + "/annotations/*.gbk"


if namePrefix == "K12RefGenome":
	pathOfAnnotatedScaffolds = "/Users/cazcullimore/dev/data/k-12RefGenomeAnalysisFiles/AllAssemblies/allBovineScaffolds/*.gbk"
	pathToAssemblies = DATA_DIR + "k-12RefGenomeAnalysisFilesWithM12/AllAssemblies/allBovineScaffolds/*.fasta"
	# pathToAssemblies = TEST_DATA_DIR + "tenAssembliesFromEachCategory/*.fasta"
	pathToRefGenomeFasta = DATA_DIR + "refGenomes/k-12.fasta"#M12 closed genome (1).fasta"#
	pathToRefGenomeGb = DATA_DIR + "refGenomes/k-12.gbff" #M12.gbk"#
elif namePrefix == "M12RefGenome":
	pathOfAnnotatedScaffolds = "/Users/cazcullimore/dev/data/k-12RefGenomeAnalysisFiles/AllAssemblies/allBovineScaffolds/*.gbk"
	pathToAssemblies = DATA_DIR + "m-12RefGenomeAnalysisFilesWithK12/allBovineScaffolds/*.fasta"
	
	pathToRefGenomeFasta = DATA_DIR + "refGenomes/M12_closed_genome.fasta"
	pathToRefGenomeGb = DATA_DIR + "refGenomes/M12.gbk"

assemblyDir = "/".join(pathToAssemblies.split("/")[:-1])

pathsToGsAlignVCFs = assemblyDir + "/ragtagOutputs/longestScaffoldFiles/gsAlignOutputs/*.vcf"

snpAlignPath = DATA_DIR + "RedoingEverything/SnpAlign/" + namePrefix +".afa"
snpIndexPath = DATA_DIR + "RedoingEverything/SnpAlign/" + namePrefix +"Indexes.txt"

genomeNameToPhylogroupPath = "Phylogroups.txt"
phylogroupSnpAlignPrefix = DATA_DIR + "RedoingEverything/SnpAlign/" + namePrefix

plinkOutputPath = DATA_DIR + "RedoingEverything/plinkFiles/"
metadataFilePath = DATA_DIR + "metaDataForMetaCatsWithExtraMastitis.tsv"
megaCatStatsFilePath = DATA_DIR + "RedoingEverything/" + namePrefix + "megaCatsStats.tsv"

normalAlignPath = DATA_DIR + "RedoingEverything/" + namePrefix + "normalAlign.afa"

reRunScaffoldAndGsAlign = True
reAlignSnps = True
reRunEzclermont = False
divideByPhylogroup = True
reDoPlinkFiles = False
runPlink = False
runMegaCats = True
reDoMegaCatsStats = True
runNormalAlignment = False
runMLAnalysis = False

typesOfSnpsToSkipDuringAlignment = []#["INSERT", "DELETE"]

if reRunScaffoldAndGsAlign:
	# os.system("cd " + assemblyDir + "; mkdir ragtagOutputs;" +
	# " for fileName in *.fasta; do conda run ragtag.py scaffold '" + pathToRefGenomeFasta +
	# "' $fileName -t 16 -o ./ragtagOutputs/$fileName/; done; " +
	#           "python /Users/cazcullimore/Documents/ericksonLabCode/getLongestContig.py .; cd ragtagOutputs; cd longestScaffoldFiles; mkdir gsAlignOutputs; for fileName in *.fasta; do conda run gsAlign -r '" + pathToRefGenomeFasta + "' -q $fileName -t 16; mv *.vcf ./gsAlignOutputs/$fileName.vcf; done;")
	os.system("cd " + assemblyDir + "; mkdir gsAlignOutputs; for fileName in *.fasta; do conda run gsAlign -r '" + pathToRefGenomeFasta + "' -q $fileName -t 16; mv *.vcf ./gsAlignOutputs/$fileName.vcf; done;")
	os.system("mkdir annotations; for fileName in *.fasta; do conda run prokka $fileName --force --centre X --compliant; cp PROKKA*/*.gbk ./annotations/$fileName.gbk; rm -r PROKKA*/; done;")
	if reAlignSnps:
		alignVcfSnps(pathsToGsAlignVCFs, outFilePath=snpAlignPath, thingsToSkip=typesOfSnpsToSkipDuringAlignment, ignoreRefSeq=False, refSeqPath=pathToRefGenomeFasta, numSnpsRequired=1)
	
pedFilePath = plinkOutputPath + namePrefix + ".ped"
covFilePath = plinkOutputPath + namePrefix + ".cov"
mapFilePath = plinkOutputPath + namePrefix + ".map"

if reRunEzclermont:
	os.system("rm " + genomeNameToPhylogroupPath + "; for fileName in " + pathToAssemblies + "; do conda run ezclermont $fileName >> " + genomeNameToPhylogroupPath + "; done;")

if divideByPhylogroup:
	divideGenomesByPhylogroup(snpAlignPath, phylogroupsPath=genomeNameToPhylogroupPath, metadataPath=metadataFilePath, outDir=phylogroupSnpAlignPrefix)


if reDoPlinkFiles:
	if not os.path.exists(plinkOutputPath):
		os.mkdir(plinkOutputPath)
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
			try:
				calculateMegaCatsStats(alignedFilePath=phylogroupSnpFile, metadataFilePath=metadataFilePath,
				                       outFilePath=megaCatStatsPhylogroupPath,
				                       matchMegacatsStyle=True)
			except ZeroDivisionError:
				print("no genomes found in this phylogroup:",phylogroupSnpFile)
			
	parseMegaCatsFile(megaCatsFile=megaCatStatsFilePath, snpGenomePath=snpAlignPath, snpIndexesPath=snpIndexPath,
	                  suffix=namePrefix, metaDataFilePath=metadataFilePath, annotatedRefGenomePath=pathToRefGenomeGb,
	                  removeSparce=True, outputDirectory=DATA_DIR + "RedoingEverything/",
	                  pathOfAnnotatedScaffolds=pathOfAnnotatedScaffolds)
	
	for phylogroupSnpFile in glob(phylogroupSnpAlignPrefix + "Phylogroup*.afa"):
		fileName = re.sub("\.afa", "",phylogroupSnpFile.split("/")[-1])
		megaCatStatsPhylogroupPath = "/".join(megaCatStatsFilePath.split("/")[:-1]) + fileName + ".tsv"
		try:
			parseMegaCatsFile(megaCatsFile=megaCatStatsPhylogroupPath, snpGenomePath=phylogroupSnpFile, snpIndexesPath=snpIndexPath,
			                  suffix=fileName, metaDataFilePath=metadataFilePath,
			                  annotatedRefGenomePath=pathToRefGenomeGb,
			                  removeSparce=True, outputDirectory=DATA_DIR + "RedoingEverything/",
			                  pathOfAnnotatedScaffolds=pathOfAnnotatedScaffolds)
		except Exception:
			print("no-data-found exception handled")
	
if runNormalAlignment:
	reconstructNormalAlignment(snpGenomePath=snpAlignPath, snpIndexesPath=snpIndexPath, refGenomePath=pathToRefGenomeGb, outputPath=normalAlignPath)
	
if runMLAnalysis:
	snpAlign = readInFastaAsDict(snpAlignPath)
	metadata = readMetaDataAsDict(metadataFilePath)
	for key in metadata.keys():
		metadata[key] = metadata[key][1:]
	mlData = makeMLData(snpAlign, metadata)
	trainRFModel(mlData)