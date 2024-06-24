import sys

try:
	from .alignVcfSnps import *
	from .megaCatsPythonVersion import *
	from .parsingMegaCatsResults import *
	from .reconstructNormalAlignment import *
	from .divideGenomesByPhylogroup import *
except ImportError:
	from alignVcfSnps import *
	from megaCatsPythonVersion import *
	from parsingMegaCatsResults import *
	from reconstructNormalAlignment import *
	from divideGenomesByPhylogroup import *
import multiprocessing as mp
import subprocess
def main():
	args = parseArgs()
	runGWASOnAssemblies(args)
	
def printHelp():
	print("usage: --ref <reference genbank file> --metadata <tsv with metadata categories> --fastas <fasta genome assembly files>\n"
	      "optional: --name: add a name to the output file names\n"
	      "          --outdir: specify the output directory")
def parseArgs():
	args = {}
	def getNextArgAndPop(flagName):
		if flagName in sys.argv:
			try:
				idx = sys.argv.index(flagName) + 1
				arg = sys.argv[idx]
				
				sys.argv.pop(idx)
				sys.argv.pop(idx - 1)
				# args[re.sub("--","",flagName)] = arg
				return arg
			except IndexError:
				print(flagName + " flag given without name")
				exit(1)

	args["name"] = getNextArgAndPop("--name")
	if args["name"] is None:
		args["name"] = ""
	
	args["ref"] = getNextArgAndPop("--ref")
	if args["ref"] is None:
		print("no annotated genbank file given for the reference sequence.")
		printHelp()
		exit(1)
	
	args["outdir"] = getNextArgAndPop("--outdir")
	if args["outdir"] is None:
		args["outdir"] = "./"
		
	args["metadata"] = getNextArgAndPop("--metadata")
	if args["metadata"] is None:
		print("no metadata found")
		printHelp()
		exit(2)
	
	skipProkkaFlag = "--skipprokka" # don't check for large deletions that aren't aligned (aka missing genes)
	args["skipprokka"] = False
	if skipProkkaFlag in sys.argv:
		args["skipprokka"] = True
		sys.argv.pop(sys.argv.index(skipProkkaFlag))
	
	try:
		args["threads"] = int(getNextArgAndPop("--threads"))
	except ValueError:
		print("threads couldn't be converted to an integer")
	if args["threads"] is None:
		args["threads"] = 16
	
	firstAssembly = getNextArgAndPop("--fastas")
	if firstAssembly is None:
		args["fastas"] = sys.argv[1:] # don't include file name in sys.argv
	else:
		args["fastas"] = [firstAssembly] + sys.argv[1:]
		
	if '*' in args["fastas"][0] and len(args["fastas"]) == 1:
		args["fastas"] = glob(args["fastas"][0]) # if it was in quotes and so it didn't expand automatically
		
	print(args)
	return args
	


def runGWASOnAssemblies(args):
	pathToRefGenomeGb = args["ref"]
	pathToAssemblies = args["fastas"]
	outDir = args["outdir"]
	namePrefix = args["name"]
	skipProkka = args["skipprokka"]
	metadataFilePath = args["metadata"]
	threads = args["threads"]
	
	if args["outdir"][-1] != "/": # TODO: test for windows compatability
		args["outdir"] += "/"

	# make ref genome fasta
	pathToRefGenomeFasta = re.sub("\..+$","",pathToRefGenomeGb.split("/")[-1]) + str(time.time())[-5:] + ".fasta"
	with open(pathToRefGenomeFasta, "w") as fastaFile:
		i = 0
		for contig in getContigs(pathToRefGenomeGb):
			fastaFile.write(">refContig"+ str(i) + "\n")
			fastaFile.write(contig + "\n")
			if i == 1:
				print("WARNING: ref genome has more than one contig, some unexpected behavior might occur") #TODO: implement support for ref genome files with more than one contig
			i += 1
	
	pathOfAnnotatedScaffolds = "/".join(pathToAssemblies[0].split("/")[:-1]) + "/annotations/*.gbk" # automatically generated
	
	
	if type(pathToAssemblies) == str:
		assemblyDir = "/".join(pathToAssemblies.split("/")[:-1])
	elif type(pathToAssemblies) == list:
		assemblyDir = "/".join(pathToAssemblies[0].split("/")[:-1])
	else:
		print("your code is broken")
	pathsToGsAlignVCFs = outDir + "gsAlignOutputs/*.vcf"
	
	snpAlignPath = outDir + namePrefix +".afa"
	snpIndexPath = outDir + namePrefix +"Indexes.txt"
	
	genomeNameToPhylogroupPath = outDir + namePrefix + "Phylogroups.txt"
	phylogroupSnpAlignPrefix = outDir + namePrefix
	
	megaCatStatsFilePath = outDir + namePrefix + "megaCatsStats.tsv"
	parsedMegaCatsFilePath = outDir
	
	normalAlignPath = outDir + namePrefix + "normalAlign.afa"
	
	reRunScaffoldAndGsAlign = True
	reAlignSnps = False
	reRunEzclermont = True
	divideByPhylogroup = True
	runMegaCats = True
	reDoMegaCatsStats = True
	runNormalAlignment = False
	
	typesOfSnpsToSkipDuringAlignment = []#["INSERT", "DELETE"]
	pool = mp.Pool(3)
	
	commandsToRun = []
	if reRunScaffoldAndGsAlign:
		# os.system("cd " + assemblyDir + "; mkdir ragtagOutputs;" +
		# " for fileName in *.fasta; do conda run ragtag.py scaffold '" + pathToRefGenomeFasta +
		# "' $fileName -t 16 -o ./ragtagOutputs/$fileName/; done; " +
		#           "python /Users/cazcullimore/Documents/ericksonLabCode/getLongestContig.py .; cd ragtagOutputs; cd longestScaffoldFiles; mkdir gsAlignOutputs; for fileName in *.fasta; do conda run gsAlign -r '" + pathToRefGenomeFasta + "' -q $fileName -t 16; mv *.vcf ./gsAlignOutputs/$fileName.vcf; done;")
		# TODO: make prokka and gsalign run multiple instances parallel to speed up
		# or maybe just run prokka parallel
		# gsalign os.system command
		commandsToRun.append("mkdir " + outDir + "gsAlignOutputs; for fileName in '" + "' '".join(pathToAssemblies) + "'; do echo $fileName; conda run gsAlign -r '" + pathToRefGenomeFasta + "' -q $fileName -t " + threads + "; mv *.vcf '" + outDir + "gsAlignOutputs/$(basename ${fileName}).vcf'; done;")
		
		if not skipProkka:
			commandsToRun.append("mkdir " + outDir + "annotations; for fileName in '" + "' '".join(pathToAssemblies) + "'; do conda run prokka '$fileName' --cpus " + threads + " --force --centre X --compliant; mv PROKKA*/*.gbk '" + outDir + "annotations/$fileName.gbk'; rm -r PROKKA*/; done;")
	if reRunEzclermont:  # can do this parallel too
		commandsToRun.append("rm " + genomeNameToPhylogroupPath + "; for fileName in '" + "' '".join(pathToAssemblies) \
		                     + "'; do conda run ezclermont $fileName >> " + genomeNameToPhylogroupPath + "; done;")
	if len(commandsToRun) != 0:
		pool.map(os.system, commandsToRun)
	if reAlignSnps:
		alignVcfSnps(pathsToGsAlignVCFs, outFilePath=snpAlignPath, thingsToSkip=typesOfSnpsToSkipDuringAlignment, ignoreRefSeq=False, refSeqPath=pathToRefGenomeFasta, numSnpsRequired=1)
	
	if divideByPhylogroup:
		divideGenomesByPhylogroup(snpAlignPath, phylogroupsPath=genomeNameToPhylogroupPath, metadataPath=metadataFilePath, outDir=phylogroupSnpAlignPrefix)
	
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
					print("no genomes found in phylogroup",phylogroupSnpFile)
				
		parseMegaCatsFile(megaCatsFile=megaCatStatsFilePath, snpGenomePath=snpAlignPath, snpIndexesPath=snpIndexPath,
		                  suffix=namePrefix, metaDataFilePath=metadataFilePath, annotatedRefGenomePath=pathToRefGenomeGb,
		                  removeSparce=True, outputDirectory=parsedMegaCatsFilePath,
		                  pathOfAnnotatedScaffolds=pathOfAnnotatedScaffolds, ignoreAnnotations=skipProkka)
		
		for phylogroupSnpFile in glob(phylogroupSnpAlignPrefix + "Phylogroup*.afa"):
			fileName = re.sub("\.afa", "",phylogroupSnpFile.split("/")[-1])
			megaCatStatsPhylogroupPath = "/".join(megaCatStatsFilePath.split("/")[:-1]) + fileName + ".tsv"
			try:
				parseMegaCatsFile(megaCatsFile=megaCatStatsPhylogroupPath, snpGenomePath=phylogroupSnpFile, snpIndexesPath=snpIndexPath,
				                  suffix=namePrefix + fileName, metaDataFilePath=metadataFilePath,
				                  annotatedRefGenomePath=pathToRefGenomeGb,
				                  removeSparce=True, outputDirectory=parsedMegaCatsFilePath,
				                  pathOfAnnotatedScaffolds=pathOfAnnotatedScaffolds, ignoreAnnotations=skipProkka)
			except Exception:
				print("no genomes in phylogroup from:",fileName)
		
	if runNormalAlignment:
		reconstructNormalAlignment(snpGenomePath=snpAlignPath, snpIndexesPath=snpIndexPath, refGenomePath=pathToRefGenomeGb, outputPath=normalAlignPath)
	
if __name__ == "__main__":
	main()