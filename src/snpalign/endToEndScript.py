import os.path
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
	
skipProkkaFlag = "skipprokka"
skipStatsFlag = "skipStats"

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
	
	
	args[skipProkkaFlag] = False # don't check for large deletions that aren't aligned (aka missing genes)
	if "--" + skipProkkaFlag in sys.argv:
		args[skipProkkaFlag] = True
		sys.argv.pop(sys.argv.index("--" + skipProkkaFlag))
	
	try:
		args["threads"] = getNextArgAndPop("--threads")
	except ValueError:
		print("threads couldn't be converted to an integer")
	if args["threads"] is None:
		args["threads"] = '16'
	args["threads"] = args["threads"]
	
	
	args[skipStatsFlag] = False
	if "--" + skipStatsFlag in sys.argv:
		args[skipStatsFlag] = True
		sys.argv.pop(sys.argv.index("--" + skipStatsFlag))
	
	firstAssembly = getNextArgAndPop("--fastas")
	if firstAssembly is None:
		args["fastas"] = sys.argv[1:] # don't include file name in sys.argv
	else:
		args["fastas"] = [firstAssembly] + sys.argv[1:]
		
	if '*' in args["fastas"][0] and len(args["fastas"]) == 1:
		args["fastas"] = glob(args["fastas"][0]) # if it was in quotes and so it didn't expand automatically
		
	return args
	


def runGWASOnAssemblies(args):
	try:
		pathToRefGenomeGb = args["ref"]
		pathToAssemblies = args["fastas"]

		namePrefix = args["name"]
		skipProkka = args[skipProkkaFlag]
		metadataFilePath = args["metadata"]
		threads = args["threads"]
		
		if args["outdir"][-1] != "/": # TODO: test for windows compatability
			args["outdir"] += "/"
			
		if not os.path.exists(args["outdir"]):
			os.makedirs(args["outdir"])
	
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
			raise Exception("your code is broken")
		pathsToGsAlignVCFs = args["outdir"] + "gsAlignOutputs/*.vcf"
		
		snpAlignPath = args["outdir"] + namePrefix +".afa"
		snpIndexPath = args["outdir"] + namePrefix +"Indexes.txt"
		
		genomeNameToPhylogroupPath = args["outdir"] + namePrefix + "Phylogroups.txt"
		phylogroupSnpAlignPrefix = args["outdir"] + namePrefix
		
		megaCatStatsFilePath = args["outdir"] + namePrefix + "megaCatsStats.tsv"
		parsedMegaCatsFilePath = args["outdir"]
		
		normalAlignPath = args["outdir"] + namePrefix + "normalAlign.afa"
		
		redoSingleAlignment = True
		reAlignSnps = True
		reRunEzclermont = True
		divideByPhylogroup = True
		runMegaCats = not args[skipStatsFlag]
		reDoMegaCatsStats = True
		runNormalAlignment = False
		
		typesOfSnpsToSkipDuringAlignment = []#["INSERT", "DELETE"]
		pool = mp.Pool(3)
		
		commandsToRun = []
		if redoSingleAlignment:
			# os.system("cd " + assemblyDir + "; mkdir ragtagOutputs;" +
			# " for fileName in *.fasta; do conda run ragtag.py scaffold '" + pathToRefGenomeFasta +
			# "' $fileName -t 16 -o ./ragtagOutputs/$fileName/; done; " +
			#           "python /Users/cazcullimore/Documents/ericksonLabCode/getLongestContig.py .; cd ragtagOutputs; cd longestScaffoldFiles; mkdir gsAlignOutputs; for fileName in *.fasta; do conda run gsAlign -r '" + pathToRefGenomeFasta + "' -q $fileName -t 16; mv *.vcf ./gsAlignOutputs/$fileName.vcf; done;")
			# TODO: make prokka and gsalign run multiple instances parallel to speed up
			# or maybe just run prokka parallel
			# gsalign os.system command
			doRagtag = False
			if doRagtag:
				commandsToRun.append((#"mkdir " + args["outdir"] + "ragtagOutputs;" +
									 #" for fileName in '" + "' '".join(pathToAssemblies) +
				                     #"'; do conda run ragtag.py scaffold '" + pathToRefGenomeFasta +
									 #"' $fileName -t 16 -o " + args["outdir"] + "ragtagOutputs/$fileName/; done; " +
					"mkdir " + args["outdir"] + "gsAlignOutputs; for fileName in " + args["outdir"] + "ragtagOutputs/$fileName/*.fasta" +
				                     "; do echo $fileName; conda run gsAlign -r '" + pathToRefGenomeFasta + "' -q $fileName -t " + threads +
				                     "; mv *.vcf '" + args["outdir"] + "gsAlignOutputs/$(basename ${fileName}).vcf'; done;").split(" "))
			else:
				# subprocess.run(["mkdir", args["outdir"] + "gsAlignOutputs"])
				commandsToRun.append(
					"mkdir " + args["outdir"] + "gsAlignOutputs; for fileName in '" + "' '".join(pathToAssemblies) +
					"'; do echo $fileName; date; conda run gsAlign -r '" + pathToRefGenomeFasta + "' -q $fileName -t " + threads +
					"; mv *.vcf " + args[
						"outdir"] + "gsAlignOutputs/$(basename ${fileName}).vcf; done;")  # .split(" "))
				# for assemblyPath in pathToAssemblies:
				# 	commandsToRun.append(["conda", "run", "gsAlign", "-r", pathToRefGenomeFasta, "-q", assemblyPath, "-t", threads, "-o", os.path.basename(assemblyPath),])
			
			if not skipProkka:
				commandsToRun.append("mkdir " + args["outdir"] + "annotations; for fileName in '" + "' '".join(pathToAssemblies) + "'; do conda run prokka '$fileName' --cpus " + threads + " --force --centre X --compliant; mv PROKKA*/*.gbk '"  + args["outdir"] + "annotations/$fileName.gbk'; rm -r PROKKA*/; done;")#.split())
		if reRunEzclermont:  # can do this parallel too
			commandsToRun.append("rm " + genomeNameToPhylogroupPath + "; for fileName in '" + "' '".join(pathToAssemblies) +
			                      "'; do conda run ezclermont $fileName >> " + genomeNameToPhylogroupPath + "; done;")#.split(" "))
		if len(commandsToRun) != 0:
			pool.map(os.system, commandsToRun)
			subprocess.run(["mv", "*.vcf",  + args[
				"outdir"] + "gsAlignOutputs/$(basename ${fileName}).vcf; done;"])
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
					megaCatStatsPhylogroupPath: str = "/".join(megaCatStatsFilePath.split("/")[:-1]) + "/" + re.sub("\.afa","",fileName) + ".tsv"
					
					try:
						calculateMegaCatsStats(alignedFilePath=phylogroupSnpFile, metadataFilePath=metadataFilePath,
						                       outFilePath=megaCatStatsPhylogroupPath,
						                       matchMegacatsStyle=True)
					except ZeroDivisionError:
						print("no genomes found in phylogroup", phylogroupSnpFile)
					
			parseMegaCatsFile(megaCatsFile=megaCatStatsFilePath, snpGenomePath=snpAlignPath, snpIndexesPath=snpIndexPath,
			                  suffix=namePrefix, metaDataFilePath=metadataFilePath, annotatedRefGenomePath=pathToRefGenomeGb,
			                  removeSparce=True, outputDirectory=parsedMegaCatsFilePath,
			                  pathOfAnnotatedScaffolds=pathOfAnnotatedScaffolds, ignoreAnnotations=skipProkka)
			
			for phylogroupSnpFile in glob(phylogroupSnpAlignPrefix + "Phylogroup*.afa"):
				fileName = re.sub("\.afa", "",phylogroupSnpFile.split("/")[-1])
				megaCatStatsPhylogroupPath: str = "/".join(megaCatStatsFilePath.split("/")[:-1]) + "/" + re.sub("\.afa","",fileName) + ".tsv"
				print(megaCatStatsPhylogroupPath)
				try:
					parseMegaCatsFile(megaCatsFile=megaCatStatsPhylogroupPath, snpGenomePath=phylogroupSnpFile, snpIndexesPath=snpIndexPath,
					                  suffix=namePrefix + fileName, metaDataFilePath=metadataFilePath,
					                  annotatedRefGenomePath=pathToRefGenomeGb,
					                  removeSparce=True, outputDirectory=parsedMegaCatsFilePath,
					                  pathOfAnnotatedScaffolds=pathOfAnnotatedScaffolds, ignoreAnnotations=skipProkka)
				except Exception as e:
					print(e)
					print("no genomes in phylogroup from:", fileName)
			
		if runNormalAlignment:
			reconstructNormalAlignment(snpGenomePath=snpAlignPath, snpIndexesPath=snpIndexPath, refGenomePath=pathToRefGenomeGb, outputPath=normalAlignPath)
	except Exception as e:
		os.remove(pathToRefGenomeFasta)
		pool.close()
		pool.terminate()
		pool.join()
		raise e
if __name__ == "__main__":
	main()