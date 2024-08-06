import os.path
import re
import sys
# MIKFSATLLATLIAASVNAATVDLRIMETTDLHSNMMDFDYYKDTATEKFGLVRTASLINDARNEVKNSVLVDNGDLIQGSPLADYISAKGLKAGDVHPVYKALNTLDYTVGTLGNHEFNYGLDYLKNALAGAKFPYVNANVIDARTKQPMFTPYLIKDTEVVDKDGKKQTLKIGYIGVPPQIMGWDKANLSGKVTVNDITETVRKYVPEMREKGADVVVVLAHSGLSADPYKVMAENSVYYLSEIPGVNAIMFGHAHAVFPGKDFADIEGADIAKGTLNGVPAVMPGMWGDHLGVVDLQLSNNSGKWQVTQAKAEARPIYDIANKKSLAAEDSKLVETLKADHDATRQFVSKPIGKSADNMYSYLALVQDDPTVQVVNNAQKAYVEHYIQGDPDLAKLPVLSAAAPFKVGGRKNDPASYVEVEKGQLTFRNAADLYLYPNTLIVVKASGKEVKEWLECSAGQFNQIDPDNTKPQSLINWDGFRTYNFDVIDGVNYQIDVTQPARYDGECQMVNANAERIKNLTFNGKPIDPNAMFLVATNNYRAYGGKFAGTGDSHIAFASPDENRSVLAAWIADESKRAGEIHPAADNNWRLAPIAGDKKLDIRFETSPSDKAAAFIKEKGQYPMNKVATDDIGFAIYQVDLSK
try:
	from .alignVcfSnps import *
	from .megaCatsPythonVersion import *
	from .parsingMegaCatsResults import *
	from .reconstructNormalAlignment import *
	from .divideGenomesByPhylogroup import *
	from .mapSNPsToGenome import *
except ImportError:
	from alignVcfSnps import *
	from megaCatsPythonVersion import *
	from parsingMegaCatsResults import *
	from reconstructNormalAlignment import *
	from divideGenomesByPhylogroup import *
	from mapSNPsToGenome import *
import multiprocessing as mp
import subprocess
from time import time
def main():
	args = parseArgs()
	runGWASOnAssemblies(args)
	

skipProkkaFlag = "noprokka"
skipStatsFlag = "nostats"
runPhylogroupFlag = "phylo"
runVisualizeSnpsFlag = "mapsnps"
readAnnotationFlag = "annotations"


def printHelp():
	print(
		"usage: --ref <reference genbank file> --metadata <tsv with metadata categories> --fastas <fasta genome assembly files>\n"
		"optional: --name: add a name to the output file names\n"
		"          --outdir: specify the output directory\n"
		"          --" + skipProkkaFlag + ": don't annotate genomes. Annotation will allow snpalign to make sure  that the snp found isn't a result of a large deletion or gap in the fasta\n"
		"          --" + readAnnotationFlag + ": read annotation files, must be genbank files (i.e. .gb, .gbk, .gbff). This replaces the prokka annotation step\n"
		"          --" + skipStatsFlag + ": don't run the chi-squared stats just make the snp alignment\n"
		"          --" + runPhylogroupFlag + ": divide genomes by clermont phylogroup and run statistics for each phylogroup\n"
		"          --" + runVisualizeSnpsFlag + ": map snps to the reference genome to see where they line up relative to the genes in the reference genome\n"
	
	)

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
				
	def setFlag(flagName, newValue,default):
		if "--" + flagName in sys.argv:
			sys.argv.pop(sys.argv.index("--" + flagName))
			return newValue
		else:
			return default

	args["name"] = getNextArgAndPop("--name")
	if args["name"] is None:
		args["name"] = ""
	
	args["ref"] = getNextArgAndPop("--ref")
	if args["ref"] is None:
		print("no annotated genbank file given for the reference sequence.")
		printHelp()
		exit(1)
	
	args[readAnnotationFlag] = [getNextArgAndPop("--" + readAnnotationFlag)] # TODO: allow for just giving the directory, and also for fastas
	
	idx = 0
	if args[readAnnotationFlag][0] is None:
		pass
	elif os.path.isdir(args[readAnnotationFlag][0]): # test this
		args[readAnnotationFlag] = glob(os.path.join(args[readAnnotationFlag][0], "*.gb*"))
	else:
		for arg in sys.argv:
			if not "--" in arg:
				idx += 1
				args[readAnnotationFlag].append(arg)
			else:
				break
	sys.argv = sys.argv[idx:]
	
	args["outdir"] = getNextArgAndPop("--outdir")
	if args["outdir"] is None:
		args["outdir"] = "./"
		
	args["metadata"] = getNextArgAndPop("--metadata")
	if args["metadata"] is None:
		print("no metadata found")
		printHelp()
		exit(2)

	args[skipProkkaFlag] = setFlag(skipProkkaFlag, True, False) # don't check for large deletions that aren't aligned (aka missing genes)
	args[skipStatsFlag] = setFlag(skipStatsFlag, True, False)
	args[runPhylogroupFlag] = setFlag(runPhylogroupFlag, True, False)
	args[runVisualizeSnpsFlag] = setFlag(runVisualizeSnpsFlag, True, False)
	
	args["threads"] = getNextArgAndPop("--threads")
	try:
		int(args["threads"])
	except TypeError:
		print("threads not given, using 16 threads")
		args["threads"] = '16'
	except ValueError:
		print("threads couldn't be converted to an integer, using 16 threads")
		args["threads"] = '16'
	
	
	
	# read in fastas
	firstAssembly = getNextArgAndPop("--fastas")
	for arg in reversed(sys.argv): # remove unknown args
		if "--" == arg[:2]:
			print('unknown argument:',arg, "ignoring")
			sys.argv.remove(arg)
	
	if firstAssembly is None:
		args["fastas"] = sys.argv[1:] # don't include file name in sys.argv
	else:
		args["fastas"] = [firstAssembly] + sys.argv[1:]
		
	if len(args["fastas"]) == 1:
		if '*' in args["fastas"][0]:
			args["fastas"] = glob(args["fastas"][0]) # if it was in quotes and so it didn't expand automatically
		elif os.path.isdir(args["fastas"][0]):
			# .fasta,.fas,.fa,.fna,.ffn,.faa,.mpfa,.frn
			args["fastas"] = glob(os.path.join(args["fastas"][0], "*.fasta")) + glob(os.path.join(args["fastas"][0], "*.fa")) #+ glob(os.path.join(args["fastas"][0], "*.mpfa"))
		else:
			print("hi")
		
	return args
	


def runGWASOnAssemblies(args):
	threwError = False
	try:
		pathToRefGenomeGb = args["ref"]
		pathToAssemblies = args["fastas"]

		namePrefix = args["name"]
		skipProkka = args[skipProkkaFlag]
		metadataFilePath = args["metadata"]
		threads = int(args["threads"])
		
		if args["outdir"][-1] != "/": # TODO: test for windows compatability
			args["outdir"] += "/"
			
		if not os.path.exists(args["outdir"]):
			os.makedirs(args["outdir"])
	
		# make ref genome fasta
		pathToRefGenomeFasta = re.sub("\..+$","",pathToRefGenomeGb.split("/")[-1]) + str(time())[-5:] + ".fasta"
		with open(pathToRefGenomeFasta, "w") as fastaFile:
			i = 0
			for contig in getContigs(pathToRefGenomeGb):
				fastaFile.write(">refContig"+ str(i) + "\n")
				fastaFile.write(contig + "\n")
				if i == 1:
					print("WARNING: ref genome has more than one contig, some unexpected behavior might occur") #TODO: implement support for ref genome files with more than one contig
				i += 1
		
		pathOfAnnotatedScaffolds = args["outdir"] + "/annotations/*.gbk" # automatically generated
		if not args[readAnnotationFlag] is None:
			pathOfAnnotatedScaffolds = args[readAnnotationFlag]
		
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
		
		redoSingleAlignment = False
		reAlignSnps = False
		reRunEzclermont = args[runPhylogroupFlag]
		divideByPhylogroup = args[runPhylogroupFlag]
		runMegaCats = not args[skipStatsFlag]
		reDoMegaCatsStats = False
		runNormalAlignment = False
		
		typesOfSnpsToSkipDuringAlignment = []#["INSERT", "DELETE"]
		pool = mp.Pool(threads)
		
		commandsToRun = []
		if redoSingleAlignment:
			# os.system("cd " + assemblyDir + "; mkdir ragtagOutputs;" +
			# " for fileName in *.fasta; do conda run ragtag.py scaffold '" + pathToRefGenomeFasta +
			# "' $fileName -t 16 -o ./ragtagOutputs/$fileName/; done; " +
			#           "python /Users/cazcullimore/Documents/ericksonLabCode/getLongestContig.py .; cd ragtagOutputs; cd longestScaffoldFiles; mkdir gsAlignOutputs; for fileName in *.fasta; do conda run gsAlign -r '" + pathToRefGenomeFasta + "' -q $fileName -t 16; mv *.vcf ./gsAlignOutputs/$fileName.vcf; done;")

			# or maybe just run prokka parallel
			# gsalign os.system command
			doRagtag = False
			t1 = time()
			if doRagtag:
				scaffoldDir = args["outdir"] + "scaffolds"
				runRagtag(pathToAssemblies, pathToRefGenomeFasta, args, pool,scaffoldDir)
			
			subprocess.run(["mkdir", args["outdir"] + "gsAlignOutputs"])
			indexPrefix = re.sub("\..+","",pathToRefGenomeFasta.split("/")[-1])
			subprocess.run(["conda", "run", "bwt_index", pathToRefGenomeFasta, indexPrefix])
			print("took", time() - t1, "to run ragtag")
			if not skipProkka and args[readAnnotationFlag] is None:
				subprocess.run(["mkdir", args["outdir"] + "annotations"])
			listToAlign = pathToAssemblies
			if doRagtag:
				listToAlign = glob(os.path.join(scaffoldDir,"*fasta"))
			print("will align: ", listToAlign)
			for assemblyPath in listToAlign:
				commandArgs = ["conda", "run", "gsAlign", "-i", indexPrefix, "-q", assemblyPath, "-o", re.sub("\..+","",os.path.basename(assemblyPath))]
				commandsToRun.append(commandArgs) # removing threads for now since threads will be part of the pool "-t", threads,
				if not skipProkka and args[readAnnotationFlag] is None:
					commandsToRun.append(["conda", "run", "prokka",
					               assemblyPath,
					                "--cpus", "1", "--force", "--centre", "X", "--compliant", "--prefix",
					                re.sub("\..+", "", os.path.basename(assemblyPath)), "--outdir", args["outdir"] + "annotations/"])
		if len(commandsToRun) != 0:
			t1 = time()
			print("started gsalign and prokka")
			pool.map(silentSubprocessRun, commandsToRun)
			print("took", time()-t1, "to run prokka and/or gsalign")
			subprocess.run(["mv"] + glob("*.vcf") + [os.path.join(args["outdir"], "gsAlignOutputs")])
			
		if reRunEzclermont:  # can do this parallel too
			commandsToRun = []
			for assemblyPath in pathToAssemblies:
				commandsToRun.append(["conda", "run", "ezclermont", assemblyPath])#.split(" "))
			t1 = time()
			outputs = pool.map(silentSubprocessRun, commandsToRun)
			print("took", time() - t1, "to run ezclermont")
			with open(genomeNameToPhylogroupPath, "w") as file:
				for outputProc in outputs:
					file.write(outputProc.stdout.decode().strip() + "\n")

		if reAlignSnps:
			print("aligning snps")
			t1 = time()
			alignVcfSnps(pathsToGsAlignVCFs, outFilePath=snpAlignPath, thingsToSkip=typesOfSnpsToSkipDuringAlignment, ignoreRefSeq=False, refSeqPath=pathToRefGenomeFasta, numSnpsRequired=1)
			print("took", time() - t1, "to run align snps")
		
		if divideByPhylogroup:
			print("dividing genomes into phylogroups")
			t1 = time()
			divideGenomesByPhylogroup(snpAlignPath, phylogroupsPath=genomeNameToPhylogroupPath, metadataPath=metadataFilePath, outDir=phylogroupSnpAlignPrefix)
			print("took", time() - t1, "to run divide phylogroups")
		
		if runMegaCats:
			if reDoMegaCatsStats:
				print("calculating chi-squared stats")
				t1 = time()
				# import cProfile, pstats, io
				# from pstats import SortKey
				# pr = cProfile.Profile()
				# pr.enable()
				
				calculateMegaCatsStats(alignedFilePath=snpAlignPath, metadataFilePath=metadataFilePath, outFilePath=megaCatStatsFilePath,
				                       matchMegacatsStyle=True)
				# pr.disable()
				# s = io.StringIO()
				# sortby = SortKey.CUMULATIVE
				# ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
				# ps.print_stats()
				# print(s.getvalue())
				
				for phylogroupSnpFile in glob(phylogroupSnpAlignPrefix + "Phylogroup*.afa"):
					fileName = phylogroupSnpFile.split("/")[-1]
					megaCatStatsPhylogroupPath: str = "/".join(megaCatStatsFilePath.split("/")[:-1]) + "/" + re.sub("\.afa","",fileName) + ".tsv"
					
					try:
						calculateMegaCatsStats(alignedFilePath=phylogroupSnpFile, metadataFilePath=metadataFilePath,
						                       outFilePath=megaCatStatsPhylogroupPath,
						                       matchMegacatsStyle=True)
					except ZeroDivisionError:
						print("no genomes found in phylogroup", phylogroupSnpFile)
				print("took", time() - t1, "to run calc stats")
			print("mapping chi-squared stats to genes")
			t1 = time()
			parseMegaCatsFile(megaCatsFile=megaCatStatsFilePath, snpGenomePath=snpAlignPath, snpIndexesPath=snpIndexPath,
			                  suffix=namePrefix, metaDataFilePath=metadataFilePath, annotatedRefGenomePath=pathToRefGenomeGb,
			                  removeSparce=True, outputDirectory=parsedMegaCatsFilePath,
			                  pathOfAnnotatedScaffolds=pathOfAnnotatedScaffolds, ignoreAnnotations=skipProkka, numThreads=threads)
			
			for phylogroupSnpFile in glob(phylogroupSnpAlignPrefix + "Phylogroup*.afa"):
				fileName = re.sub("\.afa", "",phylogroupSnpFile.split("/")[-1])
				megaCatStatsPhylogroupPath: str = "/".join(megaCatStatsFilePath.split("/")[:-1]) + "/" + re.sub("\.afa","",fileName) + ".tsv"
				print(megaCatStatsPhylogroupPath)
				try:
					parseMegaCatsFile(megaCatsFile=megaCatStatsPhylogroupPath, snpGenomePath=phylogroupSnpFile, snpIndexesPath=snpIndexPath,
					                  suffix=namePrefix + fileName, metaDataFilePath=metadataFilePath,
					                  annotatedRefGenomePath=pathToRefGenomeGb,
					                  removeSparce=True, outputDirectory=parsedMegaCatsFilePath,
					                  pathOfAnnotatedScaffolds=pathOfAnnotatedScaffolds, ignoreAnnotations=skipProkka, numThreads=threads)
				except Exception as e:
					print(e)
					print("no genomes in phylogroup from:", fileName)
			print("took", time() - t1, "to run parse megacats")
		
		if args[runVisualizeSnpsFlag]:
			# namePrefix = "M12RefGenome"  # "M12RefGenome"
			# if namePrefix == "K12RefGenome":
			# 	refGenomePath = DATA_DIR + "refGenomes/k-12.gbff"
			# elif namePrefix == "M12RefGenome":
			# 	refGenomePath = DATA_DIR + "refGenomes/M12.gbk"
			# snpAlignPath = DATA_DIR + "RedoingEverything/SnpAlign/" + namePrefix + ".afa"
			# snpIndexPath = DATA_DIR + "RedoingEverything/SnpAlign/" + namePrefix + "Indexes.txt"
			# sigSnpsPath = DATA_DIR + "RedoingEverything/snpsSortedBySignificanceWithGenesContainingThem" + namePrefix + "Pathogenicity.tsv"
			# if namePrefix == "K12RefGenome":
			# 	outFilePath = "k12RefGenomeMappedSNPs.txt"
			# elif namePrefix == "M12RefGenome":
			# 	outFilePath = "M12RefGenomeMappedSNPsAllGenomes.txt"
			
			writeMappedSnps(namePrefix=namePrefix, refGenomePath=pathToRefGenomeGb, snpAlignPath=snpAlignPath, snpIndexPath=snpIndexPath, sigSnpsPath=parsedMegaCatsFilePath,metadataPath=metadataFilePath,outDir=args["outdir"], debug=False,numThreads=threads)
		
		if runNormalAlignment:
			reconstructNormalAlignment(snpGenomePath=snpAlignPath, snpIndexesPath=snpIndexPath, refGenomePath=pathToRefGenomeGb, outputPath=normalAlignPath,numThreads=threads)
	except KeyboardInterrupt or Exception as e:
		# global e
		threwError = True
	# cleanup
	pool.close()
	pool.terminate()
	pool.join()
	# remove extra gsalign stuff
	for assembly in pathToAssemblies:
		prefix = re.sub("\..+", "", os.path.basename(assembly))
		for file in glob(prefix + ".*f"):
			if ".fasta" == file.split(".")[-1]:
				continue  # don't remove fastas
			os.remove(file)
	print("files to remove", glob(re.sub('\.fasta$', ".*", pathToRefGenomeFasta)))
	for file in glob(re.sub('\.fasta$', ".*", pathToRefGenomeFasta)):
		print(file)
		os.remove(file)
	if threwError: # pass along error after clean up
		raise e
def silentSubprocessRun(argsToRun):
	print("started:", argsToRun)
	return subprocess.run(argsToRun, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)

def runRagtag(pathToAssemblies, pathToRefGenomeFasta, args, pool, scaffoldDir):
	ragtagCommands = []
	for assemblyPath in pathToAssemblies:
		ragtagCommands.append(["conda", "run", "ragtag.py", "scaffold", pathToRefGenomeFasta, assemblyPath, "-o",
		                       args["outdir"] + re.sub("\..+", "", os.path.basename(assemblyPath))])
	
	pool.map(silentSubprocessRun, ragtagCommands)
	
	subprocess.run(["mkdir", scaffoldDir])
	for walkTuple in os.walk(args["outdir"]):
		dir = walkTuple[0]
		if dir in {args["outdir"], args["outdir"] + "scaffolds", args["outdir"] + "gsAlignOutputs"}:
			continue
		os.rename(os.path.join(dir, "ragtag.scaffold.fasta"),
		          os.path.join(scaffoldDir, os.path.basename(dir) + ".fasta"))
		for file in walkTuple[-1]:
			if file.split(".")[-1] == "fasta":
				continue
			os.remove(os.path.join(dir, file))
		os.rmdir(dir)

if __name__ == "__main__":
	main()