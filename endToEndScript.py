from secondaryPythonScripts.functions import *
from alignVcfSnps import *
pathToAssemblies = "./AllAssembliesInOneFolder/*.fasta"
pathToRefGenome = "/Users/cazcullimore/Documents/ericksonLabCode/refGenomes/k-12.fasta"

reRunScaffoldAndGsAlign = False
reAlignSnps = False

if reRunScaffoldAndGsAlign:
	os.system("cd " + pathToAssemblies + "; conda activate ragtagEnv; mkdir ragtagOutputs;" +
	" for fileName in *.fasta; do ragtag.py scaffold " + pathToRefGenome +
	" $fileName -t 16 -o ./ragtagOutputs/$fileName/; done;)" +
	"python /Users/cazcullimore/Documents/ericksonLabCode/getLongestContig.py .; cd ragtagOutputs; cd longestScaffoldFiles; conda activate gsAlign; mkdir gsAlignOutputs; for fileName in *.fasta; do gsAlign -r   /Users/cazcullimore/Documents/ericksonLabCode/refGenomes/k-12.fasta -q $fileName -o ./gsAlignOutputs/$fileName -t 16; done;" +
	"cd /Users/cazcullimore/Documents/ericksonLabCode/")
# pathToGsAlignVCFs = "/".join(pathToAssemblies.split("/")[:-1]) + "/ragtagOutputs/longestScaffoldFiles/gsAlignOutputs/*.vcf"
pathToGsAlignVCFs = "./AllAssembliesInOneFolder/carefulAssemblyGsAlign/*.vcf"
snpAlignPath = "/Users/cazcullimore/Documents/ericksonLabCode/RedoingEverything/SnpAlign/all"
if reAlignSnps:
	alignVcfSnps(pathToGsAlignVCFs,outFilePath=snpAlignPath)



os.system("cd RedoingEverything/plinkStats/; ./plink --assoc")