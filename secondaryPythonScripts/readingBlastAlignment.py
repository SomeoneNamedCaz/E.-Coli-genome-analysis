# reading a blast alignment to see if my aligner works
import re

from functions import *
from megaCatsPythonVersion import *
metadataFilePath = DATA_DIR + "metaDataForMetaCatsWithExtraMastitis.tsv"
metadata = readMetaDataAsDict(metadataFilePath)
outFilePath = "alignedQorBs.afa"
qorBMegacatsPath = "qorBMegacats.txt"

with open(DATA_DIR + "qorBalignment.txt") as file, open(outFilePath,"w") as outFile:
	nameOfGenome = "no name"
	subjectSeq = ""
	querySeq = ""
	for line in file:
		if line[0] == ">":
			if nameOfGenome != "no name":
				outFile.write(">" + nameOfGenome + "\n")
				outFile.write(subjectSeq + "\n")
				subjectSeq = ""
				querySeq = ""
			nameOfGenome = line[1:].strip()
			
		
		elif len(re.findall("Query\s+\d+", line))==1:
			querySeq += re.findall("\d+\s+(\w+)\s+\d+", line)[0]
		elif len(re.findall("Sbjct\s+\d+", line)) == 1:
			subjectSeq += re.findall("\d+\s+(\w+)\s+\d+", line)[0]
	
			

calculateMegaCatsStats(outFilePath, metadataFilePath, qorBMegacatsPath)