from parsingMegaCatsResults import *

# runs megaCats and determines if you have changed the output
megaCatsFile = "./megaCatsTestFile.tsv"
snpGenomePath = "./testAlignment.afa"
snpIndexesPath = "./testAlignmentIndexes.txt"
refGenome = "/Users/cazcullimore/Documents/ericksonLabCode/refGenomes/k-12.gbff"
metaDataFilePath = "./testMetadata.tsv"


parseMegaCatsFile(megaCatsFile, snpGenomePath, snpIndexesPath, "Test", metaDataFilePath, annotatedRefGenomePath)

print("it worked")
