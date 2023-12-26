from functions import *


def makeMLData(snpAlign, metadata):
	mlData = []
	for genomeName in snpAlign.keys():
		if genomeName == "genomeWithoutAnySnps":
			continue
		mlData.append(",".join([genomeName] + list(snpAlign[genomeName]) + metadata[genomeName]))
	return mlData


if __name__ == "__main__":
	metadataPath = DATA_DIR + "metaDataForMetaCatsWithExtraMastitis.tsv"
	snpAlignPath = DATA_DIR + "RedoingEverything/SnpAlign/50jM.afa"
	outFilePath = "/Users/cazcullimore/Documents/mlData.csv"
	snpAlign = readInFastaAsDict(snpAlignPath)
	metadata = readMetaDataAsDict(metadataPath)
	for key in metadata.keys():
		metadata[key] = metadata[key][1:]
	data = makeMLData(snpAlign, metadata)
	with open(outFilePath, "w") as outFile:
		for item in data:
			outFile.write(item + "\n")
		