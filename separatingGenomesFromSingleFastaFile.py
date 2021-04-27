from functions import *

def SeparateByLocus(filePath):
        pathLocation = '/'.join(filePath.split("/")[:-1])
        with open(filePath) as fileToSort:
            onFirstLine = True
            lastLocusName = ""
            lastStrainName = ""
            lastContigSeq = ""
            lastHeaderLine = ""
            lastStrainName = ""
            for line in fileToSort:
                if line.strip() == "":
                    continue
                if line[0] == ">":
                    if not onFirstLine:
                        genomeFilePath = pathLocation + "/" + lastStrainName + ".fasta"
                        fileData = []
                        try:
                            with open(genomeFilePath) as fileToWrite:
                                for data in fileToWrite:
                                    fileData.append(data)
                        except FileNotFoundError:
                            0#print("make newFile")
                        with open(genomeFilePath, "w") as fileToWrite:
                            for data in fileData:
                                fileToWrite.write(data)
                            fileToWrite.write(lastHeaderLine)
                            fileToWrite.write(lastContigSeq)
                            lastContigSeq = ""
                            lastHeaderLine = ""
                    lastLocusName = line[1:].split()[0] # first word without ">"
                    if "complete genome" in line:
                        lastStrainName = lastLocusName
                    elif "whole genome shotgun sequence" in line:
                        lastStrainName = lastLocusName[:-5]
                    else:
                        print("error not shotgun or complete")
                        lastStrainName = "weird_" + lastLocusName[:-5]

                    lastHeaderLine = line
                    onFirstLine = False
                else:
                    lastContigSeq += line # no strip to keep formatting



SeparateByLocus("/DownloadingFilesFromNCBI/AllCommensalStrains/assemblies/AllAssembliesInOneFile/sequence (1).fasta")