"""Because the megacats has a third group for pathogenicity that shouldn't exist,
 I'm going to combine that file with the pathogenicity only file"""

pathogenicityOnlyFilePath = "/Users/cazcullimore/Documents/ericksonLabCode/megaCATS-main/firstMegaCatsRunOnData/megaCatsInsertAndDelete/pathogenicity-rMsaInput.txt-rResultMGCStat.txt"
# can copy into all file (I'm not sure if this is true yet
fileToFixPath = "/Users/cazcullimore/Documents/ericksonLabCode/megaCATS-main/firstMegaCatsRunOnData/megaCatsInsertAndDelete/pathogenicity-rMsaInput.txt-rResultChisqTest.txt" # modified

readData = []
with open(fileToFixPath) as file:
    for line in file:
        readData.append(line)
print(len(readData))
currIndex = 0
with open(pathogenicityOnlyFilePath) as pathOnlyFile, open(fileToFixPath + "Test", "w") as fixFile:
    for pathOnlyLine in pathOnlyFile:
        pathOnlyLine = pathOnlyLine.strip()
        cols = pathOnlyLine.split("\t")
        if len(cols) <= 2 or pathOnlyLine[0] == '"' or cols[2].strip() != "2" or cols[3].strip() != "3":
            # print("skipped")
            continue
        pathOnlySnpIndex = cols[0] # string
        correctPval = cols[1]
        # print("parth index", pathOnlySnpIndex)
        # try:
        for fixLine in readData[currIndex:]:
            fixCols = fixLine.split("\t")
            # print("other inde", fixCols[0])
            currIndex += 1
            if fixCols[0].strip() == pathOnlySnpIndex.strip():
                fixFile.write(pathOnlySnpIndex + "\t" + fixCols[1] + "\t"+ correctPval + "\t" + "\t".join(fixCols[3:]))
                break

        # except:
        #     pass

