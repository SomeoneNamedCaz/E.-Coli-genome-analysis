"""Because the megacats has a third group for pathogenicity that shouldn't exist,
 I'm going to combine that file with the pathogenicity only file"""
import sys
import cProfile
from time import time

prof = cProfile.Profile()
prof.enable()

if len(sys.argv) < 3:
    print("please provide the pathogenicity only file path (should have MGCStat in the name) and the path of all of the combined stats (1AA-Summary ...)")


pathogenicityOnlyFilePath = sys.argv[1]#"/Users/cazcullimore/Documents/ericksonLabCode/megaCATS-main/firstMegaCatsRunOnData/megaCatsInsertAndDelete/pathogenicity-rMsaInput.txt-rResultMGCStat.txt"
# can copy into all file (I'm not sure if this is true yet
fileToFixPath = sys.argv[2] #"/Users/cazcullimore/Documents/ericksonLabCode/megaCATS-main/firstMegaCatsRunOnData/megaCatsInsertAndDelete/pathogenicity-rMsaInput.txt-rResultChisqTest.txt" # modified
metaDataToReplace = pathogenicityOnlyFilePath.split("/")[-1].split("-")[0] # first segment of name of file from path
addMetaData = True

readData = []
with open(fileToFixPath) as file:
    for line in file:
        readData.append(line)
print(len(readData))
currIndex = 0
t1 = time()
with open(pathogenicityOnlyFilePath) as pathOnlyFile, open(fileToFixPath + "Test", "w") as fixFile:
    for pathOnlyLine in pathOnlyFile:
        # if time() > t1 + 2*60:
        #     break
        pathOnlyLine = pathOnlyLine.strip()
        cols = pathOnlyLine.split("\t")
        if len(cols) <= 2 or pathOnlyLine[0] == '"' or cols[2].strip() != "2" or cols[3].strip() != "3":
            # try:
            #     print(len(cols) <= 2)
            # except:
            #     pass
            # try:
            #     print(pathOnlyLine[0] == '"')
            # except:
            #     pass
            # try:
            #     print(cols[2].strip() != "2")
            # except:
            #     pass
            # try:
            #     print(cols[3].strip() != "3")
            # except:
            #     pass
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
                stringToWrite = pathOnlySnpIndex + "\t" + fixCols[1] + "\t"+ correctPval + "\t" + "\t".join(fixCols[3:])[:-1]
                if addMetaData:
                    stringToWrite += "\t" + metaDataToReplace
                stringToWrite += "\n"
                fixFile.write(stringToWrite)
                # readData = readData[currIndex:]
                break
            # if currIndex > 20:
            #     print("oops")
            #     break
            # currIndex += 1
        # print(currIndex)
        # except:
        #     pass


prof.disable()
prof.print_stats(sort=1)

"""        for fixLine in readData[currIndex:]:
            fixCols = fixLine.split("\t")
            # print("other inde", fixCols[0])
            if fixCols[0].strip() == pathOnlySnpIndex.strip():
                if fixCols[-1].split("-")[0].strip() == metaDataToReplace:
                    fixFile.write(pathOnlySnpIndex + "\t" + fixCols[1] + "\t" + correctPval + "\t" + "\t".join(fixCols[3:]))
                    currIndex += 1
                else:
                    fixFile.write(fixLine)
                break"""