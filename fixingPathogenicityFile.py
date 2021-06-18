"""Because the megacats has a third group for pathogenicity that shouldn't exist,
 I'm going to combine that file with the pathogenicity only file"""

pathogenicityOnlyFilePath = "/Users/cazcullimore/Documents/ericksonLabCode/megaCATS-main/firstMegaCatsRunOnData/pathogenicity-rMsaInput.txt-rResultMGCStat.txt"
# can copy into all file
fileToFixPath = "Users/cazcullimore/Documents/ericksonLabCode/megaCATS-main/firstMegaCatsRunOnData/pathogenicity-rMsaInput.txt-rResultChisqTest.txt" # modified

readData = []
with open(fileToFixPath) as file:
    for line in file:
        readData.append(line)

with open(pathogenicityOnlyFilePath) as pathOnlyFile, open(fileToFixPath, "w") as fixFile:
    for line in pathOnlyFile:
