"""needs index path"""

from functions import *

if len(sys.argv) < 5:
    print("provide a path of the indexes of the combined snps, the indexes of the indels, the combinedStats,"
          "and the output file's suffix")
    exit(1)

allIndexesPath = sys.argv[1]
indelIndexesPath = sys.argv[2]
allStatsPath = sys.argv[3]
indelOutFilePath = "1AA_Summary_Results_Indel" + sys.argv[4] + ".tsv.txt"
# substOutFilePath = ""


allIndexes = []
indelIndexes = set()
with open(allIndexesPath) as file:
    for line in file:
        allIndexes.append(int(line.strip()))

with open(indelIndexesPath) as file:
    for line in file: # special stuff here to adjust for float method that was added
        line = float(line.strip())
        line = round((line - int(line))*1000 + int(line))
        indelIndexes.add(line)

print(len(allIndexes))

lineNumInOutFile = 1
with open(allStatsPath) as file, open(outFilePath, "w") as outFile:
    for line in file:
        cols = line.split("\t")
        if line == "" or line[0] == '"':
            continue
        # print(int(cols[0]))
        snpLocation = allIndexes[int(cols[0]) - 1]
        if snpLocation in indelIndexes:
            outFile.write(str(lineNumInOutFile) + "\t" + "\t".join(cols[1:]))
            lineNumInOutFile += 1
            indelIndexes.remove(snpLocation) # so no duplicates

