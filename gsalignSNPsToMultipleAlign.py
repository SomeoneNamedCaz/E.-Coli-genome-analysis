"""
NOTE: this isn't working yet
"""

from functions import *
gsAlignPaths = "./allMastitisAssemblies/ragtagOutputs/longestScaffoldFiles/gsAlignOutputs/*.vcf"
snps = [] # list of elements like this: [fileName, location, oldNuc, NewNuc]
snpPositions = set()
print(snpPositions)
numFiles = 0
for filePath in glob(gsAlignPaths):
    numFiles += 1
    with open(filePath) as file:
        for line in file:
            line = line.strip()
            cols = line.split("\t")
            if line == "" or line[0] == "#":
                continue
            snps.append([filePath, cols[1], cols[3], cols[4]])
            snpPositions.add(cols[1])

# for snp in snps:
#     print(snp)
print(len(snpPositions))

snps.sort(key=lambda a: a[1])
indexOfLastCheckedThing = 0
combinedSnps = {} # position -> [list of snp entries]
for position in snpPositions:

    for snp in snps[indexOfLastCheckedThing:]:
        if snp[1] == position:
            if position in combinedSnps.keys():
                combinedSnps[position].append(snp)
                combinedSnps[position].sort(key= lambda a: a[0])
            else:
                combinedSnps[position] = [snp]
        elif snp[1] > position:
            break
        indexOfLastCheckedThing += 1


for key in combinedSnps.keys():
    value = combinedSnps[key]
    oldNuc = value[0][2]
    for filePath in glob(gsAlignPaths):
        foundRecord = False
        for snpRecord in value:
            if snpRecord[0] == filePath:
                foundRecord = True
                break
            if snpRecord[0] > filePath:
                break
        if not foundRecord:
            combinedSnps[key].append([filePath, key, oldNuc, oldNuc])
            combinedSnps[key].sort(key=lambda a: a[0])

print(combinedSnps)