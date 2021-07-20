from functions import *

# with open("./substSNPcombinedGenomes/allSubstSNPsMoreThan9ActuallyWorkingMethodIndexes.txt") as file:
#      numConsecutive = 0
#      lastIndex = 0
#      for line in file:
#              line = int(line.strip())
#              if line - 1 == lastIndex:
#                      numConsecutive += 1
#              else:
#                      print(numConsecutive)
#                      numConsecutive = 0
#              lastIndex = line

 
# it looks good

# with open("/Users/cazcullimore/Documents/ericksonLabCode/megaCATS-main/metaDataForMetaCats.tsv") as file:
#      for line in file:
#              name = line.strip().split("\t")[0]
#              print(name)
#              with open("/Users/cazcullimore/Documents/ericksonLabCode/megaCATS-main/allSubstSNPsMoreThan9ActuallyWorkingMethod.afa") as otherFile:
#                      for line2 in otherFile:
#                              if name == line2.strip()[1:]:
#                                      print("found")
#              print("done")

fileNamesInAfa = set()
with open("/Users/cazcullimore/Documents/ericksonLabCode/k-12RefGenomeAnalysisFiles/snpGenomeFiles/allSnps.afa") as otherFile:

                     for line2 in otherFile:
                         if line2[0] == ">":
                            fileNamesInAfa.add(line2.strip()[1:])
                         else:
                            print(len(line2))

                            # print(line2.strip())

fileNamesInMeta = set()
with open("/Users/cazcullimore/Documents/ericksonLabCode/metaDataForMetaCatsBackUp.tsv") as file:
     for line in file:
         fileNamesInMeta.add(line.strip().split("\t")[0])
         # print(line.strip().split("\t")[0])

print(fileNamesInMeta.difference(fileNamesInAfa))
print(fileNamesInAfa.difference(fileNamesInMeta))