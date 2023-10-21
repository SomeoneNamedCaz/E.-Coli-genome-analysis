from functions import *
# counts = {} # bioprojectID
# to count
# for fileName in glob("./DownloadingFilesFromNCBI/AllCommensalStrains/*.gb"):
#     with open(fileName) as file:
#         locusName = ""
#         for line in file:
#             cols = line.split()
#             if len(cols) <= 1:
#                 continue
#             if cols[0] == "LOCUS":
#                 locusName = cols[1]
#             elif cols[1] == "BioProject:":
#                 if cols[2] in counts.keys():
#                     print("project", cols[2],"locus", locusName)
#                     counts[cols[2]] += 1
#                 else:
#                     counts[cols[2]] = 1
# print(counts)

def SeparateByLocus():
    commensalStrainPath = "./DownloadingFilesFromNCBI/AllCommensalStrains/"
    for filePath in glob(commensalStrainPath + "*.gb"):
        fileName = filePath.split("/")[-1]
        with open(filePath) as file:
            locusName = ""
            for line in file:
                cols = line.split()
                if len(cols) <= 1:
                    continue
                if cols[0] == "LOCUS":
                    locusName = cols[1]
            strainName = locusName[:-5]
            strainDir = commensalStrainPath + strainName + "/"
            if not os.path.isdir(strainDir):
                os.mkdir(strainDir)
            os.rename(filePath, strainDir + fileName) # move file to folder for strain
SeparateByLocus()