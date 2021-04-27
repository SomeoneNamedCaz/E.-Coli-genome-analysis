# sort assemblies by pathogenicity or not
from functions import *
pathOfFilesToSort = "./natureAssemblies/*/"
pathogenDestination = "./pathogenNatureAssemblies/"
commensalDestination = "./commensalNatureAssemblies/"
pathogenIDs = []
commensalIDs = []
colIndexThatLooksAtPathogenicity = 9
with open("ifstrainIsPathogenicTableNature.tsv") as pathogenicityFile:
    for line in pathogenicityFile:
        cols = line.split("\t")
        ID = cols[0]
        pathogenicityMarker = cols[colIndexThatLooksAtPathogenicity]
        if pathogenicityMarker == "Asymptomatic carriage":
            commensalIDs.append(ID)
        elif pathogenicityMarker == "Infection ":
            pathogenIDs.append(ID)

for filePath in glob(pathOfFilesToSort + "*"):
    fileName = filePath.split("/")[-1]
    identifier = fileName[:4]
    if identifier in pathogenIDs: # if pathogen
        os.rename(filePath, pathogenDestination + fileName)
    elif identifier in commensalIDs: # if commensal
        os.rename(filePath, commensalDestination + fileName)

