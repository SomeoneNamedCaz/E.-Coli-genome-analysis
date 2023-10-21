from functions import *
import sys
# read presence absence csv file

# pathOfPresenceAbsenceFile = "./reannotatedMastitisGenomes/gene_presence_absence.Rtab"
# pathOfPresenceAbsenceFile = "./reannotatedMastitisGenomes/roary3NewVersion/gene_presence_absence.Rtab"
# pathOfPresenceAbsenceFile = "./natureStuff/naturePathogenRoary/gene_presence_absence.Rtab"
# nameOfOutFile = "AvianCommensalCounts.tsv"
# pathOfPresenceAbsenceFile = "./bovineCommensalRoary/gene_presence_absence.Rtab"
# nameOfOutFile = "APECCounts.tsv"
# nameOfOutFile = "CowCommensalCounts.tsv"
# nameOfOutFile = "mastitisCountsWith79Genomes.tsv"
# nameOfOutFile = "mastitisCountsTakeFour.tsv"
# pathOfPresenceAbsenceFile = sys.argv[1]
# nameOfOutFile = sys.argv[2]

# pathOfPresenceAbsenceFile = "/Users/cazcullimore/Documents/ericksonLabCode/AllAssemblies/allMastitisAssemblies/ProkkaOutFailedBecauseTooLongContigs/gffs/roaryOutputs/gene_presence_absence.Rtab"

# pathOfPresenceAbsenceFile = "/Users/cazcullimore/Documents/ericksonLabCode/AllAssemblies/AllCommensalBovineAssemblies/ProkkaOutFailedBecauseTooLongContigs/gffs/roaryOutputs/gene_presence_absence.Rtab"
# pathOfPresenceAbsenceFile = "/Users/cazcullimore/Documents/ericksonLabCode/AllAssemblies/APEC_assemblies/ProkkaOutFailedBecauseTooLongContigs/gffs/roaryOutputs/gene_presence_absence.Rtab"
# pathOfPresenceAbsenceFile = "/Users/cazcullimore/Documents/ericksonLabCode/AllAssemblies/Avian_commensal/ProkkaOutFailedBecauseTooLongContigs/gffs/roaryOutputs/gene_presence_absence.Rtab"
pathOfPresenceAbsenceFile = "/Users/cazcullimore/Documents/ericksonLabCode/mastitisRoary/ROARY file.tsv"
# pathOfPresenceAbsenceFile = ""

nameOfOutFile = pathOfPresenceAbsenceFile.split("/")[6] + "MichaelsRoaryGeneCounts.tsv"



with open(pathOfPresenceAbsenceFile) as infile, open(nameOfOutFile,"w") as outFile:
    onFirstLine = True
    indexOfgenomesToExclude = []
    for line in infile:
        cols = line.split("\t")
        # cols = line.split(",")
        geneName = cols[0]
        # print(geneName)
        geneCount = 0
        # print(line)
        if not onFirstLine:
            colIndex = 0
            for presenceNum in cols[8:94]:
                if not colIndex in indexOfgenomesToExclude:
                    geneCount += int(presenceNum)
                colIndex += 1
            outFile.write(geneName + "\t" + str(geneCount) + "\n")
        else: # for michaels roary
            print("header", cols[8:94])
            genomes = cols[8:94]
            index = 0
            # for genome in genomes:
            #     if genome in ["M6", "M9", "M11", "M12", "M50", "M37", "M109"]:
            #         indexOfgenomesToExclude.append(index)
            #     index += 1


        onFirstLine = False

