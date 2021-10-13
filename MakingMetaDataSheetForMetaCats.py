from functions import *

metadataPath = "./metaDataForMetaCatsWithExtraMastitis.tsv"
catagoryPaths = []#["./AllAssemblies/allMastitisAssemblies/*.fasta"]
catagoryPaths.append("./AllAssemblies/AllCommensalBovineAssemblies/*.fasta")
catagoryPaths.append("./AllAssemblies/APEC_assemblies/*.fasta")
catagoryPaths.append("./AllAssemblies/Avian_commensal/*.fasta")
catagoryPaths.append("./DownloadingFilesFromNCBI/MissedMastitisCombinedFastas/*.fasta")

with open(metadataPath, "w") as file:
    # for catagoryPath in catagoryPaths:
    #     file.write(catagoryPath.split("/")[-2] + "\t")

    file.write("strainName\tanimal\tpathogenicity\n")
    for catagoryPath in catagoryPaths:
        catagory = catagoryPath.split("/")[-2]
        for filePath in glob(catagoryPath):
            fileName = filePath.split("/")[-1]
            file.write("scaffold_" + fileName + ".vcf\t")
            if catagory == "allMastitisAssemblies" or catagory == "MissedMastitisCombinedFastas":
                file.write("cow\tpathogen")
            elif catagory == "AllCommensalBovineAssemblies":
                file.write("cow\tcommensal")
            elif catagory == "APEC_assemblies":
                file.write("chicken\tpathogen")
            elif catagory == "Avian_commensal":
                file.write("chicken\tcommensal")
            file.write("\n")
    file.write("\n")