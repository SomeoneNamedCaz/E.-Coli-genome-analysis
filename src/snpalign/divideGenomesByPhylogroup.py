try:
    from .functions import *
except ImportError:
    from functions import *

def divideGenomesByPhylogroup(alignedSnpsPath, phylogroupsPath, metadataPath, outDir="./", debug=False):
    genomeNameToSnpGenome = {}
    phylogroupToGenomes = {}
    with open(alignedSnpsPath) as file:
        nextKey = ""
        for line in file:
            line = line.strip()
            if line[0] == ">":
                nextKey = re.sub("\..+","",line[1:])
            else:
                genomeNameToSnpGenome[nextKey] = line
    # if debug:
    #     print(genomeNameToSnpGenome.keys())
    with open(phylogroupsPath) as file:
        for line in file:
            cols = line.strip().split()
            if len(cols) == 0:
                continue
                
            # populate phylogroupToGenomes
            cols[0] = re.sub("\..+","",cols[0]) # remove any extensions
            if not cols[1] in phylogroupToGenomes.keys():
                phylogroupToGenomes[cols[1]] = {cols[0]}
            else:
                phylogroupToGenomes[cols[1]].add(cols[0])

    
    genomeToMetadata = readMetaDataAsDict(metadataPath)
    if debug:
        print("keys_", set(genomeToMetadata.keys()))#.difference(phylogroupToGenomes['G']))
        print("dif2",phylogroupToGenomes['G'])#.difference(set(genomeToMetadata.keys())))
    print("phylogroup, mastCount, bovCom, avianCom, APEC")
    for phylogroup in phylogroupToGenomes.keys():
        mastCount = 0
        bovCom = 0
        avianCom = 0
        APEC = 0
        for genomeName in phylogroupToGenomes[phylogroup]:
            if debug:
                print(genomeName)
            if genomeToMetadata[genomeName][0] == "cow":
                if genomeToMetadata[genomeName][1] == "pathogen":
                    mastCount += 1
                else:
                    bovCom += 1
            else:
                if genomeToMetadata[genomeName][1] == "pathogen":
                    APEC += 1
                else:
                    avianCom += 1
        print(phylogroup, mastCount, bovCom, avianCom, APEC)
        if "/" in phylogroup:
            continue
        with open(outDir + "Phylogroup" + phylogroup +".afa", "w") as outFile:
            for genomeName in phylogroupToGenomes[phylogroup]:
                try:
                    outFile.write(">" + genomeName.strip() + "\n" + genomeNameToSnpGenome[genomeName].strip() + "\n")
                except KeyError:
                    print(genomeName,"genome failed to align will be ignored")

    with open(outDir + "phylogroupDistribution.tsv", "w") as file:
        file.write("Phylogroup\tnumberOfGenomesInThePhylogroup\tgenomeNames\n")
        for phylogroup, genomes in phylogroupToGenomes.items():
            file.write(phylogroup + "\t" + str(len(genomes)) + "\t" + " ".join(genomes) + "\n")
            
if __name__ == "__main__":
    
    if len(sys.argv) < 3:
        print("please provide the phylogroupSnpFile with the aligned snps, and the phylogroup phylogroupSnpFile")
        exit(1)
    
    alignedSnpsPath = sys.argv[1]
    
    phylogroupsPath = sys.argv[2]
    divideGenomesByPhylogroup(alignedSnpsPath, phylogroupsPath)