from scr.snpalign.functions import *

def divideGenomesByPhylogroup(alignedSnpsPath, phylogroupsPath, metadataPath="./metaDataForMetaCatsWithExtraMastitis.tsv", outDir="./"):
    genomeNameToSnpGenome = {}
    phylogroupToGenomes = {}
    with open(alignedSnpsPath) as file:
        nextKey = ""
        for line in file:
            line = line.strip()
            if line[0] == ">":
                nextKey = line[1:]
            else:
                genomeNameToSnpGenome[nextKey] = line
    print(genomeNameToSnpGenome.keys())
    with open(phylogroupsPath) as file:
        for line in file:
            cols = line.strip().split()
            if len(cols) == 0:
                continue
            # populate phylogroupToGenomes
            if not cols[1] in phylogroupToGenomes.keys():
                phylogroupToGenomes[cols[1]] = {"scaffold_" + re.sub("scaffold_|.fasta.vcf","", cols[0]) + ".fasta.vcf"}
            else:
                phylogroupToGenomes[cols[1]].add("scaffold_" + re.sub("scaffold_|.fasta.vcf","", cols[0]) + ".fasta.vcf")
    
    genomeToMetadata = {}
    with open(metadataPath) as file:
        for line in file:
            line = line.strip()
            cols = line.split("\t")
            if len(cols) < 2:
                continue
            genomeToMetadata[cols[0]] = cols[1:]
    print("phylogroup, mastCount, bovCom, avianCom, APEC")
    for phylogroup in phylogroupToGenomes.keys():
        mastCount = 0
        bovCom = 0
        avianCom = 0
        APEC = 0
        for genomeName in phylogroupToGenomes[phylogroup]:
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
        try:
            with open(outDir + "Phylogroup" + phylogroup +".afa", "w") as outFile:
                for genomeName in phylogroupToGenomes[phylogroup]:
                    try:
                        outFile.write(">" + genomeName.strip() + "\n" + genomeNameToSnpGenome[genomeName].strip() + "\n")
                    except KeyError:
                        print(genomeName)
                        print("genome not found")
        except: # if the phylogroup is unsure
            pass

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