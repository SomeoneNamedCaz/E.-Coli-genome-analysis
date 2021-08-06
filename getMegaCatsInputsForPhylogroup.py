from functions import *

if len(sys.argv) < 3:
    print("please provide the file with the aligned snps, and the phylogroup file")
    exit(1)

alignedSnpsPath = sys.argv[1]
genomeNameToSnpGenome = {}
phylogroupsPath = sys.argv[2]
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
        if not cols[1] in phylogroupToGenomes.keys():
            phylogroupToGenomes[cols[1]] = ["scaffold_" + cols[0] + ".fasta.vcf"]
        else:
            phylogroupToGenomes[cols[1]].append("scaffold_" + cols[0] + ".fasta.vcf")



for phylogroup in phylogroupToGenomes.keys():
    with open("allSnpsForPhylogroup" + phylogroup.split("/")[0] +".afa", "w") as outFile:
        for genomeName in phylogroupToGenomes[phylogroup]:
            try:
                outFile.write(">" + genomeName + "\n" + genomeNameToSnpGenome[genomeName] + "\n")
            except KeyError:
                print("genome not found")



    # to check if phylogroup is restricted to a particular metadata value

    # if phylogroup != "G":
    #     continue
    # for genomeName in phylogroupToGenomes[phylogroup]:
    #     with open("/Users/cazcullimore/Documents/ericksonLabCode/megaCATS-main/metaDataForMetaCats.tsv") as metadataFile:
    #         for line in metadataFile:
    #             cols = line.strip().split("\t")
    #             if genomeName == cols[0]:
    #                 print(cols)