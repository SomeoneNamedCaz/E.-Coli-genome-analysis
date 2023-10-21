from functions import *
# # mastitisStrainsDirectory = "./ourMastitisAnnotatedGenomes/gbks/"
apecStrainsDirectory = "/Users/cazcullimore/Documents/ericksonLabCode/natureStuff/secondCommensalAnnotation/gbks/"
version1 = "GTGATTATTGATGAAAGTGCCGGTGAGGTTGTTATCGGCGCGAATACCCGTATTTGCCACGGTGCCGTTATTCAGGGGCCGGTAGTGATTGGCGCAAACTGCCTGATAGGTAATTATGCGTTTATTCGTCCTGGCACAATAATCAGCAATGGCGTAAAAATTGGTTTTGCCACCGAAATTAAAAATGCGGTTATTGAAGCGGAAGCAACGATTGGTCCGCAATGCTTTATTGCCGACTCGGTAGTTGCAAACCAGGCATATTTGGGCGCACAAGTACGTACCAGTAATCATCGTCTGGATGAACAACCCGTGTCTGTTCGAACTCCAGAGGGAATTATCGCTACCGGATGCGATAAATTAGGTTGTTATATCGGGAAGCGTTCACGCCTTGGTGTACAAGTTATTATTTTGCCTGGGCGAATTATTTCTCCGAACACACAACTTGGCCCGCGCGTGATTGTAGAACGTAATTTACCTAGTGGAACTTACTCACTCCGACAAGAACTTATCCGTACAGGAGATTAA"
version2 = "ATGTTGAATAATGCTATGAGCGTAGTGATCCTTGCCGCAGGCAAAGGCACGCGCATGTATTCCGATCTTCCGAAAGTGCTGCATACCCTTGCCGGGAAAGCGATGGTTCAGCATGTCATTGATGCTGCGAATGAATTAGGCGCAGCGCACGTTCACCTGGTTTACGGTCACGGCGGCGATCTGCTTAAACAGGCGCTGAAAGATGACAACCTGAACTGGGTGCTTCAGGCAGAACAGCTGGGTACAGGTCATGCGATGCAGCAGGCCGCACCTTTCTTTGCCGATGATGAAGACATTTTAATGCTCTACGGCGACGTGCCGCTGATCTCTGTCGAAACACTCCAGCGCCTGCGTGATGCTAAACCGCAGGGTGGCATTGGTCTGCTGACGGTAAAACTGGATGATCCGACCGGTTATGGACGTATCACCCGTGAAAACGGCAAAGTTACCGGCATTGTTGAGCACAAAGACGCCACCGACGAGCAGCGTCAGATTCAGGAGATCAACACCGGCATTCTGATTGCCAACGGCGCAGATATGAAACGCTGGCTGGCGAAGCTGACCAACAATAACGCTCAGGGTGAATACTACATCACCGACATTATTGCGCTGGCGTATCAGGAAGGGCGTGAAATCGTCGCCGTTCATCCGCAACGTTTAAGCGAAGTAGAAGGCGTGAATAACCGCCTGCAACTCTCCCGACTGGAGCGCGTTTACCAGTCCGAACAGGCTGAAAAACTGCTGTTAGCAGGCGTTATGCTGCGCGATCCGGCGCGTTTTGATCTGCGCGGTACGCTTACTCACGGGCGCGATGTTGAAATTGATACTAACGTTATCATCGAGGGCAACGTGACTCTCGGTCATCGCGTGAAAATCGGCACCGGTTGCGTGATTAAAAACAGCGTGATTGGCGATGATTGCGAAATCAGTCCGTATACCGTTGTGGAAGATGCGAATCTGGCAGCGGCCTGTACCATTGGCCCGTTTGCCCGTCTGCGTCCTGGTGCTGAGTTGCTGGAAGGTGCACACGTCGGTAACTTCGTTGAGATGAAAAAAGCGCGTCTGGGTAAAGGCTCGAAAGCTGGTCATCTGACTTACCTGGGCGATGCGGAAATTGGCGATAACGTTAACATCGGCGCGGGAACCATTACCTGCAACTACGATGGTGCGAATAAATTTAAGACCATTATCGGTGACGATGTGTTTGTCGGTTCCGACACTCAGCTGGTGGCCCCGGTAACAGTAGGCAAAGGCGCGACCATTGCTGCGGGTACAACTGTGACGCGTAATGTCGGCGAAAACGCCCTGGCGATCAGCCGTGTGCCGCAGACTCAAAAAGAAGGCTGGCGTCGTCCGGTAAAGAAAAAGTAA"
num_duplicates = 0
num_oneNormal = 0
for filePath in glob(apecStrainsDirectory + "*"):
    print(filePath)
    contigs = getContigs(filePath)
    genes = getGenesOnContigs(filePath, contigs)
    # seqToSearch = ""
    # if "glmU" in genes.keys():
    #     # seqToSearch = genes["glmU"][1]
    #     # print(seqToSearch)
    #     for geneName in genes.keys():
    #         if stringIsOnlyATCG(geneName):
    #             # print(geneName)
    #             if genesAreWithinPercentIdentical(version2, geneName, 0.8):# or genesAreWithinPercentIdentical(version2, geneName, 0.5):
    #                 print("found second glmU")

    print(filePath)
    if "glmU_2" in genes.keys():
        num_duplicates += 1
        glm1 = genes["glmU_1"][1]
        glm2 = genes["glmU_2"][1]
        print("two glmU's found")
        if geneSimilarity(glm1,glm2) > 0.8:
            print("glmU's simlar to each other")
        if geneSimilarity(glm1,version2) > 0.8:
            print("glm1 is good")
            num_oneNormal += 1
        if geneSimilarity(glm2,version2) > 0.8:
            print("glm2 is good")
            num_oneNormal += 1
        if geneSimilarity(glm1,version2) > 0.8 and geneSimilarity(glm2,version2) > 0.8:
            print("real duplicate")
        if aGeneIsAPartofAnotherGene(glm1, glm2):
            print("one is subset of the other")
        if  (genesAreWithinPercentIdentical(glm1, glm2) == True) and (False == (geneSimilarity(glm1, glm2) > 0.8)):
            print("percent id doesnt work")
        elif  (genesAreWithinPercentIdentical(glm1, glm2) == False) and (True == (geneSimilarity(glm1, glm2) > 0.8)):
            print("percent id doesnt work")


    # elif "glmU_2" in genes.keys():
    #     print(genes["glmU_2"][1])

print(" duplicates", num_duplicates) # 267 in first annotation for APEC
print("num normal ", num_oneNormal)

# conservedDomainsOfGlm = "mlnnamsvvilaagkgtrmysdlpkvlhtlagkamvqhvidaanelgaahvhlvyghggdllkqalkddnlnwvlqaeqlgtghamqqaapffaddedilmlygdvplisvetlqrlrdakpqggiglltvklddptgygritrengkvtgivehkdatdeqrqiqeintgiliangadmkrwlakltnnnaqgeyyitdiialayqegreivavhpqrlsevegvnnrlqlsrlervyqseqaeklllagvmlrdparfdlrgtlthgrdveidtnviiegnvtlghrvkigtgcviknsvigddceispytvvedanlaaactigpfarlrpgaellegahvgnfvemkkarlgkgskaghltylgdaeigdnvnigagtitcnydgankfktiigddvfvgsdtqlvapvtvgkgatiaagttvtrnvgenalaisrvpqtqkegwrrpvkkk"
# print(conservedDomainsOfGlm.upper())
# print(translate(version2))
# print(geneSimilarity(conservedDomainsOfGlm, translate(version2)))
#
# otherDomain = "asrlervyqseqaeklllagvmlrdparfdlrgtlthgrdveidtnviiegnvtlghrvkigtgcviknsvigddceispytvvedanlaaactigpfarlrpgaellegahvgnfvemkkarlgkgskaghltylgdaeigdnvnigagtitcnydgankfktiigddvfvgsdtqlvapvtvgkgatiaagttvtrnvgenalaisrvpqtqkegwrrpa"
# print(otherDomain.upper())
# print(translate(version1))