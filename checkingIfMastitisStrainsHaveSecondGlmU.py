from functions import *
mastitisStrainsDirectory = "./ourMastitisAnnotatedGenomes/gbks/"
version1 = "GTGATTATTGATGAAAGTGCCGGTGAGGTTGTTATCGGCGCGAATACCCGTATTTGCCACGGTGCCGTTATTCAGGGGCCGGTAGTGATTGGCGCAAACTGCCTGATAGGTAATTATGCGTTTATTCGTCCTGGCACAATAATCAGCAATGGCGTAAAAATTGGTTTTGCCACCGAAATTAAAAATGCGGTTATTGAAGCGGAAGCAACGATTGGTCCGCAATGCTTTATTGCCGACTCGGTAGTTGCAAACCAGGCATATTTGGGCGCACAAGTACGTACCAGTAATCATCGTCTGGATGAACAACCCGTGTCTGTTCGAACTCCAGAGGGAATTATCGCTACCGGATGCGATAAATTAGGTTGTTATATCGGGAAGCGTTCACGCCTTGGTGTACAAGTTATTATTTTGCCTGGGCGAATTATTTCTCCGAACACACAACTTGGCCCGCGCGTGATTGTAGAACGTAATTTACCTAGTGGAACTTACTCACTCCGACAAGAACTTATCCGTACAGGAGATTAA"
version2 = "ATGTTGAATAATGCTATGAGCGTAGTGATCCTTGCCGCAGGCAAAGGCACGCGCATGTATTCCGATCTTCCGAAAGTGCTGCATACCCTTGCCGGGAAAGCGATGGTTCAGCATGTCATTGATGCTGCGAATGAATTAGGCGCAGCGCACGTTCACCTGGTTTACGGTCACGGCGGCGATCTGCTTAAACAGGCGCTGAAAGATGACAACCTGAACTGGGTGCTTCAGGCAGAACAGCTGGGTACAGGTCATGCGATGCAGCAGGCCGCACCTTTCTTTGCCGATGATGAAGACATTTTAATGCTCTACGGCGACGTGCCGCTGATCTCTGTCGAAACACTCCAGCGCCTGCGTGATGCTAAACCGCAGGGTGGCATTGGTCTGCTGACGGTAAAACTGGATGATCCGACCGGTTATGGACGTATCACCCGTGAAAACGGCAAAGTTACCGGCATTGTTGAGCACAAAGACGCCACCGACGAGCAGCGTCAGATTCAGGAGATCAACACCGGCATTCTGATTGCCAACGGCGCAGATATGAAACGCTGGCTGGCGAAGCTGACCAACAATAACGCTCAGGGTGAATACTACATCACCGACATTATTGCGCTGGCGTATCAGGAAGGGCGTGAAATCGTCGCCGTTCATCCGCAACGTTTAAGCGAAGTAGAAGGCGTGAATAACCGCCTGCAACTCTCCCGACTGGAGCGCGTTTACCAGTCCGAACAGGCTGAAAAACTGCTGTTAGCAGGCGTTATGCTGCGCGATCCGGCGCGTTTTGATCTGCGCGGTACGCTTACTCACGGGCGCGATGTTGAAATTGATACTAACGTTATCATCGAGGGCAACGTGACTCTCGGTCATCGCGTGAAAATCGGCACCGGTTGCGTGATTAAAAACAGCGTGATTGGCGATGATTGCGAAATCAGTCCGTATACCGTTGTGGAAGATGCGAATCTGGCAGCGGCCTGTACCATTGGCCCGTTTGCCCGTCTGCGTCCTGGTGCTGAGTTGCTGGAAGGTGCACACGTCGGTAACTTCGTTGAGATGAAAAAAGCGCGTCTGGGTAAAGGCTCGAAAGCTGGTCATCTGACTTACCTGGGCGATGCGGAAATTGGCGATAACGTTAACATCGGCGCGGGAACCATTACCTGCAACTACGATGGTGCGAATAAATTTAAGACCATTATCGGTGACGATGTGTTTGTCGGTTCCGACACTCAGCTGGTGGCCCCGGTAACAGTAGGCAAAGGCGCGACCATTGCTGCGGGTACAACTGTGACGCGTAATGTCGGCGAAAACGCCCTGGCGATCAGCCGTGTGCCGCAGACTCAAAAAGAAGGCTGGCGTCGTCCGGTAAAGAAAAAGTAA"
for filePath in glob(mastitisStrainsDirectory + "*"):
    print(filePath)
    contigs = GetContigs(filePath)
    genes = GetGenesOnContigs(filePath,contigs)
    # seqToSearch = ""
    if "glmU" in genes.keys():
        # seqToSearch = genes["glmU"][1]
        # print(seqToSearch)
        for geneName in genes.keys():
            if StringIsOnlyATCG(geneName):
                # print(geneName)
                if GenesAreWithinPercentIdentical(version1,geneName, 0.5) or GenesAreWithinPercentIdentical(version2,geneName, 0.5):
                    print("found second glmU")
    # elif "glmU_2" in genes.keys():
    #     print(genes["glmU_2"][1])

