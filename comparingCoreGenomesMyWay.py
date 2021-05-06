from functions import *
import random
# mastitisPrevalenceFileName = "genePrevalencesMastitis.tsv" # structured as geneName + \t + count
# commensalBovinePrevalenceFileName = "genePrevalancesCowCommensalOrignalMethod.tsv"
# commensalAvianPrevalenceFileName = "genePrevalancesAcommensal.tsv"
# APECPrevalenceFileName = "genePrevalancesAPEC.tsv"

mastitisPrevalenceFileName = "mastitisCountsTakeFour.tsv"
commensalAvianPrevalenceFileName = "AvianCommensalCounts.tsv"
APECPrevalenceFileName = "APECCounts.tsv"
commensalBovinePrevalenceFileName = "CowCommensalCounts.tsv"
def findCoreGenome(prevalenceFileName): # extra are M6, M9, M11, M12,M50,M37,M109
    geneCounts = []
    with open(prevalenceFileName) as prevalenceFile:
        for line in prevalenceFile:
            line = line.strip()
            cols = line.split("\t")
            geneName = cols[0]
            geneCount = int(cols[1]) #+ random.gauss(0, 2) # very good match with 1 somwhat bad match (50%) with 2, no core genome with 3
            geneCounts.append(geneCount)
    geneCounts.sort()
    assumedNumberOfStrains = geneCounts[-100]
    print(assumedNumberOfStrains)
    # read through file to get percents
    genes99_100 = []
    genes95_99 = []
    genes15_95 = []
    genes0_15 = []
    totalGenes = 0
    with open(prevalenceFileName) as prevalenceFile:
        for line in prevalenceFile:
            line = line.strip()
            cols = line.split("\t")
            geneName = cols[0]
            geneCount = int(cols[1])
            if geneCount/assumedNumberOfStrains >= 0.95:
                genes99_100.append(geneName)
            elif geneCount/assumedNumberOfStrains >= 0.80:
                genes95_99.append(geneName)
            elif geneCount/assumedNumberOfStrains > 0.15:
                genes15_95.append(geneName)
            else:
                genes0_15.append(geneName)
            totalGenes += 1
    return [genes0_15,genes15_95,genes95_99,genes99_100]
mastGenes = findCoreGenome(mastitisPrevalenceFileName)
lowerMastGenes = mastGenes[0] + mastGenes[1]
mastGenes = mastGenes[-1]# + mastGenes[-2]
comensalBovineGenes = findCoreGenome(commensalBovinePrevalenceFileName)
lowerComensalBovineGenes = comensalBovineGenes[0] + comensalBovineGenes[1]
comensalBovineGenes = comensalBovineGenes[-1] #+ comensalBovineGenes[-2]
APECGenes = findCoreGenome(APECPrevalenceFileName)
lowerAPECGenes = APECGenes[0] + APECGenes[1]
APECGenes = APECGenes[-1] #+ APECGenes[-2]
avianGenes = findCoreGenome(commensalAvianPrevalenceFileName)
lowerAvianGenes = avianGenes[0] + avianGenes[1]
avianGenes = avianGenes[-1]# + avianGenes[-2]

sharedByAll = []
sharedByBovine = []
sharedByChicken = []
sharedByDiseased = []
sharedMastAndAvian = []
sharedcomensalBovineAndAPEC = []
sharedByComensal = []

# for mastGene in mastGenes:
#     # if gene isn't hypothetical
#     if mastGene.count("A") + mastGene.count("C") + mastGene.count("T") + mastGene.count("G") != len(mastGene):
#         if mastGene in comensalBovineGenes and mastGene in APECGenes and mastGene in avianGenes:
#             sharedByAll.append(mastGene)
#         elif mastGene in comensalBovineGenes and not mastGene in APECGenes and not mastGene in avianGenes:
#             sharedByBovine.append(mastGene)
#         elif not mastGene in comensalBovineGenes and mastGene in APECGenes and not mastGene in avianGenes:
#             sharedByDiseased.append(mastGene)
#     else:
#         if HypotheticalProteinSeqInList(mastGene,comensalBovineGenes, cutoff=0.9) and HypotheticalProteinSeqInList(mastGene,APECGenes, cutoff=0.9) and HypotheticalProteinSeqInList(mastGene,avianGenes, cutoff=0.9):
#             sharedByAll.append(mastGene)
#         elif HypotheticalProteinSeqInList(mastGene,comensalBovineGenes, cutoff=0.9) and not HypotheticalProteinSeqInList(mastGene,APECGenes, cutoff=0.9) and not HypotheticalProteinSeqInList(mastGene,avianGenes, cutoff=0.9):
#             sharedByBovine.append(mastGene)
#         elif not HypotheticalProteinSeqInList(mastGene,comensalBovineGenes, cutoff=0.9) and HypotheticalProteinSeqInList(mastGene,APECGenes, cutoff=0.9) and not HypotheticalProteinSeqInList(mastGene,avianGenes, cutoff=0.9):
#             sharedByDiseased.append(mastGene)
#
# for commensalGene in comensalBovineGenes:
#     if commensalGene.count("A") + commensalGene.count("C") + commensalGene.count("T") + commensalGene.count("G") != len(commensalGene):
#         if commensalGene in APECGenes and not commensalGene in avianGenes and not commensalGene in mastGenes:
#             sharedcomensalBovineAndAPEC.append(commensalGene)
#         elif not commensalGene in APECGenes and commensalGene in avianGenes and not commensalGene in mastGenes:
#             sharedByComensal.append(commensalGene)
#     else:
#         if HypotheticalProteinSeqInList(commensalGene,APECGenes, cutoff=0.9) and not HypotheticalProteinSeqInList(commensalGene,avianGenes, cutoff=0.9) and not HypotheticalProteinSeqInList(commensalGene,mastGenes, cutoff=0.9):
#             sharedcomensalBovineAndAPEC.append(commensalGene)
#         elif not HypotheticalProteinSeqInList(commensalGene,APECGenes, cutoff=0.9) and HypotheticalProteinSeqInList(commensalGene,avianGenes, cutoff=0.9) and not HypotheticalProteinSeqInList(commensalGene,mastGenes, cutoff=0.9):
#             sharedByComensal.append(commensalGene)
# for APECGene in APECGenes:
#     if APECGene.count("A") + APECGene.count("C") + APECGene.count("T") + APECGene.count("G") != len(APECGene):
#         if not APECGene in comensalBovineGenes and APECGene in avianGenes and not APECGene in mastGenes:
#             sharedByChicken.append(APECGene)
#     else:
#         if not HypotheticalProteinSeqInList(APECGene,comensalBovineGenes, cutoff=0.9) and HypotheticalProteinSeqInList(APECGene,avianGenes, cutoff=0.9) and not HypotheticalProteinSeqInList(APECGene,mastGenes, cutoff=0.9):
#             sharedByChicken.append(APECGene)
#
# for avianGene in avianGenes:
#     if avianGene.count("A") + avianGene.count("C") + avianGene.count("T") + avianGene.count("G") != len(avianGene):
#         if not avianGene in comensalBovineGenes and not avianGene in APECGenes and avianGene in mastGenes:
#             sharedMastAndAvian.append(avianGene)
#     else:
#         if not HypotheticalProteinSeqInList(avianGene,comensalBovineGenes, cutoff=0.9) and not HypotheticalProteinSeqInList(avianGene,APECGenes, cutoff=0.9) and HypotheticalProteinSeqInList(avianGene,mastGenes, cutoff=0.9):
#             sharedMastAndAvian.append(avianGene)

# # old comparisions
for mastGene in mastGenes:
    if mastGene in comensalBovineGenes and mastGene in APECGenes and mastGene in avianGenes:
        sharedByAll.append(mastGene)
    elif mastGene in comensalBovineGenes and mastGene in lowerAPECGenes and mastGene in lowerAvianGenes:
        sharedByBovine.append(mastGene)
    elif mastGene in lowerComensalBovineGenes and mastGene in APECGenes and mastGene in lowerAvianGenes:
        sharedByDiseased.append(mastGene)

for commensalGene in comensalBovineGenes:
    if commensalGene in APECGenes and commensalGene in lowerAvianGenes and commensalGene in lowerMastGenes:
        sharedcomensalBovineAndAPEC.append(commensalGene)
    elif commensalGene in lowerAPECGenes and commensalGene in avianGenes and commensalGene in lowerMastGenes:
        sharedByComensal.append(commensalGene)

for APECGene in APECGenes:
    if APECGene in lowerComensalBovineGenes and APECGene in avianGenes and APECGene in lowerMastGenes:
        sharedByChicken.append(APECGene)

for avianGene in avianGenes:
    if avianGene in lowerComensalBovineGenes and avianGene in lowerAPECGenes and avianGene in mastGenes:
        sharedMastAndAvian.append(avianGene)

largeDifInPrevalence = []
for mastGene in mastGenes:
    if (mastGene in lowerComensalBovineGenes and mastGene in APECGenes and mastGene in lowerAvianGenes):
        largeDifInPrevalence.append(mastGene)

print("sharedByAll")
print(sharedByAll)
print("\n\nsharedByBovine")
print(sharedByBovine)
print("\n\nsharedByChicken")
print(sharedByChicken)
print("\n\nsharedByDiseased")
print(sharedByDiseased)
print("\n\nlargeDifInPrevalence")
print(largeDifInPrevalence)
print("\n\nsharedMastAndAvian")
print(sharedMastAndAvian)
print("\n\nsharedcomensalBovineAndAPEC")
print(sharedcomensalBovineAndAPEC)
print("\n\nsharedByComensal")
print(sharedByComensal)
print('\n\nmastitisCoregenomeSize')
print(len(mastGenes))
#
# print("\n\n\n\n\n\n\n\n\n\n\n") # mast good and apec good, avian bad
# print(avianGenes)
# c95 = ['pdeI', 'zinT', 'ygjI', 'yhdJ', 'ulaE_2', 'qorB', 'tauA', 'clsB', 'recF_1', 'dbpA_1', 'cbpA', 'fhlA', 'nrdE', 'pdeF', 'yjfF', 'pagP', 'rfaG', 'trg', 'maeA', 'marR', 'glmU_2', 'xerC_1', 'ygjK', 'yafP', 'bcsC', 'nagB_2', 'fbpC2', 'rhaR_2', 'rihC', 'fdrA_1', 'uidA', 'ahr', 'group_2417', 'ompC', 'ygfK', 'ybbY', 'lsrK_2', 'flhA_2', 'bioF_1', 'ymgD', 'caiA_2', 'ygcS_2', 'fimH_2', 'fixX_3', 'speF', 'ltnD', 'glcR', 'ecpR', 'fimG_2', 'caiT_1', 'fecA']
# c98 = ['ptrB', 'garP', 'agaC_2', 'yhdJ', 'qorB', 'tauA', 'hha', 'srlR_2', 'clsB', 'recF_1', 'rbsA_1', 'gsk', 'gadE', 'dbpA_1', 'cbpA', 'ycaO', 'fhlA', 'bhsA_3', 'gss', 'nrdE', 'chbF', 'rlmA', 'xseA', 'gadW', 'xylB', 'yjfF', 'phoE', 'aqpZ', 'torA', 'por_2', 'sutR', 'rtcA', 'leuO', 'dsbG', 'ltaE', 'efeB', 'yedK', 'creD', 'maeA', 'marR', 'fbaB', 'nudI', 'glmU_2', 'xerC_1', 'comEC', 'nagB_2', 'fbpC2', 'rfaC', 'rhaR_2', 'fdrA_1', 'fucI', 'uidA', 'yjfC', 'group_2417', 'ygfK', 'ghxQ', 'fliK', 'ybbY', 'mocA', 'bioF_1', 'ymgD']
# c100 = ['yjjQ', 'malZ', 'entA', 'dhaM', 'pmrD', 'yhdY', 'frmB', 'copA', 'glsA1', 'ybcI', 'rpsU', 'fur', 'cspD', 'ygdR_2', 'ygaM', 'entS', 'nagA', 'rsfS', 'exbB', 'narH', 'rlmM', 'glgB', 'clpA', 'entH', 'rnd', 'purU', 'gltK', 'corC', 'malK', 'cysB', 'hybO', 'recF_1', 'trpS', 'dnaG', 'rfaF', 'topA_2', 'rbsA_1', 'ftsP', 'exbD', 'uppP', 'ttdB', 'plsC', 'rlmB', 'rlmH', 'ycgM', 'yeiP', 'rffH', 'yccU', 'kdpB', 'matP', 'yegT', 'gltI', 'bdm', 'truD', 'holE', 'aas', 'yebV', 'recX', 'sbp', 'gshA', 'moeA', 'ydiB', 'eptC', 'cca', 'amiB', 'ybjG', 'bcsZ', 'serC', 'mzrA', 'chbF', 'ogt', 'cmoB', 'yfcF', 'allS', 'helD', 'ugpA', 'dgcM', 'rlmA', 'arnF', 'fepE', 'xseA', 'eamB', 'ygdG', 'fucK', 'loiP', 'ygiD', 'yiaB', 'phnO', 'yjfF', 'gcd', 'frmA', 'ybdL', 'lnt', 'chiQ', 'pflB', 'ycaL', 'torA', 'por_2', 'dauA', 'narJ', 'uspF', 'manY', 'yegS', 'mglA_1', 'inaA', 'srlB', 'lplT', 'dsbC', 'pepP', 'nupG', 'qseC', 'ygiF', 'yhjB', 'lldD', 'rffG', 'wecD', 'glnA', 'thiE', 'ulaA', 'prpC', 'ssuC_2', 'lrp', 'rpsA', 'mrdA', 'creD', 'ghrA', 'cvrA', 'hslJ', 'marA', 'yebF', 'fliS', 'mtfA', 'fbaB', 'arnE', 'flk', 'pdxK', 'amiA', 'endA', 'glmU_2', 'ttdA', 'argE_1', 'argC', 'ghxP', 'alsT', 'allR', 'citX', 'mntS', 'artQ', 'dmsB_1', 'msbA_2', 'ompA', 'chbG', 'cysH_1', 'astD', 'yfcG', 'tmcA', 'nagB_2', 'torD', 'tilS', 'purT', 'ypdF', 'ydjZ_1', 'fucI', 'gcvP']
# large diff : 'ulaE_2', 'glmU_2'

with open("genesShared.tsv", "w") as sharedGenesFile:


    def writeOutSharedGeneList(sharedGeneList):
        for gene in sharedGeneList:
            sharedGenesFile.write("\t"+gene)
        sharedGenesFile.write("\n")

    sharedGenesFile.write("sharedByAll")
    writeOutSharedGeneList(sharedByAll)

    sharedGenesFile.write("sharedByBovine")
    writeOutSharedGeneList(sharedByBovine)

    sharedGenesFile.write("sharedByChicken")
    writeOutSharedGeneList(sharedByChicken)

    sharedGenesFile.write("sharedByDiseased")
    writeOutSharedGeneList(sharedByDiseased)

    sharedGenesFile.write("sharedMastAndAvian")
    writeOutSharedGeneList(sharedMastAndAvian)

    sharedGenesFile.write("sharedcomensalBovineAndAPEC")
    writeOutSharedGeneList(sharedcomensalBovineAndAPEC)

    sharedGenesFile.write("sharedByComensal")
    writeOutSharedGeneList(sharedByComensal)