import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 2:
    print('please provide one or more comma-separated file paths which contains the p-values'
          ' (you can provide the index of the column that contains'
          ' the pvalues if you want)')
    exit(1)
""" last command:
time python pValueHistogram.py /Users/cazcullimore/Documents/ericksonLabCode/megaCATS-main/firstMegaCatsRunOnData/animal-rMsaInput.txt-rResultMGCStat.txt,/Users/cazcullimore/Documents/ericksonLabCode/megaCATS-main/firstMegaCatsRunOnData/pathogenicity-rMsaInput.txt-rResultMGCStat.txt,/Users/cazcullimore/Documents/ericksonLabCode/megaCATS-main/firstMegaCatsRunOnData/megaCatsInsertAndDelete/animal-rMsaInput.txt-rResultMGCStat.txt,/Users/cazcullimore/Documents/ericksonLabCode/megaCATS-main/firstMegaCatsRunOnData/megaCatsInsertAndDelete/pathogenicity-rMsaInput.txt-rResultMGCStat.txt"""


filePaths = sys.argv[1].split(",")
try:
    pvalCol = sys.argv[2]
except IndexError:
    pvalCol = 1


pvalues = []
for path in filePaths:
    with open(path) as file:
        for line in file:
            try:
                pVal = float(line.split("\t")[pvalCol])
                # if pVal < (0.05 / 605_117):
                pvalues.append(pVal)
            except IndexError:
                pass
            except ValueError:
                pass

numComparisons = len(pvalues)
# for pValIndex in range(numComparisons - 1, -1, -1):
#     if pvalues[pValIndex] >= (0.05 / numComparisons):
#         del pvalues[pValIndex]
pvalues.sort()
print("making hist")
# pvalues = pvalues[:int(numComparisons/100)]
fig = plt.hist(pvalues, bins=100)
plt.title('Pvalue distribution')
plt.xlabel("p-value")
plt.ylabel("num p-values in bucket")
# plt.savefig("allSnpsSortedAndFloatedLessThan5.png")
plt.show()
