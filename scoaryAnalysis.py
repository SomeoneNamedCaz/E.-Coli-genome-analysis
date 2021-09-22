from functions import *

try:
    scoaryFilePath = sys.argv[1]
except IndexError:
    print("include scoary file path")
    exit(1)


with open(scoaryFilePath) as file:
    for line in file:
        cols = line.split('"')[1::2]
        if cols[11] == 'Bonferroni_p':
            print(cols[11])
            continue

        if float(cols[10]) < 0.05:
            print(cols[0:2], cols[10])