from functions import *

if (len(sys.argv) < 4):
    print("please provide the the gene count files to compare and the distance between most of the genes")
    exit(1)

distanceBetweenThem = sys.argv[3]

indexOfGeneCount = 3
genesRead = set()


with open(sys.argv[1]) as file1:
    for line1 in file1:
        cols1 = line1.split('"')[1::2] #line1.split("\t")
        # print(cols1[20] == "")
        if "group" in cols1[0] or "Gene" in cols1[0] or cols1[10].strip() != "":
            continue
        if not genesRead.isdisjoint({cols1[0]}):
            print("mutiple annotations for same gene")
        with open(sys.argv[2]) as file2:
            for line2 in file2:
                # cols2 = line2.split("\t")
                cols2 = line2.split('"')[1::2]
                # print(cols2[:30])
                if cols2[10].strip() != "":
                    continue
                if cols1[0] == cols2[0] and cols1[2].strip() == cols2[2].strip():
                    genesRead.add(cols1[0])
                    # print(cols1[:10])
                    # print(cols2[:10])
                    if abs(int(cols2[indexOfGeneCount])-int(cols1[indexOfGeneCount])) > 70: # /int(cols2[indexOfGeneCount])
                        print(cols1[0], int(cols1[indexOfGeneCount]), int(cols2[indexOfGeneCount]))#/int(cols1[1]))
                    # print("__")
                    break
