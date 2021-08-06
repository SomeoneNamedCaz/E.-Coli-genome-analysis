from functions import *

if (len(sys.argv) < 4):
    print("please provide the the gene count files to compare and the distance between most of the genes")
    exit(1)

distanceBetweenThem = sys.argv[3]

with open(sys.argv[1]) as file1:
    for line1 in file1:
        cols1 = line1.split("\t")
        if "group" in cols1[0]:
            continue
        with open(sys.argv[2]) as file2:
            for line2 in file2:
                cols2 = line2.split("\t")
                if cols1[0] == cols2[0]:
                    # print(cols1)
                    # print(cols2)
                    if abs(int(cols2[1])-int(cols1[1]))/int(cols2[1]) > 0.05:
                        print(cols1[0], int(cols1[1]), int(cols2[1]))#/int(cols1[1]))
                    # print("__")
                    break
