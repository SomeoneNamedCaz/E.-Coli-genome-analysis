from secondaryPythonScripts.functions import *
def findNamesOfGroups(snpStatPath):
    groupToName = {}
    with open(snpStatPath) as snpStatFile:
        for line in snpStatFile:
            line = line.strip()
            cols = line.split("\t")
            if line == "" or cols[0].strip() == "Position" or cols[-1] in groupToName.keys():
                continue
            groups = [cols[-2].split("(")[0], cols[-2].split("(")[1].split("|")[1]]
            print(groups)
            groupToName[cols[-1]] = groups
    return groupToName


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(
            """please give arguments: combined megaCats output file""")
        exit(1)

    snpStatPath = sys.argv[1]
    print(findNamesOfGroups(snpStatPath))