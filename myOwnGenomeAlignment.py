from functions import *

with open("./allMastitisAssemblies/ragtagOutputs/longestScaffoldFiles/AllFilesInOneTry2ForSuper.txt") as file:
    genomeShifts = []
    for line in file:
        if line[0] != ">":
