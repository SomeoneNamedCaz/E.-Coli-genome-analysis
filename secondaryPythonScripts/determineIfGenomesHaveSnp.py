from secondaryPythonScripts.functions import *

m12Path = "/Users/cazcullimore/dev/data/M12.gbk"
k12Path = "/Users/cazcullimore/dev/data/refGenomes/k-12.gbff"

m12Genes = getGenesOnContigs(m12Path,getContigs(m12Path))
k12Genes = getGenesOnContigs(k12Path, getContigs(k12Path))

geneName = "qorB"

print(m12Genes[geneName].isForward)
print(m12Genes[geneName].sequence)

print(k12Genes[geneName].isForward)
print(k12Genes[geneName].sequence)

