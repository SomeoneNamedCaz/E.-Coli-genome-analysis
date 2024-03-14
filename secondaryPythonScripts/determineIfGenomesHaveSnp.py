try:
	from secondaryPythonScripts.functions import *
except ImportError:
	from functions import *
m12Path = "/Users/cazcullimore/dev/data/refGenomes/M12.gbk"
k12Path = "/Users/cazcullimore/dev/data/refGenomes/k-12.gbff"
# k12Path = "/Users/cazcullimore/dev/data/DH5alpha.gbk"

m12Genes = getGenesOnContigs(m12Path,getContigs(m12Path))
k12Genes = getGenesOnContigs(k12Path, getContigs(k12Path))

# geneName = "ytfT"
positions = []#[672,48,276,325,326,645,647,646,774,58]
printOnlySnpsThatBothGenomesSatisfy = False

pathogenicNucs = []
commensalNucs = []
snpLines = []
refGenes = m12Genes
sortedSnpsFilePath = DATA_DIR + "RedoingEverything/snpsSortedBySignificanceWithGenesContainingThemM12RefGenomePathogenicity.tsv"

with open(sortedSnpsFilePath) as file:
	indexOfFirstGroup = 5
	indexOfSecondGroup = 8
	indexOfPosition = 2
	for line in file:
		line = line.strip()
		cols = line.split("\t")
		
		if len(cols) < 1 or cols[0] == "SNP pValue":
			continue
		# if cols[4] == geneName:
		if cols[1] == "pathogen":
			firstGroupIsPathogen = True
		else:
			firstGroupIsPathogen = False
		if firstGroupIsPathogen:
			pathogenicNucs.append(cols[indexOfFirstGroup]) # first group
			commensalNucs.append(cols[indexOfSecondGroup])
		else:
			pathogenicNucs.append(cols[indexOfSecondGroup])  # first group
			commensalNucs.append(cols[indexOfFirstGroup])
		positions.append(int(cols[indexOfPosition]))
		snpLines.append(line)

index = -1
for position in positions:
	index += 1
	printSnp = not printOnlySnpsThatBothGenomesSatisfy
	geneName = snpLines[index].split("\t")[4]
	if "NearestGeneIs:" in geneName or geneName != "cpdB":
		continue
	refSeq = snpLines[index].split("\t")[-1]
	
	try:
		if "unnamed" in geneName:
			
			m12Nuc = "not found"
			k12Nuc = "not found"
			for gene in getGenesOnContigsByPosition(m12Path,getContigs(m12Path)):
				if geneSimilarity(gene.sequence, refSeq,quiet=True) > 0.9:
					m12Nuc = gene.sequence[position]
			
			for gene in getGenesOnContigsByPosition(k12Path, getContigs(k12Path)):
				if geneSimilarity(gene.sequence, refSeq,quiet=True) > 0.9:
					k12Nuc = gene.sequence[position]
		else:
			if refGenes[geneName].isForward:
				try:
					m12Nuc = m12Genes[geneName].sequence[position]
				except KeyError:
					m12Nuc = "no gene"
				try:
					k12Nuc = k12Genes[geneName].sequence[position]
				except KeyError:
					k12Nuc = "no gene"
			else:
				try:
					m12Nuc = reverseComplement(m12Genes[geneName].sequence)[position]
				except KeyError:
					m12Nuc = "no gene"
				try:
					k12Nuc = reverseComplement(k12Genes[geneName].sequence)[position]
				except KeyError:
					k12Nuc = "no gene"
			if m12Nuc == pathogenicNucs[index] and k12Nuc == commensalNucs[index]:
				printSnp = True
				
		if printSnp:
			print(snpLines[index])
			print(position)
			if m12Nuc == pathogenicNucs[index]:
				print("m12 has pathogenic snp")
			else:
				print("m12 doesn't have pathogenic snp")
			if k12Nuc == commensalNucs[index]:
				print("k12 has commensal snp")
			else:
				print("k12 doesn't have commensal snp")
			
			print("m12 nuc", m12Nuc, "pathogenic nuc", pathogenicNucs[index])
			print("k12 nuc", k12Nuc, "commensal nuc", commensalNucs[index])
	except KeyError as e:
		if not printOnlySnpsThatBothGenomesSatisfy:
			print(e)
	except IndexError as e:
		try:
			print(m12Genes[geneName].sequence)
			print(k12Genes[geneName].sequence)
		except KeyError:
			print("failed index and key")
	if printSnp:
		print("-"*50)