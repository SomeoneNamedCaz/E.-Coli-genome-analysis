import unittest
import sys

from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
# from Bio
from Bio.Entrez import *
import re
sys.path.insert(1, '/Users/cazcullimore/Documents/ericksonLabCode/secondaryPythonScripts')
try:
    from secondaryPythonScripts.functions import *
except ImportError:
    pass

class MyTestCase(unittest.TestCase):

    def testReadIn1ContigFasta(self):
        fastaPath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/unitTests/k-12.fasta"
        fileData = readInFastaAsList(fastaPath)
        print(len(fileData))
        with open(fastaPath) as file:
            entryData = ""
            for line in file:
                if line[0] == ">":
                    self.assertEqual(fileData[0], line.strip())
                    fileData.pop(0)
                    if entryData != "":
                        self.assertEqual(fileData[0], entryData)
                        fileData.pop(0)
                        entryData = ""
                else:
                    entryData += line.strip()
        if entryData != "":
            self.assertEqual(fileData[0], entryData)
            fileData.pop(0)
        self.assertEqual(fileData, [])# should be empty

    def testReadInAssemblyFasta(self):
        fastaPath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/tenAssembliesFromEachCategory/1465_SS_220.fasta"
        fileData = readInFastaAsList(fastaPath)
        print(len(fileData))
        with open(fastaPath) as file:
            entryData = ""
            for line in file:
                if line[0] == ">":
                    if entryData != "":
                        self.assertEqual(fileData.pop(0), entryData)
                        entryData = ""
                    self.assertEqual(fileData.pop(0), line.strip())
                else:
                    entryData += line.strip()
        if entryData != "":
            self.assertEqual(fileData[0], entryData)
            fileData.pop(0)
        self.assertEqual(fileData, [])# should be empty


    def testGetAssemblyContigs(self):
        gbPath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/unitTests/k-12.gbff"
        fastaPath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/unitTests/k-12.fasta"
        index = 0
        gbContigs = getContigs(gbPath)
        fastaContigs = readInFastaAsList(fastaPath)[1:][::2]
        self.assertEqual(len(gbContigs), len(fastaContigs))
        for gbContig in gbContigs:
            foundMatch = False
            for fastaContig in fastaContigs:
                if gbContig == fastaContig:
                    foundMatch = True
            if not foundMatch:
                raise Exception()

    def testGetFullGenomeContig(self):
        gbPath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/unitTests/k-12.gbff"
        fastaPath = "/Users/cazcullimore/Documents/ericksonLabCode/tests/unitTests/k-12.fasta"
        index = 0
        self.assertEqual(getContigs(gbPath), readInFastaAsList(fastaPath)[1:][::2])

    def testGetGenes(self): # test
        # this looks at if the genes are similar to what they should look like
        testGb = "/Users/cazcullimore/Documents/ericksonLabCode/tests/testFiles/gbks/1465_SS_220.fasta.gbk"
        
        genes = getGenesOnContigs(testGb, getContigs(testGb))
        for geneName, geneSeq in genes.items():
            blastResult = NCBIWWW.qblast("blastn", "nt", geneSeq, alignments=10)
            parsedResult = NCBIXML.parse(blastResult)
            item = next(parsedResult)
            seq_record = SeqIO.read(item, "gb")
            nucleotide_accession = seq_record.annotations["db_source"]

            nucl_id = nucleotide_accession.split()[-1]

            for alignment in item.alignments:
                for hsp in alignment.hsps:
                    print(hsp)
            # blastResult
            handle = efetch(db="nucleotide", id=item.alignments[0].hit_id, rettype="gb",retmode="text", query_key="aroC", email="cazvash9.student.byu.edu")
            # for xmlLine in blastResult:
            #     if "<Hsp_hseq>" in xmlLine:
            #         re.findall(r"<Hsp_hseq>(.+)</Hsp_hseq>")
            # print("".join(blastResult))
        print("done")

    # def testTwoFilesHaveSameInfo(self):
    #     filesContainTheSameInformation()


if __name__ == '__main__':
    unittest.main()
