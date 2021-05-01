# this adds seq data to annotation files that don't have it
import sys
from glob import glob
import os
import re
import threading
import urllib.request
import time
import random

pathToFiles = sys.argv[1]
extension = sys.argv[2]

class AddSeqDataThread(threading.Thread):

    def __init__(self, name, query, path, extension, fileName):
        threading.Thread.__init__(self)
        self.name = name
        self.query = query
        self.path = path
        self.encounteredError = False
        self.fileName = fileName
        self.extention = extension
    def run(self):
        fileToDownload = ""
        fileData = [] # list of strings each string is a line
        with open(self.fileName) as annotationFile:
            for line in annotationFile:
                fileData.append(line)
                line = line.strip()
                cols = line.split()
                if len(cols) == 0:
                    continue
                if cols[0] == "CONTIG": # only appears once in file
                    # CONTIG   join(fileName:startLocationInContigs..EndLocation)
                    indexOfFileNameEnd = cols[1].find(":")
                    fileToDownload = cols[1][5:indexOfFileNameEnd]
        with open(self.fileName, "w") as annotationFile:
            for line in fileData:
                # if line.split()[0] != "CONTIG": # this assumes that the loction to download doesn't cover multiple lines
                annotationFile.write(line)
            fullQuery = self.query + fileToDownload + "&retType=" + self.extention[1:] #, self.path + fileToDownload + self.extention
            try:
                # print(str(urllib.request.urlopen(fullQuery).decode('ascii')))
                annotationFile.write(urllib.request.urlopen(fullQuery).read().decode('ascii'))
            except Exception as error:
                self.encounteredError = True
                print("error caught",error)
                time.sleep(1/3 + random.uniform(0,2))
                annotationFile.write(urllib.request.urlopen(fullQuery).read().decode('ascii'))
                print("downloaded after error")

query = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id="
for fileName in glob(os.path.join(pathToFiles, "*" + extension)):
    print(fileName)
    thread = AddSeqDataThread(fileName,query, pathToFiles, extension, fileName)
    thread.start()
    break
    # I think this program isn't done because the other files have both.