import urllib.request
import re
import time
import threading
import random
import glob
import os

#
class DownloadThread(threading.Thread):

    def __init__(self, name, ID, query, pathAndFileName):
        threading.Thread.__init__(self)
        self.name = name
        self.ID = ID
        self.query = query
        self.pathAndFileName = pathAndFileName
        self.encounteredError = False
    def run(self):
        try:
            urllib.request.urlretrieve(self.query + self.ID + "&retType=gb", self.pathAndFileName)
        except Exception as error:
            self.encounteredError = True
            print("error caught",error)
            time.sleep(2/3 + random.uniform(0,2))
            urllib.request.urlretrieve(self.query + self.ID + "&retType=gb", self.pathAndFileName)
            print("downloaded after error")
#         # self.name.exit()
# searchBeginning = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id="
# downloadPath = "./AllCommensalStrains/"
# last5Times = []
# threads = []
# with open("allCommensalEColiIDs.xml") as IDFile:
#     inIDList = False
#     IDnotDownloaded = True
#     readyToStart = False
#     for line in IDFile:
#         line = line.strip()
#         if line == '<IdList>' or line == '</IdList>':
#             inIDList = not inIDList
#         elif inIDList:
#
#             ID = re.sub(r"</?Id>", "", line)
#             if ID == "":
#                 readyToStart = True
#             if readyToStart:
#                 try:
#                     with open(downloadPath+ID+".gb") as annotationFile:
#                         IDnotDownloaded = False
#                         for line in annotationFile:
#                             if line.split()[0] == "CONTIG":
#
#                 except:
#                     IDnotDownloaded = True
#
#                 if IDnotDownloaded:
#                     t1 = time.time()
#                     try:
#                         downloadThread = DownloadThread(ID, ID, searchBeginning,downloadPath)
#                         # thread.append(downloadThread)
#                         downloadThread.start()
#                         threads.append(downloadThread)
#                     except Exception as exception:
#                         print("thread error")
#                     print("downloaded", ID)
#                     for threadIndex in range(len(threads)-1,-1,-1):
#                         thread = threads[threadIndex]
#                         if thread.encounteredError:
#                             # print(thread.is_alive())
#                             time.sleep(10)
#                         if not thread.is_alive():
#                             threads.pop(threadIndex)
#                     print(threading.activeCount())
#                     # while threading.activeCount() >= 5:
#                     #     time.sleep(0.1)
#                     time.sleep(1.3 / 3)
threads = []
with open("WGSProjectNamesWithExtraStuffRemoved.tsv") as wgsFile:
    readyToStart = False
    for line in wgsFile:
        cols = line.split("\t")
        # print(cols)
        # try:
        #     print(cols[18])
        if len(cols[-4].split("-")) == 2:
            startFileName, stopFileName = cols[-4].split("-")[0],cols[-4].split("-")[1]
            numberOfFiles = int(stopFileName[-5:])
            # print(stopFileName[-5])
            # if startFileName == "NZ_LOOA01000001":
            readyToStart = True
            if readyToStart:
                for i in range(numberOfFiles):
                    indexAsString = str(i + 1)
                    # add to 5 len
                    while len(indexAsString) < 5:
                        indexAsString = "0" + indexAsString
                    contigName =  startFileName[3:-5] + indexAsString + ".1"
                    t1 = time.time()
                    print(stopFileName)
                    print(contigName)
                    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id="# + contigName + "&retType=gb"
                    downloadThread = DownloadThread(contigName, contigName, url, "./DownloadedFromSSB/" + contigName + ".gb")
                    threads.append(downloadThread)
                    downloadThread.start()
                    for threadIndex in range(len(threads) - 1, -1, -1):
                        thread = threads[threadIndex]
                        if thread.encounteredError:
                            # print(thread.is_alive())
                            time.sleep(12)
                        if not thread.is_alive():
                            threads.pop(threadIndex)
                    # print(url)
                    # try:
                    #     urllib.request.urlretrieve(url, "./DownloadedFromSSB/"+ contigName + ".gb")
                    # except:
                    #     time.sleep(3)
                    #     urllib.request.urlretrieve(url, "./DownloadedFromSSB/" + contigName + ".gb")
                    # print("downloaded")
                    # # time.sleep(1)
                    # print("sleep time",max(2/3 - (time.time() - t1),0))
                    # print("time to download", time.time() - t1)
                    # time.sleep(max(1/3 - (time.time() - t1),0))
                    time.sleep(2/3)
        # except IndexError:
        #     0
        # elif :
        # if cols[0] == "WGS":
        #     file cols[1]
        #     "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id="
#
# def DownloadFromProjectName(concatenatedFileName, extraWaitTime= 0):
#     startFileName = concatenatedFileName.split("-")[0]
#     stopFileName = concatenatedFileName.split("-")[1]
#     numberOfFiles = int(stopFileName[-5:])
#     # print(stopFileName[-5])
#     for i in range(numberOfFiles):
#         indexAsString = str(i + 1)
#         # add to 5 len
#         while len(indexAsString) < 5:
#             indexAsString = "0" + indexAsString
#         contigName = startFileName[:-5] + indexAsString
#         t1 = time.time()
#         url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=" + contigName + "&retType=gb"
#         urllib.request.urlretrieve(url, "./DownloadedFromSSB/" + contigName + ".gb")
#         time.sleep(max(1 / 3 + extraWaitTime - (time.time() - t1), 0))
#
# def DownloadingSeqsFromNuc(nucIDfile, whereStartID=-1): # could do if master then don't read it or could read the corresponding files WGS fileStart-fileStop
#     with open(nucIDfile) as idFile:
#
#             inIDList = False
#             for line in IDFile:
#                 line = line.strip()
#                 if line == '<IdList>':
#                     inIDList = not inIDList
#                 elif inIDList:
#                     ID = re.sub(r"</?Id>", "", line)
#                     t1 = time.time()
#                     urllib.request.urlretrieve(searchBeginning + ID + "&retType=gb", downloadPath + ID + ".gb")
#                     time.sleep(max(1 - (time.time() - t1), 1))
#
#
# def DownloadingSequenceFilesFromMasterFiles(pathToSeqAndMasters):
#     onFirstLine = True
#     isMasterRecord = False
#     for line in idFile:
#         line = line.strip()
#         cols = line.split()
#         if onFirstLine and cols[1][-1] == "0":  # if master record
#             isMasterRecord = True
#         if isMasterRecord and cols[0] == "WGS":
#             elif not isMasterRecord:
#
#             onFirstLine = False