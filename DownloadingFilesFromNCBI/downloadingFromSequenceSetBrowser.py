from downloadThread import *
import sys
import os

threads = []
examplesOfLengthsRequiredForUsingFirstCol = {"DA":"DADRKS010000001", "MI":"MIWM01000001", "LO":"LOOH01000129"} # without .1

outputPathPrefix = sys.argv[1]#"./AllMastitis/"

if not os.path.exists(outputPathPrefix):
    os.mkdir(outputPathPrefix)

with open(sys.argv[2]) as wgsFile:
    readyToStart = False
    for line in wgsFile:
        cols = line.split(",")
        if cols[0] == "prefix_s":
            continue
        # print(cols)
        # try:
        #     print(cols[18])
        if len(cols[-4].split("-")) == 2:
            startFileName, stopFileName = cols[-4].split("-")[0],cols[-4].split("-")[1]
        else:
            startFileName = cols[0]
            stopFileName = startFileName
        try:
            numberOfFiles = int(stopFileName[-5:])
        except ValueError:
            numberOfFiles = int(1e5)
        # print(stopFileName[-5])
        # if startFileName == "NZ_LOOA01000001":
        readyToStart = True
        if readyToStart:
            for i in range(numberOfFiles):
                # get file info for downloading
                if numberOfFiles != int(1e5):
                    indexAsString = str(i + 1)
                    # add to 5 len
                    while len(indexAsString) < 5:
                        indexAsString = "0" + indexAsString
                    contigName =  startFileName[3:-5] + indexAsString + ".1"
                else:
                    indexAsString = str(i + 1)
                    #      add while the length of the current end file isn't as long as it should be
                    while len(indexAsString) + len(startFileName) < len(
                            examplesOfLengthsRequiredForUsingFirstCol[startFileName[:2]]):
                        indexAsString = "0" + indexAsString
                    contigName = startFileName + indexAsString + ".1"

                t1 = time.time()
                print(stopFileName)
                print(contigName)
                # check if downloaded
                outputPath = outputPathPrefix + contigName + ".fasta"
                try:
                    with open(outputPath) as annotationFile:
                        IDnotDownloaded = False
                        # for line in annotationFile:
                        #     if line.split()[0] == "CONTIG":
                except:
                    IDnotDownloaded = True

                if IDnotDownloaded:
                    print("downloading")
                    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id="# + contigName + "&retType=gb"
                    downloadThread = DownloadThread(contigName, contigName, url, outputPath, "&retType=fasta")
                    threads.append(downloadThread)
                    downloadThread.start()
                    numEncounteredError = 0
                    for threadIndex in range(len(threads) - 1, -1, -1):
                        thread = threads[threadIndex]
                        if thread.encounteredError:
                            # this is kind of a wierd way to deal with this problem. It might be better to just break at any HTTP
                            # bad access error, but this allows for some failures if the internet or ncbi is cutting out at times,
                            # which was happening earlier
                            numEncounteredError += not thread.secondTryWorked # one if didn't work
                            time.sleep(12) # so it ncbi doesn't get mad
                        if not thread.is_alive():
                            threads.pop(threadIndex)
                    if numEncounteredError >= 1:
                        while len(threads) > 0:
                            for threadIndex in range(len(threads) - 1, -1, -1):
                                thread = threads[threadIndex]
                                if not thread.is_alive():
                                    threads.pop(threadIndex)
                            time.sleep(12)
                        break
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
                    if numberOfFiles == int(1e5):
                        time.sleep(8 / 3)
                    time.sleep(2/3)
    # except IndexError:
    #     0
    # elif :
    # if cols[0] == "WGS":
    #     file cols[1]
    #     "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id="