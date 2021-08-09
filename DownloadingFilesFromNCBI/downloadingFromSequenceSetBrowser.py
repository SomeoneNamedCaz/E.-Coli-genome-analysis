from downloadThread import *

threads = []
with open("WGSProjectNamesWithExtraStuffRemoved.tsv") as wgsFile:
    readyToStart = False
    for line in wgsFile:
        cols = line.split("\t")
        if cols[0][0] != "D" and cols[0][0] != "M":
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
            numberOfFiles = 1
        # print(stopFileName[-5])
        # if startFileName == "NZ_LOOA01000001":
        readyToStart = True
        if readyToStart:
            for i in range(numberOfFiles):
                if numberOfFiles != 1:
                    indexAsString = str(i + 1)
                    # add to 5 len
                    while len(indexAsString) < 5:
                        indexAsString = "0" + indexAsString
                    contigName =  startFileName[3:-5] + indexAsString + ".1"
                else:
                    contigName = startFileName

                t1 = time.time()
                print(stopFileName)
                print(contigName)
                url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&term="# + contigName + "&retType=gb"
                downloadThread = DownloadThread(contigName, contigName, url, "./MissedMastitis/" + contigName + ".xml", "")
                threads.append(downloadThread)
                downloadThread.start()
                for threadIndex in range(len(threads) - 1, -1, -1):
                    thread = threads[threadIndex]
                    if thread.encounteredError:
                        # print(thread.is_alive())
                        time.sleep(12) # so it ncbi doesn't get mad
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
                time.sleep(1.5/3)
    # except IndexError:
    #     0
    # elif :
    # if cols[0] == "WGS":
    #     file cols[1]
    #     "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id="