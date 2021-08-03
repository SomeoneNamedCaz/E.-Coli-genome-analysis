from downloadThread import *

exit(0)
searchBeginning = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id="
downloadPath = "./AllCommensalStrains/"
last5Times = []
threads = []
with open("allCommensalEColiIDs.xml") as IDFile:
    inIDList = False
    IDnotDownloaded = True
    readyToStart = False
    for line in IDFile:
        line = line.strip()
        if line == '<IdList>' or line == '</IdList>':
            inIDList = not inIDList
        elif inIDList:

            ID = re.sub(r"</?Id>", "", line)
            if ID == "":
                readyToStart = True
            if readyToStart:
                try:
                    with open(downloadPath+ID+".gb") as annotationFile:
                        IDnotDownloaded = False
                        # for line in annotationFile:
                        #     if line.split()[0] == "CONTIG":
                except:
                    IDnotDownloaded = True

                if IDnotDownloaded:
                    t1 = time.time()
                    try:
                        downloadThread = DownloadThread(ID, ID, searchBeginning,downloadPath)
                        # thread.append(downloadThread)
                        downloadThread.start()
                        threads.append(downloadThread)
                    except Exception as exception:
                        print("thread error")
                    print("downloaded", ID)
                    for threadIndex in range(len(threads)-1,-1,-1):
                        thread = threads[threadIndex]
                        if thread.encounteredError:
                            # print(thread.is_alive())
                            time.sleep(10)
                        if not thread.is_alive():
                            threads.pop(threadIndex)
                    print(threading.activeCount())
                    # while threading.activeCount() >= 5:
                    #     time.sleep(0.1)
                    time.sleep(1.3 / 3)

