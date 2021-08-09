import urllib.request
import re
import time
import threading
import random
import glob
import os

class DownloadThread(threading.Thread):

    def __init__(self, name, ID, query, pathAndFileName, retType= "&retType=gb"):
        threading.Thread.__init__(self)
        self.name = name
        self.ID = ID
        self.query = query
        self.pathAndFileName = pathAndFileName
        self.encounteredError = False
        self.retType = retType
    def run(self):
        try:
            urllib.request.urlretrieve(self.query + self.ID + self.retType, self.pathAndFileName)
        except Exception as error:
            self.encounteredError = True
            print("error caught",error)
            time.sleep(2/3 + random.uniform(0,2))
            urllib.request.urlretrieve(self.query + self.ID + self.retType, self.pathAndFileName)
            print("downloaded after error")