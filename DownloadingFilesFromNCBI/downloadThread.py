import urllib.request
import urllib.error
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
            # if error == urllib.error.HTTPError("HTTP Error 400: Bad Request"):
            #     self.encounteredError = True
            #     self.secondTryWorked = False
            #     print("recognized bad request")
            #     return

            print("error caught",error)
            time.sleep(5/3 + random.uniform(0,2))
            self.encounteredError = True
            self.secondTryWorked = False
            urllib.request.urlretrieve(self.query + self.ID + self.retType, self.pathAndFileName)
            print("downloaded after error")
            self.secondTryWorked = True