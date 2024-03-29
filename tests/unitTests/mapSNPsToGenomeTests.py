import sys
sys.path.insert(1, '/Users/cazcullimore/dev/ericksonLabCode/')
sys.path.insert(1, '/Users/cazcullimore/dev/ericksonLabCode/secondaryPythonScripts')
import unittest
from secondaryPythonScripts.mapSNPsToGenome import *
from secondaryPythonScripts.functions import *


class MyTestCase(unittest.TestCase):
	def test_something(self):
		writeMappedSnps()
		

if __name__ == '__main__':
	unittest.main()
