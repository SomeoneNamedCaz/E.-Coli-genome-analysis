import sys
sys.path.insert(1, '/Users/cazcullimore/dev/ericksonLabCode/')
sys.path.insert(1, '/scr/snpalign/secondaryPythonScripts')
import unittest
from scr.snpalign.secondaryPythonScripts.mapSNPsToGenome import *


class MyTestCase(unittest.TestCase):
	def test_something(self):
		writeMappedSnps()
		

if __name__ == '__main__':
	unittest.main()
