"""use a random forest to get the most significant groups of snps"""
import sys
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier

from secondaryPythonScripts.functions import *

def trainRFModel(data):
	numGenomes = 0
	percentValidation = 0.1
	numIter = 1
	
	
	x = []
	y = []
	nucToNum = {"A":0,"C":1,"T":2,"G":3, "-":4, "N":5}
	metadataToNum = {"commensal":0,"pathogen": 1}
	
	# with open(sys.argv[1]) as dataFile:
	for line in data:
		cols = line.split(",")
		numGenomes += 1
		x.append([nucToNum[elt] for elt in cols[1:-1]])
		y.append(metadataToNum[cols[-1].strip()])
	
	x = np.array(x)
	y = np.array(y)
	
	
	x_train = x[:int((1-percentValidation)*numGenomes)]
	x_test = x[int((1-percentValidation)*numGenomes):]
	
	
	y_train = y[:int((1-percentValidation)*numGenomes)]
	y_test = y[int((1-percentValidation)*numGenomes):]
	
	rf = MLPClassifier(hidden_layer_sizes=(100,100))
	for i in range(numIter):
		rf.fit(x_train, y_train)
		print("train acc", sum([y_hat == y_real for y_hat, y_real in zip(rf.predict(x_train), y_train)])/len(y_train))
		print("test acc", sum([y_hat == y_real for y_hat, y_real in zip(rf.predict(x_test), y_test)])/len(y_test))
		

# if __name
# if len(sys.argv) < 2:
# 	print("please provide the phylogroupSnpFile from makeMLFiles.py")