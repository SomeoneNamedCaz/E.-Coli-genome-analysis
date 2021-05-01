import sys
import urllib
import webbrowser

# nuclotide search doesn't return enough data as it should
# bioProjectsFileName = sys.argv[1] # as csv
nucleotideSearchQuery = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=(bovine%20mastitis)%20AND%20%22Escherichia%20coli%22[porgn:__txid562]"

# with open(nucleotideSearchQuery) as bioProjectsFile:
#     for line in bioProjectsFile:
#         cols = line.split(",")