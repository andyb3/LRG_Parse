import xml.etree.ElementTree as etree
import os
import pickle

LRG_dir = "/Users/Paul/STP/LRG/LRG_Data/"

out_dir = "/Users/Paul/STP/LRG/LRG_Parse/"

# genedict = {}

# for files in os.listdir(LRG_dir):
# 	print(files)
# 	LRG_ID = files.replace(".xml", "")
# 	tree = etree.parse(LRG_dir + files)
# 	for annotation_set in tree.findall('updatable_annotation/'):
# 		if annotation_set.attrib['type'] == 'lrg':
# 			genename = annotation_set.find('lrg_locus').text
# 			genedict[LRG_ID] = genename

# with open(out_dir + 'genedict.pickle', 'wb') as handle:
# 	pickle.dump(genedict, handle)

# print("Done")

with open('genedict.pickle', 'rb') as handle:
	foo = pickle.load(handle)

for k, v in foo.iteritems():
	print(k + ": " + v)