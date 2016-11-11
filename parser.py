import xml.etree.ElementTree as etree
import os
import csv
import sys
import pickle

def checkGeneDict(usergene):
    
    with open('genedict.pickle', 'rb') as handle:
        genedict = pickle.load(handle)

    for LRG, genename in genedict.iteritems():
        if usergene == genename:
            return LRG

def readLRG(lrg):
    
    filepath = '../LRG_Data/' + lrg + '.xml'
    if os.path.isfile(filepath) == False:
         print("Please specify a valid LRG number")
         return
    else:
        tree = etree.parse(filepath)
        if tree.getroot().tag == 'lrg' and tree.getroot().attrib['schema_version'] == "1.9":
            return tree
        else:
            print("File does not conform to lrg version 1.9 specification.")
            return

def getGeneLevData(tree):
   
    '''Retrieves data that applies to whole LRG file (as opposed to individual exons) and returns in a dictionary'''

    geneData = {} # Initialise empty dictionary  
    for annotation_set in tree.findall('updatable_annotation/'):
        if annotation_set.attrib['type'] == 'lrg':
            for mapping in annotation_set.findall('mapping'):
                if mapping.attrib['type'] == 'main_assembly':
                    geneData['build'] = mapping.attrib['coord_system']
                    geneData['chrom'] = mapping.attrib['other_name']
                    geneData['other_start'] = int(mapping.find('mapping_span').attrib['other_start'])
                    geneData['other_end'] = int(mapping.find('mapping_span').attrib['other_end'])
                    geneData['strand'] = mapping.find('mapping_span').attrib['strand']
                    return geneData

def getExons(tree, dnaSeq, geneData):
    
    exonList = []
    for transcript in tree.findall('./fixed_annotation/transcript'):
        tx = transcript.attrib['name']
        for exon in transcript.findall('./exon'):
            exonnumber = exon.attrib['label']
            coordinates = exon.find('coordinates')
            coord_sys = coordinates.attrib['coord_system']
            start = coordinates.attrib['start']
            end = coordinates.attrib['end']
            if geneData['strand'] == "1":
                gDNA_start = (int(start) + int(geneData['other_start'])) - 1
                gDNA_end = (int(end) + int(geneData['other_start'])) - 1
            elif geneData['strand'] == "-1":
                gDNA_start = (int(geneData['other_end']) - int(start)) + 1
                gDNA_end = (int(geneData['other_end']) - int(end)) + 1
            seq = dnaSeq[int(start)-1:int(end)]
            tmp = [coord_sys, geneData['build'], tx, exonnumber, geneData['chrom'], gDNA_start, gDNA_end, seq]
            exonList.append(tmp)
    return exonList

def writeCSV(headers, exonList):
    
    exonList.insert(0, headers)
    outputFile = os.getcwd() + "/output.csv"
    
    # Write the output file
    with open(outputFile, 'wb') as f:
        writer = csv.writer(f)
        writer.writerows(exonList)

###############################################################################

# Check that the user has provided the correct number of arguments
assert (len(sys.argv) == 3), "Incorrect number of arguments"

# Check that the user has provided one of the two correct flags
assert (sys.argv[1] in ['-g', '-l']), "Wrong flag"

# Depending on whether the user has used the '-g' or '-l' flag...
if sys.argv[1] == '-g':
    print(checkGeneDict(sys.argv[2]))
elif sys.argv[1] == '-l':
    print("You used the -l flag!")

# tree = readLRG(sys.argv[1])
# if tree:
#     dnaSeq = tree.find('./fixed_annotation/sequence').text
#     geneData = getGeneLevData(tree)
#     exonList = getExons(tree, dnaSeq, geneData)
#     headers = ['LRG_Number', 'Build', 'Transcript_ID', 'Exon_no', 'Chrom', 'Start', 'End', 'Sequence']
#     # writeCSV(headers, exonList)
#     for row in exonList:
#         print(row)


# for i in range(250):
#     lrg = 'LRG_' + str(i)
#     print(lrg)
#     tree = readLRG(lrg)
#     if tree:
#         dnaSeq = tree.find('./fixed_annotation/sequence').text
#         build, chrom, conv_offset = getBuildChrom(tree)
#         exonList = getExons(tree, dnaSeq, build, chrom, conv_offset)
#         headers = ['LRG_Number', 'Build', 'Transcript_ID', 'Exon_no', 'Chrom', 'Start', 'End', 'Sequence']
#         # writeCSV(headers, exonList)
#         for row in exonList:
#             print(row)

# dnaSeq = tree.find('./fixed_annotation/sequence').text

# mylist = []

# build, chrom, conv_offset = getBuildChrom()

# for row in mylist:
# 	print(row)

# Create the headers
# headers = ['Transcript_ID', 'Exon_no', 'Coordinate_system', 'Start', 'End', 'Sequence']

# Add headers to the output
# mylist.insert(0, headers)

# Define the output file
# myoutput = os.getcwd() + "/output.csv"

# Write the output file
# with open(myoutput, 'wb') as f:
# 	writer = csv.writer(f)
# 	writer.writerows(mylist)

# print("Done")
