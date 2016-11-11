import urllib2
import xml.etree.ElementTree as etree
import os
import csv
import sys
import pickle


def checkGeneDict(usergene):

    '''Opens a pickle dictionary containing all of the current (Nov 2016) LRG codes and their corresponding genes
    and returns the correct LRG code for the user specified gene'''

    with open('genedict.pickle', 'rb') as handle:
        genedict = pickle.load(handle)

    for LRG, genename in genedict.iteritems():
        if usergene == genename:
            return LRG

def checkValidLRG(userLRG):

    '''Opens a pickle dictionary containing all of the current (Nov 2016) LRG codes and their corresponding genes
    and checks that the specified LRG code exists in the dict, then returns the specified LRG code'''

    with open('genedict.pickle', 'rb') as handle:
        genedict = pickle.load(handle)

    if userLRG in genedict:
        return userLRG
    else:
        print("Invalid LRG code")
        exit()

def readLRG(lrg):

    '''Takes the LRG code and fetches the relevant XML file from the EBI ftp site'''

    lrgURL = 'http://ftp.ebi.ac.uk/pub/databases/lrgex/' + lrg + '.xml'
    try:
        tree = etree.parse(urllib2.urlopen(lrgURL))
    #Catch any HTTP errors and exit program with error message if site can't be reached.
    except urllib2.HTTPError, e:
        print("Unable to retrieve LRG file from web")
        return #return a Null value so that program does not continue
    else:
        #Check that the file contains the lrg tag and conforms to version 1.9 standards.
        #An error message is displayed and program exits if these conditions are not met.
        if tree.getroot().tag == 'lrg' and tree.getroot().attrib['schema_version'] == "1.9":
            return tree
        else:
            print("File does not conform to lrg version 1.9 specification.")
            return #return a Null value so that program does not continue

def getGeneLevData(tree):

    '''Retrieves data that applies to whole LRG file (as opposed to individual exons) and returns data in a dictionary'''

    geneData = {} # Initialise empty dictionary
    geneData['hgncID'] = tree.find('./fixed_annotation/hgnc_id').text # HGNCID
    for annotation_set in tree.findall('updatable_annotation/'):
        if annotation_set.attrib['type'] == 'lrg':
            geneData['geneSymbol'] = annotation_set.find('lrg_locus').text # Gene Symbol
            for mapping in annotation_set.findall('mapping'):
                if mapping.attrib['type'] == 'main_assembly':
                    geneData['build'] = mapping.attrib['coord_system'] # Genome build
                    geneData['chrom'] = mapping.attrib['other_name'] # Chromosome
                    geneData['other_start'] = int(mapping.find('mapping_span').attrib['other_start']) # genomic DNA start pos
                    geneData['other_end'] = int(mapping.find('mapping_span').attrib['other_end']) # genomic DNA end pos
                    geneData['strand'] = mapping.find('mapping_span').attrib['strand'] # Starnd orientation (1 = forward, -1 = reverse)
    return geneData

def getExons(tree, dnaSeq, geneData):

    '''
    Retrieves exon specific data and stores data for every exon in a list of lists, along with gene level data,
    so that it can easily be converted to csv file
    '''

    exonList = [] # Initialise empty list to hold exon data
    for transcript in tree.findall('./fixed_annotation/transcript'):
        tx = transcript.attrib['name'] # Transcript
        # Loop through exons in each transcript and get details
        for exon in transcript.findall('./exon'):
            exonnumber = exon.attrib['label'] # Exon Number
            coordinates = exon.find('coordinates')
            coord_sys = coordinates.attrib['coord_system'] # LRG number
            start = coordinates.attrib['start'] # LRG start position
            end = coordinates.attrib['end'] # LRG End position
            # Convert LRG coordinates to genomic coordinates
            # If the gene is located on the forward strand...
            if geneData['strand'] == "1":
                gDNA_start = (int(start) + int(geneData['other_start'])) - 1
                gDNA_end = (int(end) + int(geneData['other_start'])) - 1
            # If the gene is located on the reverse strand...
            elif geneData['strand'] == "-1":
                gDNA_start = (int(geneData['other_end']) - int(start)) + 1
                gDNA_end = (int(geneData['other_end']) - int(end)) + 1
            # Use the LRG coordinates to take a slice the DNA and retrieve exon sequence
            seq = dnaSeq[int(start)-1:int(end)]
            # Create list containing all details to be included in CSV row for that exon
            tmp = [coord_sys, geneData['build'], tx, geneData['hgncID'], geneData['geneSymbol'], exonnumber, geneData['chrom'], gDNA_start, gDNA_end, seq]
            # Append list to the main exon list.
            exonList.append(tmp)
    return exonList

def writeCSV(headers, exonList, LRG):

    '''Writes the output to a CSV file'''

    exonList.insert(0, headers)
    outputFile = os.getcwd() + "/" + LRG + "_output.csv"

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
    LRG = checkGeneDict(sys.argv[2])
elif sys.argv[1] == '-l':
    LRG = checkValidLRG(sys.argv[2])

tree = readLRG(LRG)
if tree:
    dnaSeq = tree.find('./fixed_annotation/sequence').text
    geneData = getGeneLevData(tree)
    exonList = getExons(tree, dnaSeq, geneData)
    headers = ['LRG_Number', 'Build', 'HGNCID', 'Gene', 'Transcript_ID', 'Exon_no', 'Chrom', 'Start', 'End', 'Sequence']
    writeCSV(headers, exonList, LRG)
    #for row in exonList:
    #    print(row)

#print("Done")
