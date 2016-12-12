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

    '''Takes the LRG code and fetches the relevant XML file from the EBI ftp server'''

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
    
    # Initialise empty dictionary
    geneData = {} 
    # Add HGNCID to geneData dictionary
    geneData['hgncID'] = tree.find('./fixed_annotation/hgnc_id').text 
    # Get data under LRG annotation set node
    annotation_set = tree.find("updatable_annotation/annotation_set[@type='lrg']")
    # Add Gene Symbol to geneData dictionary
    geneData['geneSymbol'] = annotation_set.find('lrg_locus').text
    # Get data under main assembly mapping node
    mapping = annotation_set.find("mapping[@type='main_assembly']")
    # Add genome build to geneData dictionary
    geneData['build'] = mapping.attrib['coord_system']
    # Add Chromosome name to geneData dictionary
    geneData['chrom'] = mapping.attrib['other_name']
    # Add genomic start position to geneData dictionary
    geneData['other_start'] = int(mapping.find('mapping_span').attrib['other_start']) # genomic DNA start pos
    # Add genomic end position to geneData dictionary
    geneData['other_end'] = int(mapping.find('mapping_span').attrib['other_end']) # genomic DNA end pos
    # Add strand orientation to geneData dictionary
    geneData['strand'] = mapping.find('mapping_span').attrib['strand'] # Starnd orientation (1 = forward, -1 = reverse)
    # return geneData dictionary
    return geneData

def getSeq(tree):
    '''
    Maps LRG sequence to main assembly sequence.
    Accounts for all differences (mismatch, insertions and deletions) between LRG and main assembly sequence.
    Returns main assembly DNA sequence and lookup dictionary to convert LRG sequence positions to main assembly sequence positions. 
    '''
    #Get complete LRG sequence and create copy that will be converted to main assembly sequence 
    lrgSeqFull = tree.find('./fixed_annotation/sequence').text
    mainAssemSeq = lrgSeqFull
    #Get end position of LRG sequence and check LRG sequence is expected length
    mainAssemData = tree.find("updatable_annotation/annotation_set[@type='lrg']/mapping[@type='main_assembly']/mapping_span")
    lrgEnd = int(mainAssemData.attrib['lrg_end'])
    assert (len(lrgSeqFull) == lrgEnd), "Discrepency between stated length of LRG sequennce and sequence provided"
    #Initialise lookup dictionary  with the LRG positions
    convertPosDict = {}
    for i in range(lrgEnd):
        convertPosDict[i+1] = i+1
    #Loop through differences between LRG and main assembly
    #Modify mainAssemSeq and converPosDict to account for each difference.
    for diff in mainAssemData.findall("diff"):
        #Convert LRG start and end positions of the variant to the current main assembly start and end positions
        startPos = convertPosDict[int(diff.attrib['lrg_start'])]
        endPos = convertPosDict[int(diff.attrib['lrg_end'])]
        #Get the LRG and main assembly sequence at the location of variant 
        lrgSeq = diff.attrib['lrg_sequence']
        otherSeq = diff.attrib['other_sequence']
        #If it's a mismatch...
        if diff.attrib['type'] == 'mismatch':
            #Check the stated LRG sequence matches the actual sequence at that position
            assert (mainAssemSeq[startPos-1:endPos] == lrgSeq), "Error when converting to genomic sequence. LRG sequence not as expected."
            #Add the change into the main assembly sequence
            mainAssemSeq = mainAssemSeq[:startPos-1] + otherSeq + mainAssemSeq[endPos:]
        #If it's an insertion into the main assembly...
        elif diff.attrib['type'] == 'other_ins':
            #Check there is no LRG sequence
            assert (lrgSeq == "-"), "Error when converting to genomic sequence. LRG sequence not as expected."
            #Add the insertion into the main assembly sequence
            mainAssemSeq = mainAssemSeq[:startPos] + otherSeq + mainAssemSeq[endPos-1:]
            #Update the lookup dictionary so all positions after the insertion are shifted by the insertion length 
            insLength = len(otherSeq)
            for i in range(endPos, lrgEnd+1):
                convertPosDict[i] += insLength 
        #If it's an insertion into the LRG sequence (effectively a deletion in main assembly)...
        elif diff.attrib['type'] == 'lrg_ins':
            #Check the stated LRG sequence matches the actual sequence at that position            
            assert (mainAssemSeq[startPos-1:endPos] == lrgSeq), "Error when converting to genomic sequence. LRG sequence not as expected."
            #Check there is no main assembly sequence 
            assert (otherSeq == "-"), "Error when converting to genomic sequence. LRG sequence not as expected."
            #Remove the deleted sequence from the main assembly sequence    
            mainAssemSeq = mainAssemSeq[:startPos-1] + mainAssemSeq[endPos:]
            delLength = len(lrgSeq)
            #For LRG positions that lie within deleted region, 
            #update the lookup dictionary so they map to the first base following the deletion in main assembly
            posConvDel = endPos - (delLength-1)
            for i in range(startPos, endPos+1):
                convertPosDict[i] = posConvDel
            #Update the lookup dictionary so all positions after the deletion are shifted by the deletion length 
            for i in range(endPos+1, lrgEnd+1):
                convertPosDict[i] -= delLength 
        else:
            raise ValueError("An unexpected diff type was encountered when converting to genomic sequence. Couldn't process LRG file.")
    return mainAssemSeq, convertPosDict


def getExons(tree, geneData, mainAssemSeq, convertPosDict):

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
            start = convertPosDict[int(start)]
            end = coordinates.attrib['end'] # LRG End position
            end = convertPosDict[int(end)]
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
            seq = mainAssemSeq[int(start)-1:int(end)]
            # Create list containing all details to be included in CSV row for that exon
            tmp = [coord_sys, geneData['build'], geneData['hgncID'], geneData['geneSymbol'], tx, exonnumber, geneData['chrom'], gDNA_start, gDNA_end, seq]
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
    mainAssemSeq, convertPosDict = getSeq(tree)
    geneData = getGeneLevData(tree)
    exonList = getExons(tree, geneData, mainAssemSeq, convertPosDict)
    headers = ['LRG_Number', 'Build', 'HGNCID', 'Gene', 'Transcript_ID', 'Exon_no', 'Chrom', 'Start', 'End', 'Sequence']
    writeCSV(headers, exonList, LRG)
    #for row in exonList:
    #    print(row)

#print("Done")
