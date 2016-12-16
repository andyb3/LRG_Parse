'''
Authors:        Paul Acklam and Andy Bond
Created:        December 2016
Description:    Unit testing suite for LRGparser.py
Usage:          Clone entire repo from https://github.com/andyb3/LRG_Parse and run the tests.py script.
                See README for more information
'''


import unittest
import LRGparser
import xml.etree.ElementTree as etree
import pickle

class TestStringMethods(unittest.TestCase):

    def test_checkGeneDict(self):
        '''Tests that checkGeneDict() function produces expected output'''
        # With valid gene symbol...
        self.assertEqual(LRGparser.checkGeneDict('BRCA1'), "LRG_292", "checkGeneDict() did not return correct LRG for BRCA1")
        # With invalid gene symbol...
        self.assertEqual(LRGparser.checkGeneDict('xxxxxxx'), None, "checkGeneDict() did not return Null value when invalid gene symbol provided")

    def test_checkValidLRG(self):
        '''Tests that checkValidLRG() function produces expected output'''
        # With valid LRG...
        self.assertEqual(LRGparser.checkValidLRG('LRG_1'), 'LRG_1', "valid LRG did not pass checkValidLRG() validation")
        # With invalid LRG...
        self.assertEqual(LRGparser.checkValidLRG('xxxxxxx'), None, "invalid LRG was not caught by checkValidLRG() validation")

    def test_readLRG(self):
        '''Tests that readLRG() function produces expected output'''
        # With valid URL for LRG FTP server...
        tree = LRGparser.readLRG('LRG_1')
        self.assertTrue(tree, "No element tree object returned from readLRG() when requested")
        self.assertEqual(tree.find('./fixed_annotation/hgnc_id').text, '2197', "Element tree object does not represent correct LRG file")
        # With invalid URL for LRG FTP server
        tree = LRGparser.readLRG('xxxxxxx')
        self.assertFalse(tree, "readLRG() did not return Null value when unable to download LRG file")

    def test_getGeneLevData(self):
        '''Tests that getGeneLevData() function produces expected output'''
        tree = etree.parse('./Test_Files/LRG_1.xml')
        geneData = LRGparser.getGeneLevData(tree)
        # Check first and last pieces of data extracted by function are correct
        self.assertEqual(geneData['hgncID'], '2197')
        self.assertEqual(geneData['strand'], '-1')

    def test_getGenomicSeq(self):
        '''Tests that getGenomicSeq() function produces expected output'''
        # BRCA1 file contains multiple mismatches, insertions and deletions between LRG and genomic main assembly sequence
        # Test that getGenomicSeq can convert from LRG to genomic main assembly correctly
        tree = etree.parse('./Test_Files/LRG_292.xml')
        mainAssemSeq, convertPosDict = LRGparser.getGenomicSeq(tree)
        handle = open('./Test_Files/BRCA_38_Seq.txt', 'r')
        BRCA_38_Seq = handle.read()
        handle.close()
        # Check that the sequence returned matches the reference sequence
        self.assertEqual(mainAssemSeq, BRCA_38_Seq, "Sequence output from getGenomicSeq() is incorrect")
        # Check that the position conversion dictionary is correct before and after each indel
        self.assertEqual(convertPosDict[13168], 13168, "Error in convertPosDict")
        self.assertEqual(convertPosDict[13171], 13169, "Error in convertPosDict")
        self.assertEqual(convertPosDict[15459], 15457, "Error in convertPosDict")
        self.assertEqual(convertPosDict[15460], 15459, "Error in convertPosDict")
        self.assertEqual(convertPosDict[15460], 15459, "Error in convertPosDict")

    def test_getExons(self):
        '''Tests that getGenomicSeq() function produces expected output'''
        # Test with the following LRG files to test full range of program features
        test_files = {'COL1A1':'./Test_Files/LRG_1.xml', #Reverse strand, 1 transcript, No differences between LRG/main assembly.
                      'BRCA1':'./Test_Files/LRG_292.xml', #Reverse strand. 1 transcript, numerous mismatch and indels between LRG/main assembly
                      'NF1':'./Test_Files/LRG_214.xml',#Forward strand, 2 transcripts, No differences between LRG/main assembly
                      'IndelTest':'./Test_Files/LRG_IndelTest.xml'#Modified to include exon start, exon stop and whole exon deletions in main assembly (compared to LRG)
                      }
        for gene, filepath in test_files.iteritems():
            tree = etree.parse(filepath)
            mainAssemSeq, convertPosDict = LRGparser.getGenomicSeq(tree)
            geneData = LRGparser.getGeneLevData(tree)
            exonList = LRGparser.getExons(tree, geneData, mainAssemSeq, convertPosDict)
            handle = open('./Test_Files/'+ gene + "_exonList.pickle", 'rb')
            refExonList = pickle.load(handle)
            handle.close()
            # Test the output lists match the expected lists
            # The genomic coordinates and sequences have all been checked against Ensembl genome browser
            self.assertEqual(exonList, refExonList, "getExons() output doesn't match reference for gene: " + gene)

if __name__ == '__main__':
    print "\nNOTE: IGNORE ANY ERROR MESSAGES ABOVE THE FIRST DOTTED LINE\n"
    unittest.main()
