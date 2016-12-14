import unittest
import parser
import xml.etree.ElementTree as etree
import pickle

class TestStringMethods(unittest.TestCase):

    def test_checkGeneDict(self):        
        self.assertEqual(parser.checkGeneDict('BRCA1'), "LRG_292", "checkGeneDict() did not return correct LRG for BRCA1")
        self.assertEqual(parser.checkGeneDict('xxxxxxx'), None, "checkGeneDict() did not return Null value when invalid gene symbol provided")
        
    def test_checkValidLRG(self):
        self.assertEqual(parser.checkValidLRG('LRG_1'), 'LRG_1', "valid LRG did not pass checkValidLRG() validation")
        self.assertEqual(parser.checkValidLRG('xxxxxxx'), None, "invalid LRG was not caught by checkValidLRG() validation")
    
    def test_readLRG(self):
        tree = parser.readLRG('LRG_1')
        self.assertTrue(tree, "No element tree object returned from readLRG() when requested")
        self.assertEqual(tree.find('./fixed_annotation/hgnc_id').text, '2197', "Element tree object does not represent correct LRG file")
        tree = parser.readLRG('xxxxxxx')
        self.assertFalse(tree, "readLRG() did not return Null value when unable to download LRG file")

    def test_getGeneLevData(self):
        tree = etree.parse('./Test_Files/LRG_1.xml')
        geneData = parser.getGeneLevData(tree)
        self.assertEqual(geneData['hgncID'], '2197')
        self.assertEqual(geneData['strand'], '-1')
        
    def test_getGenomicSeq(self):
        tree = etree.parse('./Test_Files/LRG_292.xml')
        mainAssemSeq, convertPosDict = parser.getGenomicSeq(tree)
        handle = open('./Test_Files/BRCA_38_Seq.txt', 'r')
        BRCA_38_Seq = handle.read()
        handle.close()
        self.assertEqual(mainAssemSeq, BRCA_38_Seq, "Sequence output from getGenomicSeq() is incorrect")
        self.assertEqual(convertPosDict[13168], 13168, "Error in convertPosDict")
        self.assertEqual(convertPosDict[13171], 13169, "Error in convertPosDict")
        self.assertEqual(convertPosDict[15459], 15457, "Error in convertPosDict")
        self.assertEqual(convertPosDict[15460], 15459, "Error in convertPosDict")
        self.assertEqual(convertPosDict[15460], 15459, "Error in convertPosDict")

        
        
if __name__ == '__main__':
    unittest.main()