'''
Created on 22/09/2020

@author: mmp
'''
import unittest, os
from msa_masker import import_seqs, import_depth, process_data

class Test(unittest.TestCase):


	def setUp(self):
		pass


	def tearDown(self):
		pass


	def test_msa_masker(self):
		
		file_aligned = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tests/msa_masker/aligned_sequences.fasta")
		self.assertTrue(os.path.exists(file_aligned))
		file_deep = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tests/msa_masker/hit.depth.gz")
		self.assertTrue(os.path.exists(file_deep))
		
		msa_seqs = import_seqs(file_aligned)
		sample_map = import_depth(file_deep)
		
		reference = msa_seqs[0]
		reference_id = reference[0]
		samples = msa_seqs[1:]
		cutoff = 190
		mask_gaps = False
		process_data(reference[1], samples, sample_map, cutoff, mask_gaps)
		self.assertEqual("hit", samples[0][0])
		self.assertEqual("NNNAANTTTTNCCCCC-GGNNN", samples[0][1])

	def test_msa_masker_2(self):
		
		file_aligned = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tests/msa_masker/aligned_sequences.fasta")
		self.assertTrue(os.path.exists(file_aligned))
		file_deep = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tests/msa_masker/hit.depth")
		self.assertTrue(os.path.exists(file_deep))
		
		msa_seqs = import_seqs(file_aligned)
		sample_map = import_depth(file_deep)
		
		reference = msa_seqs[0]
		reference_id = reference[0]
		samples = msa_seqs[1:]
		cutoff = 80
		mask_gaps = False
		process_data(reference[1], samples, sample_map, cutoff, mask_gaps)
		self.assertEqual("hit", samples[0][0])
		self.assertEqual("AAAAATTTTTNCCCCC-GGGGA", samples[0][1])

	def test_msa_masker_3(self):
		
		file_aligned = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tests/msa_masker/aligned_sequences.fasta")
		self.assertTrue(os.path.exists(file_aligned))
		file_deep = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tests/msa_masker/hit.depth")
		self.assertTrue(os.path.exists(file_deep))
		
		msa_seqs = import_seqs(file_aligned)
		sample_map = import_depth(file_deep)
		
		reference = msa_seqs[0]
		reference_id = reference[0]
		samples = msa_seqs[1:]
		cutoff = 80
		mask_gaps = True
		process_data(reference[1], samples, sample_map, cutoff, mask_gaps)
		self.assertEqual("hit", samples[0][0])
		self.assertEqual("AAAAATTTTTNCCCCC-GGGGA", samples[0][1])

	def test_msa_masker_4(self):
		
		file_aligned = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tests/msa_masker/aligned_sequences_2.fasta")
		self.assertTrue(os.path.exists(file_aligned))
		file_deep = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tests/msa_masker/hit_2.depth")
		self.assertTrue(os.path.exists(file_deep))
		
		msa_seqs = import_seqs(file_aligned)
		sample_map = import_depth(file_deep)
		
		reference = msa_seqs[0]
		reference_id = reference[0]
		samples = msa_seqs[1:]
		cutoff = 80
		mask_gaps = True
		process_data(reference[1], samples, sample_map, cutoff, mask_gaps)
		self.assertEqual("hit", samples[0][0])
		self.assertEqual("--AAATTTTTNCCCCC-GGG--", samples[0][1])

	def test_msa_masker_5(self):
		
		file_aligned = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tests/msa_masker/aligned_sequences_2.fasta")
		self.assertTrue(os.path.exists(file_aligned))
		file_deep = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tests/msa_masker/hit_2.depth")
		self.assertTrue(os.path.exists(file_deep))
		
		msa_seqs = import_seqs(file_aligned)
		sample_map = import_depth(file_deep)
		
		reference = msa_seqs[0]
		reference_id = reference[0]
		samples = msa_seqs[1:]
		cutoff = 200
		mask_gaps = True
		process_data(reference[1], samples, sample_map, cutoff, mask_gaps)
		self.assertEqual("hit", samples[0][0])
		self.assertEqual("--NNANNNTTNCCCCC-GGG--", samples[0][1])
		
	def test_import_seqs(self):
		
		file_aligned = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tests/msa_masker/aligned_sequences.fasta")
		self.assertTrue(os.path.exists(file_aligned))
		vect_seqs = import_seqs(file_aligned)
		self.assertEqual(2, len(vect_seqs))
		self.assertEqual("reference", vect_seqs[0][0])
		self.assertEqual("AAAA-TTTTTTCCCCCCGGGGA", vect_seqs[0][1])
		self.assertEqual("hit", vect_seqs[1][0])
		self.assertEqual("AAAAATTTTTTCCCCC-GGGGA", vect_seqs[1][1])
		
		
	def test_import_depth(self):
		
		file_deep = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tests/msa_masker/hit.depth.gz")
		self.assertTrue(os.path.exists(file_deep))
		sample_map = import_depth(file_deep)
		self.assertTrue("hit" in sample_map)
		self.assertEqual(1, len(sample_map))
		self.assertEqual(21, len(sample_map["hit"]))
		self.assertEqual(3, len(sample_map["hit"][0]))
		self.assertEqual('hit', sample_map["hit"][0][0])
		self.assertEqual('1', sample_map["hit"][0][1])
		self.assertEqual('186', sample_map["hit"][0][2])
		self.assertEqual('hit', sample_map["hit"][20][0])
		self.assertEqual('21', sample_map["hit"][20][1])
		self.assertEqual('136', sample_map["hit"][20][2])
		
	def test_import_depth_2(self):
		
		file_deep = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tests/msa_masker/hit.depth")
		self.assertTrue(os.path.exists(file_deep))
		sample_map = import_depth(file_deep)
		self.assertTrue("hit" in sample_map)
		self.assertTrue("hit2" in sample_map)
		self.assertEqual(2, len(sample_map))
		self.assertEqual(21, len(sample_map["hit"]))
		self.assertEqual(5, len(sample_map["hit2"]))
		self.assertEqual(3, len(sample_map["hit"][0]))
		self.assertEqual('hit', sample_map["hit"][0][0])
		self.assertEqual('1', sample_map["hit"][0][1])
		self.assertEqual('186', sample_map["hit"][0][2])
		self.assertEqual('hit', sample_map["hit"][20][0])
		self.assertEqual('21', sample_map["hit"][20][1])
		self.assertEqual('136', sample_map["hit"][20][2])

	def test_import_depth_folder(self):
		
		path_deep = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tests/msa_masker/depth_files")
		self.assertTrue(os.path.exists(path_deep))
		sample_map = import_depth(path_deep, ["hit", "hit2"])
		self.assertTrue("hit" in sample_map)
		self.assertTrue("hit2" in sample_map)
		self.assertEqual(2, len(sample_map))
		self.assertEqual(21, len(sample_map["hit"]))
		self.assertEqual(5, len(sample_map["hit2"]))
		self.assertEqual(3, len(sample_map["hit"][0]))
		self.assertEqual('hit', sample_map["hit"][0][0])
		self.assertEqual('1', sample_map["hit"][0][1])
		self.assertEqual('186', sample_map["hit"][0][2])
		self.assertEqual('hit', sample_map["hit"][20][0])
		self.assertEqual('21', sample_map["hit"][20][1])
		self.assertEqual('136', sample_map["hit"][20][2])
if __name__ == "__main__":
	#import sys;sys.argv = ['', 'Test.test_msa_masker']
	unittest.main()