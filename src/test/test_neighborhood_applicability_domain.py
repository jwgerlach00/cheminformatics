import unittest
from rdkit import Chem, DataStructs
from sklearn.ensemble import RandomForestClassifier
import numpy as np

import naclo

from chemi import NeighborhoodApplicabilityDomain


class TestNeighborhoodApplicabilityDomain(unittest.TestCase):

    def test_mols_2_fps(self):
        test_mol = Chem.MolFromSmiles('CCO')
        
        # ECFP
        ecfps = NeighborhoodApplicabilityDomain.mols_2_fps([test_mol], fp_type='ecfp')
        self.assertIsInstance(ecfps[0], DataStructs.cDataStructs.ExplicitBitVect)
        self.assertEqual(len(ecfps[0]), NeighborhoodApplicabilityDomain.ecfp_n_bits)
        
        # MACCS
        maccs = NeighborhoodApplicabilityDomain.mols_2_fps([test_mol], fp_type='maccs')
        self.assertIsInstance(maccs[0], DataStructs.cDataStructs.ExplicitBitVect)
        
        # Invalid fp_type
        with self.assertRaises(ValueError):
            NeighborhoodApplicabilityDomain.mols_2_fps([test_mol], fp_type='')

    def test_calculate_bias(self):
        predictions = [0.5, 0.3, 0.1, 0.9]
        true = [1, 0, 0, 1]
        expected = [0.5, 0.7, 0.9, 0.9]
        self.assertEqual(
            NeighborhoodApplicabilityDomain.calculate_bias(predictions, true),
            expected
        )

    def test_calculate_std(self):
        test_mols = naclo.smiles_2_mols([
            'CCO',
            'CCC',
            'CCCCCCCCCN',
            'O=C1C=CC=CC1',
        ])
        test_labels = [1, 1, 0, 0]
        
        rf_model = RandomForestClassifier(class_weight="balanced")
        
        # Test all fingerprints
        for fp_type in ['ecfp', 'maccs']:
            fps = NeighborhoodApplicabilityDomain.mols_2_fps(test_mols, fp_type=fp_type)
            rf_model.fit(fps, test_labels)
            
            stdevs = NeighborhoodApplicabilityDomain.calculate_std(fps, rf_model)
            self.assertEqual(len(stdevs), len(test_mols))  # Only test length, outputs vary (RANDOM)
        

    def test_train_surrogate_model(self):
        test_mols = naclo.smiles_2_mols([
            'CCO',
            'CCC',
            'CCCCCCCCCN',
            'O=C1C=CC=CC1',
        ])
        test_labels = [1, 1, 0, 0]
        
        # Test all fingerprints
        for fp_type in ['ecfp', 'maccs']:
            ad = NeighborhoodApplicabilityDomain()
            
            self.assertFalse(ad.trained)
            ad.train_surrogate_model(domain_mols=test_mols, domain_bin_labels=test_labels, fp_type=fp_type)
            self.assertTrue(ad.trained)
            
            self.assertEqual(len(ad.biases), len(test_mols))
            self.assertEqual(len(ad.stdevs), len(test_mols))

    def test_calculate_applicability(self):
        test_mols = naclo.smiles_2_mols([
            'CCO',
            'CCC',
            'CCCCCCCCCN',
            'O=C1C=CC=CC1',
        ])
        test_labels = [1, 1, 0, 0]
        
        test_query_mols = naclo.smiles_2_mols([
            'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
            'O=C1C=CC=CC1',
            'CC(CC1=CC=CC=C1)NC'
        ])
        
        ad = NeighborhoodApplicabilityDomain()
        
        with self.assertRaises(RuntimeError):  # Not trained
            ad.calculate_applicability(query_mols=test_query_mols)

        for fp_type in ['ecfp', 'maccs']:
            ad.train_surrogate_model(domain_mols=test_mols, domain_bin_labels=test_labels, fp_type=fp_type)
            scores = ad.calculate_applicability(query_mols=test_query_mols)
            self.assertEqual(len(scores), len(test_query_mols))
            self.assertEqual(scores[1], 1)  # Test query_mol 2 is in domain
            
            # Test that staticmethod obtains the same results
            query_fps = NeighborhoodApplicabilityDomain.mols_2_fps(test_query_mols, fp_type=fp_type)
            static_scores = NeighborhoodApplicabilityDomain.calculate_applicability(query_fps, ad.domain_fps,
                                                                                    biases=ad.biases, stdevs=ad.stdevs)
            self.assertTrue(
                np.array_equal(
                    scores,
                    static_scores
                )
            )
