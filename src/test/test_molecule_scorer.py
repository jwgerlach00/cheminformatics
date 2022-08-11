import unittest
from rdkit import Chem

from chemi import MoleculeScorer


class TestMoleculeScorer(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls) -> None:
        cls.smiles = 'CC(=O)Oc1ccccc1C(=O)O'
        cls.mol = Chem.MolFromSmiles(cls.smiles)
        return super().setUpClass()
    
    def test_generate_properties(self):
        properties = MoleculeScorer.generate_properties(self.mol)
        self.assertIsInstance(properties, dict)
        self.assertEqual(
            list(properties.keys()),
            ['QED', 'MW', 'ALOGP', 'HBA', 'HBD']
        )
        
        for key in properties.keys():
            self.assertTrue(
                type(properties[key]) in [float, int]
            )
            self.assertTrue(
                properties[key] > 0
            )
            
    def test_score_similarity(self):
        same_smiles_score = MoleculeScorer.score_similarity(self.smiles, self.smiles)
        self.assertEqual(
            same_smiles_score,
            (1, 1)
        )
        
        unrelated_score = MoleculeScorer.score_similarity(self.smiles, 'S')
        self.assertEqual(
            unrelated_score,
            (0, 0)
        )
        
        somewhat_related_score = MoleculeScorer.score_similarity(self.smiles, 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C')
        self.assertTrue(
            somewhat_related_score[0] > 0 and somewhat_related_score[0] < 1
        )
        self.assertTrue(
            somewhat_related_score[1] > 0 and somewhat_related_score[1] < 1
        )
        
    def test_weighted_prediction_score(self):
        score = MoleculeScorer.weighted_prediction_score(
            model_scores=[0.5, 0.5],
            qed=0.5,
            similarity=0.5,
            similarity_weight=0.7,
            qed_weight=0.3
        )
        self.assertEqual(score, 1.5)
        
        # No similarity
        score = MoleculeScorer.weighted_prediction_score(
            model_scores=[0.5, 0.5],
            qed=0.5,
            similarity=None,
            similarity_weight=0.7,
            qed_weight=0.3
        )
        self.assertEqual(score, 1.15)
        
    def test_score_on_models(self):
        self.fail()
        
    def test_score_on_model(self):
        self.fail()


if __name__ == '__main__':
    unittest.main()
