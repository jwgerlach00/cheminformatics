import unittest
from rdkit import Chem
from rdkit.Chem import Recap

from chemi import ChemFrag


class TestChemFrag(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls) -> None:
        cls.smiles = 'CC(=O)Oc1ccccc1C(=O)O'
        cls.expected_murcko_scaffold = 'CC(C)CC1CCCCC1C(C)C'
        cls.expected_recap = [
            'CC=O',
            'CCC',
            'O=C(O)c1ccccc1O',
            'CC(C)C1CCCCC1C',
            'O',
            'C',
            'CC(=O)Oc1ccccc1C=O',
            'CCC1CCCCC1CC(C)C',
            'O=C(O)c1ccccc1',
            'CC(C)C1CCCCC1'
        ]
        cls.expected_fragments = [
            cls.smiles,
            cls.expected_murcko_scaffold
        ]
        cls.expected_fragments += cls.expected_recap
        return super().setUpClass()
    
    def test_murcko_scaffold(self):
        self.assertEqual(
            ChemFrag.murcko_scaffold(self.smiles),
            self.expected_murcko_scaffold
        )
    
    def test_recap_decomp(self):
        self.assertEqual(
            ChemFrag.recap_decomp(Chem.MolFromSmiles(self.smiles)),
            self.expected_recap
        )
        
    def test_find_fragments(self):
        self.assertEqual(
            ChemFrag.find_fragments(self.smiles),
            self.expected_fragments
        )
    
    def test_fragments(self):
        CF = ChemFrag(self.smiles)
        self.assertEqual(
            CF.fragments,
            self.expected_fragments
        )
        

if __name__ == '__main__':
    unittest.main()
