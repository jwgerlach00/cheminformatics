import unittest
import numpy as np

from chemi import SmilesVocab


class TestSmilesVocab(unittest.TestCase):
    
    def test_vocab_smiles(self):
        self.assertEqual(
            SmilesVocab(['CCO', 'CC[Br]O', 'CC[Br]O[Cl]']).vocab_smiles,
            ['CCO', 'CC[Br]O', 'CC[Br]O[Cl]']
        )
        
    
    def test_vocab_table(self):
        self.assertEqual(
            SmilesVocab(['CCO', 'CC[Br]O', 'CC[Br]O[Cl]']).vocab_table,
            {'<B>': 0, '<EOS>': 1, 'C': 2, 'O': 3, '[L]': 4, '[R]': 5}
        )

    def test_unique_tokens(self):
        self.assertEqual(
            SmilesVocab.unique_tokens(['CCO', 'CC[Br]O', 'CC[Br]O[Cl]']),
            ['C', 'O', '[L]', '[R]']  # sorted
        )

    def test_create_vocab_table(self):
        self.assertEqual(
            SmilesVocab.create_vocab_table(['CCO', 'CC[Br]O', 'CC[Br]O[Cl]']),
            {'<B>': 0, '<EOS>': 1, 'C': 2, 'O': 3, '[L]': 4, '[R]': 5}
        )

    def test_encode(self):
        self.assertTrue(
            np.array_equal(
                SmilesVocab.encode('CCO', {'<B>': 0, '<EOS>': 1, 'C': 2, 'O': 3, '[L]': 4, '[R]': 5}),
                np.array([0, 2, 2, 3, 1])
            )
        )

    def test_decode(self):
        self.assertEqual(
            SmilesVocab.decode(np.array([0, 2, 2, 3, 1]), {'<B>': 0, '<EOS>': 1, 'C': 2, 'O': 3, '[L]': 4, '[R]': 5}),
            'CCO'
        )


if __name__ == '__main__':
    unittest.main()
