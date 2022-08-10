import unittest

from chemi import TokenizeSmiles


class TestTokenizeSmiles(unittest.TestCase):

    def test_sub_halogens(self):
        self.assertEqual(TokenizeSmiles.sub_halogens('BrCl'), 'RL')

    def test_unsub_halogens(self):
        self.assertEqual(TokenizeSmiles.unsub_halogens('RL'), 'BrCl')

    def test_regex_split(self):
        self.assertEqual(
            TokenizeSmiles.regex_split('CCO'),
            ['C', 'C', 'O']
        )
        self.assertEqual(
            TokenizeSmiles.regex_split('CC[Br]O'),
            ['C', 'C', '[R]', 'O']
        )
        self.assertEqual(
            TokenizeSmiles.regex_split('CC[Br]O[Cl]'),
            ['C', 'C', '[R]', 'O', '[L]']
        )

    def test_tokenize(self):
        self.assertEqual(
            TokenizeSmiles.tokenize('CCO'),
            ['<B>', 'C', 'C', 'O', '<EOS>']
        )
        self.assertEqual(
            TokenizeSmiles.tokenize('CC[Br]O'),
            ['<B>', 'C', 'C', '[R]', 'O', '<EOS>']
        )
        self.assertEqual(
            TokenizeSmiles.tokenize('CC[Br]O[Cl]'),
            ['<B>', 'C', 'C', '[R]', 'O', '[L]', '<EOS>']
        )


if __name__ == '__main__':
    unittest.main()
