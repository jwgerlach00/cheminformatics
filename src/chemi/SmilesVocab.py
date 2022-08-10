import copy
import re
from typing import Iterable, List


class SmilesVocab:
    __smiles_regex = '(\[[^\[\]]{1,6}\])'
    __start_token = '<B>'
    __end_token = '<EOS>'
    
    def __init__(self, smiles_list:Iterable[str]) -> None:
        self.__smiles_list = copy.deepcopy(smiles_list)
        self.__vocab_table = self._instance_create_vocab_table()
        
    def __str__(self) -> str:
        return self.__vocab_table.__str__()
        
    def __len__(self) -> int:
        return len(self.__vocab_table)
        
    @property
    def smiles_list(self) -> List[str]:
        return self.__smiles_list
    
    @property
    def vocab_table(self) -> dict:
        return self.__vocab_table
    
    @staticmethod
    def replace_halogens(smiles:str) -> str:
        '''Regex to replace Br and Cl with single letters'''
        halogens = ['Br', 'Cl']
        substitues = ['R', 'L']
        for h, s in zip(halogens, substitues):
            smiles = smiles.replace(h, s)
        return smiles
    
    @staticmethod
    def regex_split(smiles:str) -> List[str]:
        smiles = SmilesVocab.replace_halogens(smiles)
        frags = re.split(SmilesVocab.__smiles_regex, smiles)
        
        out = []
        for frag in frags:
            if frag.startswith('['):
                out.append(frag)
            else:
                out += list(frag)
                
        return out
    
    @staticmethod
    def tokenize_smiles(smiles:str) -> List[str]:
        tokenized = SmilesVocab.regex_split(smiles)
        tokenized.insert(0, SmilesVocab.__start_token)
        tokenized.append(SmilesVocab.__end_token)
        return tokenized
    
    @staticmethod
    def unique_tokens(smiles_list:Iterable[str]) -> List[str]:
        '''Returns a list of unique tokens in a list of SMILES strings'''
        return list(set([SmilesVocab.regex_split(smiles) for smiles in smiles_list]))
    
    @staticmethod
    def create_vocab_table(tokens:Iterable[str]) -> dict:
        '''Assigns integer values to each token in the vocabulary'''
        vocab_table = {}
        for token in tokens:
            if token not in vocab_table:
                vocab_table[token] = len(vocab_table)
        return vocab_table
    
    def _instance_create_vocab_table(self) -> dict:
        tokens = self.unique_tokens(self.__smiles_list)
        return self.create_vocab_table(tokens)
