import copy
from typing import Dict, Iterable, List
import numpy as np

from chemi import TokenizeSmiles


class SmilesVocab:
    __start_token_index = 0
    __end_token_index = 1
    
    def __init__(self, vocab_smiles_list:Iterable[str]) -> None:
        self.__vocab_smiles_list = copy.deepcopy(vocab_smiles_list)
        self.__vocab_table = self._instance_create_vocab_table()
        
        self.encode = self._instance_encode
        self.decode = self._instance_decode
        
    def __str__(self) -> str:
        return self.__vocab_table.__str__()
        
    def __len__(self) -> int:
        return len(self.__vocab_table)
        
    @property
    def vocab_smiles(self) -> List[str]:
        return self.__vocab_smiles_list
    
    @property
    def vocab_table(self) -> dict:
        return self.__vocab_table

    @staticmethod
    def unique_tokens(smiles_list:Iterable[str]) -> List[str]:
        all_tokens = []
        for smiles in smiles_list:
            all_tokens += TokenizeSmiles.tokenize(smiles, add_start_end_tokens=False)
        return sorted((set(all_tokens)))
        
    @staticmethod
    def create_vocab_table(smiles_list:Iterable[str]) -> Dict[str, int]:
        '''Assigns integer values to each token in the vocabulary'''
        tokens = SmilesVocab.unique_tokens(smiles_list)
        
        # Add start and end tokens
        tokens.insert(SmilesVocab.__start_token_index, TokenizeSmiles.start_token)
        tokens.insert(SmilesVocab.__end_token_index, TokenizeSmiles.end_token)
        
        vocab_table = {}
        for token in tokens:
            vocab_table[token] = len(vocab_table)
            
        return vocab_table
        
    def _instance_create_vocab_table(self) -> dict:
        tokens = self.unique_tokens(self.__vocab_smiles_list)
        return self.create_vocab_table(tokens)
    
    @staticmethod
    def encode(smiles:str, vocab_table:Dict[str, int]) -> np.ndarray:
        '''Converts a list of tokens to an array of integers as defined in the vocab table'''
        tokenized = TokenizeSmiles.tokenize(smiles)
        return np.array([vocab_table[token] for token in tokenized])
    
    def _instance_encode(self, smiles:str) -> np.ndarray:
        return SmilesVocab.encode(smiles, self.__vocab_table)
    
    @staticmethod
    def decode(smiles_arr:np.ndarray, vocab_table:Dict[str, int]) -> str:
        reversed_vocab = {v: k for k, v in vocab_table.items()}
        
        error_flag = False
        decoded = []
        for n in smiles_arr:
            if n == SmilesVocab.__start_token_index:
                pass
            else:
                if n == SmilesVocab.__end_token_index:
                    break
                try:
                    token = reversed_vocab[n]
                    decoded.append(token)
                except:
                    error_flag = True
                    break
            
        return None if error_flag else TokenizeSmiles.unsub_halogens(''.join(decoded))

    def _instance_decode(self, smiles_arr:np.ndarray) -> str:
        SmilesVocab.decode(smiles_arr, self.__vocab_table)
