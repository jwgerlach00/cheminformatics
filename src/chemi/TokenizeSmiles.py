from typing import List
import re


class TokenizeSmiles:
    __smiles_regex = '(\[[^\[\]]{1,6}\])'
    __start_token = '<B>'
    __end_token = '<EOS>'
    __halogen_substitutes = {  # Single letter substitutes for halogens
        'Br': 'R',
        'Cl': 'L',
    }
    
    @staticmethod
    def start_token():
        return TokenizeSmiles.__start_token
    
    @staticmethod
    def end_token():
        return TokenizeSmiles.__end_token
    
    @staticmethod
    def sub_halogens(smiles:str) -> str:
        '''Regex to replace Br and Cl with single letters'''
        for h, s in TokenizeSmiles.__halogen_substitutes.items():
            smiles = smiles.replace(h, s)
        return smiles
    
    @staticmethod
    def unsub_halogens(subbed_smiles:str) -> str:
        for h, s in TokenizeSmiles.__halogen_substitutes.items():
            smiles = smiles.replace(s, h)
        return smiles
    
    @staticmethod
    def regex_split(smiles:str) -> List[str]:
        smiles = TokenizeSmiles.sub_halogens(smiles)
        frags = re.split(TokenizeSmiles.__smiles_regex, smiles)
        
        out = []
        for frag in frags:
            if frag.startswith('['):
                out.append(frag)
            else:
                out += list(frag)
                
        return out
    
    @staticmethod
    def tokenize(smiles:str, add_start_end_tokens:bool=True) -> List[str]:
        tokenized = TokenizeSmiles.regex_split(smiles)
        
        if add_start_end_tokens:
            tokenized.insert(0, TokenizeSmiles.__start_token)
            tokenized.append(TokenizeSmiles.__end_token)
            
        return tokenized
