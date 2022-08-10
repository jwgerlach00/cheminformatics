from typing import List
import re


class TokenizeSmiles:
    start_token = '<B>'
    end_token = '<EOS>'
    
    __smiles_regex = '(\[[^\[\]]{1,6}\])'
    __halogen_substitutes = {  # Single letter substitutes for halogens
        'Br': 'R',
        'Cl': 'L',
    }
    
    @staticmethod
    def sub_halogens(smiles:str) -> str:
        '''Regex to replace Br and Cl with single letters'''
        for h, s in TokenizeSmiles.__halogen_substitutes.items():
            smiles = smiles.replace(h, s)
        return smiles
    
    @staticmethod
    def unsub_halogens(subbed_smiles:str) -> str:
        for h, s in TokenizeSmiles.__halogen_substitutes.items():
            subbed_smiles = subbed_smiles.replace(s, h)
        return subbed_smiles
    
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
            tokenized.insert(0, TokenizeSmiles.start_token)
            tokenized.append(TokenizeSmiles.end_token)
            
        return tokenized
