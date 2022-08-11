from rdkit import Chem
from rdkit.Chem import Recap
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem
from typing import List


class ChemFrag:
    
    def __init__(self, smiles:str) -> None:
        self.__smiles = smiles
        self.__mol = Chem.MolFromSmiles(smiles)
        
        self.murcko_scaffold = self._instance_murcko_scaffold
        self.recap_decomp = self._instance_recap_decomp
        
        self.__fragments = self.find_fragments(self.__smiles)
        
    @property
    def fragments(self) -> List[str]:
        return self.__fragments
        
    @staticmethod
    def murcko_scaffold(smiles:str) -> str:
        return Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(Chem.MolFromSmiles(smiles)))
    
    def _instance_murcko_scaffold(self) -> str:
        return ChemFrag.murcko_scaffold(self.__smiles)
    
    @staticmethod
    def recap_decomp(mol:Chem.rdchem.Mol) -> List[str]:
        rec_list = list(Recap.RecapDecompose(mol).children.keys())
        
        out = []
        for rec in rec_list:
            mol = Chem.MolFromSmiles(rec)
            
            # Cannonicalize the smiles, removing dummy atom
            pattern = Chem.MolFromSmiles('*')
            mol = AllChem.ReplaceSubstructs(mol, pattern, Chem.MolFromSmiles('[H]'), True)[0]
            mol = Chem.RemoveHs(mol)
            
            scaff = Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(mol))
            smiles = Chem.MolToSmiles(mol)
            
            out.append(smiles)
            if scaff != smiles:
                out.append(scaff)
        
        return out
    
    def _instance_recap_decomp(self) -> List[str]:
        return ChemFrag.recap_decomp(self.__mol)
    
    @staticmethod
    def find_fragments(smiles:str) -> List[str]:
        out = [smiles]
        out.append(ChemFrag.murcko_scaffold(smiles))
        out += ChemFrag.recap_decomp(Chem.MolFromSmiles(smiles))
        return out
