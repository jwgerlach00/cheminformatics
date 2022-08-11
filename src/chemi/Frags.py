from rdkit import Chem
from rdkit.Chem import Recap
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem


class FragmentSmiles:
    
    def __init__(self, smiles:str) -> None:
        self.__smiles = smiles
        self.__mol = Chem.MolFromSmiles(smiles)
        
    @staticmethod
    def murcko_scaffold(smiles:str) -> str:
        return Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(Chem.MolFromSmiles(smiles)))
    
    def _instance_murcko_scaffold(self) -> str:
        return self.murcko_scaffold(self.__smiles)
    
    @staticmethod
    def recap_decomp(mol:Chem.rdchem.Mol) -> list:
        rec_list = list(Recap.RecapDecompose(mol).children.keys())
        
        for rec in rec_list:
            mol = Chem.MolFromSmiles(rec)
            
            # Cannonicalize the smiles, removing dummy atom
            pattern = Chem.MolFromSmiles('*')
            mol = AllChem.ReplaceSubstructs(mol, pattern, Chem.MolFromSmiles('[H]'), True)[0]
            mol = Chem.RemoveHs(mol)
            
            scaff = Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(x2))
            x2 = Chem.MolToSmiles(x2)
    
    def _instance_recap_decomp(self) -> list:
        return self.recap_decomp(self.__mol)
    
    
    
    
    


def chem_frag(smiles:str):
    out = [smiles]
    
    # Append the murcko scaffold
    out.append(Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(Chem.MolFromSmiles(smiles))))
    

    # Decompose Mol using recap rules
    mol = Chem.MolFromSmiles(smiles)
    rec_mol = Recap.RecapDecompose(mol)
    rec_list = list(rec_mol.children.keys())

    # cononicalize each smiles fragment, append to list
    for frag in rec_list:
        # turn the fragment into a molecule
        frag = Chem.MolFromSmiles(frag)

        # canonicalize the smiles, removing the dummy atom
        x = Chem.MolFromSmiles('*')
        x1 = AllChem.ReplaceSubstructs(
            frag, x, Chem.MolFromSmiles('[H]'), True)[0]
        x2 = Chem.MolToSmiles(x1)
        x2 = Chem.RemoveHs(x1)

        # get the generic murcko scaffold
        scaff = Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(x2))
        x2 = Chem.MolToSmiles(x2)

        # if the scaffold is not the same as the fragment, append it
        if scaff != x2:
            fragments.append(scaff)

        fragments.append(x2)

    return fragments