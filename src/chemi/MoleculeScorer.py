import copy
from typing import Dict, Iterable, Tuple, List, Union
from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem.Descriptors import qed
from rdkit.Chem import QED
import numpy as np
from torch import nn
import pandas as pd

import naclo
from chemi import ChemFrag


class MoleculeScorer:
    '''Scores using drug-likeness, similarity to target, and model score'''
    
    def __init__(self, smiles:str, model_names:Iterable[str], desired_model_outcomes:Iterable[int],
                 model_importances:Iterable[float], models_map:Dict[str, nn.Module], smiles_query=None) -> None:
        self.__model_names = copy.deepcopy(model_names)
        self.__desired_model_outcomes = copy.deepcopy(desired_model_outcomes)
        self.__model_importances = copy.deepcopy(model_importances)
        
        self.__models_map = models_map
        
        self.__smiles_query = smiles_query
        
        self.__smiles = smiles
        self.__mol = Chem.MolFromSmiles(smiles)
        self.__fingerprint = naclo.mols_2_ecfp([self.__mol], radius=3, n_bits=1024, return_numpy=True)[0].astype(int)
        
        self.properties = self.generate_properties(self.__mol)
        
        self.score_similarity = self._instance_score_similarity
        self.score_on_models = self._instance_score_on_models
    
    @staticmethod
    def generate_properties(mol:Chem.rdchem.Mol) -> Dict[str, Union[float, int]]:
        _QED = qed(mol)
        qed_properties = QED.properties(mol)
        
        return {
            'QED': _QED,
            'MW': qed_properties.MW,
            'ALOGP': qed_properties.ALOGP,
            'HBA': qed_properties.HBA,
            'HBD': qed_properties.HBD
        }
    
    @staticmethod
    def score_similarity(smiles:str, smiles_query:str) -> Tuple[float, float]:
        query_mol = Chem.MolFromSmiles(smiles)
        reference_mol = Chem.MolFromSmiles(smiles_query)

        # Convert both the reference list and our molecule of interest into fingerprints
        query_fingerprint = naclo.mols_2_ecfp([query_mol], radius=3, n_bits=1024)[0]
        reference_fingerprint = naclo.mols_2_ecfp([reference_mol], radius=3, n_bits=1024)[0]
        total_similarity = []
        
        # Append full molecule similarity
        full_similarity = DataStructs.FingerprintSimilarity(query_fingerprint, reference_fingerprint)
        total_similarity.append(full_similarity)
        
        # Append decomposed fragments similiarity
        fragments = ChemFrag(smiles_query).fragments
        for frag_smiles in fragments:
            frag_mol = Chem.MolFromSmiles(frag_smiles)
            frag_fingerprint = naclo.mols_2_ecfp([frag_mol], radius=3, n_bits=1024, return_numpy=False)[0]
            total_similarity.append(DataStructs.FingerprintSimilarity(query_fingerprint, frag_fingerprint))

        # Return similarity against full refernce molecule and the max of all similarities
        return full_similarity, max(total_similarity)
    
    def _instance_score_similarity(self) -> Tuple[float, float]:
        return MoleculeScorer.score_similarity(self.__smiles, self.__smiles_query)
    
    @staticmethod
    def weighted_prediction_score(model_scores, qed, similarity=None, similarity_weight=0.7, qed_weight=0.3):
        if similarity is not None:
            score = sum(model_scores) + (similarity*similarity_weight) + (qed*qed_weight)
        else:
            score = sum(model_scores) + (qed*qed_weight)
        return score
    
    @staticmethod
    def score_on_models(fingerprint:np.ndarray, model_names:Iterable[str], desired_model_outcomes:Iterable[int],
                        model_importances:Iterable[float], models_map:Dict[str, nn.Module]) -> List[float]:
        predictions = []
        model_scores = []
        for model_name, desired_outcome, importance in zip(model_names, desired_model_outcomes, model_importances):
            prediction = models_map[model_name].predict_proba(fingerprint.reshape(1, -1))[:, 1]
            predictions.append(prediction[0])
            
            if desired_outcome == 0:
                prediction = 1 - prediction
                
            model_scores.append(prediction*importance)
            
        return model_scores
    
    def _instance_score_on_models(self) -> List[float]:
        return self.score_on_models(self.__fingerprint, self.__model_names, self.__desired_model_outcomes,
                                    self.__model_importances, self.__models_map)
            
    def score(self):
        model_scores = self.score_on_models()
        similarity = self.score_similarity() if self.__smiles_query else None
        return self.weighted_prediction_score(model_scores, self.properties['QED'], similarity=similarity)
