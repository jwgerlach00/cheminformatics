import rdkit.Chem as Chem
from sklearn.ensemble import RandomForestClassifier
from rdkit.Chem import DataStructs
from typing import Iterable, Union, Optional, List
import numpy as np

import naclo

Fingerprint = Union[np.ndarray, DataStructs.cDataStructs.ExplicitBitVect]


class NeighborhoodApplicabilityDomain:
    ecfp_radius = 3
    ecfp_n_bits = 1024
    
    def __init__(self) -> None:
        self.surrogate_random_forest = RandomForestClassifier(class_weight="balanced")
        
        # Calculated during training
        self.domain_fps = None
        self.biases = []
        self.stdevs = []
        
        # Flags
        self.trained = False
        self.train_fp_type = None
        self.calculate_applicability = self.__calculate_applicability
    
    @staticmethod
    def mols_2_fps(mols:Iterable[Chem.rdchem.Mol], fp_type:str='ecfp') -> List[Fingerprint]:
        '''Convert a list of RDKit molecules to a list of fingerprints'''
        if fp_type == 'ecfp':
            fps = naclo.mols_2_ecfp(mols, radius=NeighborhoodApplicabilityDomain.ecfp_radius,
                                    n_bits=NeighborhoodApplicabilityDomain.ecfp_n_bits)
        elif fp_type == 'maccs':
            fps = naclo.mols_2_maccs(mols)
        else:
            raise ValueError('fp_type must be "ecfp" or "maccs"')
        
        return fps
    
    @staticmethod
    def calculate_bias(predictions:Iterable[float], true_values:Iterable[int]) -> List[float]:
        '''Calculate the bias for each molecule in the domain'''
        return [pred if true == 1 else 1 - pred for pred, true in zip(predictions, true_values)]
    
    @staticmethod
    def calculate_std(fps:Iterable[Fingerprint], trained_rf:RandomForestClassifier) -> List[float]:
        '''Calculate the standard deviation of the predictions for each molecule in the domain'''
        predictions = [estimator.predict_proba(fps)[:, 1] for estimator in trained_rf.estimators_]
        return np.std(predictions, axis=0)
    
    def train_surrogate_model(self, domain_mols:Iterable[Chem.rdchem.Mol], domain_bin_labels:Iterable[int],
                              fp_type:str='ecfp', random_seed:Optional[int]=42) -> RandomForestClassifier:
        '''Train a surrogate random forest model to predict applicability of molecules to a domain'''
        if random_seed:
            np.random.seed(random_seed)
        
        # Fit model
        self.domain_fps = self.mols_2_fps(domain_mols, fp_type=fp_type)
        self.surrogate_random_forest.fit(self.domain_fps, domain_bin_labels)
        
        # Get biases
        predictions = self.surrogate_random_forest.predict_proba(self.domain_fps)
        
        if predictions.shape != (len(domain_mols), 2):
            raise ValueError('Too many classes')
        
        self.biases = self.calculate_bias(predictions[:, 1], domain_bin_labels)  # Only use class 1
        
        # Get standard deviations
        self.stdevs = self.calculate_std(self.domain_fps, self.surrogate_random_forest)
        
        # Set trained flag
        self.trained = True
        self.train_fp_type = fp_type
        
        return self.surrogate_random_forest
    
    @staticmethod
    def calculate_applicability(query_fps, domain_fps, biases, stdevs) -> List[float]:
        '''Calculate the applicability of a list of molecules to a domain'''
        applicability_scores = []
        for fp in query_fps:
            tanimoto_scores = DataStructs.BulkTanimotoSimilarity(fp, domain_fps)
            max_index = tanimoto_scores.index(max(tanimoto_scores))
            
            if tanimoto_scores[max_index] == 1:
                applicability_scores.append(1)
            else:
                W = biases[max_index]*(1 - stdevs[max_index])
                applicability_scores.append(tanimoto_scores[max_index]*W)
                
        return applicability_scores
        
    def __calculate_applicability(self, query_mols:Iterable[Chem.rdchem.Mol]) -> List[float]:
        '''Calculate the applicability of a list of molecules to a domain'''
        if not self.trained:
            raise RuntimeError('Must train model before calculating applicability')
        else:
            fps = self.mols_2_fps(query_mols, fp_type=self.train_fp_type)
            return NeighborhoodApplicabilityDomain.calculate_applicability(query_fps=fps, domain_fps=self.domain_fps,
                                                                           biases=self.biases, stdevs=self.stdevs)
