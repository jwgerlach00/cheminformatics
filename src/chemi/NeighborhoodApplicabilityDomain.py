import rdkit.Chem as Chem
from sklearn.ensemble import RandomForestClassifier
from rdkit.Chem import DataStructs
from typing import Iterable, Union
import numpy as np

import naclo

Fingerprint = Union[np.array, DataStructs.cDataStructs.ExplicitBitVect]


class NeighborhoodApplicabilityDomain:
    ecfp_radius = 3
    ecfp_n_bits = 1024
    
    def __init__(self):
        self.surrogate_random_forest = RandomForestClassifier(class_weight="balanced")
        
        # Calculated during training
        self.trained = False
        self.domain_fps = None
        self.biases = []
        self.stdevs = []
    
    @staticmethod
    def mols_2_fps(mols:Iterable[Chem.rdchem.Mol], fp_type:str='ecfp'):
        if fp_type == 'ecfp':
            fps = naclo.mols_2_ecfp(mols, radius=NeighborhoodApplicabilityDomain.ecfp_radius,
                                    n_bits=NeighborhoodApplicabilityDomain.ecfp_n_bits, return_numpy=True)
        elif fp_type == 'maccs':
            fps = naclo.mols_2_maccs(mols)
        else:
            raise(ValueError, 'fp_type must be "ecfp" or "maccs"')
        
        return fps
    
    @staticmethod
    def train_surrogate_rf(fps:Iterable[Fingerprint], bin_labels:Iterable[int]):
        return NeighborhoodApplicabilityDomain.surrogate_random_forest.fit(fps, bin_labels)
    
    @staticmethod
    def calculate_bias(predictions:Iterable[float], true_values:Iterable[int]):
        return [pred if true == 1 else 1 - pred for pred, true in zip(predictions, true_values)]
    
    @staticmethod
    def calculate_std(fps:Iterable[Fingerprint], trained_rf:RandomForestClassifier):
        predictions = [estimator.predict_proba(fps)[:, 1] for estimator in trained_rf.estimators_]
        return np.std(predictions)
    
    def train_surrogate_model(self, domain_mols:Iterable[Chem.rdchem.Mol], domain_bin_labels:Iterable[int],
                              fp_type:str='ecfp'):
        # Fit model
        self.domain_fps = self.mols_2_fps(domain_mols, fp_type=fp_type)
        self.surrogate_random_forest.fit(self.domain_fps, domain_bin_labels)
        
        # Get biases
        self.surrogate_random_forest.predict_proba(self.domain_fps)
        self.biases = self.calculate_bias(self.domain_fps, domain_bin_labels)
        
        # Get standard deviations
        self.stdevs = self.calculate_std(self.domain_fps, self.surrogate_random_forest)
        
        # Set trained flag
        self.trained = True
        
    def calculate_applicability(self, query_mols:Iterable[Chem.rdchem.Mol], fp_type:str='ecfp'):
        if not self.trained:
            raise(ValueError, 'Must train model before calculating applicability')
        else:
            fps = self.mols_2_fps(query_mols, fp_type=fp_type)
            
            applicability_scores = []
            for fp in fps:
                tanimoto_scores = DataStructs.BulkTanimotoSimilarity(fp, self.surrogate_fps)
                max_index = tanimoto_scores.index(max(tanimoto_scores))
                
                if tanimoto_scores[max_index] == 1:
                    applicability_scores.append(1)
                else:
                    W = self.biases[max_index]*(1 - self.stdevs[max_index])
                    applicability_scores.append(tanimoto_scores[max_index]*W)
                    
            return applicability_scores
