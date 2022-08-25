from typing import Iterable
import pandas as pd
import glob
import os
import rdkit.Chem as Chem
import numpy as np
from scipy.stats import gaussian_kde


class AutoThreshold:

    @staticmethod
    def numerical_distance(values:Iterable[float]) -> float:
        # get track the index of the greatest numerical distance
        greatest_distance = 0
        threshold = 0
        # for each value in the list, starting from 1 to the list len-1
        for j in range(len(values) - 2):
            # get the distance between the previous and next point
            dist_1 = abs(values[j + 1] - values[j])
            dist_2 = abs(values[j + 1] - values[j + 2])

            # if the distance is greater than the current threshold, replace the current value
            if dist_1 + dist_2 > threshold:
                greatest_distance = values[j + 1]
                threshold = dist_1 + dist_2

        return greatest_distance

def gaussian_kde_threshold(values:Iterable[float]) -> float:
    '''Returns minima of gaussian kde'''
    kde = gaussian_kde(values)

    samples = np.linspace(min(values), max(values), len(values))
    probs = kde.evaluate(samples)
    return values[probs.argmin()]

def auto_threshold(df:pd.DataFrame, value='converted_units', min='auto', max='auto', method='numerical'):

    # >/< can look to weird results. Remove those for the results of threshold finding sounds like the right thing.
    df = df[df['Standard Relation'] == "'='"]

    # if min and max not supplied, auto-find the threshold
    if min == 'auto' and max == 'auto':
        # first, get the values for the 20% and 40%
        quant_40 = np.percentile(np.array(df[value].to_list()), 85)
        quant_20 = np.percentile(np.array(df[value].to_list()), 50)
    # else, use the provided range for look for a suitable threshold
    else:
        quant_40 = min
        quant_20 = max

    print('Auto-range: {0} \t{1}'.format(quant_20, quant_40))
    # next, subset the DF by those values
    thresh_df = df[(df[value] <= quant_40) & (df[value] >= quant_20)]
    thresh_df.sort_values(by=[value], inplace=True)

    # get the threshold to split on
    if method == "numerical":
        threshold = numerical_distance(thresh_df[value].to_list())
    else:
        threshold = gaussian_kde_threshold(thresh_df[value].to_list())
    # return the max value
    return threshold