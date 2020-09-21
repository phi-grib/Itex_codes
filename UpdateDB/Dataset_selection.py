"""
    Created by: Eric March Vila (eric.march@upf.edu)
    On: 21/09/2020, 13:45 PM
"""

import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import PandasTools
from sklearn.model_selection import train_test_split
from typing import Union, Optional, Tuple

from .imbalance_correction import imb_cor_functions as imbcor

class Selection(object):
    """
        This object aims to authomatise the dataset selection from CII database to build the models.
        It also includes the imbalance correction, which is applied if the user needs to.
    """

    def __init__(self, dataframe: pd.DataFrame, train_prop: float, test_prop: float, imbalance_algorithm: str = None):
        """
            Initializes class

            :param dataframe:
            :param train_prop:
            :param test_prop:
            :param imbalance_algorithm: default is None. Otherwise, it can be Oversampling, Subsampling or SMOTEEN
        """

        self.main_data = dataframe
        self.train_prop = train_prop
        self.test_prop = test_prop
        self.imbalance_algorithm = imbalance_algorithm
    
    ### Selection main function

    def split_main_dataset(self) -> pd.DataFrame:
        """
            This is the main function that returns the training set and the test set
            after applying the different proportions and the imbalance correction to the main set

            :return train_set, test_set:
        """

        train_set, test_set = self.get_sets(self.main_data, self.train_prop, self.test_prop)
        
        if self.imbalance_algorithm.lower() == 'oversampling':
            train_set, test_set = imbcor.random_oversampler(train_set,test_set)
        elif self.imbalance_algorithm.lower() == 'subsampling':
            train_set, test_set = imbcor.random_subsampler(train_set, test_set)
        elif self.imbalance_algorithm.lower == 'smoteen':
            train_set, test_set = imbcor.smoteen_resample_sets(train_set, test_set)

        return train_set, test_set

    ### Set selection

    def get_sets(self, df: pd.DataFrame, train_prop: float, test_prop: float) -> pd.DataFrame:
        """
            This function performs the train and test set selection as explained in the link below:
            https://stats.stackexchange.com/questions/394056/splitting-into-train-and-test-sets-keeping-class-proportions

            :param df:
            :param train_prop:
            :param test_prop:

            :return train_set, test_set:
        """

        y = df['activity']
        x = df.drop(columns='activity')
        
        X_train, X_test, y_train, y_test = train_test_split(x, y, train_size=train_prop, test_size=test_prop, stratify = y, random_state=42)
        
        train_set = pd.concat([X_train, y_train], axis=1).reindex(X_train.index)
        test_set = pd.concat([X_test, y_test], axis=1).reindex(X_test.index)
        
        return train_set, test_set