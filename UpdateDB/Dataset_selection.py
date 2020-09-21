"""
    Created by: Eric March Vila (eric.march@upf.edu)
    On: 21/09/2020, 13:45 PM
"""

import numpy as np
import pandas as pd

from imblearn.combine import SMOTEENN, SMOTETomek
from imblearn.over_sampling import RandomOverSampler
from imblearn.under_sampling import RandomUnderSampler
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import PandasTools
from sklearn.model_selection import train_test_split
from typing import Union, Optional, Tuple

class Selection(object):
    """
        This object aims to authomatise the dataset selection from CII database to build the models.
        It also includes the imbalance correction, which is applied if the user needs to

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

    ### Imbalance correction functions

    def random_subsampler(self, train_: pd.DataFrame, test_: pd.DataFrame) -> pd.DataFrame:
        """
            Random Undersampling
            https://imbalanced-learn.readthedocs.io/en/stable/under_sampling.html

            :param train_:
            :param test_: 

            :return resampled_train_set, resampled_test_set:
        """

        rus = RandomUnderSampler(random_state=0)
        
        y_train = train_['activity']
        x_train = train_.drop(columns='activity')

        y_test = test_['activity']
        x_test = test_.drop(columns='activity')
        
        regular_cols = x_train.columns
        
        x_train_resampled, y_train_resampled = rus.fit_resample(x_train, y_train)
        x_test_resampled, y_test_resampled = rus.fit_resample(x_test, y_test)
        
        resampled_train_set = pd.DataFrame(data=x_train_resampled, columns=regular_cols)
        resampled_train_set.loc[:,'activity'] = y_train_resampled
        
        resampled_test_set = pd.DataFrame(data=x_test_resampled, columns=regular_cols)
        resampled_test_set.loc[:,'activity'] = y_test_resampled
        
        return resampled_train_set, resampled_test_set
    
    def random_oversampler(self, train_: pd.DataFrame, test_: pd.DataFrame) -> pd.DataFrame:
        """
            Random Oversampling
            https://imbalanced-learn.readthedocs.io/en/stable/over_sampling.html

            :param train_:
            :param test_:

            :return resampled_train_set, resampled_test_set:
        """

        ros = RandomOverSampler(random_state=0)
        
        y_train = train_['activity']
        x_train = train_.drop(columns='activity')

        y_test = test_['activity']
        x_test = test_.drop(columns='activity')
        
        regular_cols = x_train.columns
        
        x_train_resampled, y_train_resampled = ros.fit_resample(x_train, y_train)
        x_test_resampled, y_test_resampled = ros.fit_resample(x_test, y_test)
        
        resampled_train_set = pd.DataFrame(data=x_train_resampled, columns=regular_cols)
        resampled_train_set.loc[:,'activity'] = y_train_resampled
        
        resampled_test_set = pd.DataFrame(data=x_test_resampled, columns=regular_cols)
        resampled_test_set.loc[:,'activity'] = y_test_resampled
        
        return resampled_train_set, resampled_test_set
    
    def process_datasets(self, dataset: pd.DataFrame) -> Tuple[np.ndarray, pd.DataFrame]:
        """
            this function is used to reshape the dataset in order to prepare it for SMOTEEN.

            :param dataset:

            :return x_reshaped, y_:
        """

        y_ = dataset['activity']
        x_ = dataset.drop(columns='activity')
        x_reshaped = self.reshape_datasets(x_)
        
        return x_reshaped, y_

    def reshape_datasets(self, dataset: pd.DataFrame) -> np.ndarray:
        """
            This functions applies numpy's reshape to the dataset

            :param dataset:

            :return reshaped_dataset:
        """

        reshaped_dataset = np.reshape(dataset.index, (len(dataset.index), 1))
        
        return reshaped_dataset

    def get_proper_datasets(self, dataset: pd.DataFrame, indexes: np.ndarray) -> pd.DataFrame:
        """
            Removes duplicate registers after SMOTEEN algorithm is applied.
            Some of the registers that were being created had the same id in both datasets, thus
            giving errors.
            We've decided to remove those registers in order to keep datasets cleaned of noise.

            :param dataset:

            :return proper_dataset:
        """

        proper_dataset = pd.DataFrame(columns=dataset.columns,index=indexes).reset_index()
        index_to_drop = set()
        
        for i, row in proper_dataset.iterrows():
            idx = row['index']
            if idx in dataset.index:
                proper_dataset.loc[i] = dataset.loc[idx]
            else:
                index_to_drop.add(idx)
        
        proper_dataset.drop(proper_dataset.loc[proper_dataset['index'].isin(index_to_drop)].index, inplace=True)
        
        return proper_dataset
                
    def smoteen_resample_sets(self, train_: pd.DataFrame, test_: pd.DataFrame) -> pd.DataFrame:
        """
            This function applies SMOTEEN once train and test set are split.

            :param train_:
            :param test_:

            :return resampled_train_set, resampled_test_set:
        """

        smote_enn = SMOTEENN(random_state=42)
        
        x_train_index_reshape, y_train = self.process_datasets(train_)
        x_test_index_reshape, y_test = self.process_datasets(test_)

        x_train_resampled, y_train_resampled = smote_enn.fit_resample(x_train_index_reshape, y_train)
        x_test_resampled, y_test_resampled = smote_enn.fit_resample(x_test_index_reshape, y_test)
        
        x_train_res_flatten = x_train_resampled.flatten()
        x_test_res_flatten = x_test_resampled.flatten()
        
        resampled_train_set = self.get_proper_datasets(train_, x_train_res_flatten)
        resampled_test_set = self.get_proper_datasets(test_, x_test_res_flatten)
        
        return resampled_train_set, resampled_test_set

    def resample_main_set(self, main_set:pd.DataFrame) -> pd.DataFrame:
        """
            This function applies SMOTEEN to the whole dataset before splitting into train and test sets.
            It is deprecated since we found that is better to apply it after splitting, but we keep the function.

            :param main_set:

            :return resampled_set:
        """

        smote_enn = SMOTEENN(random_state=42)
        
        x, y = self.process_datasets(main_set)
        
        x_resampled, y_resampled = smote_enn.fit_resample(x,y)
        
        x_flatten = x_resampled.flatten()
        
        resampled_set = self.get_proper_datasets(main_set, x_flatten)
        
        return resampled_set
