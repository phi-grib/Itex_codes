"""
    Created by: Eric March Vila (eric.march@upf.edu)
    On: 21/09/2020, 14:33 PM
"""

import numpy as np
import pandas as pd

from imblearn.combine import SMOTEENN, SMOTETomek
from imblearn.over_sampling import RandomOverSampler
from imblearn.under_sampling import RandomUnderSampler
from typing import Union, Optional, Tuple

### Random subsampler

def random_subsampler(train_: pd.DataFrame, test_: pd.DataFrame) -> pd.DataFrame:
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

### Random oversampler

def random_oversampler(train_: pd.DataFrame, test_: pd.DataFrame) -> pd.DataFrame:
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

#### SMOTEEN part

def process_datasets(dataset: pd.DataFrame) -> Tuple[np.ndarray, pd.DataFrame]:
    """
        this function is used to reshape the dataset in order to prepare it for SMOTEEN.

        :param dataset:

        :return x_reshaped, y_:
    """

    y_ = dataset['activity']
    x_ = dataset.drop(columns='activity')
    x_reshaped = reshape_datasets(x_)
    
    return x_reshaped, y_

def reshape_datasets(dataset: pd.DataFrame) -> np.ndarray:
    """
        This functions applies numpy's reshape to the dataset

        :param dataset:

        :return reshaped_dataset:
    """

    reshaped_dataset = np.reshape(dataset.index, (len(dataset.index), 1))
    
    return reshaped_dataset

def get_proper_datasets(dataset: pd.DataFrame, indexes: np.ndarray) -> pd.DataFrame:
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
            
def smoteen_resample_sets(train_: pd.DataFrame, test_: pd.DataFrame) -> pd.DataFrame:
    """
        This function applies SMOTEEN once train and test set are split.

        :param train_:
        :param test_:

        :return resampled_train_set, resampled_test_set:
    """

    smote_enn = SMOTEENN(random_state=42)
    
    x_train_index_reshape, y_train = process_datasets(train_)
    x_test_index_reshape, y_test = process_datasets(test_)

    x_train_resampled, y_train_resampled = smote_enn.fit_resample(x_train_index_reshape, y_train)
    x_test_resampled, y_test_resampled = smote_enn.fit_resample(x_test_index_reshape, y_test)
    
    x_train_res_flatten = x_train_resampled.flatten()
    x_test_res_flatten = x_test_resampled.flatten()
    
    resampled_train_set = get_proper_datasets(train_, x_train_res_flatten)
    resampled_test_set = get_proper_datasets(test_, x_test_res_flatten)
    
    return resampled_train_set, resampled_test_set

def resample_main_set(main_set:pd.DataFrame) -> pd.DataFrame:
    """
        This function applies SMOTEEN to the whole dataset before splitting into train and test sets.
        It is deprecated since we found that is better to apply it after splitting, but we keep the function.

        :param main_set:

        :return resampled_set:
    """

    smote_enn = SMOTEENN(random_state=42)
    
    x, y = process_datasets(main_set)
    
    x_resampled, y_resampled = smote_enn.fit_resample(x,y)
    
    x_flatten = x_resampled.flatten()
    
    resampled_set = get_proper_datasets(main_set, x_flatten)
    
    return resampled_set
