import functions_to_process as fp
import numpy as np
import pandas as pd
import re

class IDB(object):
    def __init__(self):
        self._main_dataframe = None
        self._column = None

    def define_dataframe(self, dataframe: pd.DataFrame, column: str):
        self._main_dataframe = dataframe
        self._column = column

    def process_dataframe(self):
        self.set_, self.cas_, self.index_, self.ec_, self.additional_ = fp.clean_df(self._main_dataframe, self._column)

    def get_col_name(self):
        return self._column

    def get_main_set(self):
        return self.set_

    def get_cas_set(self):
        return self.cas_

    def get_index_set(self):
        return self.index_
    
    def get_ec_set(self):
        return self.ec_
    
    def get_additional_info_set(self):
        return self.additional_