import numpy as np
import pandas as pd
import processing_columns_functions as pcf
import re


from typing import Optional, Union, Tuple, Callable

def extend_class_name(dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    This function extends the class name to empty cells when there's a non-empty cell at its right,
    meaning there's a preferred name in our dataframe but the class name has not been read well
    """
    cl_name = False
    for index, row in dataframe.iterrows():
        if type(row[0]) == str and type(row[1]) == str:
            cl_name = row[0]
        elif type(row[0]) == float:
            if cl_name:
                dataframe.iloc[index, 0] = cl_name

def get_set_or_list(row: list, func: Callable[[str],str], list_ = None) -> Union[list,set]:
    if list_:
        new_list = []
    else:
        _set = set()
    for element in row:
        new_element = func(element)
        if list_:
            if isinstance(new_element, list):
                new_list.extend(new_element)
            else:
                new_list.append(new_element)
        else:
            if isinstance(new_element, list):
                _set.update(set(new_element))
            else:
                _set.add(new_element)
    
    if list_:
        return new_list
    else:
        return _set

def change_nan_to_none(value: Union[float,str,int]) -> Optional[str]:
    """
    Function to get rid of nan values from the dataframe before adding them into the database.
    None will be read as NULL in Postgres.
    """
    if type(value) is float:
        value_to_return = None
        return value_to_return
    if value == 'nan':
        value_to_return = None
        return value_to_return
    else:
        return value

def check_cas_ec_index(cas: str, ec: str, index:str) -> Optional[str]:
    if index:
        if cas:
            if cas[0] in index[0]:
                cas = None
            else:
                pass
        if ec:
            if ec[0] in index[0]:
                ec = None
            else:
                pass
        return cas,ec,index
    else:
        return cas,ec,index

#
# First functions process CAS, EC and Index numbers
#

def get_cas(element: str) -> Optional[str]:
    """
    Uses the pattern to extract a CAS registry number from the input. Returns the CAS if found, otherwise, returns None
    """
    cas_pattern = r"[0-9]{1,6}-\d{2}-\d\.?"
    casrn = re.findall(cas_pattern, element)
    if casrn:
        return casrn
    else:
        return None

def get_ec(element: str) -> Optional[str]:
    """
    Uses the pattern to extract an EC registry number from the input. Returns the EC if found, otherwise, returns None
    """
    ec_pattern = r"\.?[0-9]{1,3}-\d{3}-\d\.?"
    ec = re.findall(ec_pattern, element)
    if ec:
        return ec
    else:
        return None
    
def get_index_number(element: str) -> Optional[str]:
    """
    Uses the pattern to extract an Index registry number from the input. Returns the index if found, otherwise, returns None
    """
    index_pattern = r"[0-9]{1,3}-[0-9]{1,3}-\d{2}-[0-9A-Z]"
    index_ = re.findall(index_pattern, element)
    if index_:
        return index_
    else:
        return None
    
def get_cas_numbers(row: list) -> set:
    """
    Gets the row from the DataFrame as a list and for each element it checks if there's a CAS there.
    Stores each CAS in a set, to avoid duplications. Returns the set.
    """
    cas_set = set()
    for element in row:
        casrn = get_cas(element)
        if casrn:
            for number in casrn:
                cas_set.add(number.strip())
                
    return cas_set

def get_index_numbers(row: list) -> set:
    """
    Gets the row from the DataFrame as a list and for each element it checks if there's an EC there.
    Stores each EC in a set, to avoid duplications. Returns the set.
    """
    index_set = set()
    for element in row:
        index_ = get_index_number(element)
        if index_:
            for number in index_:
                index_set.add(number.strip())
    
    return index_set

def get_ec_numbers(row: list) -> set:
    """
    Gets the row from the DataFrame as a list and for each element it checks if there's an Index number there.
    Stores each Index number in a set, to avoid duplications. Returns the set.
    """
    ec_set = set()
    for element in row:
        ec = get_ec(element)
        if ec:
            for number in ec:
                ec_set.add(number.strip())
    
    return ec_set


#
# START Common block to treat strings
#

def get_additional_info(element: str, col_name: str) -> Optional[str]:
    """
    Uses col_name to decide which functions to apply to select the additional information, specific to each column.
    """
    if col_name == 'eu_subreg_reach_subeval_annexiii':
        add_info = pcf.get_suspected(element)
    elif col_name == 'eu_subreg_reach_subeval_general_eval':
        add_info = pcf.get_list_cleaned(element)
    elif col_name == 'eu_subreg_reach_svhc_cl_preliminary':
        add_info = pcf.get_add_info_svhc(element)
    elif col_name == 'eu_subreg_reach_svhc_auth_annexiv':
        add_info = pcf.get_date_string(element)
    elif col_name == 'eu_subreg_reach_svhc_cl_definitive' or col_name == 'eu_subreg_reach_svhc_auth_prelim_inc_annexiv'\
        or col_name == 'eu_subreg_reach_svhc_auth_approved' or col_name == 'eu_subreg_reach_svhc_auth_application' \
        or col_name == 'eu_subreg_reach_svhc_auth_rejected' or col_name == 'eu_subreg_clp_harmonized_cl_spec_conc_limits'\
        or col_name == 'eu_subreg_pop_prohib_annexi_part_a' or col_name == 'eu_subreg_pop_prohib_annexi_part_b'\
        or col_name == 'eu_subreg_pop_sublist_restrict' or col_name == 'eu_subreg_pop_sublist_release_reduct'\
        or col_name == 'eu_subreg_pop_sublist_waste_management' or col_name == 'eu_subreg_pop_waste_management'\
        or col_name == 'eu_subreg_wfd_sublist_water_policy' or col_name == 'eu_subreg_wfd_env_quality' or col_name == 'eu_subreg_wfd_env_quality'\
        or col_name == 'eu_subreg_pbt_vpvb_assess_reach' or col_name == 'eu_subreg_pbt_vpvb_reach_annexxiii_reg_info'\
        or col_name == 'eu_subreg_pbt_vpvb_reach_annexxiii_cand_list' or col_name == 'eu_subreg_pbt_vpvb_other'\
        or col_name == 'eu_subreg_endocrine_dis_biocidal_cmr_subs' or col_name == 'eu_subreg_endocrine_dis_eu_com'\
        or col_name == 'eu_label_clp_harmonized_cl_pictogram' or col_name == 'eu_label_clp_harmonized_cl_signal_word'\
        or col_name == 'eu_label_clp_harmonized_cl_hazard' or col_name == 'eu_label_clp_harmonized_cl_sup_hazard' \
        or col_name == 'sw_subclass_kemi_env_hazard_score' or col_name == 'sw_subclass_kemi_env_monit_data' \
        or col_name == 'gr_subclass_com_ord_restrict_subs_articles':
        add_info = None
    elif col_name == 'eu_subreg_reach_restrictedsub_annexvii_preliminary' or col_name == 'eu_subreg_reach_restrictedsub_annexvii_restrict_submix'\
        or col_name == 'eu_subreg_reach_restrictedsub_annexvii_restrict_articles':
        add_info = pcf.clean_add_info_axvii_prelim(element)
    elif col_name == 'eu_subreg_clp_preliminary_cli_reg_dossier' or col_name == 'eu_subreg_clp_preliminary_cli_notification':
        add_info = pcf.clean_add_info_registration_dossier_clp(element)
    elif col_name == 'eu_subreg_clp_preliminary_harmonised_cl':
        add_info = pcf.clean_add_info_harm_cl(element)[0]
    elif col_name == 'eu_subreg_clp_harmonized_cl_index_num' or col_name == 'eu_subreg_clp_harmonized_cl_hazard_class'\
        or col_name == 'eu_subreg_clp_harmonized_cl_hazard':
        add_info = pcf.clean_add_info_index(element)
    elif col_name == 'eu_subreg_clp_harmonized_cl_notes':
        add_info = pcf.clean_add_info_notes(element)
    elif col_name == 'sw_subclass_chemprod_ord_restrict_subs_submix' or col_name == 'sw_subclass_chemprod_ord_restrict_subs_articles':
        add_info = pcf.clean_add_info_sw_restrict_submix(element)
    elif col_name == 'gr_subclass_chemicals_act_restrict_submix':
        add_info = pcf.clean_gr_chemicals_act(element)

    return add_info

def check_exceptions_by_column(row: list, col_name: str, list_ = None) -> Union[Tuple[list, set], Optional[Tuple[list, set]]]:
    """
    Uses col_name to decide which functions to apply in order to process the elements in the input row.
    """
    new_row = row
    
    if col_name == 'eu_subreg_reach_subeval_annexiii':
        new_row = pcf.clean_asbestos_row(new_row)
        row_set = pcf.clean_elements_set(new_row, list_)
        additional_info_set = pcf.get_suspect_set(row_set)

    elif col_name == 'eu_subreg_reach_subeval_general_eval':
        raw_row_set = pcf.clean_elements_set(new_row)
        raw_additional_info_set = pcf.get_list_set(raw_row_set, raw = True)
        row_set = raw_row_set - raw_additional_info_set
        additional_info_set = pcf.get_list_set(raw_row_set)
    
    elif col_name == 'eu_subreg_reach_svhc_cl_preliminary':
        row_set = get_set_or_list(new_row,pcf.clean_svhc_cl_preliminary,list_)
        additional_info_set = pcf.get_svhc_add_info_set(row_set)
    
    elif col_name == 'eu_subreg_reach_svhc_cl_definitive':
        row_set = get_set_or_list(new_row, pcf.clean_elements_svhc_def, list_)
        additional_info_set = set()
    
    elif col_name == 'eu_subreg_reach_svhc_auth_prelim_inc_annexiv':
        row_set = get_set_or_list(new_row, pcf.clean_svhc_prelim_axiv, list_)
        additional_info_set = set()
    
    elif col_name == 'eu_subreg_reach_svhc_auth_annexiv':
        row_set = get_set_or_list(new_row, pcf.clean_svhc_subauth_axiv, list_)
        additional_info_set = pcf.get_add_info_svhc_subauth(row_set)

    elif col_name == 'eu_subreg_reach_svhc_auth_approved' or col_name == 'eu_subreg_reach_svhc_auth_application'\
        or col_name == 'eu_subreg_reach_svhc_auth_rejected':
        row_set = get_set_or_list(new_row, pcf.clean_svhc_approved_auth_axiv, list_) 
        additional_info_set = set()
    
    elif col_name == 'eu_subreg_reach_restrictedsub_annexvii_preliminary' or col_name == 'eu_subreg_reach_restrictedsub_annexvii_restrict_submix'\
        or col_name == 'eu_subreg_reach_restrictedsub_annexvii_restrict_articles':
        raw_row_set = get_set_or_list(new_row,pcf.clean_axvii_preliminary, list_)
        additional_info_set = pcf.get_add_info_axvii_prelim_set(raw_row_set)
        if list_:
            row_set = raw_row_set
        else:
            raw_additional_info_set = pcf.get_add_info_axvii_prelim_set(raw_row_set, raw = True)
            row_set = raw_row_set - raw_additional_info_set

    elif col_name == 'eu_subreg_clp_preliminary_cli_reg_dossier' or col_name == 'eu_subreg_clp_preliminary_cli_notification':
        raw_row_set = get_set_or_list(new_row,pcf.clean_registration_dossier_clp,list_)
        additional_info_set = pcf.get_add_info_registration_dossier_clp(raw_row_set)
        if list_:
            row_set = raw_row_set
        else:
            raw_additional_info_set = pcf.get_add_info_registration_dossier_clp(raw_row_set, raw = True)
            row_set = raw_row_set - raw_additional_info_set

    elif col_name == 'eu_subreg_clp_preliminary_harmonised_cl':
        raw_row_set = get_set_or_list(new_row,pcf.clean_harmonized_cl_clp,list_)
        additional_info_set, endpoints_set = pcf.get_add_info_harm_cl(raw_row_set)
        if list_:
            endpoints_list = list(endpoints_set)
            raw_row_set.extend(endpoints_list)
            row_set = raw_row_set
        else:
            raw_additional_info_set = pcf.get_add_info_harm_cl(raw_row_set, raw=True)[0]
            row_set = (raw_row_set - raw_additional_info_set) | endpoints_set

    elif col_name == 'eu_subreg_clp_harmonized_cl_index_num':
        row_set = get_set_or_list(new_row,pcf.clean_harmonized_index_number_clp,list_)
        additional_info_set = pcf.get_add_info_index(row_set)

    elif col_name == 'eu_subreg_clp_harmonized_cl_hazard_class':
        row_set = get_set_or_list(new_row,pcf.clean_hazard_class_clp,list_)
        additional_info_set = pcf.get_add_info_index(row_set)

    elif col_name == 'eu_subreg_clp_harmonized_cl_hazard':
        row_set = get_set_or_list(new_row,pcf.clean_hazard_clp,list_)
        additional_info_set = pcf.get_add_info_index(row_set)

    elif col_name == 'eu_subreg_clp_harmonized_cl_spec_conc_limits':
        row_set = get_set_or_list(new_row, pcf.clean_conc_limits_clp,list_)
        additional_info_set = set()

    elif col_name == 'eu_subreg_clp_harmonized_cl_notes':
        row_set = get_set_or_list(new_row, pcf.clean_notes_clp,list_)
        additional_info_set = pcf.get_add_info_notes(row_set)

    elif col_name == 'eu_subreg_pop_prohib_annexi_part_a' or col_name == 'eu_subreg_pop_prohib_annexi_part_b'\
        or col_name == 'eu_subreg_pop_sublist_restrict' or col_name == 'eu_subreg_pop_sublist_release_reduct'\
        or col_name == 'eu_subreg_pop_sublist_waste_management' or col_name == 'eu_subreg_pop_waste_management'\
        or col_name == 'eu_subreg_pbt_vpvb_assess_reach' or col_name == 'sw_subclass_kemi_env_hazard_score'\
        or col_name == 'sw_subclass_kemi_env_monit_data':
        row_set = get_set_or_list(new_row,pcf.clean_pops_regulations, list_)
        additional_info_set = set()

    elif col_name == 'eu_subreg_wfd_sublist_water_policy' or col_name == 'eu_subreg_wfd_env_quality':
        row_set = get_set_or_list(new_row,pcf.clean_WFD_regulations,list_)
        additional_info_set = set()

    elif col_name == 'eu_subreg_pbt_vpvb_reach_annexxiii_reg_info' or col_name == 'eu_subreg_pbt_vpvb_reach_annexxiii_cand_list'\
        or col_name == 'eu_subreg_pbt_vpvb_other':
        row_set = get_set_or_list(new_row,pcf.clean_pbt_vpvb_reg_info,list_)
        additional_info_set = set()

    elif col_name == 'eu_subreg_endocrine_dis_biocidal_cmr_subs' or col_name == 'eu_subreg_endocrine_dis_eu_com':
        row_set = get_set_or_list(new_row,pcf.clean_endoc_dis_info,list_)
        additional_info_set = set()

    elif col_name == 'eu_label_clp_harmonized_cl_pictogram' or col_name == 'eu_label_clp_harmonized_cl_signal_word':
        row_set = get_set_or_list(new_row,pcf.clean_lab_clp_pictogram,list_)
        additional_info_set = set()

    elif col_name == 'eu_label_clp_harmonized_cl_hazard':
        row_set = get_set_or_list(new_row,pcf.clean_lab_clp_hazard_statement,list_)
        additional_info_set = set()

    elif col_name == 'eu_label_clp_harmonized_cl_sup_hazard':
        row_set = get_set_or_list(new_row,pcf.clean_lab_clp_hazard_supplementary, list_)
        additional_info_set = set()

    elif col_name == 'sw_subclass_chemprod_ord_restrict_subs_submix' or col_name == 'sw_subclass_chemprod_ord_restrict_subs_articles':
        row_set = get_set_or_list(new_row,pcf.clean_sw_classif_restric_submix, list_)
        additional_info_set = pcf.get_add_info_sw_restrict_submix(row_set)

    elif col_name == 'gr_subclass_com_ord_restrict_subs_articles':
        row_set = get_set_or_list(new_row,pcf.clean_gr_gco_restric_articles,list_)
        additional_info_set = set()
    
    elif col_name == 'gr_subclass_chemicals_act_restrict_submix':
        row_set = get_set_or_list(new_row,pcf.clean_gr_gco_restric_articles,list_)
        additional_info_set = pcf.get_add_info_gr_chemicals_act(row_set)

    return new_row, row_set, additional_info_set

def string_comparison(s1: str, s2: str) -> bool:
    s1_trimed = s1.replace(" ", "").lower()
    s2_trimed = s2.replace(" ", "").lower()
    
    return s1_trimed == s2_trimed

def count_mayus(string_list: list) -> Tuple[str,int]:
    for str_ in string_list:
        yield str_, sum(1 for letter in str_ if letter.isupper())

def check_set_case_elements(main_set: set) -> set:
    lower_set = set(map(lambda x: x.lower(), main_set))
    copy_set = main_set.copy()

    for element in lower_set:
        i = 0
        main_el_collect = []
        for main_el in main_set:
            if element == main_el.lower():
                i += 1
                main_el_collect.append(main_el)
        if i > 1:
            selected_string = max(count_mayus(main_el_collect), key=lambda x: x[1])[0]
            main_el_collect.remove(selected_string)
            for element in main_el_collect:
                copy_set.remove(element)
    
    return copy_set

def clean_row_df(row: list, col: str) -> set:
    """
    General function to obtain a big set with the names, cas set, ec set, index set and additional information set. It uses the functions above.
    Col allows to decide which specific column has to be treated.
    """
    new_row, row_set, additional_info = check_exceptions_by_column(row, col)
    cas_set = get_cas_numbers(new_row)
    index_set = get_index_numbers(new_row)
    ec_set = get_ec_numbers(new_row)
    no_cas_index_row_set = row_set - cas_set - index_set - ec_set - additional_info
    checked_no_cas_index_row_set = check_set_case_elements(no_cas_index_row_set)
    
    return checked_no_cas_index_row_set, cas_set, index_set, ec_set, additional_info

def clean_df(df: pd.DataFrame, col: str) -> set:
    """
    Main function. Takes the original dataframe, splits it and passes it to a list. Then calls clean_row_df to finally get the sets.
    """
    _clean = df.loc[:,col].apply(lambda x: str(x).strip().split('\n'))
    _clean_list = [sublist for li in _clean for sublist in li]
    _set, cas_set, index_set, ec_set, additional_info = clean_row_df(_clean_list, col)
    
    return _set, cas_set, index_set, ec_set, additional_info

def get_proper_value(variable: str, col_name: Optional[str] = None) -> Optional[list]:
    """
    This function splits the strings from the dataframe columns. Is it used when we want to finally add the elements to the big dictionary,
    before obtaining the big dataframe and insert it into the database.
    """
    new_var = None
    try:
        new_var = variable.strip().split('\n')
        if col_name:
            new_var = check_exceptions_by_column(new_var, col_name, list_ = True)[1]
        else:
            new_var = pcf.clean_elements_set(new_var, list_=True)
    except AttributeError:
        pass
    
    return new_var

#
# END block to treat strings
#
#
# START block to add data to the dictionary
#

def df_to_concat_df(df_: pd.DataFrame, set_: set, id_key: str, value_name: str) -> pd.DataFrame:
    """
    Concatenates updated set with new values to old ones. Returns updated dataframe
    """
    
    diff_set = set_.difference(set(df_[value_name].values))
    max_id = df_[id_key].values.max()
    dif_dict = {id_key:[], value_name:[]}

    for i,elements in enumerate(diff_set):
        id_ = max_id + i + 1
        dif_dict[id_key].append(id_)
        dif_dict[value_name].append(elements)

    dif_df = pd.DataFrame(data=dif_dict)
    new_df = pd.concat([df_, dif_df], sort=False).reset_index().drop(['index'], axis=1)

    return new_df


def dict_to_df(set_: set, id_key: str, value_name: str) -> pd.DataFrame:
    """
    Function to obtain a dataframe from differents sets. We get a specific id for each element in the set.
    """
    dict_ = {id_key:[], value_name:[]}

    for i_, name in enumerate(set_):
        if name:
            i = i_ + 1
            dict_[id_key].append(i)
            dict_[value_name].append(name)

    df_ = pd.DataFrame(data=dict_)
    
    return df_

def add_to_dict(dictionary: dict, id_: int, subs_id: int, reg_country_id: int, reg_type_id: int, gen_reg_id: int, spec_reg_id: int,
                subspec_reg_id: int, special_cases_id: Optional[int], additional_classification_id: float, additional_classificaton_name: Optional[str],
                chem_id: Optional[str], chem_type_id: Union[int,float], regulaton_id: int, regulation_explained: str):
    """
    This function adds to the dictionary all the variables passed as input.
    """
    dictionary['id'].append(id_)
    dictionary['subs_id'].append(subs_id)
    dictionary['reg_country_id'].append(reg_country_id)
    dictionary['reg_type_id'].append(reg_type_id)
    dictionary['gen_reg_id'].append(gen_reg_id)
    dictionary['spec_reg_id'].append(spec_reg_id)
    dictionary['subspec_reg_id'].append(subspec_reg_id)
    dictionary['special_cases_id'].append(special_cases_id)
    dictionary['additional_classification_id'].append(additional_classification_id)
    dictionary['additional_classificaton_name'].append(additional_classificaton_name)
    dictionary['chem_id'].append(chem_id)
    dictionary['chem_type_id'].append(chem_type_id)
    dictionary['regulation_id'].append(regulaton_id)
    dictionary['regulation_explained'].append(regulation_explained)

def add_element_to_dictionary(regulation: Optional[list], reg_country_id: int, reg_type_id: int, gen_reg_id: int, spec_reg_id: int,
                                subspec_reg_id: int, special_cases_id: Optional[int], id_: int, subs_id: int, name_dataframe: pd.DataFrame,
                                classification_dataframe: pd.DataFrame, dictionary: dict, col_name: str) -> int:
    """
    This function compares the cleaned names on the set previosly generated to the raw dataframe. Depending on which names detects (CAS, EC...)
    it will pass different arguments to the add_to_dict function. Returns the updated id_ for the following row of the dataframe
    """
    for i, names in name_dataframe.iterrows():
        id_name = names[0]
        name = names[1]
        flag_cas = 0
        flag_index = 0
        flag_ec = 0
        flag_classif = 0
        if regulation:
            for element in regulation:
                casrn = get_cas(element)
                ec = get_ec(element)
                index_num = get_index_number(element)
                casrn, ec, index_num = check_cas_ec_index(casrn, ec, index_num)
                classif = get_additional_info(element,col_name)
                if classif:
                    flag_classif = 1
                    classif_id = classification_dataframe.loc[classification_dataframe['additional_info'] == classif, 'id'].iat[0]
                    classif_name = classif
                if casrn:
                    flag_cas = 1
                    cas_to_add = casrn[0]
                if ec:
                    flag_ec = 1
                    ec_to_add = ec[0]
                if index_num:
                    flag_index = 1
                    index_to_add = index_num[0]
                if name.lower() == element.lower():
                    id_ += 1
                    if flag_cas == 1 and flag_index == 1 and flag_ec == 1 and flag_classif == 1:
                        add_to_dict(dictionary, id_, subs_id, reg_country_id, reg_type_id, gen_reg_id, spec_reg_id, subspec_reg_id, special_cases_id, 
                                    classif_id, classif_name, cas_to_add, 1, id_name, name)
                        id_ += 1
                        add_to_dict(dictionary, id_, subs_id, reg_country_id, reg_type_id, gen_reg_id, spec_reg_id, subspec_reg_id, special_cases_id, 
                                    classif_id, classif_name, ec_to_add, 2, id_name, name)
                        id_ += 1
                        add_to_dict(dictionary, id_, subs_id, reg_country_id, reg_type_id, gen_reg_id, spec_reg_id, subspec_reg_id, special_cases_id, 
                                    classif_id, classif_name, index_to_add, 3, id_name, name)
                    elif flag_cas == 1 and flag_ec == 1 and flag_classif == 1:
                        add_to_dict(dictionary, id_, subs_id, reg_country_id, reg_type_id, gen_reg_id, spec_reg_id, subspec_reg_id, special_cases_id, 
                                    classif_id, classif_name, cas_to_add, 1, id_name, name)
                        id_ += 1
                        add_to_dict(dictionary, id_, subs_id, reg_country_id, reg_type_id, gen_reg_id, spec_reg_id, subspec_reg_id, special_cases_id, 
                                    classif_id, classif_name, ec_to_add, 2, id_name, name)
                    elif flag_index == 1 and flag_ec == 1 and flag_classif == 1:
                        add_to_dict(dictionary, id_, subs_id, reg_country_id, reg_type_id, gen_reg_id, spec_reg_id, subspec_reg_id, special_cases_id, 
                                    classif_id, classif_name, ec_to_add, 2, id_name, name)
                        id_ += 1
                        add_to_dict(dictionary, id_, subs_id, reg_country_id, reg_type_id, gen_reg_id, spec_reg_id, subspec_reg_id, special_cases_id, 
                                    classif_id, classif_name, index_to_add, 3, id_name, name)
                    elif flag_cas == 1 and flag_index == 1 and flag_ec == 1:
                        add_to_dict(dictionary, id_, subs_id, reg_country_id, reg_type_id, gen_reg_id, spec_reg_id, subspec_reg_id, special_cases_id, 
                                    classif, classif, cas_to_add, 1, id_name, name)
                        id_ += 1
                        add_to_dict(dictionary, id_, subs_id, reg_country_id, reg_type_id, gen_reg_id, spec_reg_id, subspec_reg_id, special_cases_id, 
                                    classif, classif, ec_to_add, 2, id_name, name)
                        id_ += 1
                        add_to_dict(dictionary, id_, subs_id, reg_country_id, reg_type_id, gen_reg_id, spec_reg_id, subspec_reg_id, special_cases_id, 
                                    classif, classif, index_to_add, 3, id_name, name)
                    elif flag_cas == 1 and flag_classif == 1:
                        add_to_dict(dictionary, id_, subs_id, reg_country_id, reg_type_id, gen_reg_id, spec_reg_id, subspec_reg_id, special_cases_id, 
                                    classif_id, classif_name, cas_to_add, 1, id_name, name)
                    elif flag_index == 1 and flag_classif == 1:
                        add_to_dict(dictionary, id_, subs_id, reg_country_id, reg_type_id, gen_reg_id, spec_reg_id, subspec_reg_id, special_cases_id, 
                                    classif_id, classif_name, index_to_add, 3, id_name, name)
                    elif flag_ec == 1 and flag_classif == 1:
                        add_to_dict(dictionary, id_, subs_id, reg_country_id, reg_type_id, gen_reg_id, spec_reg_id, subspec_reg_id, special_cases_id, 
                                    classif_id, classif_name, ec_to_add, 2, id_name, name)
                    elif flag_cas == 1 and flag_ec == 1:
                        add_to_dict(dictionary, id_, subs_id, reg_country_id, reg_type_id, gen_reg_id, spec_reg_id, subspec_reg_id, special_cases_id, 
                                    classif, classif, cas_to_add, 1, id_name, name)
                        id_ += 1
                        add_to_dict(dictionary, id_, subs_id, reg_country_id, reg_type_id, gen_reg_id, spec_reg_id, subspec_reg_id, special_cases_id, 
                                    classif, classif, ec_to_add, 2, id_name, name)
                    elif flag_cas == 1 and flag_index == 1:
                        add_to_dict(dictionary, id_, subs_id, reg_country_id, reg_type_id, gen_reg_id, spec_reg_id, subspec_reg_id, special_cases_id, 
                                    classif, classif, cas_to_add, 1, id_name, name)
                        id_ += 1
                        add_to_dict(dictionary, id_, subs_id, reg_country_id, reg_type_id, gen_reg_id, spec_reg_id, subspec_reg_id, special_cases_id, 
                                    classif, classif, index_to_add, 3, id_name, name)
                    elif flag_ec == 1 and flag_index == 1:
                        add_to_dict(dictionary, id_, subs_id, reg_country_id, reg_type_id, gen_reg_id, spec_reg_id, subspec_reg_id, special_cases_id, 
                                    classif, classif, ec_to_add, 2, id_name, name)
                        id_ += 1
                        add_to_dict(dictionary, id_, subs_id, reg_country_id, reg_type_id, gen_reg_id, spec_reg_id, subspec_reg_id, special_cases_id, 
                                    classif, classif, index_to_add, 3, id_name, name)
                    elif flag_classif == 1:
                        add_to_dict(dictionary, id_, subs_id, reg_country_id, reg_type_id, gen_reg_id, spec_reg_id, subspec_reg_id, special_cases_id, 
                                    classif_id, classif_name, casrn, ec, id_name, name)
                    elif flag_cas == 1:
                        add_to_dict(dictionary, id_, subs_id, reg_country_id, reg_type_id, gen_reg_id, spec_reg_id, subspec_reg_id, special_cases_id, 
                                    classif, classif, cas_to_add, 1, id_name, name)
                    elif flag_index == 1:
                        add_to_dict(dictionary, id_, subs_id, reg_country_id, reg_type_id, gen_reg_id, spec_reg_id, subspec_reg_id, special_cases_id, 
                                    classif, classif, index_to_add, 3, id_name, name)
                    elif flag_ec == 1:
                        add_to_dict(dictionary, id_, subs_id, reg_country_id, reg_type_id, gen_reg_id, spec_reg_id, subspec_reg_id, special_cases_id, 
                                    classif, classif, ec_to_add, 2, id_name, name)
                    else:
                        add_to_dict(dictionary, id_, subs_id, reg_country_id, reg_type_id, gen_reg_id, spec_reg_id, subspec_reg_id, special_cases_id, 
                                    classif, classif, casrn, ec, id_name, name)
    return id_