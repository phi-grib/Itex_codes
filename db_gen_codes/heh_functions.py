"""
functions to handle Heh dataframes. Similar to functions_to_process.py
I just don't want to mess that package and I store these functions here

Created by: Eric March

"""

import numpy as np
import pandas as pd
import psycopg2
import re


def get_cas_dict(df):
    casdict = {}
    for i, row in df.iterrows():
        if type(row[2]) is not float:
            sub_id = i + 1
            cas = row[2].split()
            casdict.update({sub_id:cas})

    return casdict

def get_cas(element):
    cas_pattern = r"^[0-9]{1,6}-\d{2}-\d\.?"
    casrn = re.findall(cas_pattern, element)
    if casrn:
        return casrn
    else:
        return None

def get_ec(element):
    ec_pattern = r"^[0-9]{1,3}-\d{3}-\d$"
    ec = re.findall(ec_pattern, element)
    if ec:
        return ec
    else:
        return None
    
def get_index_number(element):
    index_pattern = r"[0-9]{1,3}-[0-9]{1,3}-\d{2}-[0-9A-Z]"
    index_ = re.findall(index_pattern, element)
    if index_:
        return index_
    else:
        return None
    
def get_cas_numbers(row):
    cas_set = set()
    for element in row:
        casrn = get_cas(element)
        if casrn:
            for number in casrn:
                cas_set.add(number.strip())
                
    return cas_set

def get_index_numbers(row):
    index_set = set()
    for element in row:
        index_ = get_index_number(element)
        if index_:
            for number in index_:
                index_set.add(number.strip())
    
    return index_set

def get_ec_numbers(row):
    ec_set = set()
    for element in row:
        ec = get_ec(element)
        if ec:
            for number in ec:
                ec_set.add(number.strip())
    
    return ec_set

def get_classification(element):
    if 'Self-classification' in element or 'Harmonised C&L' in element:
        return element
    else:
        return None

def get_classif_set(row):
    classif_set = set()
    for element in row:
        classif = get_classification(element)
        if classif:
            classif_set.add(classif)
    return classif_set

def clean_elements_set(row, list_ = None):
    element_set = set()
    if list_:
        element_list = []
    for element in row:
        clean_element = element.strip().strip(':').split(':')
        for el in clean_element:
            new_element = el.split(';')
            for e in new_element:
                final_element = e.lstrip('-').lstrip()
                final_element = re.sub(r"\sas\s?"," ", final_element)
                if final_element:
                    element_set.add(final_element.lstrip().strip())
                    if list_:
                        element_list.append(final_element)
    
    if list_:
        return element_list
    else:
        return element_set

def clean_row_df(row):
    cas_set = get_cas_numbers(row)
    index_set = get_index_numbers(row)
    ec_set = get_ec_numbers(row)
    row_set = clean_elements_set(row)
    classification_set = get_classif_set(row_set)
    no_cas_index_row_set = row_set - cas_set - index_set - ec_set - classification_set

    return no_cas_index_row_set, cas_set, index_set, ec_set, classification_set

def clean_heh_df(df, col):
    heh_clean = df.loc[:,col].apply(lambda x: str(x).strip().split('\n'))
    heh_clean_list = [sublist for li in heh_clean for sublist in li]
    heh_set, cas_set, index_set, ec_set, classification_set = clean_row_df(heh_clean_list)
    
    return heh_set, cas_set, index_set, ec_set, classification_set

def get_proper_value(variable):
    new_var = None
    try:
        new_var = variable.strip().split('\n')
        new_var = clean_elements_set(new_var, list_=True)
    except AttributeError:
        pass
    
    return new_var

def add_to_dict(dictionary, id_, subs_id, cas, index, ec, endpoint_id, classif_id, classif_name, id_name, name):
    dictionary['id'].append(id_)
    dictionary['subs_id'].append(subs_id)
    dictionary['cas_number'].append(cas)
    dictionary['index_number'].append(index)
    dictionary['ec_number'].append(ec)
    dictionary['heh_type_id'].append(endpoint_id)
    dictionary['heh_classification_id'].append(classif_id)
    dictionary['heh_classification_name'].append(classif_name)
    dictionary['heh_name_id'].append(id_name)
    dictionary['heh_name'].append(name)
                        
def add_element_to_dictionary(heh_general_dict, endpoint, endpoint_id, id_, subs_id, name_dataframe, classification_dataframe):
    for i, names in name_dataframe.iterrows():
        id_name = names[0]
        name = names[1]
        flag_cas = 0
        flag_index = 0
        flag_ec = 0
        flag_classif = 0
        if endpoint:
            for element in endpoint:
                casrn = get_cas(element)
                index_num = get_index_number(element)
                ec = get_ec(element)
                classif = get_classification(element)
                if classif:
                    flag_classif = 1
                    classif_id = classification_dataframe.loc[classification_dataframe['classification'] == classif, 'id'].iat[0]
                    classif_name = classif
                if casrn:
                    flag_cas = 1
                    cas_to_add = casrn[0]
                if index_num:
                    flag_index = 1
                    index_to_add = index_num[0]
                if ec:
                    flag_ec = 1
                    ec_to_add = ec[0]
                if name == element:
                    id_ += 1
                    if flag_cas == 1 and flag_index == 1 and flag_ec == 1 and flag_classif == 1:
                        add_to_dict(heh_general_dict, id_, subs_id, cas_to_add, index_to_add, ec_to_add, endpoint_id, classif_id,classif_name, id_name, name)
                    elif flag_cas == 1 and flag_ec == 1 and flag_classif == 1:
                        add_to_dict(heh_general_dict, id_, subs_id, cas_to_add, index_num, ec_to_add, endpoint_id, classif_id,classif_name, id_name, name)
                    elif flag_index == 1 and flag_ec == 1 and flag_classif == 1:
                        add_to_dict(heh_general_dict, id_, subs_id, casrn, index_to_add, ec_to_add, endpoint_id, classif_id, classif_name, id_name, name)
                    elif flag_cas == 1 and flag_index == 1 and flag_ec == 1:
                        add_to_dict(heh_general_dict, id_, subs_id, cas_to_add, index_to_add, ec_to_add, endpoint_id, classif, classif, id_name, name)
                    elif flag_cas == 1 and flag_classif == 1:
                        add_to_dict(heh_general_dict, id_, subs_id, cas_to_add, index_num, ec, endpoint_id, classif_id,classif_name, id_name, name)
                    elif flag_index == 1 and flag_classif == 1:
                        add_to_dict(heh_general_dict, id_, subs_id, casrn, index_to_add, ec, endpoint_id, classif_id, classif_name, id_name, name)
                    elif flag_ec == 1 and flag_classif == 1:
                        add_to_dict(heh_general_dict, id_, subs_id, casrn, index_num, ec_to_add, endpoint_id, classif_id, classif_name, id_name, name)
                    elif flag_classif == 1:
                        add_to_dict(heh_general_dict, id_, subs_id, casrn, index_num, ec, endpoint_id, classif_id, classif_name, id_name, name)
                    elif flag_cas == 1:
                        add_to_dict(heh_general_dict, id_, subs_id, cas_to_add, index_num, ec, endpoint_id, classif, classif, id_name, name)
                    elif flag_index == 1:
                        add_to_dict(heh_general_dict, id_, subs_id, casrn, index_to_add, ec, endpoint_id, classif, classif, id_name, name)
                    elif flag_ec == 1:
                        add_to_dict(heh_general_dict, id_, subs_id, casrn, index_num, ec_to_add, endpoint_id, classif, classif, id_name, name)
                    else:
                        add_to_dict(heh_general_dict, id_, subs_id, casrn, index_num, ec, endpoint_id, classif, classif, id_name, name)
        
    return id_