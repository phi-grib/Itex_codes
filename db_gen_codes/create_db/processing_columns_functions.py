"""
Functions that process the different columns of the Inditex Excel.

Created by: Eric March Vila (eric.march@upf.edu)
On: 02/04/2019, 10:00 AM
"""

import numpy as np
import pandas as pd
import re

from typing import Optional, Union, Tuple, Callable

#
# START Annex 3 spec functions.
#

def get_suspected(element: str) -> Optional[str]:
    """ 
    Gets the element with Suspect within, otherwise None
    """
    if 'Suspected' in element:
        return element
    else:
        return None

def get_suspect_set(row: list) -> set:
    """
    Uses the funciton above to check over the elements of the list (df row). Stores found elements in set to avoid duplicates
    """
    suspected_set = set()
    for element in row:
        suspect = get_suspected(element)
        if suspect:
            suspected_set.add(suspect)
    return suspected_set

def clean_asbestos_row(row: list) -> list:
    """
    It treates the asbestos string specifically since it has its own format. Splits everything in the proper way and removes it from the previous row.
    Extends new_row with the new elements and returns new_row.
    """
    new_row = row
    for i, el in enumerate(row):
        if 'asbestos' in el:
            index_to_remove = i
            row_to_format = el.replace('#', ' ').strip().split()
            ec_num, cas_num = row_to_format[1:3]
            classif_1 = ' '.join(row_to_format[3:7])
            classif_2 = ' '.join(row_to_format[7:])
            new_row.extend([ec_num, classif_1, classif_2, cas_num, classif_1, classif_2])
            del new_row[index_to_remove]
            break

    return new_row

def clean_elements_set(row: list, list_: Optional[list] = None) -> Union[list,set]:
    """
    General function to split the input row into different readable elements that will be stored either in a set or in a list,
    depending on which function is calling this one.
    """
    start_rgx = re.compile(r'^\s?-?\s?\(?')
    end_rgx = re.compile(r'#?\)?-?:?$')
    element_set = set()
    
    if list_:
        element_list = []
    for element in row:
        clean_element = start_rgx.sub('',element)
        clean_element = end_rgx.sub('',clean_element)
        clean_element = clean_element.replace('#', ', ').lstrip().strip().strip(':').split(', ')
        if clean_element:
            for el in clean_element:
                if 'asbestos' in el:
                    asbestos_row = el.split(' ', 3)[1:]
                    for asb_el in asbestos_row:
                        element_set.add(asb_el)
                        if list_:
                            element_list.append(asb_el)
                else:
                    element_set.add(el)
                    if list_:
                        element_list.append(el)
    
    if list_:
        return element_list
    else:
        return element_set

#
# END Annex 3 spec functions.
#
#
# START general eval functions
#

def get_list_cleaned(element: str, raw = None) -> Optional[str]:
    """
    Check if word 'list' is in the string and returns it. Otherwise None
    """
    if raw:
        if 'list' in element or 'List' in element:
            return element
        else:
            return None
    else:
        if 'list' in element:
            clean_list_string = element.replace(':','. ').split('. ')[1]
            return clean_list_string
        elif 'List' in element:
            clean_list_string = element.split('. ')[1]
            return clean_list_string
        else:
            return None

def get_list_set(row: list, raw = None) -> set:
    """
    Gets the strings with the word list within it and stores them into a set.
    """
    list_set = set()
    for element in row:
        list_ = get_list_cleaned(element, raw)
        if list_:
            list_set.add(list_)
    return list_set

#
# END general eval functions
#
#
# START SVHC prelim functions
#

def get_add_info_svhc(element: str, raw = None) -> Optional[str]:
    """
    Use the pattern to find all the additional information on this column
    """
    if 'Submitted' in element or 'Withdrawn' in element:
        return element
    else:
        return None

def get_svhc_add_info_set(row: list, raw = None) -> set:
    """
    Using the above function it generates a set to avoid duplicates.
    """
    new_set = set()
    for element in row:
        new_element = get_add_info_svhc(element, raw)
        if new_element:
            new_set.add(new_element)

    return new_set

def clean_svhc_cl_preliminary(element: str) -> str:
    clean_element = element.lstrip('1.').lstrip('2.').strip()
    return clean_element

#
# END SVHC prelim functions
#
#
# START SVHC def functions
#

def clean_elements_svhc_def(element: str) -> str:
    """
    Instead of using general cleaning function, we use this one to clean the elements of this column
    """
    info_pattern = re.compile(r"^\d\.\s?")
    clean_info = info_pattern.sub('', element)
    clean_info = clean_info.lstrip('- ')
    return clean_info

#
# END SVHC def functions
#
#
# START SVHC prelim annex xiv functions
#

def clean_svhc_prelim_axiv(element: str) -> str:
    start_pattern = re.compile(r"^\d\.\s")
    clean_element = start_pattern.sub('', element)
    clean_element = clean_element.replace('(Recommended', 'Recommended').replace('July 2015)', 'July 2015').strip()
    return clean_element

#
# END SVHC prelim annex xiv functions
#
#
# START SVHC subject authorisation annex xiv functions
#

def get_date_string(element: str) -> Optional[str]:
    if 'date' in element or 'Exempted' in element:
        return element
    else:
        return None

def get_add_info_svhc_subauth(row: set) -> set:
    _set = set()
    for element in row:
        new_el = get_date_string(element)
        if new_el:
            _set.add(new_el)
    return _set

def clean_svhc_subauth_axiv(element: str) -> str:
    clean_element = element.replace(':','').strip()
    return clean_element

#
# END SVHC subject authorisation annex xiv functions
#
#
# START SVHC approved/rejected/application authorisation annex xiv functions
#

def clean_svhc_approved_auth_axiv(element: str) -> str:
    clean_element = element.replace('-','').strip()
    if 'No revelant authorisation' in clean_element:
        clean_element = 'No relevant authorisation'
    return clean_element

#
# END SVHC approved/rejected/application authorisation annex xiv functions
#
#
# START annex XVII preliminary/substance-mixture functions
#

def clean_add_info_axvii_prelim(element: str, raw = None) -> Optional[str]:
    start_pattern = r"^\d{1}\.\s.+"
    add_info = re.findall(start_pattern, element)
    if add_info:
        if raw:
            return add_info[0]
        else:
            sub_pattern = re.compile(r"^\d{1}\.\s")
            clean_add_info = sub_pattern.sub('',add_info[0]).strip(':').strip(': ').strip()
            return clean_add_info
    else:
        return None

def get_add_info_axvii_prelim_set(row: set, raw = None) -> set:
    new_set = set()
    for element in row:
        new_element = clean_add_info_axvii_prelim(element, raw)
        if new_element:
            new_set.add(new_element)

    return new_set

def clean_axvii_preliminary(element: str) -> str:
    start_pattern = r"^\d{1}\.\d\s.+"
    found_element = re.findall(start_pattern, element)
    if found_element:
        sub_pattern = re.compile(r"^\d{1}\.\d\s")
        clean_element = sub_pattern.sub('',found_element[0]).strip(':').strip(': ').strip()
        return clean_element
    else:
        replaced_element = element.replace("‘","").replace("’","").replace('revelant','relevant').replace("   "," ")
        clean_element = replaced_element.strip().lstrip('- ').strip(':').strip(';').strip("'").rstrip(',').rstrip('.')
        
        return clean_element

#
# END annex XVII preliminary/substance-mixture functions
#
#
# START registration dossier/notification CLP functions
#

def clean_add_info_registration_dossier_clp(element: str, raw = None) -> Optional[str]:
    start_pattern = r"^\w{1}\)\s.+"
    add_info = re.findall(start_pattern, element)
    if add_info:
        if raw:
            return add_info[0]
        else:
            sub_pat = re.compile(r"^\w{1}\)\s")
            clean_add_info = sub_pat.sub('',add_info[0])
            return clean_add_info
    else:
        return None

def get_add_info_registration_dossier_clp(row: set, raw = None) -> set:
    new_set = set()
    for element in row:
        new_element = clean_add_info_registration_dossier_clp(element, raw)
        if new_element:
            new_set.add(new_element)
    
    return new_set

def clean_registration_dossier_clp(element: str) -> str:
    clean_element = element.replace('( ','(').replace('.2','. 2').replace('(',' (').replace('   ',' ').replace('  ',' ').strip(': ').strip(':').strip()
    cas_backslash_pattern = re.compile(r"\d\/.+")
    rep_pattern = re.compile(r"^.+\d{1}\s[A-Z]{1}$")
    clean_element_replaced = re.findall(rep_pattern,clean_element)
    clean_element_backslash_rep = re.findall(cas_backslash_pattern,clean_element)
    if clean_element_replaced:
        clean_element = rep_pattern.sub('', clean_element_replaced[0])
        return clean_element
    elif clean_element_backslash_rep:
        clean_element = cas_backslash_pattern.sub('', clean_element_backslash_rep[0])
        return clean_element
    else:
        return clean_element

#
# END registration dossier/notification CLP functions
#
#
# START harmonised C&L CLP functions
#

def clean_add_info_harm_cl(element: str, raw = None) -> Optional[Union[str, Tuple[str,set]]]:
    if element.startswith('Entry'):
        if raw:
            return element
        else:
            add_info_list = element.split('as')
            entry_string = add_info_list[0].strip()
            if len(add_info_list) > 1:
                endpoint_string = add_info_list[1].lstrip().split(';')
                endpoint_string = set([el.lstrip() for el in endpoint_string])
            else:
                endpoint_string = set()
            return entry_string, endpoint_string
    else:
        return None, None

def get_add_info_harm_cl(row: set, raw = None) -> set:
    info_set = set()
    endpoint_set = set()
    for element in row:
        if raw:
            raw_element = clean_add_info_harm_cl(element, raw)
            if raw_element:
                info_set.add(raw_element)
        else:
            clean_element, endpoints = clean_add_info_harm_cl(element)
            if clean_element:
                info_set.add(clean_element)
            if endpoints:
                endpoint_set.update(endpoints)
    return info_set, endpoint_set

def clean_harmonized_cl_clp(element: str) -> list:
    clean_element = element.replace('1.','1,').replace(' –',',').strip(':').strip('.').strip().strip(',').split(';')
    clean_element = [c_el.lstrip().rstrip() for c_el in clean_element]
    return clean_element

#
# END harmonised C&L CLP functions
#
#
# START harmonised index number CLP functions
#

def clean_add_info_index(element: str) -> Optional[str]:
    if element.startswith('Amendment'):
        return element
    else:
        return None

def get_add_info_index(row: set) -> set:
    new_set = set()
    for element in row:
        info_element = clean_add_info_index(element)
        if info_element:
            new_set.add(info_element)
    return new_set

def clean_harmonized_index_number_clp(element: str) -> Union[str,list]:
    index_pattern = r"[0-9]{1,3}-[0-9]{1,3}-\d{2}-[0-9A-Z]"
    clean_element = element.replace('Containing', 'containing').replace('w) particles','w) of particles').strip(': ').strip(':').lstrip('-').lstrip('‘').strip()
    index_ = re.findall(index_pattern, clean_element)
    if index_:
        list_element = clean_element.split(' ',1)
        if len(list_element) > 1:
            clean_element = list_element
    else:
        clean_element = clean_element.split('/ ')

    return clean_element

#
# END harmonised index number CLP functions
#
#
# START hazard class CLP functions
#

def clean_hazard_class_clp(element: str) -> Union[str,list]:
    stripped_element = element.lstrip('-').strip(': ').strip()
    replaced_element = stripped_element.replace('AcuteTox.','Acute Tox. ').replace('STOTRE','STOT RE').\
        replace('RE2', 'RE 2').replace('   (',' (').replace('  (',' (').replace('( ','(')
    clean_element = replaced_element.split('/ ')

    return clean_element

#
# END hazard class CLP functions
#
#
# START hazard CLP functions
#

def clean_hazard_clp(element: str) -> Union[str,list]:
    number_pattern = re.compile(r'^([0-9]{3}.*?[.:])')
    subs_pattern = re.compile(r'\d\:\w')
    stripped_element = element.lstrip('-').strip(': ').strip('.').replace(':  ',': ').replace('( ','(').replace('.',':',1).\
        replace('Suspect ed', 'Suspected').replace(' Df:','Df:').replace('inhalated','inhaled').replace('inhale','inhaled').replace('dd','d').\
            rstrip('H').rstrip(' 2').strip()
    sub_el = re.findall(subs_pattern, stripped_element)
    numb_el = re.findall(number_pattern, stripped_element)
    if sub_el:
        stripped_element = stripped_element.replace(':',': ')
    if numb_el:
        stripped_element = ''.join(['H',stripped_element])
        
    clean_element = stripped_element.split('/ ')

    return clean_element

#
# END hazard CLP functions
#
#
# START specific concentration limits CLP functions
#

def clean_conc_limits_clp(element: str) -> str:
    clean_element = element.lstrip('-').strip(': ').strip(';').strip('*').strip('’').strip()
    replaced_element = clean_element.replace('C', ' C ',1).replace('C  ', 'C ',1).replace('  C', ' C',1).replace(' C arc','Carc').replace('M=','M =').\
        replace(' %','%').replace('m=','M =').replace('=1', '= 1').replace('. Liq','.Liq').replace('C orr', 'Corr').replace('<','< ').replace('<  ','< ').\
            replace('  ≥',' ≥').replace('≥','≥ ').replace('≥  ','≥ ').replace('.1','. 1')
    if 'STOTSE3' in replaced_element:
        replaced_element = 'STOT SE 3'
    return replaced_element

#
# END specific concentration limits CLP functions
#
#
# START notes CLP functions
#

def clean_add_info_notes(element: str) -> Optional[str]:
    if element.startswith('Note'):
        return element
    else:
        return None

def get_add_info_notes(row: set) -> set:
    new_set = set()
    for element in row:
        info_element = clean_add_info_notes(element)
        if info_element:
            new_set.add(info_element)
    return new_set

def clean_notes_clp(element: str) -> str:
    clean_element = element.replace('NOTE','Note').replace('*','').strip(':').strip('.').strip()
    return clean_element

#
# END notes CLP functions
#
#
# START POPs functions
#

def clean_pops_regulations(element: str) -> str:
    clean_element = element.strip()
    if 'limit :' in clean_element:
        clean_element = clean_element.replace('limit :','limit:')
    return clean_element

#
# END POPs functions
#
#
# START WFD functions
#

def clean_WFD_regulations(element: str) -> str:
    clean_element = element.replace(' :',':').replace(',','.').replace('Σ = ','Σ=').strip(':').strip()
    return clean_element

#
# END WFD functions
#
#
# START pbt/vpvb functions
#

def clean_pbt_vpvb_reg_info(element: str) -> str:
    clean_element = element.replace(' /','/').replace('  ',' ').strip(':').strip()
    return clean_element

#
# END pbt/vpvb assessment reach functions
#
#
# START endocrine disruption functions
#

def clean_endoc_dis_info(element: str) -> str:
    clean_element = element.replace('CAT','CAT ').replace('  ',' ').strip(':').strip(': ').strip()
    return clean_element

#
# END endocrine disruption functions
#
#
# START Labelling CLP pictogram/signal word functions
#

def clean_lab_clp_pictogram(element: str) -> str:
    clean_element = element.replace(':  ',': ').strip(':').strip('-').strip()
    if ',' in clean_element:
        clean_element = [c_el.strip() for c_el in clean_element.split(',')]
    elif ';' in clean_element:
        clean_element = [c_el.strip() for c_el in clean_element.split(';')]
    elif '/ ' in clean_element:
        clean_element = [c_el.strip() for c_el in clean_element.split('/')]

    return clean_element

#
# END Labelling CLP pictogram/signal word functions
#
#
# START Labelling CLP hazard statement functions
#

def clean_lab_clp_hazard_statement(element: str) -> str:
    number_pattern = re.compile(r'^([0-9]{3}.*?[.:])')
    clean_element = element.replace(':M',': M').replace(':  ',': ').replace(' : ',': ').replace(' – ',': ').replace('( ','(').replace("'","").\
        replace('inhalated','inhaled').replace('inhale','inhaled').replace('dd','d').replace('Suspect ed','Suspected').replace('H302.','H302:').\
        replace('lasting effects','lasting effect').strip(':').strip('-').strip('.').strip()
    sub_el = re.findall(number_pattern, clean_element)
    if sub_el:
        clean_element = ''.join(['H',clean_element])
    if ';' in clean_element:
        clean_element = [c_el.strip() for c_el in clean_element.split(';')]
    elif '/ ' in clean_element:
        clean_element = [c_el.strip() for c_el in clean_element.split('/')]

    return clean_element

#
# END Labelling CLP hazard statement functions
#
#
# START Labelling CLP hazard supplementary functions
#

def clean_lab_clp_hazard_supplementary(element: str) -> str:
    clean_element = element.lstrip('-').strip(':').strip()

    return clean_element

#
# END Labelling CLP hazard supplementary functions
#
#
# START Swedish regulations chemical products restriction submix functions
#

def clean_add_info_sw_restrict_submix(element: str) -> Optional[str]:
    if 'The prohibition shall not apply' in element:
        return element
    else:
        return None

def get_add_info_sw_restrict_submix(row: set) -> set:
    new_set = set()
    for element in row:
        info_element = clean_add_info_sw_restrict_submix(element)
        if info_element:
            new_set.add(info_element)
    return new_set

def clean_sw_classif_restric_submix(element: str) -> str:
    clean_element = element.lstrip('-').strip(': ').strip('.').strip()

    return clean_element

#
# END Swedish regulations chemical products restriction submix functions
#
#
# START German substance classification GCO restricted substances in articles functions
#

def clean_gr_chemicals_act(element: str) -> Optional[str]:
    if 'Substances, mixtures and products' in  element or 'If condition ' in element or 'If conditions ' in element:
        return element
    else:
        return None

def get_add_info_gr_chemicals_act(row: set) -> set:
    new_set = set()
    for element in row:
        info_element = clean_gr_chemicals_act(element)
        if info_element:
            new_set.add(info_element)
    return new_set

def clean_gr_gco_restric_articles(element: str) -> str:
    clean_element = element.replace('menas','means').replace('and  ','and ').replace('1. ','').replace('2. ','').replace('3. ','').strip('.').strip(':').strip()

    return clean_element

#
# END German substance classification GCO restricted substances in articles functions
#
