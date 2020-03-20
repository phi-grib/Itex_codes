import numpy as np
import pandas as pd

from typing import Optional, Union, Tuple, Callable


def add_info_dict(empty_dict: dict, id_:int, endpoint: str, regulation: Optional[str], element: str):
    """
        This functions adds the different information into the dictionary

        :param empty_dict: dictionary to store the annotations. Will be converted into a dataframe
        :param id_: substance id
        :param endpoint: endpoint we're working on (CMR, PBT, vPvB, Endocrine disruptor, Sensitiser, Other)
        :param element: regulation obtained from the database
    """

    empty_dict['subs_id'].append(id_)
    empty_dict['endpoint_type'].append(endpoint)
    empty_dict['regulation'].append(regulation)
    empty_dict['name'].append(element)

def fixed_generate_dataframe(substance_id_list: np.ndarray, endpoint: str, regulations_df: pd.DataFrame, endpoint_annotations: pd.DataFrame) -> dict:
    """
        General function that replicate HEH annotations from Excel

        :param substance_id_list: list with substances id
        :param endpoint: endpoint we're working on (CMR, PBT, vPvB, Endocrine disruptor, Sensitiser, Other)
        :param regulation_df: datafame with regulations from the dataframe
        :param endpoint_annotations: dataframe with annotations for each endpoint coming from HEH in the database

        :returns captured_ann_dict: dataframe containing annotations that replicates original HEH
    """

    captured_ann_dict = {'subs_id':[],'endpoint_type':[],'regulation':[],'name':[]}

    for substance_id in substance_id_list:
        cmr_flag = 0
        no_flag = 0
        
        # Generate dataframes from different regulation sources
        
        svhc_cand_list = regulations_df.loc[(regulations_df['subs_id'] == substance_id) &
                        (regulations_df['specific_regulation_name'] == 'svhc') &
                        (regulations_df['subspecific_regulation_name'] == 'candidate_list'), ['specific_regulation_name','names']]
        
        pbt_vpvb_endoc = regulations_df.loc[(regulations_df['subs_id'] == substance_id) &
                        (regulations_df['general_regulation_name'].isin(['pbt_vpvb','endocrine_disruptors'])), ['specific_regulation_name','names']]

        reg_dos_not = regulations_df.loc[(regulations_df['subs_id'] == substance_id) &
                        (regulations_df['special_cases_name'].isin(['registration_dossier'])), ['additional_information_name','names']]
        
        not_df = regulations_df.loc[(regulations_df['subs_id'] == substance_id) &
                (regulations_df['special_cases_name'].isin(['notification'])), ['additional_information_name','names']]

        draft_df = regulations_df.loc[(regulations_df['subs_id'] == substance_id) &
                  (regulations_df['additional_information_name'].isin(['Submitted SVHC intentions',
                  'Withdrawn SVHC intentions and submissions','Amendment 2016/1179'])), ['additional_information_name','names']]
        
        # First step: check in SVHC candidate list if it's CMR
        
        for element in endpoint_annotations:
            if element in svhc_cand_list.names.values:
                regulation = svhc_cand_list.loc[svhc_cand_list['names'] == element, 'specific_regulation_name'].iat[0]
                cmr_flag = 1
                add_info_dict(captured_ann_dict, substance_id, endpoint, regulation, element)
            
            if cmr_flag == 1:
                
                if element in reg_dos_not.names.values:
                    reg_rd = reg_dos_not.loc[reg_dos_not['names'] == element, 'additional_information_name'].iat[0]
                    cmr_flag = 1
                    format_reg_rd = ' '.join([reg_rd,'(RD)'])
                    add_info_dict(captured_ann_dict, substance_id, endpoint, format_reg_rd, element)
                    
                elif element in not_df.names.values :
                    reg_n = not_df.loc[not_df['names'] == element, 'additional_information_name'].iat[0]
                    cmr_flag = 1
                    format_reg_n = ' '.join([reg_n,'(N)'])
                    add_info_dict(captured_ann_dict, substance_id, endpoint, format_reg_n, element)
                    
        if cmr_flag == 1:
            add_info_dict(captured_ann_dict, substance_id, endpoint, None, 'YES')
            
        # Second step: check in Harmonized C&L if it's CMR

        if cmr_flag == 0:
            
            clp_reg = regulations_df.loc[(regulations_df['subs_id'] == substance_id) &
                    (regulations_df['specific_regulation_name'] == 'harmonised_c&l') &
                    (regulations_df['subspecific_regulation_name'].isin(['hazard_class'])), ['specific_regulation_name','names']]
            
            for element in endpoint_annotations:
                if element in clp_reg.names.values:
                    regulation = clp_reg.loc[clp_reg['names'] == element, 'specific_regulation_name'].iat[0]
                    cmr_flag = 1
                    add_info_dict(captured_ann_dict, substance_id, endpoint, regulation, element)
                
                if cmr_flag == 1:
                    
                    if element in reg_dos_not.names.values:
                        reg_rd = reg_dos_not.loc[reg_dos_not['names'] == element, 'additional_information_name'].iat[0]
                        cmr_flag = 1
                        format_reg_rd = ' '.join([reg_rd,'(RD)'])
                        add_info_dict(captured_ann_dict, substance_id, endpoint, format_reg_rd, element)

                    elif element in not_df.names.values:
                        reg_n = not_df.loc[not_df['names'] == element, 'additional_information_name'].iat[0]
                        cmr_flag = 1
                        format_reg_n = ' '.join([reg_n,'(N)'])
                        add_info_dict(captured_ann_dict, substance_id, endpoint, format_reg_n, element)
                        
            if cmr_flag == 1:
                add_info_dict(captured_ann_dict, substance_id, endpoint, None, 'YES')
                
        # Third step: check Registration Dossier (REACH) and Notification (CLP) for CMR. First of all, we check
        # if there is ANY annotation within these regulations and then check if this annotations are related to any
        # endpoint of interest.

        if cmr_flag == 0:

            if 'No Notification' not in not_df.names.values and 'No Registration Dossier' not in reg_dos_not.names.values:
                pending_flag = '' 
                for element in endpoint_annotations:
                    if element in not_df.names.values and element in reg_dos_not.names.values:
                        reg_rd = reg_dos_not.loc[reg_dos_not['names'] == element, 'additional_information_name'].iat[0]
                        reg_n = not_df.loc[not_df['names'] == element, 'additional_information_name'].iat[0]
                        format_reg_rd = ' '.join([reg_rd,'(RD)'])
                        format_reg_n = ' '.join([reg_n,'(N)'])
                        add_info_dict(captured_ann_dict, substance_id, endpoint, format_reg_rd, element)
                        add_info_dict(captured_ann_dict, substance_id, endpoint, format_reg_n, element)
                        pending_flag = 'Pending (3b)'
                        
                    elif element in reg_dos_not.names.values and 'Not included' in not_df.names.values and 'No information' in not_df.names.values\
                    and 'Not information available' in not_df.names.values:
                        regulation = reg_dos_not.loc[reg_dos_not['names'] == element, 'additional_information_name'].iat[0]
                        cmr_flag = 1
                        format_reg_rd = ' '.join([regulation,'(RD)'])
                        add_info_dict(captured_ann_dict, substance_id, endpoint, format_reg_rd, element)
                        pending_flag = 'Pending (2)'
                    
                    elif element in not_df.names.values and 'Not included' in reg_dos_not.names.values and 'No information' in reg_dos_not.names.values\
                    and 'Not information available' in reg_dos_not.names.values:
                        regulation = not_df.loc[not_df['names'] == element, 'additional_information_name'].iat[0]
                        cmr_flag = 1
                        format_reg_n = ' '.join([regulation,'(N)'])
                        add_info_dict(captured_ann_dict, substance_id, endpoint, format_reg_n, element)
                        pending_flag = 'Pending (4)'
                        
                    elif element in reg_dos_not.names.values and element not in not_df.names.values:
                        reg_rd = reg_dos_not.loc[reg_dos_not['names'] == element, 'additional_information_name'].iat[0]
                        cmr_flag = 1
                        format_reg_rd = ' '.join([reg_rd,'(RD)'])
                        add_info_dict(captured_ann_dict, substance_id, endpoint, format_reg_rd, element)
                        pending_flag = 'Pending (3a)'
                        
                    elif element in not_df.names.values and element not in reg_dos_not.names.values:
                        reg_n = not_df.loc[not_df['names'] == element, 'additional_information_name'].iat[0]
                        cmr_flag = 1
                        format_reg_n = ' '.join([reg_n,'(N)'])
                        add_info_dict(captured_ann_dict, substance_id, endpoint, format_reg_n, element)
                        pending_flag = 'Pending (3c)'

                    elif element not in not_df.names.values and element not in reg_dos_not.names.values:
                        no_flag = 1

                if cmr_flag == 1:
                    add_info_dict(captured_ann_dict, substance_id, endpoint, None, pending_flag)

        # Fourth step: check in SVHC-draft and CLP-draft (Amendment) if CMR
        
        if cmr_flag == 0:
            pending_flag = '' 
            for element in endpoint_annotations:
                if element in draft_df.names.values:
                    reg_n = draft_df.loc[draft_df['names'] == element, 'additional_information_name'].iat[0]
                    cmr_flag = 1
                    add_info_dict(captured_ann_dict, substance_id, endpoint, reg_n, element)
                    pending_flag = 'Pending (1)'
            
            if cmr_flag == 1:
                add_info_dict(captured_ann_dict, substance_id, endpoint, None, pending_flag)
        
        # Fifth step: check for PBT-vPvB-Endocrine disruptors regulations

        if cmr_flag == 0:
            for element in endpoint_annotations:
                if element in pbt_vpvb_endoc.names.values:
                    reg_n = pbt_vpvb_endoc.loc[pbt_vpvb_endoc['names'] == element, 'specific_regulation_name'].iat[0]
                    cmr_flag = 1
                    add_info_dict(captured_ann_dict, substance_id, endpoint, reg_n, element)
                
            if cmr_flag == 1:
                add_info_dict(captured_ann_dict, substance_id, endpoint, None, 'YES')
        
        if cmr_flag == 0:
            if no_flag == 0:
                add_info_dict(captured_ann_dict, substance_id, endpoint, None, 'No information')
            elif no_flag == 1:
                add_info_dict(captured_ann_dict, substance_id, endpoint, None, 'NO')

    return captured_ann_dict

def compare_heh_dataframes(original_heh_df: pd.DataFrame, created_heh_df: pd.DataFrame, substance_id_list: np.ndarray, endpoint: Optional[str] = None) -> dict:
    """
        function that compares HEH dataframes annotations and created a dict with the results of the comparison

        :param original_heh_df: HEH dataframe from Inditex Excel/CII database
        :param created_heh_df: HEH dataframe obtained from implementing manually USC workflow
        :param substance_id_list: list with substance id's
        :param endpoint: option to obtain a comparison dict by specific enpodint. Endpoints available are:
                        - CMR
                        - PBT
                        - vPvB
                        - Endocrine_disruptor
                        - Sensitiser 
                        - Other

        :returns comparison_dict: dictionary containing comparison informaton. To be converted into pandas DataFrame
    """

    comparison_dict = {'subs_id':[],'heh_excel_annotation_number':[], 'heh_generated_annotation_number':[], 'shared_annotations':[],'%_recuperation':[],
                       'non_shared_old_annotations':[],'old_annotations_not_captured':[],'non_shared_new_annotations':[],'new_annotations_from_generated_heh':[]}

    
    for subs_id in substance_id_list:
        if endpoint:
            org_names = set(original_heh_df[(original_heh_df['subs_id'] == subs_id) &
                                            (original_heh_df['type'] == endpoint)].name.unique())
            new_names = set(created_heh_df[(created_heh_df['subs_id'] == subs_id) &
                                           (created_heh_df['endpoint_type'] == endpoint)].name.unique())
        else:
            org_names = set(original_heh_df[original_heh_df['subs_id'] == subs_id].name.unique())
            new_names = set(created_heh_df[created_heh_df['subs_id'] == subs_id].name.unique())

        intersec_ = org_names.intersection(new_names)
        comparison_dict['subs_id'].append(subs_id)
        comparison_dict['heh_excel_annotation_number'].append(len(org_names))
        comparison_dict['heh_generated_annotation_number'].append(len(new_names))
        comparison_dict['shared_annotations'].append(len(intersec_))
        
        if len(org_names) == 0:
            comparison_dict['%_recuperation'].append(100)
            comparison_dict['non_shared_old_annotations'].append(None)
            comparison_dict['old_annotations_not_captured'].append(None)
            comparison_dict['non_shared_new_annotations'].append(None)
            comparison_dict['new_annotations_from_generated_heh'].append(None)
        
        elif len(new_names) > len(org_names):
            diff_ann = list(new_names.difference(org_names))
            diff_new_str = str(diff_ann).strip('[').strip(']').replace("'","")
            comparison_dict['%_recuperation'].append(100)
            comparison_dict['non_shared_old_annotations'].append(None)
            comparison_dict['old_annotations_not_captured'].append(None)
            comparison_dict['non_shared_new_annotations'].append(len(diff_ann))
            comparison_dict['new_annotations_from_generated_heh'].append(diff_new_str)

        else:
            if len(intersec_) == 0:
                comparison_dict['%_recuperation'].append(0)
            else:
                recuperation_ = (len(new_names)/len(org_names))*100
                comparison_dict['%_recuperation'].append(recuperation_)

            old_diff = list(org_names.difference(new_names))
            if old_diff:
                comparison_dict['non_shared_old_annotations'].append(len(old_diff))
                diff_old_str = str(old_diff).strip('[').strip(']').replace("'","")
                comparison_dict['old_annotations_not_captured'].append(diff_old_str)
            else:
                comparison_dict['non_shared_old_annotations'].append(None)
                comparison_dict['old_annotations_not_captured'].append(None)

            new_diff = list(new_names.difference(org_names))
            if new_diff:
                comparison_dict['non_shared_new_annotations'].append(len(new_diff))
                diff_new_str = str(new_diff).strip('[').strip(']').replace("'","")
                comparison_dict['new_annotations_from_generated_heh'].append(diff_new_str)
            else:
                comparison_dict['non_shared_new_annotations'].append(None)
                comparison_dict['new_annotations_from_generated_heh'].append(None)


        # if len(new_names) > len(org_names):
        #     diff_ann = list(new_names.difference(org_names))
        #     comparison_dict['%_recuperation'].append(100)
        #     comparison_dict['non_shared_old_annotations'].append(None)
        #     comparison_dict['old_annotations_not_captured'].append(None)
        #     comparison_dict['non_shared_new_annotations'].append(len(diff_ann))
        #     comparison_dict['new_annotations_from_generated_heh'].append(diff_ann)
        # elif len(new_names) == len(org_names) == 0:
        #     comparison_dict['%_recuperation'].append(100)
        #     comparison_dict['non_shared_old_annotations'].append(None)
        #     comparison_dict['old_annotations_not_captured'].append(None)
        #     comparison_dict['non_shared_new_annotations'].append(None)
        #     comparison_dict['new_annotations_from_generated_heh'].append(None)
        # else:
        #     recuperation_ = (len(new_names)/len(org_names))*100
        #     comparison_dict['%_recuperation'].append(recuperation_)
        #     comparison_dict['non_shared_old_annotations'].append(None)
        #     comparison_dict['old_annotations_not_captured'].append(None)
        #     comparison_dict['non_shared_new_annotations'].append(None)
        #     comparison_dict['new_annotations_from_generated_heh'].append(None)
    
    return comparison_dict

def compare_unique_annotations_dataframes(original_heh_df: pd.DataFrame, created_heh_df: pd.DataFrame, endpoint: Optional[str] = None) -> dict:
    """
        Compares which unique annotations are shared between full dataframes instead of compound per compound

        :param original_heh_df: HEH dataframe from Inditex Excel/CII database
        :param created_heh_df: HEH dataframe obtained from implementing manually USC workflow
        :param endpoint: option to obtain a comparison dict by specific enpodint. Endpoints available are:
                        - CMR
                        - PBT
                        - vPvB
                        - Endocrine_disruptor
                        - Sensitiser 
                        - Other

        :returns comparison_dict: dictionary containing comparison informaton. To be converted into pandas DataFrame
    """

    comparison_dict = {'dataframe':[],'heh_excel_annotation_number':[], 'heh_generated_annotation_number':[], 'shared_annotations':[], '%_recuperation':[],
                       'non_shared_old_annotations':[], 'old_annotations_not_captured':[],'non_shared_new_annotations':[], 'new_annotations':[]}
    
    
    if endpoint:
        org_names = set(original_heh_df[original_heh_df['type'] == endpoint].name.unique())
        new_names = set(created_heh_df[created_heh_df['endpoint_type'] == endpoint].name.unique())
        dataframe_name = '_'.join([endpoint,'dataframe'])
    else:
        org_names = set(original_heh_df.name.unique())
        new_names = set(created_heh_df.name.unique())
        dataframe_name = 'whole_dataframe'

    intersec_ = org_names.intersection(new_names)
    difference_ = org_names.difference(new_names)
    comparison_dict['dataframe'].append(dataframe_name)
    comparison_dict['heh_excel_annotation_number'].append(len(org_names))
    comparison_dict['heh_generated_annotation_number'].append(len(new_names))
    comparison_dict['shared_annotations'].append(len(intersec_))

    recuperation_ = (len(new_names)/len(org_names))*100
    comparison_dict['%_recuperation'].append(recuperation_)

    comparison_dict['non_shared_old_annotations'].append(len(difference_))

    diff_ann_old = list(difference_)
    if diff_ann_old:
        diff_old_str = str(diff_ann_old).strip('[').strip(']').replace("'","")
        comparison_dict['old_annotations_not_captured'].append(diff_old_str)
    else:
        comparison_dict['old_annotations_not_captured'].append(None)

    diff_ann_new = list(new_names.difference(org_names))
    comparison_dict['non_shared_new_annotations'].append(len(diff_ann_new))
    if diff_ann_new:
        diff_new_str = str(diff_ann_new).strip('[').strip(']').replace("'","")
        comparison_dict['new_annotations'].append(diff_new_str)
    else:
        comparison_dict['new_annotations'].append(None)

    return comparison_dict

def confusion_matrix_subs(original_heh_df: pd.DataFrame, created_heh_df: pd.DataFrame, substance_id_list: Optional[np.ndarray] = None,
                            endpoint: Optional[str] = None) -> dict:
    """
        function that compares HEH dataframes annotations and created a dict with the results of the comparison

        :param original_heh_df: HEH dataframe from Inditex Excel/CII database
        :param created_heh_df: HEH dataframe obtained from implementing manually USC workflow
        :param substance_id_list: list with substance id's
        :param endpoint: option to obtain a comparison dict by specific enpodint. Endpoints available are:
                        - CMR
                        - PBT
                        - vPvB
                        - Endocrine_disruptor
                        - Sensitiser 
                        - Other

        :returns comparison_dict: dictionary containing comparison informaton. To be converted into pandas DataFrame
    """

    if substance_id_list:
        conf_mat_dict = {'subs_id':[],'Anotacion_en_inventario':[], 'Anotacion_generada':[], 'Anotacion_no_presente_en_inventario':[],
                        'Anotacion_no_generada':[],'Anotaciones_compartidas':[]}

        for subs_id in substance_id_list:
            org_names = set(original_heh_df[(original_heh_df['subs_id'] == subs_id) &
                                            (original_heh_df['type'] == endpoint)].name.unique())
            new_names = set(created_heh_df[(created_heh_df['subs_id'] == subs_id) &
                                            (created_heh_df['endpoint_type'] == endpoint)].name.unique())
            
            conf_mat_dict['subs_id'].append(subs_id)
            conf_mat_dict['Anotacion_en_inventario'].append(len(org_names))
            conf_mat_dict['Anotacion_generada'].append(len(new_names))

            negativo_inventario = list(new_names.difference(org_names))
            negativo_pregenerado = list(org_names.difference(new_names))

            conf_mat_dict['Anotacion_no_presente_en_inventario'].append(len(negativo_inventario))
            conf_mat_dict['Anotacion_no_generada'].append(len(negativo_pregenerado))

            intersec_ = org_names.intersection(new_names)
            conf_mat_dict['Anotaciones_compartidas'].append(len(intersec_))

    else:
        conf_mat_dict = {'Endpoint':[],'Anotaciones_en_inventario':[], 'Anotaciones_generadas':[], 'Anotaciones_no_presente_en_inventario':[],
                        'Anotaciones_no_generadas':[],'Anotaciones_compartidas':[]}

        endpoint_list = ['CMR','PBT','vPvB','Endocrine_disruptor','Sensitiser','Other']

        for ep in endpoint_list:
            org_names = set(original_heh_df[original_heh_df['type'] == ep].name.unique())
            new_names = set(created_heh_df[created_heh_df['endpoint_type'] == ep].name.unique())
            
            conf_mat_dict['Endpoint'].append(ep)
            conf_mat_dict['Anotaciones_en_inventario'].append(len(org_names))
            conf_mat_dict['Anotaciones_generadas'].append(len(new_names))

            negativo_inventario = list(new_names.difference(org_names))
            negativo_pregenerado = list(org_names.difference(new_names))

            conf_mat_dict['Anotaciones_no_presente_en_inventario'].append(len(negativo_inventario))
            conf_mat_dict['Anotaciones_no_generadas'].append(len(negativo_pregenerado))

            intersec_ = org_names.intersection(new_names)
            conf_mat_dict['Anotaciones_compartidas'].append(len(intersec_))

    return conf_mat_dict

def confusion_matrix_pos_neg(original_heh_df: pd.DataFrame, created_heh_df: pd.DataFrame, substance_id_list: Optional[np.ndarray] = None,
                        endpoint: Optional[str] = None) -> dict:
    """
        function that compares HEH dataframes annotations and created a dict with the results of the comparison

        :param original_heh_df: HEH dataframe from Inditex Excel/CII database
        :param created_heh_df: HEH dataframe obtained from implementing manually USC workflow
        :param substance_id_list: list with substance id's
        :param endpoint: option to obtain a comparison dict by specific enpodint. Endpoints available are:
                        - CMR
                        - PBT
                        - vPvB
                        - Endocrine_disruptor
                        - Sensitiser 
                        - Other

        :returns comparison_dict: dictionary containing comparison informaton. To be converted into pandas DataFrame
    """
    ann_info_list = ['YES','NO','No information','Pending', 'Pending (1)','Pending (2)',
                    'Pending (3)', 'Pending (3a)','Pending (3b)','Pending (3c)','Pending (4)']

    conf_mat_dict = {'subs_id':[],'Anotacion_en_inventario':[], 'Anotacion_generada':[], 'Anotacion_no_presente_en_inventario':[],
                    'Anotacion_no_generada':[]}

    for subs_id in substance_id_list:
        flag_org = 0
        flag_new = 0
        flag_ni = 0
        flag_np = 0
        org_names = set(original_heh_df[(original_heh_df['subs_id'] == subs_id) &
                                        (original_heh_df['type'] == endpoint)].name.unique())
        new_names = set(created_heh_df[(created_heh_df['subs_id'] == subs_id) &
                                        (created_heh_df['endpoint_type'] == endpoint)].name.unique())
        
        conf_mat_dict['subs_id'].append(subs_id)

        negativo_inventario = list(new_names.difference(org_names))
        negativo_pregenerado = list(org_names.difference(new_names))

        for pending in ann_info_list:
            if flag_org == 1:
                pass
            else:  
                if pending in org_names:
                    conf_mat_dict['Anotacion_en_inventario'].append(pending)
                    flag_org = 1
            if flag_new == 1:
                pass
            else:
                if pending in new_names:
                    conf_mat_dict['Anotacion_generada'].append(pending)
                    flag_new = 1
            if pending in negativo_inventario:
                conf_mat_dict['Anotacion_no_presente_en_inventario'].append(pending)
                flag_ni = 1
            if pending in negativo_pregenerado:
                flag_np = 1
                conf_mat_dict['Anotacion_no_generada'].append(pending)
        
        if flag_ni == 0:
            conf_mat_dict['Anotacion_no_presente_en_inventario'].append(None)
        if flag_np == 0:
            conf_mat_dict['Anotacion_no_generada'].append(None)
    
    return conf_mat_dict