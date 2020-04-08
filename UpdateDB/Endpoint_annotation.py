"""
    Created by: Eric March Vila (eric.march@upf.edu)
    On: 08/04/2020, 12:04 AM
"""

import numpy as np
import pandas as pd

from typing import Union
from UpdateDB.Connect_CII import Connector

class Endpoint(Connector):
    """
        Child class of Connector. Uses Conenctor functions to interact with CII.
        This class aims to check in CII the presence of hazard annotations for each substance and generate 
        a new annotation for a given endpoint related to a certain hazard.
        The theoretical background comes from the USC Workflow to generate the endpoint annotations.

        Example:
        Formaldehyde has Carc. 2 hazard annotation. This means is positive for CMR, then we will have a YES annotation
        for formaldehyde in CMR.
        If no annotation is found, then No information is given.
    """

    def __init__(self, host: str, dbname: str, user: str, password: str):
        """
            Initializes class with main arguments for psycopg2 connection to the database

            :param host:
            :param dbname:
            :param user:
            :param password:
        """

        super().__init__(host, dbname, user, password)
        super().open_connection()
    
    def get_annotations_for_substance(self, substances_id: np.ndarray, endpoint_annotations: dict) -> pd.DataFrame:
        """
            For each substance id in the list, checks regulations table if there are certain annotations,
            which are the values of the dict endpoint_annotations. 

            :param substances_id:
            :param endpoint_annotations: dictionary which keys are endpoints (CMR, PBT...) and values are
                                        the hazard annotations associated to that endpoints

            :return substance_endpoint_annotations:
        """

        substance_endpoint_annotations = pd.DataFrame(index=range(len(substances_id)))
       
        for i, id_ in enumerate(substances_id):
            substance_endpoint_annotations.loc[substance_endpoint_annotations.index == i, 'subs_id'] = id_
            for endpoint in endpoint_annotations.keys():
                annotations = endpoint_annotations[endpoint]
                final_annotation = self.get_annotation_per_endpoint(id_, annotations)
                substance_endpoint_annotations.loc[substance_endpoint_annotations.index == i, endpoint] = final_annotation
        
        return substance_endpoint_annotations

    def get_annotation_per_endpoint(self, subs_id: int, annotations: list) -> str:
        """
            Checks if there are annotations for that substance id

            :param subs_id:
            :param annotations:

            :return final_annotation:
        """

        substance_annotations = self.check_presence_in_table(subs_id, annotations)

        if substance_annotations.empty:
            final_annotation = 'No information'
        else:
            final_annotation = self.check_source_of_annotation(substance_annotations)

        return final_annotation

    def check_presence_in_table(self, subs_id: int, annotations: str) -> pd.DataFrame:
        """
            Ask CII if there are annotations for the input substance

            :param subs_id:
            :param annotations:

            :return substance_annotations:
        """

        query_ = """SELECT reg.id, reg.subs_id, rco.country, rt."type", rg.general_regulation_name, 
                    rspec.specific_regulation_name, rsub.subspecific_regulation_name, 
                    rsc.special_cases_name, addr.additional_information_name, cid."name" as chemical_identifier, 
                    ct."type" as type_of_identifier, regn.names
                    FROM regulations reg
                    LEFT JOIN substance sub ON sub.id = reg.subs_id
                    left join regulation_country rco on rco.id = reg.reg_country_id
                    left join regulation_type rt on rt.id = reg.reg_type_id
                    left join general_regulation rg on rg.id = reg.gen_reg_id
                    left join specific_regulation rspec on rspec.id = reg.spec_reg_id
                    LEFT JOIN subspecific_regulation rsub ON rsub.id = reg.subspec_reg_id
                    left join special_cases_regulation rsc on rsc.id = reg.special_cases_id
                    left join additional_information_regulation addr on addr.id = reg.additional_information_id
                    LEFT JOIN chem_id cid ON cid.chem_type_id = reg.chem_type_id
                    LEFT JOIN chem_type ct ON ct.id = cid.chem_type_id
                    LEFT JOIN regulation_names regn ON regn.id = reg.regulation_id 
                    WHERE reg.subs_id = {} and regn.names in {}""".format(subs_id, tuple(annotations))
        
        substance_annotations = pd.read_sql_query(query_, self.conn)

        return substance_annotations
    
    def check_source_of_annotation(self, substance_annotations: pd.DataFrame) -> str:
        """
            Checks which source the annotation comes from and assign either a YES or a Pending annotation
            for the endpoint of interest.

            :param substance_annotations:

            :return final_annotation
        """
        
        sources = substance_annotations[['general_regulation_name','specific_regulation_name','subspecific_regulation_name',
        'special_cases_name','additional_information_name']].drop_duplicates()
        
        # if 'svhc' in sources['specific_regulation_name'].values:
        #     final_annotation = 'YES'
        # elif 'clp' in sources['general_regulation_name'].values and 'harmonised_c&l' in sources['specific_regulation_name'].values:
        #     final_annotation = 'YES'
        # elif 'pbt_vpvb' in sources['general_regulation_name'].values or 'endocrine_disruptors' in sources['general_regulation_name'].values

        if not sources['general_regulation_name'].isin(['clp', 'pbt_vpvb', 'endocrine_disruptors']).empty:
            final_annotation = 'YES'
        elif not sources['specific_regulation_name'].isin(['svhc', 'harmonised_C&L']).empty:
            final_annotation = 'YES'
        else:
            final_annotation = 'Pending'

        return final_annotation

