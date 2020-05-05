"""
    Created by: Eric March Vila (eric.march@upf.edu)
    On: 14/02/2020, 10:00 AM
"""

import numpy as np
import pandas as pd
import psycopg2

from compoundDB import inputtools as it
from psycopg2.extensions import register_adapter, AsIs, Float
from typing import Union

class Connector():
    """
       Main class to connect to CII and CR database.
    """

    def __init__(self, host: str = None, dbname: str = None, user: str = None, password: str = None):
        """
            Initializes class with main arguments for psycopg2 connection to the database.
            
            :param host:
            :param dbname:
            :param user:
            :param password:
        """

        self.host = host
        self.dbname = dbname
        self.user = user
        self.password = password

        register_adapter(float, self.nan_to_null)

    def open_connection(self) -> psycopg2.extensions.connection:
        """
            Opens a psycopg2 connection to the PostgreSQL DB and returns the connection object.

            :return self.conn:
        """

        self.conn = psycopg2.connect(host=self.host, dbname=self.dbname, user=self.user, password=self.password)
        self.curs = self.conn.cursor()

        return self.conn

    def compoundDB_connection(self) -> psycopg2.extensions.connection:
        """
            Connects to compoundDB (CR)

            :return self.compounddb_conn:
        """

        self.compounddb_conn = it.openconnection(host='gea', password='DBAdmin')

    def nan_to_null(self, f, _NULL: str =AsIs('NULL'), _NaN: np.nan =np.NaN, _Float: float =Float) -> Union[float,str]:
        """
            Special function to handle NaN values from numpy when passed as part of a query in postgres.
            NaN should be interpreted as NULL.

            :param _NULL:
            :param _NaN:
            :param _Float:

            :returns _Float(f)
            :returns _NULL

            BUG: not sure if this is working
        """

        if f is not _NaN:
            return _Float(f)
        return _NULL

    #### Substances block
    #### Functions that interact with the substances tables and extract them as dataframes

    def get_substances(self) -> pd.DataFrame:
        """
            Get subsances from database

            :return substances_df
        """

        substances_df = pd.read_sql_query("""SELECT id, class_name, preferred_name, mol_formula FROM public.substance;""",self.conn)

        return substances_df
        
    def get_substances_with_chemical_identifiers(self) -> pd.DataFrame:
        """
            Get substances with chemical identifiers (CAS, EC, Index) from database

            :return substances_chem_id:
        """
        
        substances_chem_id = pd.read_sql_query("""SELECT sub.id as substance_id, class_name, preferred_name, mol_formula, 
                                                cid."name" as Chemical_identifier, ct."type" as Type_of_identifier
                                                FROM substance sub
                                                left join chem_id cid on cid.subs_id = sub.id
                                                left join chem_type ct on ct.id = cid.chem_type_id""",self.conn)

        return substances_chem_id
    
    def get_substances_with_structure(self) -> pd.DataFrame:
        """
            Get substances with SMILES from database

            :return substance_structures:
        """

        substance_structures = pd.read_sql_query("""SELECT s.class_name, s.preferred_name, cid."name" , struc."structure"
                                                FROM substance s
                                                left join substance_structure struc on struc.subs_id = s.id
                                                left join chem_id cid on cid.subs_id  = s.id
                                                where cid.chem_type_id = 1
                                                order by s.id ASC""", self.conn)

        return substance_structures
    
    def get_endpoints(self) -> pd.DataFrame:
        """
            Gets all the content from endpoint_annotation table

            :return endpoint_annotation:
        """

        endpoint_annotation = pd.read_sql_query("""SELECT *
                                                FROM endpoint_annotation  
                                                order by id ASC""", self.conn)
        
        return endpoint_annotation
        
    def get_substances_with_endpoint_annotations_and_structure(self) -> pd.DataFrame:
        """
            Get substances with SMILES and endpoint annotations. The aim is to generate an sdf from 
            the resulting dataframe

            :return sub_ann_struc:
        """

        sub_ann_struc = pd.read_sql_query("""SELECT class_name, preferred_name, mol_formula, 
                                            str."structure", ep.cmr, ep.pbt, ep.vpvb, 
                                            ep.sensitiser, ep.endocrine_disruptor
                                            FROM substance sub
                                            left join substance_structure str on str.subs_id = sub.id
                                            left join endpoint_annotation ep on ep.subs_id = sub.id
                                            where str."structure" is not null  
                                            order by sub.id ASC""", self.conn)
        
        return sub_ann_struc

    #### Chemical Identifier block (CAS/EC/Index)
    #### Functions that interact with chem id tables and extract them as dataframes

    def get_chem_id_dataframe(self) -> pd.DataFrame:
        """
            Get chem id dataframe from database

            :return chem_id_df
        """
        
        chem_id_df = pd.read_sql_query("""SELECT id, "name", chem_type_id, subs_id FROM public.chem_id;""",self.conn)

        return chem_id_df

    def get_chem_id_type_dataframe(self) -> pd.DataFrame:
        """
            Get chemical identifier type from database

            :return chem_type_df:
        """
        
        chem_type_df = pd.read_sql_query("""SELECT id, "type" FROM public.chem_type;""",self.conn)

        return chem_type_df

    #### Regulations block
    #### Functions that interact with the regulations tables and extract them as dataframes

    def get_regulation_country(self) -> pd.DataFrame:
        """
            Get the regulation country from the database

            :return regulation_country:
        """

        regulation_country = pd.read_sql_query("""SELECT id, country
                                                FROM public.regulation_country;""", self.conn)
        
        return regulation_country

    def get_regulation_type(self) -> pd.DataFrame:
        """
            Get the regulation type from the database

            :return regulation_type:
        """

        regulation_type = pd.read_sql_query("""SELECT id, "type"
                                            FROM public.regulation_type;""", self.conn)

        return regulation_type

    def get_general_regulations(self) -> pd.DataFrame:
        """
            Get the general regulations from the database

            :return general_regulations:
        """

        general_regulations = pd.read_sql_query("""SELECT id, general_regulation_name
                                                FROM public.general_regulation;""", self.conn)
        
        return general_regulations
    
    def get_specific_regulations(self) -> pd.DataFrame:
        """
            Get the specific regulations from the database

            :return specific_regulations:
        """

        specific_regulations = pd.read_sql_query("""SELECT id, specific_regulation_name
                                                FROM public.specific_regulation;""", self.conn)
        
        return specific_regulations

    def get_subspecific_regulations(self) -> pd.DataFrame:
        """
            Get subspecific regulations from the database

            :return subspecific_regulations:
        """

        subspecific_regulations = pd.read_sql_query("""SELECT id, subspecific_regulation_name
                                                    FROM public.subspecific_regulation;""",self.conn)
        
        return subspecific_regulations

    def get_special_cases(self) -> pd.DataFrame:
        """
            Get special cases names from the database

            :return special_cases:
        """

        special_cases = pd.read_sql_query("""SELECT id, special_cases_name
                                            FROM public.special_cases_regulation;""",self.conn)
        
        return special_cases

    def get_additional_information(self) -> pd.DataFrame:
        """
            Get additional information names from the database

            :return additional_information:
        """

        additional_information = pd.read_sql_query("""SELECT id, additional_information_name
                                                    FROM public.additional_information_regulation;""",self.conn)
        
        return additional_information
    
    def get_regulation_names(self) -> pd.DataFrame:
        """
            Get regulation names from the database

            :return regulation_names:
        """

        regulation_names = pd.read_sql_query("""SELECT id, names
                                            FROM public.regulation_names;""",self.conn)

        return regulation_names
    
    def get_big_regulations_table(self) -> pd.DataFrame:
        """
            Gets the big regulations table wiht only the id's

            :return big_regulations_table:
        """

        big_regulations_table = pd.read_sql_query("""SELECT id, subs_id, reg_country_id, reg_type_id, gen_reg_id, 
                                                    spec_reg_id, subspec_reg_id, special_cases_id, additional_information_id, 
                                                    chem_id_name, chem_type_id, regulation_id FROM public.regulations;""",self.conn)
        
        return big_regulations_table
        
    def get_regulations_per_substance(self) -> pd.DataFrame:
        """
            Get regulations for each substance from the database

            :return regulations_per_substance:
        """

        regulations_per_substance = pd.read_sql_query("""SELECT reg.id, reg.subs_id, sub.class_name, sub.preferred_name, 
                                                    rco.country, rt."type", rg.general_regulation_name, 
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
                                                    LEFT JOIN regulation_names regn ON regn.id = reg.regulation_id""",self.conn)

        return regulations_per_substance
    
    def get_regulations_with_id(self) -> pd.DataFrame:
        """
            Get regulations with id from the database

            :return regulations_with_id:
        """

        regulations_with_id = pd.read_sql_query("""SELECT DISTINCT rco.country, rco.id as "country_id", rt."type", rt.id as "reg_type_id",
                                                    rg.general_regulation_name, rg.id as "general_reg_id", rspec.specific_regulation_name,
                                                    rspec.id as "spec_reg_id", rsub.subspecific_regulation_name, rsub.id as "subspec_reg_id",
                                                    rsc.special_cases_name, rsc.id as "special_cases_id"
                                                    FROM regulations reg
                                                    left join regulation_country rco on rco.id = reg.reg_country_id
                                                    left join regulation_type rt on rt.id = reg.reg_type_id
                                                    left join general_regulation rg on rg.id = reg.gen_reg_id
                                                    left join specific_regulation rspec on rspec.id = reg.spec_reg_id
                                                    LEFT JOIN subspecific_regulation rsub ON rsub.id = reg.subspec_reg_id
                                                    left join special_cases_regulation rsc on rsc.id = reg.special_cases_id""",self.conn)

        return regulations_with_id
    
    # List of regulation dataframes to iterate over

    def get_dict_of_regulation_dataframes(self) -> dict:
        """
            Puts in a list the main regulation dataframe so one can iterate over them

            :return self.df_dict:
        """

        gen_reg = self.get_general_regulations()
        spec_reg = self.get_specific_regulations()
        subspec_reg = self.get_subspecific_regulations()
        special_reg = self.get_special_cases()

        self.df_dict = {'general_regulation':gen_reg, 'specific_regulation':spec_reg, 
                        'subspecific_regulation':subspec_reg, 'special_cases_names':special_reg}
        
        return self.df_dict