"""
    Created by: Eric March Vila (eric.march@upf.edu)
    On: 14/02/2020, 10:00 AM
"""

import numpy as np
import pandas as pd
import pubchempy as pcp

from CreateDB.functions_to_process import get_cas, get_ec, get_index_number
from phitools import moleculeHelper as mh
from typing import Union, Optional
from UpdateDB.Connect_CII import Connector

class UpdateDB(Connector):
    """
        Child class of Connector. Allows to update the database with new information.

        Main functions iterate through the input dataframe
        and updates CII records depending on what we want to add:
            - all the information (done)
            - only the substances (TODO)
            - only chemical identifiers (TODO)
            - only structures (TODO)
            - only sources (TODO)
            - only annotations (TODO)
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
        super().compoundDB_connection()

    #### Main functions

    def add_all_information_from_dataframe(self, dataframe: pd.DataFrame, preferred_name_field: str, chem_id_field: str, 
                                    chem_id_type: str, sourceName_field: str, regulation_field: str, class_name_field: Optional[str] = None):
        """
            Extracts the information from a Pandas Dataframe and adds it into the database.

            :param dataframe:
            :param class_name_field:
            :param preferred_name_field:
            :param chem_id_field:
            :param chem_id_type:
            :param sourceName_field:
            :param regulation_field:
        """

        for idx, row in dataframe.iterrows():

            sourceName = row[sourceName_field]
            if 'Inditex' in str(sourceName):
                continue

            # Substance names and CAS/EC/Index numbers processing

            if class_name_field:
                class_name = row[class_name_field]
            else:
                class_name = None

            preferred_name = row[preferred_name_field]

            chem_identifier = self.process_string_or_list_from_pandas(row[chem_id_field])
            
            if self.match_chemical_identifier_pattern(chem_identifier, chem_id_type):
                chem_id, substance_id = self.add_substance_chemical_identifier(class_name,preferred_name,chem_identifier, chem_id_type)
            else:
                chem_id = substance_id = None

            if not substance_id:
                substance_id = self.add_substance(class_name,preferred_name)

            if substance_id and chem_id:
                self.add_structure(substance_id, chem_id, chem_identifier)

            # Source Name processing

            if isinstance(sourceName,float):
                sourceID = None
            else:
                sourceID = self.add_regulation(sourceName)
            
            # Annotation processing

            try:
                annotation = self.process_string_or_list_from_pandas(row[regulation_field])
            except AttributeError:
                annotation = None

            if isinstance(annotation, list):
                for ann in annotation:
                    annotationID = self.add_annotation(ann)
                    # Big regulations table updating
                    big_reg_id = self.add_regulations_and_annotaitons(substance_id, sourceID, annotationID)
            elif isinstance(annotation, str):
                annotationID = self.add_annotation(annotation)
                # Big regulations table updating
                big_reg_id = self.add_regulations_and_annotaitons(substance_id, sourceID, annotationID)
    
    def add_substances_from_dataframe(self, dataframe: pd.DataFrame, preferred_name_field: str, class_name_field: Optional[str] = None):
        """
            Adds new substances from dataframe

            :param dataframe:
            :param preferred_name_field:
            :param class_name_field:
        """

        for idx, row in dataframe.iterrows():

            if class_name_field:
                class_name = row[class_name_field]
            else:
                class_name = None

            preferred_name = row[preferred_name_field]

            self.add_substance(class_name,preferred_name)
    
    def add_chemical_identifier_from_dataframe(self, dataframe: pd.DataFrame, chem_id_field: str, chem_id_type: str, preferred_name_field: str,
                                                class_name_field: Optional[str] = None):
        """
            Adds new chemical identifiers (CAS/EC/Index) from dataframe

            :param dataframe:
            :param chem_id_field:
            :param chem_id_type:
        """

        for idx, row in dataframe.iterrows():

            # Substance names and CAS/EC/Index numbers processing

            if class_name_field:
                class_name = row[class_name_field]
            else:
                class_name = None

            preferred_name = row[preferred_name_field]

            chem_identifier = self.process_string_or_list_from_pandas(row[chem_id_field])
            
            if self.match_chemical_identifier_pattern(chem_identifier, chem_id_type):
                chem_id, substance_id = self.add_substance_chemical_identifier(class_name,preferred_name,chem_identifier, chem_id_type)
            else:
                chem_id = substance_id = None

            if not substance_id:
                substance_id = self.add_substance(class_name,preferred_name)

    #### Input string processing

    def process_string_or_list_from_pandas(self, str_: str) -> Union[str,list]:
        """
            Process fields from dataframe and returns either a string or a list

            :param str_: content from field in the dataframe
            
            :return str_proc: is either a string or a list containing the processed field related to the substance
        """

        preproc_str = str_.split(',')

        if len(preproc_str) > 1:
            str_proc = [str_element.replace("'","").strip() for str_element in preproc_str]
            if '-' in str_proc:
                str_proc.remove('-')
        else:
            str_proc = str(preproc_str).strip('[').strip(']').replace("'","").strip()

        return str_proc

    #### CAS/EC/Index table

    def add_substance_chemical_identifier(self, class_name:str, preferred_name:str, chem_id_to_add:str, chem_id_type:str) -> int:
        """
            Adds new CAS/EC/Index number in the database.
            If it's already there, just gets chemical identifier id.

            :param chem_id_to_add:
            :param chem_id_type:

            :return chem_identifier_id:
        """

        chemid_cmd = "SELECT id FROM chem_id WHERE name = '{}';".format(chem_id_to_add)
        chem_identifier_id = self.check_presence_or_absence(chemid_cmd)
        
        if not chem_identifier_id:
            substanceID = self.add_substance(class_name, preferred_name)
            max_id_cmd = """SELECT max(id) from chem_id"""
            insert_cmd = """INSERT INTO public.chem_id (id, name, chem_type_id, subs_id) VALUES({}, '{}', {}, {});"""
            chem_type_id = self.get_chem_type_id(chem_id_type)
            chem_identifier_id = self.insert_in_database(max_id_cmd, insert_cmd, chem_id_to_add, chem_type_id, substanceID)
        else:
            sub_cmd = "SELECT subs_id FROM chem_id WHERE name = '{}';".format(chem_id_to_add)
            substanceID = self.check_presence_or_absence(sub_cmd)
        
        return chem_identifier_id, substanceID
    
    def match_chemical_identifier_pattern(self, chemical_identifier: Union[str,float], chem_id_type: str) -> Optional[str]:
        """
            Uses the functions in functions_to_process.py to check whether the chem id is CAS, EC or Index given the chem_id_type
            Returns it if True, otherwise False

            :param chemical_identifier:
            :param chem_id_type:

            :return chemical_identifier_match:
        """

        if 'cas' in chem_id_type.lower():
            chemical_identifier_match = get_cas(chemical_identifier)
        elif 'ec' in chem_id_type.lower():
            chemical_identifier_match = get_ec(chemical_identifier)
        elif 'index' in chem_id_type.lower():
            chemical_identifier_match = get_index_number(chemical_identifier)
        
        return chemical_identifier_match

    def get_chem_type_id(self, chem_type:str) -> int:
        """
            Gets the chem type id using the query to ask the database

            :param chem_type:

            :return chem_type_id:
        """

        cmd = "SELECT id FROM chem_type where type like '%{}%'".format(chem_type)
        chem_type_id = self.check_presence_or_absence(cmd)

        return chem_type_id

    #### Substance table

    def add_substance(self, class_name: str, preferred_name: str) -> int:
        """
            Adds a new substance to the database if it's not present.
            If it's already there, just retrieves the associated substance id.

            :param class_name:
            :param preferred_name:

            :return substanceID:
        """

        cmd = "SELECT id FROM substance WHERE preferred_name = '{}';".format(preferred_name)
        substanceID = self.check_presence_or_absence(cmd)

        if not substanceID:
            # check on substance table using name
            max_id_cmd = """select max(id) from substance"""
            insert_cmd = """INSERT INTO public.substance (id, class_name, preferred_name, mol_formula)
                    VALUES ('{}','{}','{}','{}')"""
            substanceID = self.insert_in_database(max_id_cmd, insert_cmd, class_name, preferred_name, None)

        return substanceID

    #### Structure table

    def add_structure(self, subs_id: int, chem_id: int, chemical_identifier: str):
        """
            Checks if substance has structure in database. If not, tries to calculate SMILES from chem_id
            and adds the susbtance if found

            :param subs_id:
            :param chem_id: chemical identifier id in database
            :param chemical_identifier: CAS number
        """

        cmd = """SELECT id FROM substance_structure 
                WHERE subs_id = {} AND chem_id = {}""".format(subs_id, chem_id)
        structureID = self.check_presence_or_absence(cmd)
        
        if not structureID:
            structure = self.calculate_structure(chemical_identifier)
            if structure:
                max_id_cmd = """SELECT max(id) FROM substance_structure;"""
                insert_cmd = """INSERT INTO public.substance_structure (id, subs_id, chem_type_id, chem_id, structure)
                                VALUES({}, {}, {}, {}, '{}');"""
                structureID = self.insert_in_database(max_id_cmd, insert_cmd, subs_id, 1, chem_id, structure)

    def calculate_structure(self, chemical_identifier:str) -> str:
        """
            Gets SMILES from CAS number looking in three different sources:
                - Standardiser
                - CACTUS server
                - PubChem
            
            :param chemical_identifier:

            :return strucure:
        """

        structure = mh.resolveCAS(chemical_identifier, conn=self.compounddb_conn)

        if not structure:
            struc_list = []
            try:
                results = pcp.get_compounds(chemical_identifier, 'name')
            except:
                results = None
            if results:
                for compound in results:
                    smiles = compound.canonical_smiles
                    struc_list.append(smiles)
                if len(struc_list) > 1:
                    struc_list = list(set(struc_list))
                structure = struc_list[0]
            else:
                structure = None

        return structure

    #### Regulation tables

    def add_regulation(self, sourceName: str) -> int:
        """
            Checks if the regulation is already in the database and if not, it adds it.

            :param sourceName:

            :return sourceID:
        """

        sourceID = self.check_presence_or_absence_of_regulation(sourceName)

        if not sourceID:
            format_sourcename = sourceName.lower()
            max_id_cmd = """select max(id) from general_regulation"""
            insert_cmd = """INSERT INTO public.general_regulation (id, general_regulation_name)
                            VALUES ('{}','{}')"""
            sourceID_number = self.insert_in_database(max_id_cmd, insert_cmd, format_sourcename)
            sourceID = [sourceID_number,'general_regulation']

        elif 'annex' in sourceName.lower():
            reach_, annex = sourceName.lower().split(' ',1)
            annex_format = '_'.join(annex.split())
            annex_source_id = self.check_presence_or_absence_of_regulation(annex_format)

            if not annex_source_id:
                max_id_cmd = """select max(id) from specific_regulation"""
                insert_cmd = """INSERT INTO public.specific_regulation
                            (id, specific_regulation_name) VALUES('{}', '{}')"""
                annex_id = self.insert_in_database(max_id_cmd, insert_cmd, annex_format)
                annex_source_id = (annex_id, 'specific_regulation')

            sourceID.extend(annex_source_id)
        
        elif 'registration' in sourceName.lower():
            # Registration in REACH refers to registration dossiers, so we change the name of this
            # regulation to registration_dossier, which is actually in the database
            reach_, reg = sourceName.lower().split()
            reg_dos = '_'.join([reg,'dossier'])
            reg_source_id = self.check_presence_or_absence_of_regulation(reg_dos)
            sourceID.extend(reg_source_id)
        
        return sourceID
    
    #### Annotation (regulation names) table

    def add_annotation(self, annotation: str) -> int:
        """
            Check if annotation is in database and if not, adds it

            :param annotation:

            :return annotationID:
        """
        
        cmd = "SELECT id FROM regulation_names WHERE names = '{}';".format(annotation)
        annotationID = self.check_presence_or_absence(cmd)

        if not annotationID:
            max_id_cmd = """SELECT max(id) FROM regulation_names;"""
            insert_cmd = """INSERT INTO public.regulation_names (id, names) VALUES({}, '{}');"""
            annotationID = self.insert_in_database(max_id_cmd, insert_cmd, annotation)

        return annotationID

    #### Big table with regulations and annotations

    def add_regulations_and_annotaitons(self, subs_id: int, source_id: list, ann_id: int) -> int:
        """
            Adds annotation with regulations in big table if it's not present

            :param subs_id:
            :param source_id:
            :param ann_id:

            I put 1 in regulation country id which correspond to europe because 
            most annotations come from EU

            :return reg_id:
        """
        
        sources_dict_query = self.create_sources_dict(source_id)
       
        new_cmd = """SELECT id FROM public.regulations where subs_id = {} {} and regulation_id = {}""".format(subs_id,
                                                                                            ' '.join(sources_dict_query['check_query']), ann_id)
        
        reg_id = self.check_presence_or_absence(new_cmd)

        if not reg_id:
            max_id_cmd = """SELECT max(id) FROM regulations;"""
            #### values_to_add: I made it like this to overcome an error of tuple index out of range in insert_in_database function, since
            #### sources_dict_query['id'] is a list and when it had more than one element, the function wasn't capturing it well
            values_to_add = [subs_id, 1]
            values_to_add.extend(sources_dict_query['id'])
            values_to_add.append(ann_id)

            values_str = ''.join(["VALUES ({},",','.join(["{}".format(str(element)) for element in values_to_add]),");"])
            insert_cmd = """INSERT INTO public.regulations (id, subs_id, reg_country_id, 
                        {} regulation_id) {}""".format(' '.join(sources_dict_query['insert_query']),values_str)
           
            reg_id = self.insert_in_database(max_id_cmd, insert_cmd, values_to_add)

        return reg_id
    
    def create_sources_dict(self, sources:list) -> dict:
        """
            Creates a dictionary with different formatted strings for the queries to CII db.

            :param sources: list of sources obtained from dataframe

            :return sources_dict:
        """
        
        sources_dict = {'check_query':[],'insert_query':[],'id':[]}
        
        for i, source in enumerate(sources):
            if 'general_regulation' in source[1]:
                gen_reg_id = source[0]
                check_query = 'and gen_reg_id = {}'.format(gen_reg_id)
                insert_query = 'gen_reg_id,'
                
                sources_dict['check_query'].append(check_query)
                sources_dict['insert_query'].append(insert_query)
                sources_dict['id'].append(gen_reg_id)

            elif 'specific_regulation' in source[1]:
                spec_reg_id = source[0]
                check_query = 'and spec_reg_id = {}'.format(spec_reg_id)
                insert_query = 'spec_reg_id,'
                
                sources_dict['check_query'].append(check_query)
                sources_dict['insert_query'].append(insert_query)
                sources_dict['id'].append(spec_reg_id)

            elif 'subspecific_regulation' in source[1]:
                subspec_reg_id = source[0]
                check_query = 'and subspec_reg_id = {}'.format(subspec_reg_id)
                insert_query = 'subspec_reg_id,'
                
                sources_dict['check_query'].append(check_query)
                sources_dict['insert_query'].append(insert_query)
                sources_dict['id'].append(subspec_reg_id)
    
            elif 'special_cases_names' in source[1]:
                special_cases_id = source[0]
                check_query = 'and special_cases_id = {}'.format(special_cases_id)
                insert_query = 'special_cases_id,'
                
                sources_dict['check_query'].append(check_query)
                sources_dict['insert_query'].append(insert_query)
                sources_dict['id'].append(special_cases_id)
        
        return sources_dict

    #### Check presence in database. If ID is returned, then is used to update DB for the new compound.
    #### If not, a new one is generated while new regulation is added

    def check_presence_or_absence(self, query: str) -> int:
        """
            Checks if the desired input is already in DB and returns its corresponding id

            :param query: query to check

            :return id_:
        """

        self.curs.execute(query)
        try:
            id_ = self.curs.fetchone()[0]
            self.conn.commit()
        except TypeError:
            id_ = None
        
        return id_   
        
    def check_presence_or_absence_of_regulation(self, regulation_to_add: str) -> list:
        """
            Compares the regulation from the dataframe with the ones in the DB.
            Returns a list if there are more than one ID's for it or just a string
            it it only finds one

            :param regulation_to_add:
            
            :return source_id_list:

            TODO: change the way it's done. Instead of using dataframe, try to implement direct database ID extraction.
            The problem lies in the substring check when I have clp notification and I want it to be clp in one side and notification 
            on the other, thus having a general regulation clp and special cases regulation notification and not a new general regulation
            clp notification
        """
        
        source_id_list = []
        format_reg = regulation_to_add.lower()

        regulations_dict = self.get_dict_of_regulation_dataframes()
        
        for reg_type, dataframe in regulations_dict.items():
            for idx, row in dataframe.iterrows():
                for element in row:
                    if str(element) in format_reg:
                        source_id_list.append((row['id'],reg_type))

        # gen_reg_cmd = "SELECT id FROM general_regulation WHERE general_regulation_name in ('{}')".format(format_reg)
        # spec_reg_cmd = "SELECT id FROM specific_regulation WHERE specific_regulation_name in ('{}')".format(format_reg)
        # subspec_reg_cmd = "SELECT id FROM subspecific_regulation WHERE subspecific_regulation_name in ('{}')".format(format_reg)
        # special_cases_cmd = "SELECT id FROM special_cases_regulation WHERE special_cases_name in ('{}')".format(format_reg)

        # gen_reg_id = self.check_presence_or_absence(gen_reg_cmd)
        # spec_reg_id = self.check_presence_or_absence(spec_reg_cmd)
        # subspec_reg_id = self.check_presence_or_absence(subspec_reg_cmd)
        # special_cases_id = self.check_presence_or_absence(special_cases_cmd)

        # if gen_reg_id:
        #     source_id_list.append((gen_reg_id,'general_regulation'))
        # if spec_reg_id:
        #     source_id_list.append((spec_reg_id,'specific_regulation'))
        # if subspec_reg_id:
        #     source_id_list.append((subspec_reg_id,'subspecific_regulation'))
        # if special_cases_id:
        #     source_id_list.append((special_cases_id,'special_cases_names'))
        
        return source_id_list

    #### Insert new data into the database

    def insert_in_database(self, max_db_cmd: str, insert_cmd: str, *args: tuple) -> int:
        """
            Gets max ID from selected table in query and adds new values with new ID into that table
            
            :param max_db_cmd: gets max id from database table
            :param insert_cmd: insert into db the new register
            :param *args: arguments to pass to the query as formated string

            :return new_id: new id generated from the query
        """
        
        self.curs.execute(max_db_cmd)
        ID_number = self.curs.fetchone()[0] + 1
        self.conn.commit()
        
        self.curs.execute(insert_cmd.format(ID_number, *args))
        self.conn.commit()
        
        return ID_number
    
    #### Delete data from the database

    def delete_records_larger_than_id(self, checkpoint_df: pd.DataFrame = None, table: str = None, id_: int = None):
        """
            Deletes the records from tables that are larger than the ID provided. Can be executed either for 
            a whole dataframe with table names and ids or for a single table with a selected id

            :param checkpoint_df:
            :param table:
            :param id_:
        """

        delete_cmd = """DELETE FROM {0} WHERE id > {1}"""

        if isinstance(checkpoint_df, pd.DataFrame):
            for idx, row in checkpoint_df.iterrows():
                tab = row[0]
                max_id = row[1]
                del_cmd_format = delete_cmd.format(tab,max_id)
                self.curs.execute(del_cmd_format)
                self.conn.commit()
        elif table and id_:
            del_cmd_format = delete_cmd.format(table, id_)
            self.curs.execute(del_cmd_format)
            self.conn.commit()