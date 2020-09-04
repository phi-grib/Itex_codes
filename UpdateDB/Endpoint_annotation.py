"""
    Created by: Eric March Vila (eric.march@upf.edu)
    On: 08/04/2020, 12:04 AM
"""

import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import PandasTools
from typing import Union, Optional
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

    def __init__(self, host: str = None, dbname: str = None, user: str = None, password: str = None):
        """
            Initializes class with main arguments for psycopg2 connection to the database

            :param host:
            :param dbname:
            :param user:
            :param password:
        """

        if 'cii' in dbname.lower() or 'inventory' in dbname.lower():
            super().__init__(host, dbname, user, password)
            super().open_connection()
            self.db_tag = 'cii'
        elif dbname.lower() in ['cr']:
            super().compoundDB_connection()
            self.db_tag = 'cr'
    
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
                final_annotation = self.get_annotation_per_endpoint(id_, endpoint, annotations)
                substance_endpoint_annotations.loc[substance_endpoint_annotations.index == i, endpoint] = final_annotation
        
        return substance_endpoint_annotations

    def get_annotation_per_endpoint(self, subs_id: int, endpoint: str, annotations: list) -> str:
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
            if self.db_tag == 'cii':
                final_annotation = self.check_source_of_annotation(endpoint, substance_annotations)
            elif self.db_tag == 'cr':
                final_annotation = self.check_cr_source(endpoint, substance_annotations)

        return final_annotation

    def get_total_annotations_per_endpoint(self, substance_endpoint_annotations: pd.DataFrame) -> pd.DataFrame:
        """
            Calculates the total number of annotations per endpoint in the input dataframe

            :param substance_endpoint_annotations:

            :return total_annotations_endpoint:
        """

        endpoints = substance_endpoint_annotations.columns[1:]
        total_annotations_endpoint = pd.DataFrame(index=range(len(endpoints)))

        for i, endpoint in enumerate(endpoints):
            yes_count = len(substance_endpoint_annotations.loc[substance_endpoint_annotations[endpoint] == 'YES',endpoint])
            pen_count = len(substance_endpoint_annotations.loc[substance_endpoint_annotations[endpoint] == 'Pending',endpoint])
            no_count = len(substance_endpoint_annotations.loc[substance_endpoint_annotations[endpoint] == 'No information',endpoint])

            total_annotations_endpoint.loc[total_annotations_endpoint.index == i, 'Endpoints'] = endpoint
            total_annotations_endpoint.loc[total_annotations_endpoint.index == i, 'YES'] = yes_count
            total_annotations_endpoint.loc[total_annotations_endpoint.index == i, 'Pending'] = pen_count
            total_annotations_endpoint.loc[total_annotations_endpoint.index == i, 'No information'] = no_count

        return total_annotations_endpoint

    def check_presence_in_table(self, subs_id: int, annotations: str) -> pd.DataFrame:
        """
            Ask CII/CR if there are annotations for the input substance

            :param subs_id:
            :param annotations:

            :return substance_annotations:
        """

        if self.db_tag == 'cii':
            query_ = """SELECT reg.id, reg.subs_id, rco.country, rt."type", rg.general_regulation_name, 
                        rspec.specific_regulation_name, rsub.subspecific_regulation_name, 
                        rsc.special_cases_name, addr.additional_information_name, reg.chem_id_name as chemical_identifier, 
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
                        LEFT JOIN chem_type ct ON ct.id = reg.chem_type_id
                        LEFT JOIN regulation_names regn ON regn.id = reg.regulation_id 
                        WHERE reg.subs_id = {} and regn.names in {}""".format(subs_id, tuple(annotations))
            
            conn = self.conn
        
        elif self.db_tag == 'cr':
            query_ = """SELECT synonym.name as reg_number, source.name as source_name, 
                        subs_ann.original_annotation, annotation.general, annotation.category
                        FROM substance sub
                        left join synonym on synonym.subsid = sub.id
                        left join source on source.id = sub.sourceid
                        left join subs_ann on subs_ann.subsid = sub.id
                        left join annotation on annotation.id = subs_ann.annid
                        where synonym.name = '{}'
                        and subs_ann.original_annotation in {}""".format(subs_id,tuple(annotations))
            
            conn = self.compounddb_conn
        
        substance_annotations = pd.read_sql_query(query_, conn)
        substance_annotations.drop_duplicates(inplace=True)

        return substance_annotations
    
    def check_cr_source(self, endpoint: str, substance_annotations: pd.DataFrame) -> str:
        """
            Checks which source the annotation comes from in CR database.

            :param substance_annotations:

            :return final_annotation:
        """

        sources = substance_annotations['source_name'].drop_duplicates()

        # List of sources that should be checked. 
        # Decision taken:
        # - REACH Annex VI is considered as yes_source since is the one including Harmonised entries for CLP (Harmonised C&L)
        # - EPA Genetox is considered yes_sources since it's a compendium of 3000 substances analyzed for genetic toxicology, currently outdated
        # - CosmosDB, EFSA OpenFoodTox and Inditex related sources aren't considered since they don't have GHS hazard annotations

        not_considered = ['CosmosDB', 'EFSA OpenFoodTox', 'Inditex',
                          'Inditex MRSL for Parts Suppliers',
                          'Inditex MRSL for Wet Processing Units']

        yes_sources = ['SVHC','REACH Annex VI','EPA Genetox']
        pending_sources = ['CLP Notification', 'REACH Registration','REACH Annex III']

        yes_ann = self.check_yes(sources_df=sources, cr_source=yes_sources)
        if yes_ann:
            final_annotation = yes_ann
        else:
            pending_ann = self.check_pending(sources_df=sources,cr_source=pending_sources)
            final_annotation = pending_ann

        return final_annotation

    def check_source_of_annotation(self, endpoint: str, substance_annotations: pd.DataFrame) -> str:
        """
            Checks which source the annotation comes from and assign either a YES or a Pending annotation
            for the endpoint of interest.

            :param substance_annotations:

            :return final_annotation
        """
        
        pbt = None
        if endpoint in ['PBT', 'vPvB', 'Endocrine Disruptor']:
            pbt = True
        
        sources = substance_annotations[['general_regulation_name','specific_regulation_name','subspecific_regulation_name',
        'special_cases_name','additional_information_name','names']].drop_duplicates()
        
        # We use this lists to check the presence of annotations in these regulations,
        # which are the ones that are used in the USC Workflow. 
        # TODO: check other regulations or add new ones.
        gen_regs = ['clp', 'source1']
        pbt_endoc = ['pbt_vpvb', 'endocrine_disruptors']
        spec_regs = ['svhc', 'harmonised_C&L','annex_vi','source3']
        subspec_regs = ['candidate_list','hazard_class','source2']

        # Registration dossiers and notification. To decide whether to add Pending or not
        reg_dos_not = ['registration_dossier','notification']
        drafts = ['Submitted SVHC intentions','Withdrawn SVHC intentions and submissions','Amendment 2016/1179']

        # These list include terms that indicates no presence of annotation
        # no_presence = ['No Notification','No Registration Dossier','Not included','No information']

        if pbt:
            pbt_endoc_ann = self.check_pbt_vpvb_endoc(sources, pbt_endoc, spec_regs)
            final_annotation = pbt_endoc_ann
        else:
            
            yes_ann = self.check_yes(sources_df=sources, general_regs=gen_regs, specific_regs=spec_regs, subspec_regs=subspec_regs,
            spec_cases=reg_dos_not, drafts=drafts)
            if yes_ann:
                final_annotation = yes_ann
            else:
                pending_ann = self.check_pending(sources_df=sources, spec_cases=reg_dos_not, drafts=drafts)
                final_annotation = pending_ann

        return final_annotation
    
    def check_yes(self, sources_df: pd.DataFrame, cr_source: list = None, general_regs: list = None, 
                    specific_regs: list = None, subspec_regs: list = None, spec_cases: list = None, drafts: list = None) -> Optional[str]:
        """
            Checks regulations to get a YES

            :param sources_df:
            :param general_regs: general regulations list
            :param specific_regs: specific regulations list
            :param no_presence: annotations that indicate no presence of hazard annotation

            :return annotation:
        """

        if self.db_tag == 'cii':
            # yes_df = sources_df[((sources_df['general_regulation_name'].isin(general_regs)) &
            #                     (sources_df['subspecific_regulation_name'].isin(subspec_regs)) |
            #                     (sources_df['specific_regulation_name'].isin(specific_regs)) &
            #                     (sources_df['subspecific_regulation_name'].isin(subspec_regs)))]
            # yes_df = sources_df[((sources_df['general_regulation_name'].isin(general_regs)) &
            #                     (sources_df['subspecific_regulation_name'].isin(subspec_regs)) |
            #                     (sources_df['specific_regulation_name'].isin(specific_regs)) &
            #                     (sources_df['subspecific_regulation_name'].isin(subspec_regs)) |
            #                     (sources_df['general_regulation_name'].isin(general_regs)) | 
            #                     (sources_df['subspecific_regulation_name'].isin(subspec_regs)) |
            #                     (sources_df['specific_regulation_name'].isin(specific_regs)))]
            # yes_df = sources_df[((sources_df['general_regulation_name'].isin(general_regs)) |
            #                     (sources_df['specific_regulation_name'].isin(specific_regs)) |
            #                     (sources_df['subspecific_regulation_name'].isin(subspec_regs)))]
            yes_df = sources_df[(((sources_df['general_regulation_name'].isin(general_regs)) |
                    (sources_df['specific_regulation_name'].isin(specific_regs)) |
                    (sources_df['subspecific_regulation_name'].isin(subspec_regs))) &
                    ~(sources_df['special_cases_name'].isin(spec_cases)) &
                    ~(sources_df['additional_information_name'].isin(drafts)))]

        elif self.db_tag == 'cr':
            yes_df = sources_df.isin(cr_source)
        
        if not yes_df.empty:
            annotation = 'YES'
        else:
            annotation = None
        
        return annotation

    def check_pbt_vpvb_endoc(self, sources_df: pd.DataFrame, pbt_vpvb_endoc_regs: list, spec_regs:list) -> str:
        """
            Checks regulations to get YES for PBT, vPvB and Endocrine Disruptors, which have specific regulations containing its
            hazard annotations

            :param sources_df:
            :param pbt_vpvb_endoc_regs:
            :param spec_regs:

            :return annotation:
        """

        pbt_vpvb_endoc_df = sources_df[((sources_df['general_regulation_name'].isin(pbt_vpvb_endoc_regs)) |
                                        (sources_df['specific_regulation_name'].isin(spec_regs)))]
        
        if not pbt_vpvb_endoc_df.empty:
            annotation = 'YES'
        else:
            annotation = 'No information'
        
        return annotation

    def check_pending(self, sources_df: pd.DataFrame, cr_source: list = None, 
                        spec_cases: list = None, drafts: list = None) -> str:
        """
            Checks regulations to get Pending if YES hasn't been annotated previously.

            :param sources_df:
            :param spec_cases: general regulations list
            :param drafts: specific regulations list

            :return annotation:
        """

        if self.db_tag == 'cii':
            pending_df = sources_df[(sources_df['special_cases_name'].isin(spec_cases)) |
                                (sources_df['additional_information_name'].isin(drafts))]
        elif self.db_tag == 'cr':
            pending_df = sources_df.isin(cr_source)
        
        if not pending_df.empty:
            annotation = 'Pending'
        else:
            annotation = 'No information'
        
        return annotation
    
    def add_annotations_to_database(self, substance_endpoint_annotations:pd.DataFrame):
        """
            Adds dataframe information to CII database.

            :param substance_endpoint_annotations
        """

        for i, row in substance_endpoint_annotations.iterrows():
            idx = i
            subs_id = int(row['subs_id'])
            cmr = row['CMR']
            pbt = row['PBT']
            vpvb = row['vPvB']
            endoc = row['Endocrine Disruptor']
            sens = row['Sensitiser']
            query = """INSERT INTO public.endpoint_annotation (id, subs_id, cmr, pbt, vpvb, endocrine_disruptor, sensitiser)
                        VALUES ({},{},'{}','{}','{}','{}','{}')""".format(idx,subs_id,cmr,pbt,vpvb,endoc,sens)
            self.curs.execute(query)
            self.conn.commit()
    
    def get_sdf_files_from_endpoint_annotation(self):
        """
            Generates sdf files, both for train and test, from endpoint annotations table.

            TODO: in process
        """

        # First we get the substances with annotations and structures and set all the molecule names to one column
        # by adding the class name into the preferred name column where there are NaN

        sub_df = self.get_substances_with_endpoint_annotations_and_structure()
        sub_df.loc[sub_df['preferred_name'].isna(),'preferred_name'] = sub_df.loc[sub_df['preferred_name'].isna(),'class_name']
        sub_df.drop('class_name', axis=1, inplace=True)
        sub_df.rename(columns={'preferred_name':'name'},inplace = True)

        # Second, we generate the mol object from RDKit and then remove all the substances which
        # mol object hasn't been generated. Finally we add the Hydrogens

        PandasTools.AddMoleculeColumnToFrame(sub_df,'structure')
        no_mol = sub_df[sub_df['ROMol'].isna()]
        sub_df.drop(no_mol.index, axis=0, inplace=True)
        sub_df['ROMol'] = [Chem.AddHs(x) for x in sub_df['ROMol'].values.tolist()]

        cmr_df = sub_df[['name','mol_formula','structure','cmr','ROMol']]
        pbt_df = sub_df[['name','mol_formula','structure','pbt','ROMol']]
        vpvb_df = sub_df[['name','mol_formula','structure','vpvb','ROMol']]
        sens_df = sub_df[['name','mol_formula','structure','sensitiser','ROMol']]
        endoc_df = sub_df[['name','mol_formula','structure','endocrine_disruptor','ROMol']]

        cmr_df.loc[cmr_df['cmr'].isin(['YES','Pending']), 'activity'] = 1
        cmr_df.loc[~cmr_df['cmr'].isin(['YES','Pending']), 'activity'] = 0
        pbt_df.loc[pbt_df['pbt'].isin(['YES','Pending']), 'activity'] = 1
        pbt_df.loc[~pbt_df['pbt'].isin(['YES','Pending']), 'activity'] = 0
        vpvb_df.loc[vpvb_df['vpvb'].isin(['YES','Pending']), 'activity'] = 1
        vpvb_df.loc[~vpvb_df['vpvb'].isin(['YES','Pending']), 'activity'] = 0
        sens_df.loc[sens_df['sensitiser'].isin(['YES','Pending']), 'activity'] = 1
        sens_df.loc[~sens_df['sensitiser'].isin(['YES','Pending']), 'activity'] = 0
        endoc_df.loc[endoc_df['endocrine_disruptor'].isin(['YES','Pending']), 'activity'] = 1
        endoc_df.loc[~endoc_df['endocrine_disruptor'].isin(['YES','Pending']), 'activity'] = 0

        cmr_pos = cmr_df[cmr_df['cmr'].isin(['YES','Pending'])]
        cmr_neg = cmr_df[~cmr_df['cmr'].isin(['YES','Pending'])]
        pbt_pos = pbt_df[pbt_df['pbt'].isin(['YES','Pending'])]
        pbt_neg = pbt_df[~pbt_df['pbt'].isin(['YES','Pending'])]
        vpvb_pos = vpvb_df[vpvb_df['vpvb'].isin(['YES','Pending'])]
        vpvb_neg = vpvb_df[~vpvb_df['vpvb'].isin(['YES','Pending'])]
        sens_pos = sens_df[sens_df['sensitiser'].isin(['YES','Pending'])]
        sens_neg = sens_df[~sens_df['sensitiser'].isin(['YES','Pending'])]
        endoc_pos = endoc_df[endoc_df['endocrine_disruptor'].isin(['YES','Pending'])]
        endoc_neg = endoc_df[~endoc_df['endocrine_disruptor'].isin(['YES','Pending'])]