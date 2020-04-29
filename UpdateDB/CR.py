"""
    Created by: Eric March Vila (eric.march@upf.edu)
    On: 30/03/2020, 12:27 AM
"""

import pandas as pd

from UpdateDB.Connect_CII import Connector

class CR(Connector):
    """
        This class uses different functions to interact with CR database.
        Useful for updating CII annotations, since we create the CR object here and
        the dataframes generated from the inner queries get automatically cleaned here,
        saving some time in data cleaning in jupyter notebooks
    """

    def __init__(self):
        """
            This class gets initialized by using the self.compoundDB_connection function from Connector parent class.
            We get the self.compounddb_conn that will allow us to interact with CR database
        """

        self.compoundDB_connection()
        self.endpoint_dict = {'CMR':['Carc. 1 A',
                            'Carc. 1A',
                            'Carc. 1B',
                            'Carc. 2',
                            'Carcinogenic (Article 57a)',
                            'Carcinogenicity','Effect on or via lactation','Forward gene mutation at the HPRT locus','Genotoxic','Harmonised classification for carcinogenicity',
                            'Harmonised classification for mutagenicity',
                            'Harmonised classification for reprotoxicity','Harmonised classification for effects on or via lactation','Histidine reverse gene mutation','Lact.',
                            'Muta. 1A',
                            'Muta. 1B',
                            'Muta. 2',
                            'Mutagenic (Article 57b)',
                            'Mutagenicity','Repr. 1A',
                            'Repr. 1B',
                            'Repr. 2',
                            'Repr.1B',
                            'Reprotoxicity','Toxic for reproduction (Article 57c)','chromosome aberrations','Suspected carcinogen','Suspected mutagen','Suspected toxic for reproduction'],
                            'PBT':['PBT','PBT (Article 57d)','Suspected bioaccumulative','Suspected persistent in the environment'],
                            'vPvB':['vPvB','vPvB (Article 57e)','Suspected bioaccumulative','Suspected persistent in the environment'],
                            'Sensitiser':['Resp. Sens. 1',
                            'Resp. Sens. 1A',
                            'Resp. Sens. 1B','Sensitiser','Skin Sens. 1',
                            'Skin Sens. 1 A',
                            'Skin Sens. 1 B',
                            'Skin Sens. 1A',
                            'Skin Sens. 1A; H317',
                            'Skin Sens. 1B',
                            'Skin Sens. 1B H317','Harmonised classification for respiratory sensitisation',
                            'Harmonised classification for skin sensitisation','Respiratory sensitising properties (Article 57(f) - human health)','Suspected respiratory sensitiser',
                            'Suspected skin sensitiser'],
                            'Endocrine Disruptor':['ENDOCRINE - EFFECTS','Endocrine disrupting properties (Article 57(f) - environment)',
                            'Endocrine disrupting properties (Article 57(f) - human health)',
                            'Endocrine disruptor']}
    
    def get_annotations_per_CAS(self) -> pd.DataFrame:
        """
            Retrieve all the annotations in CR with is respective sources for each CAS.

            :returns annotation_dataframe:
        """

        annotation_dataframe = pd.read_sql_query("""SELECT synonym.name as reg_number, source.name as source_name, 
                                                subs_ann.original_annotation, annotation.annotation, annotation.general, annotation.category
                                                FROM substance sub
                                                left join synonym on synonym.subsid = sub.id
                                                left join source on source.id = sub.sourceid
                                                left join subs_ann on subs_ann.subsid = sub.id
                                                left join annotation on annotation.id = subs_ann.annid
                                                where synonym.type like '%CAS%'
                                                order by synonym.name ASC""", self.compounddb_conn)
        
        annotation_dataframe = self.clean_annotations_dataframe(annotation_dataframe)
        
        return annotation_dataframe
    
    def clean_annotations_dataframe(self, ann_df: pd.DataFrame):
        """
            Removes invalid values from the dataframe and changes column name and dataframe order.

            :param ann_df:

            :return ann_df
        """

        ann_df.drop_duplicates(inplace=True)
        ann_df.rename(columns={'reg_number':'CAS'},inplace=True)
        ann_df.drop(labels=ann_df[ann_df['CAS'].isin(['-','_','---','—','Multiple CAS No. Possible','NOCAS_CMS-88808'])].index, axis=0, inplace=True)
        ann_df.drop(labels=ann_df[ann_df['CAS'].str.contains('/')].index,axis=0, inplace=True)
        ann_df.drop(labels=ann_df[ann_df['CAS'].str.contains('nan')].index,axis=0, inplace=True)
        ann_df.loc[ann_df['CAS'].str.contains('-,'), 'CAS'] = ann_df[ann_df['CAS'].str.contains('-,')].CAS.apply(lambda x: str(x).replace('-,  ',''))
        ann_df.sort_values(by=['CAS','source_name','original_annotation','annotation'], inplace=True)

        return ann_df