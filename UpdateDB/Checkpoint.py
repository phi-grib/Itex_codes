"""
    Created by: Eric March Vila (eric.march@upf.edu)
    On: 26/03/2020, 10:55 AM
"""

import pandas as pd

from UpdateDB.Connect_CII import Connector

class Checkpoint(Connector):
    """
        This class aims to get the maximum IDs of a stable version of CII
        in order to be able to go back to the stable state once we update 
        the development version and we see any error
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
    
    def get_max_id_for_each_table(self) -> pd.DataFrame:
        """
            Get maximum IDs for each table from CII database

            :return tables_ids:
        """

        list_of_tables = self.get_list_of_tables()
        id_tab_df = pd.DataFrame(index=range(len(list_of_tables)), columns=['Table_name','max_id'])

        for i, table in enumerate(list_of_tables):
            max_id_cmd = """SELECT max(id) FROM {}""".format(table)
            max_id = self.get_id(max_id_cmd)
            id_tab_df.loc[id_tab_df.index == i, 'Table_name'] = table
            id_tab_df.loc[id_tab_df.index == i, 'max_id'] = max_id

        id_tab_df = self.sort_dataframe(id_tab_df)
        
        return id_tab_df

    def get_list_of_tables(self) -> pd.DataFrame:
        """
            Gets a list of all the tables in the database.

            :returns list_of_tables:
        """

        list_of_tables = pd.read_sql_query("""SELECT * FROM information_schema.tables WHERE table_schema = 'public'""", self.conn)

        tables = list_of_tables['table_name']

        return tables

    def get_id(self, query: str) -> int:
        """
            Gets max id from the input query

            :param query: query to check

            :return id_:
        """

        self.curs.execute(query)
        id_ = self.curs.fetchone()[0]
        self.conn.commit()
        
        return id_
    
    def sort_dataframe(self, table_dataframe: pd.DataFrame) -> pd.DataFrame:
        """
            Puts regulations table the first one because if afterwards we try to delete
            a row from that table, code will break because of Referential Integrity, since there are other tables with foreign key constraints
            that refer to table regulations. Instead of applying ON DELETE CASCADE, we prefer to just remove first the referenced rows and then the
            rows that are used as foreign keys.

            :param table_dataframe:

            :return sorted_table_dataframe:
        """
        
        reg_row = table_dataframe[table_dataframe['Table_name'] == 'regulations'].reset_index(drop=True)
        new_tab_df = table_dataframe.drop(index=table_dataframe[table_dataframe['Table_name'] == 'regulations'].index).reset_index(drop=True)
        sorted_table_dataframe = pd.concat([reg_row,new_tab_df]).reset_index(drop=True)

        return sorted_table_dataframe
