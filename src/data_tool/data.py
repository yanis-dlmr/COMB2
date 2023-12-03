import pandas as pd

from .list import List

class Data(object):
    def __init__(self, data_path: str, sheet_name: list = None):
        """Data class for loading data from data_path which is a .xlsx file"""
        self.__data_path = data_path
        self.__sheet_name = sheet_name
        self.__data = self.__load_data()

    def __load_data(self) -> pd.DataFrame:
        """Load data from data_path"""
        return pd.read_excel(self.__data_path, sheet_name=self.__sheet_name)

    def __get_data_sheet(self, sheet_name: str) -> pd.DataFrame:
        """Return data"""
        if sheet_name not in self.__data.keys():
            raise Exception(f"sheet_name '{sheet_name}' not in {self.__data.keys()}")
        return self.__data[sheet_name]
    
    def get_list(self, sheet_name: str, list_name: str) -> list:
        """
        Return a list from sheet_name with the label list_name.
        Clear the list before returning it.
        """
        page = self.__get_data_sheet(sheet_name)
        all_column = self.__get_all_column(page, list_name)
        list_column = List(all_column, list_name)
        final_list = list_column.clear()
        if final_list != []:
            return final_list
        all_raw = self.__get_all_raw(page, list_name)
        list_raw = List(all_raw, list_name)
        final_list = list_raw.clear()
        if final_list != []:
            return final_list
        return None
    
    def __get_all_column(self, page: pd.DataFrame, list_name: str) -> list:
        """Return all column as a list"""
        for i in range(page.shape[0]):
            for j in range(page.shape[1]):
                if page.iloc[i, j] == list_name:
                    return page.iloc[:, j].tolist()
        raise Exception(f"list_name '{list_name}' not in this page")
    
    def __get_all_raw(self, page: pd.DataFrame, list_name: str) -> list:
        """Return all raw as a list"""
        for i in range(page.shape[0]):
            for j in range(page.shape[1]):
                if page.iloc[i, j] == list_name:
                    return page.iloc[i, :].tolist()
        raise Exception(f"list_name '{list_name}' not in this page")