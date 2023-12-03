import pandas as pd

class List:
    
    def __init__(self, data: list, column_name: str):
        self.__data = data
        self.__column_name = column_name
        
    def __get_data(self) -> list:
        """Return data"""
        return self.__data
    
    def __set_data(self, data: list) -> None:
        """Set data"""
        self.__data = data
    
    def clear(self) -> list:
        """Return list cleared, remove values before column_name and from first NaN"""
        self.__remove_values_before_column_name()
        self.__remove_values_after_first_NaN()
        return self.__get_data()

    def __remove_values_before_column_name(self) -> None:
        """Remove values before column_name, column_name included"""
        data = self.__get_data()
        for i, value in enumerate(data):
            if value == self.__column_name:
                self.__set_data(data[i+1:])                
                return
        return
    
    def __remove_values_after_first_NaN(self) -> None:
        """Remove values from first NaN"""
        data = self.__get_data()
        for i, value in enumerate(data):
            if pd.isna(value):
                self.__set_data(data[:i])
                return
        return