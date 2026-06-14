import re

import pandas as pd
import numpy as np

from typing import Union
from difflib import SequenceMatcher
from collections import namedtuple


from difflib import SequenceMatcher
from collections import namedtuple

import pypath.share.session as session

_logger = session.Logger(name = 'inputs.expression_atlas.processor')
_log = _logger._log


baseline_experiement = namedtuple("BaselineExperiement", ["tpm_values", "experiement_design", "matching_factor"])

class BaselineExperiementDataProcessor:
    def __init__(self, 
                 tpm_file_path: str,
                 experiment_design_file_path: str,
                 matching_factor: str,
                 skip_bad_data: bool = False,
                 baseline_experiement_namedtuple: namedtuple = baseline_experiement) -> None:
        
        self.tpm_file_path = tpm_file_path
        self.experiment_design_file_path = experiment_design_file_path
        self.matching_factor = matching_factor
        self.skip_bad_data = skip_bad_data
        self.baseline_experiement_namedtuple = baseline_experiement_namedtuple
    

    def read_data(self, file_path: str) -> Union[pd.DataFrame, bool]:
        """
        Reads a tab-separated values (TSV) file into a pandas DataFrame, 
        drops columns and rows that contain only NaN values, and resets the index.

        Args:
            file_path (str): The path to the TSV file to be read.

        Returns:
            pd.DataFrame: A pandas DataFrame containing the data from the file, 
                  with columns and rows containing only NaN values removed, 
                  and the index reset.
        """
        df = pd.read_csv(file_path, sep="\t", comment='#',)

        if self.skip_bad_data and df.empty:
            return False
        elif df.empty:
            raise ValueError(f"No data found in the file: {file_path}")

        # drop columns with all NaN values except the specific column
        protected_columns = [
            "Gene ID",
            "Gene Name",
            f"Sample Characteristic Ontology Term[{self.matching_factor}]",
            f"Factor Value Ontology Term[{self.matching_factor}]",
        ]

        columns_to_drop = [
            col for col in df.columns
            if col not in protected_columns
               and df[col].isna().all()
        ]
        df.drop(columns=columns_to_drop, inplace=True)

        # drop rows with all NaN values
        df.dropna(axis=0, how='all', inplace=True)

        return df.reset_index(drop=True)

    def read_tpm_data(self) -> Union[pd.DataFrame, bool]:
        """Read data from the file path."""
        df = self.read_data(file_path=self.tpm_file_path)

        if df is False:
            return False
        
        if any([True if " and " in col else False for col in df.columns]):
            _log(f"Columns contain ' and ' in their names. Columns: {df.columns.to_list()}. File: {self.tpm_file_path}")
            df.drop(columns=[col for col in df.columns if " and " in col], inplace=True)


        df.columns = df.columns.str.strip()

        assert all(df.columns[:2] == ["Gene ID", "Gene Name"]), "Gene ID and Gene Name columns not found"
        # drop rows if Gene ID is NaN
        df.dropna(subset="Gene ID", inplace=True, axis=0)

        return df.reset_index(drop=True) 
    
    def read_experiement_design_data(self) -> Union[pd.DataFrame, bool]:
        """Read experiment design data from the file path."""
        df = self.read_data(file_path=self.experiment_design_file_path)

        if df is False:
            return False
        
        if self.skip_bad_data and f"Sample Characteristic[{self.matching_factor}]" not in df.columns:
            raise False
        elif f"Sample Characteristic[{self.matching_factor}]" not in df.columns:
            raise ValueError(f"Sample Characteristic[{self.matching_factor}] column not found")

         # drop rows if Sample Characteristic[matching_factor] is NaN
        df.dropna(subset=f"Sample Characteristic[{self.matching_factor}]", inplace=True, axis=0)

        # Drop rows if the "Sample Characteristic[matching_factor]" column contains " and "
        df = df[~df[f"Sample Characteristic[{self.matching_factor}]"].str.contains(" and ", na=False)]

        df[f"Sample Characteristic[{self.matching_factor}]"] = np.where((pd.notna(df[f"Factor Value[{self.matching_factor}]"])) & (df[f"Sample Characteristic[{self.matching_factor}]"] != df[f"Factor Value[{self.matching_factor}]"]), df[f"Factor Value[{self.matching_factor}]"], df[f"Sample Characteristic[{self.matching_factor}]"])
        
        if self.matching_factor == "cell line":
            df.loc[df[f"Sample Characteristic Ontology Term[{self.matching_factor}]"].isna(), f"Sample Characteristic Ontology Term[{self.matching_factor}]"] = df[f"Sample Characteristic[{self.matching_factor}]"]
            df.loc[df[f"Factor Value Ontology Term[{self.matching_factor}]"].isna(), f"Factor Value Ontology Term[{self.matching_factor}]"] = df[f"Factor Value[{self.matching_factor}]"]
        
        if self.skip_bad_data and not {f"Sample Characteristic Ontology Term[{self.matching_factor}]", f"Factor Value Ontology Term[{self.matching_factor}]"}.issubset(set(df.columns)):
            return False
        elif not {f"Sample Characteristic Ontology Term[{self.matching_factor}]", f"Factor Value Ontology Term[{self.matching_factor}]"}.issubset(set(df.columns)):
            raise ValueError(f"Sample Characteristic Ontology Term[{self.matching_factor}] and Factor Value Ontology Term[{self.matching_factor}] columns not found")
        
        if self.skip_bad_data and df[[f"Sample Characteristic Ontology Term[{self.matching_factor}]", f"Factor Value Ontology Term[{self.matching_factor}]"]].isna().all().all():
            return False
        elif df[[f"Sample Characteristic Ontology Term[{self.matching_factor}]", f"Factor Value Ontology Term[{self.matching_factor}]"]].isna().all().all():
            raise ValueError(f"Sample Characteristic Ontology Term[{self.matching_factor}] and Factor Value Ontology Term[{self.matching_factor}] columns are full of NaN")
        
        if df[f"Sample Characteristic Ontology Term[{self.matching_factor}]"].isna().all():
            df.loc[df[f"Sample Characteristic Ontology Term[{self.matching_factor}]"].isna(), f"Sample Characteristic Ontology Term[{self.matching_factor}]"] = df[f"Sample Characteristic[{self.matching_factor}]"]
        
        if not pd.Series.equals(df[f"Sample Characteristic Ontology Term[{self.matching_factor}]"], df[f"Factor Value Ontology Term[{self.matching_factor}]"]):
            df = df[df[f"Sample Characteristic Ontology Term[{self.matching_factor}]"] == df[f"Factor Value Ontology Term[{self.matching_factor}]"]]

            if self.skip_bad_data and df.empty == False:
                return False
            elif df.empty == False:
                raise ValueError(f"Sample Characteristic Ontology Term[{self.matching_factor}] and Factor Value Ontology Term[{self.matching_factor}] columns do not match")
        
        ontology_terms_to_extract = ["organism part", "organism", "cell type", "cell line", "compound", "infect", "disease"]
        selected_columns = []
        new_columns_names = []
        for col in df.columns:
            out = re.search(r'\[(.*?)\]', col).group(1) if re.search(r'\[(.*?)\]', col) else col
            if col.startswith("Sample Characteristic") and "Ontology Term" in col and out in ontology_terms_to_extract:
                selected_columns.append(col)
                new_columns_names.append(f'{out.replace(" ","_")}_ontology_term')
        
        self.unique_sample_characteristic_factor_terms = df[f"Sample Characteristic[{self.matching_factor}]"].unique()
        new_df = df[selected_columns]
        new_df.columns = new_columns_names

        if new_df.empty:
            print("No matching columns found")
            return False

        for col in new_df.columns:
            new_df.loc[:, col] = new_df[col].apply(lambda x: x.rsplit('/', 1)[-1].replace('_', ':') if pd.notna(x) else x)

        new_df[self.matching_factor] = df[f"Sample Characteristic[{self.matching_factor}]"]

        new_df.drop_duplicates(subset=[self.matching_factor, f'{self.matching_factor.replace(" ","_")}_ontology_term'], inplace=True, ignore_index=True)

        assert new_df.shape[0] == len(df[f"Sample Characteristic[{self.matching_factor}]"].unique()), "Number of unique sample characteristic terms do not match"
        
        return new_df.reset_index(drop=True)
    
    def create_dictionary_from_experiement_design_data(self) -> Union[dict, bool]:
        """Create a dictionary from the experiment design data."""
        df = self.read_experiement_design_data()
        if df is False:
            return False
        
        df = df.where(pd.notna(df), None)  # Replace NaNs with None
        return df.set_index(self.matching_factor).T.to_dict()
    
    def merge_tpm_and_experiment_design_data(self) -> Union[bool, 
                                                            tuple[
                                                                namedtuple, 
                                                                list[dict[str, Union[str, float, None]]], 
                                                                dict[str, dict]
                                                                ]
                                                            ]:
        """
        Merges TPM (Transcripts Per Million) data with experiment design data.
        This method performs the following steps:
        1. Reads TPM data into a DataFrame.
        2. Creates a dictionary from the experiment design data.
        3. Finds the best string matches between the TPM data columns and the experiment design keys.
        4. Drops columns from the TPM data DataFrame that do not have a match in the matches dictionary.
        5. Renames columns in the TPM data DataFrame to match the keys in the matches dictionary.
        6. Converts the TPM data DataFrame into a list of dictionaries, with specific formatting for certain keys.
        7. Creates a namedtuple using the TPM values, experiment design dictionary, and a matching factor.
        Returns:
            tuple: A tuple containing:
            - namedtuple: A namedtuple containing the merged data.
            - list[dict[str, Union[str, float, None]]]: A list of dictionaries representing the TPM values.
            - dict[str, dict]: A dictionary representing the experiment design data.
        """

        tpm_data_df = self.read_tpm_data()
        experiment_design_dict = self.create_dictionary_from_experiement_design_data()

        if tpm_data_df is False or experiment_design_dict is False:
            return False
        
        matches = self.find_best_string_matchs(tpm_data_df.columns[2:].to_list(), list(experiment_design_dict.keys()))

        # Drop columns from tpm_data_df if their names (from second to last) do not have a match in the matches dict
        tpm_data_df.drop(columns=[col for col in tpm_data_df.columns[2:] if matches.get(col) is None], inplace=True)

        # Rename columns in tpm_data_df to match the keys in the matches dict
        tpm_data_df.rename(columns=matches, inplace=True)

        tpm_values = []
        for _, row in tpm_data_df.iterrows():
            row_dict = {k.lower().replace(" ","_") if k in ["Gene ID","Gene Name"] else k: (v if v != 0 else None) for k, v in row.where(pd.notna(row), None).to_dict().items()}
            tpm_values.append(row_dict)


        nmdtuple = self.baseline_experiement_namedtuple(tpm_values, experiment_design_dict, self.matching_factor)

        return nmdtuple, tpm_values, experiment_design_dict
    
    @staticmethod
    def calculate_string_similarity(a: str, b: str) -> float:
        """
        Calculate the similarity ratio between two strings.

        Args:
            a (str): The first string to compare.
            b (str): The second string to compare.

        Returns:
            float: A similarity ratio between 0 and 1.
        """
        return SequenceMatcher(None, a, b).ratio()
    
    def find_best_string_matchs(self, a: list[str], b: list[str]) -> dict[str, str]:
        """
        Find the best string matches between two lists of strings.
        This method compares each string in the first list `a` with all strings in the second list `b`
        to find the best match based on a string similarity score. If the best match score is less than or equal
        to 0.55, it considers that NO best match is found for that string.
        Args:
            a (list[str]): The first list of strings to be matched.
            b (list[str]): The second list of strings to match against.
        Returns:
            dict[str, str]: A dictionary where keys are strings from the first list `a` and values
                         are the best matching strings from the second list `b`.
        """

        matches = {}
        for first_el in a:
            best_score = 0
            best_match = None
            for second_el in b:
                score = self.calculate_string_similarity(first_el, second_el)
                if score > best_score:
                    best_score = score
                    best_match = second_el
            
            if best_score <= 0.55:
                _log(f"No best match found for {first_el}")
            else:
                matches[first_el] = best_match

        return matches



differential_experiement = namedtuple("DifferentialExperiement", ["log2fold_change_values", "experiement_design", "matching_factor",
                                                                  "comparison_factors_raw"])

class DifferentialExperiementDataProcessor:

    def __init__(self,
                 log2foldchange_file_path: str,
                 experiment_design_file_path: str,
                 matching_factor: str,                    
                skip_bad_data: bool = False,
                differential_experiement_namedtuple: namedtuple = differential_experiement) -> None:                

        self.log2foldchange_file_path = log2foldchange_file_path
        self.experiment_design_file_path = experiment_design_file_path
        self.matching_factor = matching_factor
        self.skip_bad_data = skip_bad_data
        self.differential_experiement_namedtuple = differential_experiement_namedtuple
        
        

    def read_data(self, file_path: str) -> Union[pd.DataFrame, bool]:
        """
        Reads a tab-separated values (TSV) file into a pandas DataFrame, 
        drops columns and rows that contain only NaN values, and resets the index.

        Args:
            file_path (str): The path to the TSV file to be read.

        Returns:
            pd.DataFrame: A pandas DataFrame containing the data from the file, 
                  with columns and rows containing only NaN values removed, 
                  and the index reset.
        """
        df = pd.read_csv(file_path, sep="\t", comment='#',)

        if self.skip_bad_data and df.empty:
            return False
        elif df.empty:
            raise ValueError(f"No data found in the file: {file_path}")

        # drop columns with all NaN values
        df.dropna(axis=1, how='all', inplace=True)

        # drop rows with all NaN values
        df.dropna(axis=0, how='all', inplace=True)

        return df.reset_index(drop=True)
    
    def read_log2fold_data(self) -> Union[pd.DataFrame, bool]:
        """Read log2fold data from the file path."""
        df = self.read_data(file_path=self.log2foldchange_file_path)

        if df is False:
            return False
        
        if self.skip_bad_data and any([True if " and " in col else False for col in df.columns]):
            return False
        elif any([True if " and " in col else False for col in df.columns]):
            raise ValueError(f"Columns contain ' and ' in their names. Columns: {df.columns.to_list()}. File: {self.log2foldchange_file_path}")
        
        assert df.columns[0] == "Gene ID", "Gene ID column not found"
        # drop rows if Gene ID is NaN
        df.dropna(subset="Gene ID", inplace=True, axis=0)

        if "Design Element" in df.columns:
            df.drop(columns="Design Element", inplace=True)

        return df.reset_index(drop=True)
    
    def read_experiement_design_data(self) -> Union[pd.DataFrame, bool]:
        
        df = self.read_data(file_path=self.experiment_design_file_path)

        if df is False:
            return False
        
        if self.skip_bad_data and (f"Factor Value[{self.matching_factor}]" not in df.columns or f"Factor Value Ontology Term[{self.matching_factor}]" not in df.columns):
            return False
        elif f"Factor Value[{self.matching_factor}]" not in df.columns or f"Factor Value Ontology Term[{self.matching_factor}]" not in df.columns:
            raise ValueError(f"Factor Value[{self.matching_factor}] or Factor Value Ontology Term[{self.matching_factor}] column(s) not found")

        # drop rows if Sample Characteristic[matching_factor] is NaN
        df.dropna(subset=f"Factor Value[{self.matching_factor}]", inplace=True, axis=0)

        # Drop rows if the "Factor Value[matching_factor]" column contains " and "
        df = df[~df[f"Factor Value[{self.matching_factor}]"].str.contains(" and ", na=False)]

        if self.matching_factor == "cell line":
            df.loc[df[f"Factor Value Ontology Term[{self.matching_factor}]"].isna(), f"Factor Value Ontology Term[{self.matching_factor}]"] = df[f"Factor Value[{self.matching_factor}]"]
        

        if df[f"Factor Value Ontology Term[{self.matching_factor}]"].isna().all():
            df.loc[df[f"Factor Value Ontology Term[{self.matching_factor}]"].isna(), f"Factor Value Ontology Term[{self.matching_factor}]"] = df[f"Factor Value[{self.matching_factor}]"]
        
        ontology_terms_to_extract = ["organism part", "organism", "cell type", "cell line", "compound", "infect", "disease"]
        selected_columns = []
        new_columns_names = []
        for col in df.columns:
            out = re.search(r'\[(.*?)\]', col).group(1) if re.search(r'\[(.*?)\]', col) else col
            if (col.startswith("Sample Characteristic") or col.startswith("Factor Value")) and "Ontology Term" in col and out in ontology_terms_to_extract:
                if col == f"Factor Value Ontology Term[{self.matching_factor}]" and f"Sample Characteristic Ontology Term[{self.matching_factor}]" in selected_columns:
                    continue
                    
                selected_columns.append(col)
                new_columns_names.append(f'{out.replace(" ","_")}_ontology_term')

        self.unique_sample_characteristic_factor_terms = df[f"Factor Value[{self.matching_factor}]"].unique()
        new_df = df[selected_columns]
        new_df.columns = new_columns_names
        
        assert new_df.empty == False, "No matching columns found"

        for col in new_df.columns:
            new_df.loc[:, col] = new_df[col].apply(lambda x: x.rsplit('/', 1)[-1].replace('_', ':') if pd.notna(x) else x)

        new_df[self.matching_factor] = df[f"Factor Value[{self.matching_factor}]"]

        new_df.drop_duplicates(subset=self.matching_factor, inplace=True, ignore_index=True)

        assert new_df.shape[0] == len(self.unique_sample_characteristic_factor_terms), "Number of unique sample characteristic terms do not match"

        return new_df.reset_index(drop=True)

    def create_dictionary_from_experiement_design_data(self) -> Union[dict, bool]:
        """Create a dictionary from the experiment design data."""
        df = self.read_experiement_design_data()
        
        if df is False:
            return False
        
        df = df.where(pd.notna(df), None)  # Replace NaNs with None
        return {k:v for k, v in df.set_index(self.matching_factor).T.to_dict().items() if k}
        
    
    def merge_log2fold_and_experiment_design_data(self) -> Union[bool, namedtuple]:
        log2foldchange_df = self.read_log2fold_data()
        experiment_design_dict = self.create_dictionary_from_experiement_design_data()

        if log2foldchange_df is False or experiment_design_dict is False:
            return False
        
        if "Gene Name" not in log2foldchange_df.columns:
            comparison_factors = self.extract_comparison_factors(log2foldchange_df.columns[2:])
        else:
            comparison_factors = self.extract_comparison_factors(log2foldchange_df.columns[1:])

        if self.skip_bad_data and not comparison_factors:
            return False
        elif not comparison_factors:
            raise ValueError(f"Comparison factors not found. File: {self.log2foldchange_file_path}. ")

        matches = self.find_best_string_matchs([list(d.keys())[0] for d in comparison_factors], list(experiment_design_dict.keys()))
        
        if self.skip_bad_data and not matches:
            return False
        elif not matches:
            raise ValueError(f"No best matches found for comparison factors. File: {self.log2foldchange_file_path}")
        
        # Drop columns from log2foldchange_df if their names (from second to last) do not have a match in the matches dict
        rename_dict = {}
        self.comparison_factors_raw = []
        founds = []
        for col in log2foldchange_df.columns:
            if col in ["Gene ID", "Gene Name"]:
                continue
            
            found = False
            for k, v in matches.items():
                
                if k in col:
                    if col.endswith("p-value"):
                        rename_dict[col] = f"{v}_p_value"
                        self.comparison_factors_raw.append(col.replace("p-value", ""))
                    elif col.endswith("t-statistic"):
                        rename_dict[col] = f"{v}_t_statistic"
                    elif col.endswith("log2foldchange"):
                        rename_dict[col] = f"{v}_log2foldchange"
                    else:
                        raise ValueError(f"Column name endwith invalid string. Column name: {col}. File: {self.log2foldchange_file_path}")
                    
                    found = True
                    break
            
            founds.append(found)

        if not any(founds):
            raise ValueError(f"Matching factor not found in the File: {self.log2foldchange_file_path}")
        

        # if the rows does not in rename_dict, drop them
        log2foldchange_df.drop(columns=[col for col in log2foldchange_df.columns if col not in rename_dict.keys() and col not in ["Gene ID", "Gene Name"]], inplace=True)
        
        if self.skip_bad_data and (log2foldchange_df.empty or set(log2foldchange_df.columns) == {"Gene ID", "Gene Name"}):
            return False
        elif (log2foldchange_df.empty or set(log2foldchange_df.columns) == {"Gene ID", "Gene Name"}):
            raise ValueError(f"No matching data found in the file: {self.log2foldchange_file_path}")
        

        
        log2foldchange_df.rename(columns=rename_dict, inplace=True)

        log2foldchange_values = []
        for _, row in log2foldchange_df.iterrows():
            row_dict = {k.lower().replace(" ","_") if k in ["Gene ID","Gene Name"] else k: (v if v != 0 else None) for k, v in row.where(pd.notna(row), None).to_dict().items()}
            log2foldchange_values.append(row_dict)


        nmdtuple = self.differential_experiement_namedtuple(log2foldchange_values, experiment_design_dict, self.matching_factor, self.comparison_factors_raw)

        return nmdtuple
    
    def extract_comparison_factors(self, columns: list[str]) -> Union[list[dict[str, str]], bool]:
        factors = []
        for col in columns:
            factor = self.extract_comparison_factor(col)
            if factor and {factor[0]: factor[1]} not in factors:
                factors.append({factor[0]: factor[1]})
        
        return factors if factors else False
    
    def extract_comparison_factor(self, column: str) -> Union[tuple[str, str], bool]:

        selected_comparison_2_factors = ["normal", "none", "control", "uninfected", "healthy",
                                         "vehicle", "DMSO", "mock", "water", "placebo", "DMSO control"]

        pattern = r"'([^']*)'\s+vs\s+'([^']*)'\.\w+"
        match = re.search(pattern, column)
        if match:
            comparison_factor_1, comparison_factor_2 = match.groups()
            comparison_factor_2 = comparison_factor_2.strip().strip(";")
            if comparison_factor_2 in selected_comparison_2_factors:
                return comparison_factor_1, comparison_factor_2
            
            return False
        
        return False
    
    @staticmethod
    def calculate_string_similarity(a: str, b: str) -> float:
        """
        Calculate the similarity ratio between two strings.

        Args:
            a (str): The first string to compare.
            b (str): The second string to compare.

        Returns:
            float: A similarity ratio between 0 and 1.
        """
        return SequenceMatcher(None, a, b).ratio()
    
    def find_best_string_matchs(self, a: list[str], b: list[str]) -> dict[str, str]:
        """
        Find the best string matches between two lists of strings.
        This method compares each string in the first list `a` with all strings in the second list `b`
        to find the best match based on a string similarity score. If the best match score is less than or equal
        to 0.55, it considers that NO best match is found for that string.
        Args:
            a (list[str]): The first list of strings to be matched.
            b (list[str]): The second list of strings to match against.
        Returns:
            dict[str, str]: A dictionary where keys are strings from the first list `a` and values
                         are the best matching strings from the second list `b`.
        """

        matches = {}
        for first_el in a:
            best_score = 0
            best_match = None
            for second_el in b:
                score = self.calculate_string_similarity(first_el, second_el)
                if score > best_score:
                    best_score = score
                    best_match = second_el
            
            if best_score <= 0.49:
                _log(f"No best match found for {first_el}. Score: {best_score}")
            else:
                matches[first_el] = best_match

        return matches

        
