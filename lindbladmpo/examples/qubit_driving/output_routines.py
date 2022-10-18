# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

"""
Routines for managing the output of research projects running multiple simulations.
"""

import csv
import os.path
import pandas as pd
from pathlib import Path
from pandas import DataFrame
from lindbladmpo.plot_routines import *


def generate_paths(s_output_path: str, b_make_paths=True,
                   s_data_subdir = 'data/', s_plot_subdir = 'figures/'):
    """Concatenate a data directory and figures directory path, and optionally create the directories.

    Args:
            s_output_path: The output path, a base directory for the data and figures directories
            b_make_paths: If True the directories are created if missing.
            s_data_subdir: Subdirectory of the output folder where the data files are.
            s_plot_subdir: Subdirectory of the output folder where the figures files are.

    Returns:
            A tuple with the data and plot directory strings.
    """
    s_data_path = ""
    s_plot_path = ""
    if not os.path.exists(s_output_path):
        if b_make_paths:
            os.mkdir(s_output_path)
    s_data_path = s_output_path + s_data_subdir
    if not os.path.exists(s_data_path):
        if b_make_paths:
            os.mkdir(s_data_path)
    s_plot_path = s_output_path + s_plot_subdir
    if not os.path.exists(s_plot_path):
        if b_make_paths:
            os.mkdir(s_plot_path)
    return s_data_path, s_plot_path


def save_to_db(s_db_path: str, sim_metadata: dict):
    """Save a data line into the .csv dataframe file using pandas.

    Args:
            s_db_path: The full filename for the .csv dataframe file.
            sim_metadata: The dictionary giving one line of db data.
    """
    if not os.path.isfile(s_db_path):
        # New database, write header line based on metadata keys
        with open(s_db_path, "w") as f:
            header = sim_metadata.keys()
            writer = csv.writer(f)
            writer.writerow(header)
            f.close()

    df = pd.read_csv(s_db_path)
    db_line = {}
    for key in sim_metadata.keys():
        db_line[key] = [sim_metadata[key]]
    db_data = pd.DataFrame(db_line)
    df = df.append(db_data)
    df.to_csv(s_db_path, index=False)


def find_db_files(s_db_path: str):
    """Find all csv files in a given directory and store them in a list.

    Args:
            s_db_path: A string with the data path to the database folder.

    Returns:
        files: A list with db (.csv) files in the folder.
    """
    if os.path.isfile(s_db_path):
        files = [s_db_path]
    else:
        files = []
        s_data_path = Path(s_db_path)
        for item in s_data_path.iterdir():
            if item.suffix in [".csv"]:
                files.append(Path.resolve(item))
    return files


def query_simulations(
    files: List[str],
    s_filter_query: str,
    sort_by: Optional[Any] = None,
    ascending: Union[bool, List[bool]] = True,
    na_position="last",
):
    """Find the simulations according to the criteria query and metadata dictionaries in a list.

    Args:
            files: A list with db (.csv) files.
            s_filter_query: A string with the desired query.
            sort_by: If not None, defines sorting using the method sort_values() of the data frame.
            ascending: If sort_by is not None, defines an option for the sort
            na_position: If sort_by is not None, defines an option for the sort

    Returns:
            A list with the relevant simulation dicts.
    """
    sims = []
    for file in files:
        df = pd.read_csv(file)
        df_2 = df.query(s_filter_query)
        if sort_by is not None:
            df_3 = df_2.sort_values(
                sort_by, ascending=ascending, na_position=na_position
            )
        else:
            df_3 = df_2
        sims.extend(_take_list(df_3))
    return sims


def get_simulation_dict(s_output_path: str, s_unique_id: str):
    """Returns a dictionary corresponding to a unique simulation id from the dataframe.

    Args:
            s_output_path: A string with the data path to the database folder.
            s_unique_id: The uuid of the simulation to return.

    Returns:
            The requested simulation dict.
    """
    files = find_db_files(s_output_path)
    sim_dict = None
    for file in files:
        df = pd.read_csv(file)
        if 'unique_id' in df.keys():
            df_2 = df.query(f"unique_id == '{s_unique_id}'")
        else:
            df_2 = None
        if df_2 is not None and not df_2.empty:
            sim_dict = _take_list(df_2)[0]
            break
    return sim_dict


def _take_list(df: DataFrame):
    sims = []
    for i in range(len(df)):
        sim_dict = df.take([i]).to_dict("list")
        for key in sim_dict.keys():
            sim_dict[key] = sim_dict[key][0]
        sims.append(sim_dict)
    return sims
