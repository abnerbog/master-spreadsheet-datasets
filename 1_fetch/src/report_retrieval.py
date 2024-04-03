# -*- coding: utf-8 -*-
"""
Created on Thur Mar 28 14:01:27 2024

@author: AbnerBogan
"""

import pickle
import requests
import pandas as pd
from io import StringIO

def get_data_sharing_report(endpoint):
    """
    Get data sharing report from FastAPI endpoint.

    Parameters
    ----------
    endpoint : string
        FastAPI endpoint generating a .csv containing metadata
        of datasets shared by clusters.

    Returns
    -------
    df : pandas dataframe
        Dataframe containing metadata of datasets shared by clusters.

    """
    
    # make get request to API endpoint
    response = requests.get(endpoint)

    # check if request was successful
    if response.status_code == 200:
        # read in csv with first column (dataset number) as index
        df = pd.read_csv(StringIO(response.text),index_col=0)
    else:
        print(f"Failed with status code {response.status_code}")
        
    return df

def save_data_sharing_report(df,out_file):
    """
    Save data sharing report to file as a pickled Python object.

    Parameters
    ----------
    df : pandas dataframe
        Dataframe containing metadata of datasets shared by clusters.
    out_file : string
        Name of output file containing pickled dataframe.

    Returns
    -------
    None.

    """
    
    # pickle (serialize) the data to file
    with open(out_file, 'wb') as file:
        pickle.dump(df, file)  
        
        
def main(endpoint,out_file):
    """
    Main function.

    Parameters
    ----------
    endpoint : string
        FastAPI endpoint generating a .csv containing metadata
        of datasets shared by clusters.
    out_file : string
        Name of output file containing pickled dataframe.

    Returns
    -------
    None.

    """
    
    # retrieve data sharing report from API endpoint (note this can take a few minutes)
    df = get_data_sharing_report(endpoint)
    
    # save dataframe to file
    save_data_sharing_report(df,out_file)
    
        
if __name__ == '__main__':

    api_endpoint_csv = str(snakemake.params['api_endpoint'])

    output_filename = str(snakemake.output['out_filename'])
    
    main(api_endpoint_csv,output_filename)
    
    
    