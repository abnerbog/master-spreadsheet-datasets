# -*- coding: utf-8 -*-
"""
Created on Thurs Mar 28 15:07:13 2024

@author: AbnerBogan
"""

import pickle
import openpyxl 
from openpyxl.utils.dataframe import dataframe_to_rows

def load_cleaned_dataframe(in_file):
    """
    Load cleaned dataframe via unpickling.

    Parameters
    ----------
    in_file : string
        Name of input file containing cleaned dataframe.

    Returns
    -------
    df : pandas dataframe
        Cleaned dataframe loaded from file.

    """
    
    # load file via unpickling
    with open(in_file, 'rb') as file:
        df = pickle.load(file)
        
    return df

def split_dataframe_by_cluster(df,clusters=['Bedrock','Big Data','CINET','Coastal',
                                            'Drylands','DUST^2','Dynamic Water','GeoMicroBio','Urban']):
    """
    Split dataframe containing all cluster dataset metadata into a dictionary.

    Parameters
    ----------
    df : pandas dataframe
        Cleaned dataframe loaded from file.
    clusters : list, optional.
        List of strings containing clusters of interest for creating summary spreadsheet. The default is
        ['Bedrock','Big Data','CINET','Coastal','Drylands','DUST^2','Dynamic Water','GeoMicroBio','Urban']

    Returns
    -------
    clusters_dict : dict
        Dictionary containing keys as each cluster and values as a dataframe of associated dataset metadata.

    """
    
    # initialize dictionary
    clusters_dict = {}

    # iterate through list of clusters for populating dictionary
    for cluster in clusters:

        cluster_key = cluster
        
        # modify DUST^2 as it is written in API endpoint and add '\\' for subsequent 
        # regex pattern search with special characters
        if cluster == 'DUST^2':
            cluster = 'Dust\\^2'
        
        # filter rows of dataframe corresponding to a given cluster,
        # enter key value pair as cluster and filtered dataframe respectively
        clusters_dict[cluster_key] = df[df['clusters'].str.contains(cluster)]
        
    return clusters_dict

        
def create_spreadsheet(clusters_dict):
    """
    Make summary spreadsheet.

    Parameters
    ----------
    clusters_dict : dict
        Dictionary containing keys as each cluster and values as a dataframe of associated dataset metadata.

    Returns
    -------
    wb : openpxyl workbook
        Workbook object containing spreadsheet where each sheet contains dataset metadata of a given cluster.

    """
    
    # create loop for building spreadsheet
    for i, cluster in enumerate(clusters_dict):
        
        # retrieve dataframe for a given cluster key
        cluster_data = clusters_dict[cluster]

        # # rename Dust^2 due to special character
        # if cluster == 'Dust\\^2':
        #     cluster = 'Dust^2'
        
        # initialize workbook
        if i == 0:
            wb = openpyxl.Workbook()
            ws = wb.active
            
            # define title of sheet as cluster name
            ws.title = cluster
            
        else:
            ws = wb.create_sheet(title=cluster)
        
        # populate spreadsheet with dataframe
        for r in dataframe_to_rows(cluster_data, index=False, header=True):
            ws.append(r)        
              
    return wb

def save_temp_spreadsheet(wb,tmp_file):
    """
    Save spreadsheet as a temporary file.

    Parameters
    ----------
    wb : openpxyl workbook
        Workbook object containing spreadsheet where each sheet contains dataset metadata of a given cluster.
    tmp_file : string
        Name of temporary file containing spreadsheet of cluster dataset metadata.

    Returns
    -------
    None.

    """
    
    wb.save(tmp_file)
    
def main(in_file,tmp_file):
    """
    Main function.

    Parameters
    ----------
    in_file : string
        Name of input file containing pickled dataframe only containing rows with clusters of interest.
    tmp_file : string
        Name of temporary file containing spreadsheet of cluster dataset metadata.
    clusters : list
        List of strings containing clusters of interest for creating summary spreadsheet.

    Returns
    -------
    None.

    """
    
    # load pickled dataframe file
    df = load_cleaned_dataframe(in_file)
    
    # split dataframe into dictionary by cluster and associated dataset metadata
    clusters_dict = split_dataframe_by_cluster(df)
    
    # create spreadsheet where each sheet contains datasets shared by a cluster
    wb = create_spreadsheet(clusters_dict)
    
    # save temporary spreadsheet to file
    save_temp_spreadsheet(wb, tmp_file)
    

if __name__ == '__main__':

    input_filename = str(snakemake.input['in_filename'])

    output_filename = str(snakemake.output['out_filename'])
    
    main(input_filename,output_filename)
