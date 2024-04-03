# -*- coding: utf-8 -*-
"""
Created on Thurs Mar 28 13:08:53 2024

@author: AbnerBogan
"""

import pickle
import pandas as pd

def load_data_sharing_report(in_file):
    """
    Load in data sharing report from file.

    Parameters
    ----------
    in_file : string
        Name of input file containing data sharing report.

    Returns
    -------
    df : pandas dataframe
        Dataframe containing shared dataset metdata.

    """
    
    # load file via unpickling
    with open(in_file, 'rb') as file:
        df = pickle.load(file)
        
    return df

def remove_columns(df,columns=['legacy']):
    """
    Remove defined columns from dataframe.

    Parameters
    ----------
    df : pandas dataframe
        Dataframe containing shared dataset metdata. 
    columns : list, optional
        List of column headers to remove. The default is ['legacy'].

    Returns
    -------
    df_dropped : pandas dataframe
        Dataframe (minus dropped columns) containing shared dataset metdata.

    """
    
    df_dropped = df.drop(columns,axis=1)
    
    return df_dropped

def standardize_timestamp_format(df, date_columns=['dateCreated', 'datePublished']):
    """
    Standardize the format and type of timestamps in input dataframe to 'YYYY-MM-DD HH:MM:SS'.

    Parameters
    ----------
    df : pandas dataframe
        Dataframe containing shared dataset metdata.
    date_columns : list, optional
        Column headers containing dates. The default is ['dateCreated', 'datePublished'].

    Returns
    -------
    df : pandas dataframe
        Dateframe containing shared dataset metdata with a standardized date format.

    """
    
    # iterate over list of date column headers
    for column in date_columns:
        # some timestamps include microseconds; remove from string by truncating values after period
        df[column] = df[column].apply(lambda x: x.split('.')[0] if pd.notnull(x) and '.' in x else x)
        # convert the truncated string to a datetime object
        df[column] = pd.to_datetime(df[column])
        
    return df

def filter_datasets_by_clusters(df,clusters=['Bedrock','Big Data','CINET','Coastal',
                                             'Drylands','DUST^2','Dynamic Water','GeoMicroBio','Urban']):
    """
    Filter datasets by clusters of interest.

    Parameters
    ----------
    df : pandas dataframe
        Dataframe containing shared dataset metdata.
    clusters : list, optional
        List of strings containing clusters of interest. The default is ['Bedrock','Big Data',
        'CINET','Coastal','Drylands','DUST^2','Dynamic Water','GeoMicroBio','Urban'].

    Returns
    -------
    df_filtered : pandas dataframe
        Dataframe containing shared dataset metdata from a subset of clusters.

    """
    
    # modify DUST^2 as it is written in API endpoint and add '\\' for subsequent 
    # regex pattern search with special characters
    try:
        clusters[clusters.index('DUST^2')] = 'Dust\\^2'
    except:
        pass
    
    # create regex pattern from list of clusters for filtering
    pattern = '|'.join(clusters)
    
    # run matching based off of regex pattern and include rows where matching is observed
    df_filtered = df[df['clusters'].str.contains(pattern, na=False)]
    
    return df_filtered

def fill_missing_dates(df):
    """
    Fill in any missing 'dateCreated' values with corresponding 'datePublished' value.

    Parameters
    ----------
    df : pandas dataframe
        Dataframe containing shared dataset metdata.

    Returns
    -------
    df : pandas dataframe
        Dataframe containing shared dataset metdata with filled in values in 'dateCreated' column.

    """
    
    # find rows where 'dateCreated' column is empty
    blank_rows = df['dateCreated'].isna()
    
    # Replace empty values in 'dateCreated with values from 'datePublished'
    df.loc[blank_rows, 'dateCreated'] = df.loc[blank_rows, 'datePublished']
    
    return df

def save_cleaned_dataframe(df,out_file):
    """
    Save cleaned dataframe to file via pickling.

    Parameters
    ----------
    df : pandas dataframe
        Cleaned dataframe containing shared dataset metdata by subset of clusters.
    out_file : string
        Name of output file containing pickled dataframe.

    Returns
    -------
    None.

    """
    
    # pickle (serialize) the data to a file
    with open(out_file, 'wb') as file:
        pickle.dump(df, file)
        
def main(in_file,out_file):
    """
    Main function.

    Parameters
    ----------
    in_file : string
        Name of input file containing pickled dataframe generated through API endpoint.
    out_file : string
        Name of output file containing pickled dataframe only containing rows with clusters of interest.
    clusters : list
        List of strings containing clusters of interest for creating summary spreadsheet.

    Returns
    -------
    None.

    """
    
    # load data sharing report from file
    df = load_data_sharing_report(in_file)
    
    # drop unimportant columns from dataframe
    df_dropped = remove_columns(df)
    
    # standardize format/type of timestamps indicating time of data sharing
    df_standardized = standardize_timestamp_format(df_dropped)
    
    # filter dataframe to include datasets shared by a subset of clusters
    df_filtered = filter_datasets_by_clusters(df_standardized)
    
    # fill in any missing values in 'dateCreated' field with 'datePublished'
    df_cleaned = fill_missing_dates(df_filtered)
    
    # save dataframe to file via pickling
    save_cleaned_dataframe(df_cleaned,out_file)
    


if __name__ == '__main__':

    input_filename = str(snakemake.input['in_filename'])

    output_filename = str(snakemake.output['out_filename'])
   
    main(input_filename,output_filename)
    
    
    
    