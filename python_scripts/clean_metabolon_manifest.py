import pandas as pd
import numpy as np

'''
This script is used to help create the prop_manifest.txt and tolsurf_manifest file for the infants
'''

if __name__ == '__main__':
    tolsurf_metabolon_table = pd.read_excel("data/input_data/tolsurf_manifest_full.xlsx")
    tolsurf_metabolon_table = tolsurf_metabolon_table[tolsurf_metabolon_table['Time Point']=='2nd']
    tolsurf_metabolon_table.to_csv("data/input_data/tolsurf_manifest.txt", sep='\t', index=False)

    # prop_metabolon_table = pd.read_excel("data/input_data/prop_manifest_sample2.xlsx")
    # print(prop_metabolon_table)
    prop_metabolon_table = pd.read_excel("data/input_data/UCSF-03-19VW CDT.XLSX","OsmoNormImpData")
    prop_manifest = prop_metabolon_table.iloc[0:37,12:].copy()

    prop_manifest = prop_manifest.T
    new_header = prop_manifest.iloc[0]  # grab the first row for the header
    prop_manifest = prop_manifest[1:]  # take the data less the header row
    prop_manifest.columns = new_header


    prop_manifest.insert(loc=0, column='CLIENT IDENTIFIER', value = prop_manifest.index)
    prop_manifest = prop_manifest.rename(columns={f'{prop_manifest.columns[-1]}': 'Group'})
    prop_manifest = prop_manifest[prop_manifest['Group'].notnull()]
    #
    prop_manifest = prop_manifest[prop_manifest['TIME POINT'].isin(['2nd'])]
    prop_manifest.to_csv("data/input_data/prop_manifest.txt", sep='\t', index=False)
    # print(prop_manifest)