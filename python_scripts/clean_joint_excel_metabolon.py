import pandas as pd
import numpy as np

'''
    This file will clean the metabolite matrix produced by metabolon for the joint metabolite abundance matrix,
    along with thier meta information. 
'''

if __name__ == '__main__':
    joint_measurment_filepath = 'data/input_data/prop_merged_metabolon_manifest.XLSX'
    joint_imputed_data = pd.read_excel(joint_measurment_filepath, "MergedScaledImpData")

    chem_manifest = joint_imputed_data.iloc[:,0:12].copy()
    chem_manifest.dropna(axis=0, how='all', inplace=True)
    new_header = chem_manifest.iloc[0]  # grab the first row for the header
    chem_manifest = chem_manifest[1:]  # take the data less the header row
    chem_manifest.columns = new_header  # set the header row as the df header
    chem_manifest["CHEMICAL ID"] = "X_" + chem_manifest["CHEMICAL ID"].astype(str)


    joint_imputed_data = joint_imputed_data.iloc[:,12:]


    column_names = np.concatenate([joint_imputed_data.iloc[0:15,0].astype(str),
                                   chem_manifest["CHEMICAL ID"].astype(str)])
    column_names = np.insert(column_names, 15, 'Group')

    joint_imputed_data.drop(columns=joint_imputed_data.columns[0], axis=1, inplace=True)
    joint_imputed_data = joint_imputed_data.reindex().T
    joint_imputed_data.columns = column_names
    joint_imputed_data.insert(loc=0, column='Sample_ID', value=joint_imputed_data.index.values)
    joint_imputed_data.to_csv("data/roc_data/joint_cleaned_manifest.csv", index=False)

    cols_to_keep = chem_manifest["CHEMICAL ID"].astype(str).tolist()
    tolsurf_sample_id = "Tolsurf_" + joint_imputed_data['Sample_ID'][0:171].astype(str)
    prop_sample_id = joint_imputed_data['Sample_ID'][171:].astype(str)
    sample_id = pd.concat([tolsurf_sample_id, prop_sample_id])


    joint_sample_1 = joint_imputed_data[joint_imputed_data['TIME POINT'] == '1st'].copy()
    joint_s1_matrix = joint_sample_1[cols_to_keep]
    joint_s1_matrix.insert(loc=0, column='Sample_ID', value=sample_id)
    joint_s1_matrix.to_csv("data/roc_data/joint_metabolite_matrix_tp1.txt", sep='\t', index=False)

    print(joint_s1_matrix)

    joint_sample_2 = joint_imputed_data[joint_imputed_data['TIME POINT']=='2nd'].copy()
    joint_s2_matrix = joint_sample_2[cols_to_keep]
    joint_s2_matrix.insert(loc=0, column='Sample_ID', value=sample_id)
    joint_s2_matrix.to_csv("data/roc_data/joint_metabolite_matrix_tp2.txt", sep='\t', index=False)

    true_tpn_status = joint_sample_2[['RECENTLY ON TPN']].copy()
    true_tpn_status['TPN_status'] = joint_sample_2['RECENTLY ON TPN'].replace({'yes': 1, 'no': 0})
    true_tpn_status.drop(columns='RECENTLY ON TPN', axis=1, inplace=True)
    true_tpn_status.insert(loc=0, column='Sample_ID', value=sample_id)
    true_tpn_status.to_csv('data/roc_data/true_joint_tpn_status.txt', sep='\t', index=False)

