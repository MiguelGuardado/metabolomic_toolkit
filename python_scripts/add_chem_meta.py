import pandas as pd
import numpy as np
import argparse


'''
This script is used as a supplement to run_association_analysis.py. This script will be used to add information of the
chemical id used by metabolon. The metabolite_id we use for Prop and Tolsurf are the chemical ID provided by metabolon.
This script will take in the individual association analysis summary stats file, and add information on the list of 
candidate biochemicals, assigning their name, and the sub pathway and super pathway of the biochemical. 

'''

# def check_file_input(user_summ_stats, user_meta_file):
#     chem_id_names = ['metabolite_id', 'chemical_id']
#     chem_id_types_isec = np.intersect1d(chem_id_names, list(user_summ_stats.columns))
#     print(chem_id_types_isec)
#     if(len(chem_id_types_isec)>0):
#         print(user_summ_stats.metabolite_id.to_list())
#         print(type(user_summ_stats.metabolite_id.to_numpy()))
#         summ_stats_chem_id = [i.split('_')[0] for i in user_summ_stats.metabolite_id.to_numpy()]
#         print(list(summ_stats_chem_id))
#         exit(0)
#         chem_id_isec = np.intersect1d(summ_stats_chem_id, user_meta_file['CHEMICAL_ID'])
#         print(len(chem_id_isec))
#         print(len(user_summ_stats[chem_id_types_isec]))
#     else:
#         print("Error: Summary stats file that does not provide the column 'chemical_id','metabolite_id'")
#         exit(0)


def clean_chem_id(user_summ_stats, user_meta_file):
    user_summ_stats.dropna(inplace=True)
    summ_stats_chem_id = [int(i.split('_')[1]) for i in user_summ_stats.metabolite_id.to_numpy()]

    chem_id_isec = np.intersect1d(summ_stats_chem_id, user_meta_file.CHEMICAL_ID)
    user_summ_stats.insert(loc=0, column='Chemical_ID', value=chem_id_isec)
    user_summ_stats.drop(columns=['metabolite_id'], inplace=True)

    chem_meta_file = user_meta_file[user_meta_file['CHEMICAL_ID'].isin(chem_id_isec)]
    if (len(chem_id_isec) == len(summ_stats_chem_id)):
        user_summ_stats
        return(user_summ_stats, chem_meta_file)
    else:
        user_summ_stats = user_summ_stats[user_summ_stats['metabolite_id'].isin(chem_id_isec)]
        return(user_summ_stats, chem_meta_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--summ_stats_file')
    parser.add_argument('-c','--chem_meta_file')
    user_args = parser.parse_args()

    summ_stats_file = pd.read_csv(str(user_args.summ_stats_file), sep='\t')
    chem_meta_file = pd.read_csv(str(user_args.chem_meta_file), sep='\t')

    summ_stats_file, chem_meta_file = clean_chem_id(user_summ_stats=summ_stats_file, user_meta_file=chem_meta_file)

    summ_stats_file.insert(loc=0, column='Biochemical', value=chem_meta_file['BIOCHEMICAL'].to_numpy())
    summ_stats_file.insert(loc=1, column='Super_Pathway', value=chem_meta_file['SUPER_PATHWAY'].to_numpy())
    summ_stats_file.insert(loc=2, column='Sub_Pathway', value=chem_meta_file['SUB_PATHWAY'].to_numpy())

    output_prefix = str(user_args.summ_stats_file).split('.')[0]
    output_filepath = f'{output_prefix}_w_meta.txt'

    summ_stats_file.to_csv(output_filepath, index=False, sep='\t')


