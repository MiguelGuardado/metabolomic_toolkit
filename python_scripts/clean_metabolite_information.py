"""


"""

import pandas as pd
import numpy as np

# Load in infant manifest files, chem meta information, and metabolite abudnance values.
manifest_tolsurf = pd.read_csv("data/input_data/tolsurf_manifest.txt", sep='\t')
manifest_prop = pd.read_csv("data/input_data/prop_manifest.txt", sep='\t')

m_matrix_pre_tolsurf = pd.read_table('data/input_data/tolsurf_metabolite_matrix_pre.txt', sep='\t')
m_matrix_pre_prop = pd.read_table('data/input_data/prop_metabolite_matrix_pre.txt', sep='\t')

m_matrix_post_tolsurf = pd.read_table('data/input_data/tolsurf_metabolite_matrix_post.txt', sep='\t')
m_matrix_post_prop = pd.read_table('data/input_data/prop_metabolite_matrix_post.txt', sep='\t')

m_matrix_fc_tolsurf = pd.read_table('data/input_data/tolsurf_metabolite_matrix_fc.txt', sep='\t')
m_matrix_fc_prop = pd.read_table('data/input_data/prop_metabolite_matrix_fc.txt', sep='\t')

chem_meta_tolsurf = pd.read_table('data/input_data/metabolites_PercentMissing_tolsurf.txt', sep='\t')
chem_meta_prop = pd.read_table('data/input_data/metabolites_PercentMissing_prop_noanchors.txt', sep='\t')

##  We first will extract the metabolomic matrices we will use for our individual association analysis. We will
## include metabolites that have metabolomic information available for 30% of the samples.

chem_meta_tolsurf_m30 = chem_meta_tolsurf[chem_meta_tolsurf['PercentMissing'] < 30]
m_avil_tolsurf = ["X_" + str(s) for s in chem_meta_tolsurf_m30['CHEMICAL_ID']]

tolsurf_sample_id = "Tolsurf_" + manifest_tolsurf['Unique Sample ID'].astype(str)
prop_sample_id = manifest_prop['CLIENT IDENTIFIER'].astype(str)


tolsurf_pre_m30_matrix = m_matrix_pre_tolsurf[m_avil_tolsurf]
tolsurf_pre_m30_matrix.insert(loc=0, column='Sample_ID', value=tolsurf_sample_id)
tolsurf_pre_m30_matrix.to_csv("data/indiv_assoc_data/tolsurf_m30_pre_metabolite_matrix.txt", sep='\t', index=False)

tolsurf_post_m30_matrix = m_matrix_post_tolsurf[m_avil_tolsurf]
tolsurf_post_m30_matrix.insert(loc=0, column='Sample_ID', value=tolsurf_sample_id)
tolsurf_post_m30_matrix.to_csv("data/indiv_assoc_data/tolsurf_m30_post_metabolite_matrix.txt", sep='\t', index=False)

tolsurf_fc_m30_matrix = m_matrix_fc_tolsurf[m_avil_tolsurf]
tolsurf_fc_m30_matrix.insert(loc=0, column='Sample_ID', value=tolsurf_sample_id)
tolsurf_fc_m30_matrix.to_csv("data/indiv_assoc_data/tolsurf_m30_fc_metabolite_matrix.txt", sep='\t', index=False)

#Prop Cohort
chem_meta_prop_m30 = chem_meta_prop[chem_meta_prop['PercentMissing'] < 30]
m_avil_prop = ["X_" + str(s) for s in chem_meta_prop_m30['CHEMICAL_ID']]

prop_pre_m30_matrix = m_matrix_pre_prop[m_avil_prop]
prop_pre_m30_matrix.insert(loc=0, column='Sample_ID', value=prop_sample_id)
prop_pre_m30_matrix.to_csv("data/indiv_assoc_data/prop_m30_pre_metabolite_matrix.txt", sep='\t', index=False)


prop_post_m30_matrix = m_matrix_post_prop[m_avil_prop]
prop_post_m30_matrix.insert(loc=0, column='Sample_ID', value=prop_sample_id)
prop_post_m30_matrix.to_csv("data/indiv_assoc_data/prop_m30_post_metabolite_matrix.txt", sep='\t', index=False)


prop_fc_m30_matrix = m_matrix_fc_prop[m_avil_prop]
prop_fc_m30_matrix.insert(loc=0, column='Sample_ID', value=prop_sample_id)
prop_fc_m30_matrix.to_csv("data/indiv_assoc_data/prop_m30_fc_metabolite_matrix.txt", sep='\t', index=False)


## We additionally will attach the sample_id to the metabolite matrix, for input for run_pca.py
# tolsurf_pre_m30_matrix.insert(loc=0, column='Sample_ID', value=tolsurf_sample_id)
# tolsurf_post_m30_matrix.insert(loc=0, column='Sample_ID', value=tolsurf_sample_id)
# tolsurf_fc_m30_matrix.insert(loc=0, column='Sample_ID', value=tolsurf_sample_id)
#
# prop_pre_m30_matrix.insert(loc=0, column='Sample_ID', value=prop_sample_id)
# prop_post_m30_matrix.insert(loc=0, column='Sample_ID', value=prop_sample_id)
# prop_fc_m30_matrix.insert(loc=0, column='Sample_ID', value=prop_sample_id)
#
#
# tolsurf_pre_m30_matrix.to_csv("data/pca_data/tolsurf_pre_m30_pca_input.txt", sep='\t', index=False)
# tolsurf_post_m30_matrix.to_csv("data/pca_data/tolsurf_post_m30_pca_input.txt", sep='\t', index=False)
# tolsurf_fc_m30_matrix.to_csv("data/pca_data/tolsurf_fc_m30_pca_input.txt", sep='\t', index=False)
#
#
# prop_pre_m30_matrix.to_csv("data/pca_data/prop_pre_m30_pca_input.txt", sep='\t', index=False)
# prop_post_m30_matrix.to_csv("data/pca_data/prop_post_m30_pca_input.txt", sep='\t', index=False)
# prop_fc_m30_matrix.to_csv("data/pca_data/prop_fc_m30_pca_input.txt", sep='\t', index=False)


##  Now we will extract the phenotypic information for BPD, we will test both BPD at 36 weeks, and BPD and 40 weeks.
##  We will only extract a singular array for each phenotypes being used.
tpn_status_tolsurf = manifest_tolsurf['TPN sample 2 status'].replace({'yes': 1, 'no': 0}).to_frame()
tpn_status_tolsurf.insert(loc=0, column='Sample_ID', value=tolsurf_sample_id)
tpn_status_tolsurf.columns = ['Sample_ID','TPN_status']
tpn_status_tolsurf.to_csv('data/indiv_assoc_data/tpn_status_tolsurf.txt', index=False, sep='\t')

tpn_status_prop = manifest_prop['RECENTLY ON TPN'].replace({'yes': 1, 'no': 0}).to_frame()
tpn_status_prop.insert(loc=0, column='Sample_ID', value=prop_sample_id)
tpn_status_prop.columns = ['Sample_ID','TPN_status']
tpn_status_prop.to_csv('data/indiv_assoc_data/tpn_status_prop.txt', index=False, sep='\t')


## Finally we will obtain the covariate information we want to correct for in our individual assocation analysis.
## Covariates we will correct for:
# - Birth Weight
# - Sex
# - Maternal Race/Ethnicity
# - Infant on Corticosteroid
# - Infant on TPN


tolsurf_cov = manifest_tolsurf[['Sex','Birth Weight', 'Race / Ethnicity','Days On Corticosteroid Sample 2']].copy()
tolsurf_cov.columns = ['Sex','Birth_Weight', 'Maternal_Race_Ethnicity','Days_On_Corticosteroid_Sample_2']
tolsurf_cov['Sex'].replace({'Male': 1, 'Female': 0}, inplace=True)
tolsurf_cov.insert(loc=0, column='Sample_ID', value=tolsurf_sample_id)
tolsurf_cov['Days_On_Corticosteroid_Sample_2'] = np.where(tolsurf_cov['Days_On_Corticosteroid_Sample_2'].values > 1, 1, 0)

tolsurf_cov.to_csv("data/indiv_assoc_data/tolsurf_cov.txt", sep='\t', index=False)


prop_cov = manifest_prop[['SEX','BIRTH WEIGHT G','MOTHERS RACE', 'CORTICOSTEROID DAY 28']].copy()
prop_cov.columns = ['Sex','Birth_Weight', 'Maternal_Race_Ethnicity','Days_On_Corticosteroid_Sample_2']
prop_cov['Sex'].replace({'M': 1, 'F': 0}, inplace=True)
prop_cov.insert(loc=0, column='Sample_ID', value=prop_sample_id)
prop_cov['Maternal_Race_Ethnicity'].replace({'White': 'NHW', 'Black': 'AA', 'Hispanic':'HL'}, inplace=True)
prop_cov.to_csv("data/indiv_assoc_data/prop_cov.txt", sep='\t', index=False)

#Finally for our roc analysis, we will merge the two covariate files together.
frames = [tolsurf_cov, prop_cov]
joint_cov = pd.concat(frames)
print(joint_cov)
joint_cov.to_csv("data/roc_data/joint_cov.txt", sep='\t', index=False)