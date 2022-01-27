import pandas as pd
import numpy as np
from scipy.stats import hypergeom
from scipy.stats import fisher_exact

#Convert to user input later
meta_results = pd.read_csv("results/indiv_assoc_results/tpn_tolsurf_prop_meta.txt", sep='\t')

imp_score_threshold = 30

chem_meta_info = pd.read_csv('data/input_data/metabolites_PercentMissing_tolsurf.txt', sep='\t')

# End of user input

#Correct for bonferonni, 652x296 652 individual test in addtion to 296 done from meta

meta_bonf_df = meta_results[meta_results['p_value'] < 0.05].copy()

biochem_tested_df = chem_meta_info[chem_meta_info['PercentMissing'] < imp_score_threshold].copy()
avil_subpathway = np.unique(biochem_tested_df['SUB_PATHWAY'])

num_tested = len(biochem_tested_df['BIOCHEMICAL'])
num_sign = len(meta_bonf_df['Biochemical'])

print(num_tested,num_sign)

pathway_enrichment_results = []
for pathway in avil_subpathway:
    subpath_subset = biochem_tested_df[biochem_tested_df['SUB_PATHWAY'] == pathway].copy()
    subpath_meta_res = meta_bonf_df[meta_bonf_df['Sub_Pathway'] == pathway].copy()
    super_pathway = np.unique(subpath_subset['SUPER_PATHWAY'])[0]

    if subpath_meta_res.empty:
        print("EMPTY DATAFRAME")
        continue

    if len(subpath_meta_res.index) == 1:
        print('Singular value, wont preform meta')
        continue


    subpath_total = len(subpath_subset['BIOCHEMICAL'])
    n_sub_sign = subpath_total - len(subpath_meta_res['Biochemical'])
    print(f"subpath_total: {subpath_total}.   n_sub_sing: {n_sub_sign}")


    conf_1 = n_sub_sign
    conf_2 = num_sign - n_sub_sign
    conf_3 = subpath_total - n_sub_sign
    conf_4 = num_tested - num_sign - conf_2

    cont_table = np.array([[conf_1, conf_3], [conf_2, conf_4]])

    oddsr, p = fisher_exact(cont_table, alternative='two-sided')
    row = [super_pathway, pathway, n_sub_sign, subpath_total, oddsr, p]
    pathway_enrichment_results.append(row)
    print(oddsr, p)

pathway_enrichment_df = pd.DataFrame(pathway_enrichment_results, columns=['Super Pathway','Sub Pathway', 'Number Signifigant',
                                                                          'Subpathway Total', 'Odds Ratio', 'P_value'])

pathway_enrichment_df.to_csv('results/indiv_assoc_results/meta_ora_pathway.csv', index=False)

print(pathway_enrichment_df[pathway_enrichment_df['P_value'] < 0.05])

