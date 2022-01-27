import pandas as pd
import numpy as np
import scipy.stats
import argparse

class meta_analysis:
    def __init__(self, b1, se1, b2, se2, chemid, chem_meta, output_prefix):
        self.b1 = b1
        self.se1 = se1
        self.b2 = b2
        self.se2 = se2
        self.chemid = chemid
        self.w1 = 0
        self.w2 = 0
        self.chem_meta = chem_meta
        self.output_prefix = output_prefix

        self.run()


    def run(self):
        '''
            This function will preform a 2 way meta anlysis of the metabolomics data, preforming the meta analysis via a
            inverse variance approach. User must input 4 arrays all of the same size.
        '''
        self.w1 = 1 / (self.se1 ** 2)
        self.w2 = 1 / (self.se2 ** 2)


        se_total = np.sqrt(1 / (self.w1 + self.w2))

        beta_total = (self.w1 * self.b1 + self.w2 * self.b2) / (self.w1 + self.w2)

        z_total = beta_total / se_total

        p_total = scipy.stats.norm.sf(abs(z_total)) * 2

        ci_start = beta_total - 1.96 * se_total
        ci_end = beta_total + 1.96 * se_total



        meta_df = pd.DataFrame(
            {
             'Biochemical':self.chem_meta['Biochemical'],
             'Super_Pathway': self.chem_meta['Super_Pathway'],
             'Sub_Pathway': self.chem_meta['Sub_Pathway'],
             'chemical_id':self.chemid,
             'beta_coef': beta_total,
             'Z_score': z_total,
             'p_value': p_total,
             'ci_start':ci_start,
             'ci_end':ci_end
             })
        return (meta_df)
        # output_filepath = f'{output_prefix}_meta.txt'
        # meta_df.to_csv(output_filepath, sep='\t', index=False)





if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s1', '--summ_stats_1')
    parser.add_argument('-s2', '--summ_stats_2')
    parser.add_argument('-o', '--output_prefix')

    user_args = parser.parse_args()

    output_prefix = str(user_args.output_prefix)
    study_1 = pd.read_csv(str(user_args.summ_stats_1), sep='\t')
    study_2 = pd.read_csv(str(user_args.summ_stats_2), sep='\t')

    #Filter based of pval inclusion threshold
    # study_1 = study_1[study_1['p_value'] < float(user_args.p_meta_threshold)].copy()
    # study_2 = study_2[study_2['p_value'] < float(user_args.p_meta_threshold)].copy()

    isec_chemicals = np.intersect1d(study_1['Chemical_ID'], study_2['Chemical_ID'])

    study_1 = study_1[study_1['Chemical_ID'].isin(isec_chemicals)].copy()
    study_2 = study_2[study_2['Chemical_ID'].isin(isec_chemicals)].copy()


    meta_df = meta_analysis(b1=study_1['beta_coef'].values, se1=study_1['standard_error'].values,
                  b2=study_2['beta_coef'].values, se2=study_2['standard_error'].values,
                  chemid=study_1['Chemical_ID'].values, chem_meta=study_1, output_prefix=output_prefix)





