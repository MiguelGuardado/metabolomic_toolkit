#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Recent Update: 01/05/22
@author: Miguel Guardado Miguel.Guardado@ucsf.edu
"""
import argparse
import os
import sys
import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.preprocessing import OneHotEncoder

"""
This script will be used to run logistic regression analysis for metabolomics data, but can be used in any quantitative
genetic data. The response variable will be a ordinal binary phenotype. 

To run this script, you will need to download conda, and install the 'metab_association_analysis.yml' environment. 
If you dont have conda, you can download it here: https://conda.io/projects/conda/en/latest/user-guide/install/index.html

To install the dep of this script, you will need to download the conda environment file `association_analysis.yml`.
$  conda env create -f association_analysis.yml
$  conda activate association_analysis

Once you have the conda environment installed, you will need to provide this script 3/4 different parameters 
to run. Assuming we have n samples as rows, and m columns as metabolites, the files you need to provide are as follows.
-p phenotype array size(Nx1)
-m metabolite matrix (NxM)
-c covariates matrix (Nxk) where k represents the number of different covariates
-o output prefix of the association analysis summary statistics file.

The output of the association analysis summary statistics will be reported as a text and csv file.
ex. -o tmp will result in two files names tmp_summ_stats.csv and tmp_summ_stats.txt. Both identical, 
but one is easier for visualization, and the other is easier for data analysis input. 
"""


class metabolite_association_analysis:
    def __init__(self, m_matrix, p_arr, output_prefix, cov_matrix=None):
        self.m_matrix_filepath = m_matrix
        self.p_arr_filepath = p_arr
        self.output_prefix_filepath = output_prefix
        self.cov_matrix_filepath = cov_matrix
        self.m_matrix = []
        self.p_arr = []
        self.cov_matrix = []


        self.run()


    def load_files(self):
        self.m_matrix = pd.read_table(self.m_matrix_filepath, sep= "\t")
        self.p_arr = pd.read_table(self.p_arr_filepath, sep='\t')

        if self.cov_matrix_filepath is not None:
            self.cov_matrix = pd.read_table(self.cov_matrix_filepath, sep= "\t")


    def check_matrix(self):
        """
        This will check the dimensions of the metabolite, phenotype, and covariance matrix(if provided)
        We use the metabolite matrix to identify the number of samples and individuals recorded. We require the samples
        to be seperated by rows, and metabolites/covariates be represented by columns.
        Assuming we have N individuals over m metabolite samples, we will check the dimensions of each as follows...

        Metabolite matrix - Get N samples (number of rows) and M biochemicals (number of columns) [N x M]
        Phenotype matrix - Check if a [N x 1] matrix is provided
        Covariance matrix - Check is a [N x k] matrix is provided, where k represents how many covariates are provided
        :return: Nothing, will break from script if any dimension is not correct.
        """
        #First we check if the first column in each of the 3 inputs are equal by sample number

        matrix_id = self.m_matrix.iloc[:,0].to_numpy()
        p_id = self.p_arr.iloc[:,0].to_numpy()

        if not np.array_equal(matrix_id, p_id):
            print("Sample ID from Metabolite matrix does not equal the Sample ID of the phenotype file.")
            exit(0)

        N = list(self.m_matrix.shape)[0]
        M = list(self.m_matrix.shape)[1] - 1


        # First will check phenotype matrix
        p_num_indivs = list(self.p_arr.shape)[0]
        if (N != p_num_indivs):
            print(f" Metabolite Dimensions: {N} individuals over {M} metabolite samples")
            print(f" Phenotype Dimensions: {p_num_indivs} individuals")
            sys.exit("Dimension Error: Dimension of phenotype array does not equal number of individuals recorded in metabolie matrix, please check the dimensions")

        if self.p_arr.shape[1] != 2:
            print("No phenotype file found, expecting a file of two columns [SAMPLE ID, PHENOTYPE]")

        if self.cov_matrix_filepath is None:
            return

        cov_id = self.p_arr.iloc[:,0].to_numpy()
        if not np.array_equal(matrix_id, cov_id):
            print("Sample ID from Metabolite matrix does not equal the Sample ID of the phenotype file.")
            exit(0)

        cov_num_indivs = list(self.cov_matrix.shape)[0]
        if (N != cov_num_indivs):
            print(f" Metabolite Dimensions: {N} individuals over {M} metabolite samples")
            print(f" Covariate Dimensions: {cov_num_indivs} individuals ")
            sys.exit("Dimension Error: Dimension of covariance does not equal number of individuals recorded in metabolie matrix, please check the dimensions")

    def check_categorical_data(self):
        """
            This will check the inputted covariates matrix to determine if there are any non binary categorical inputted
            recorded. We will check each variable 1 by 1
        :return:
        """
        enc = OneHotEncoder(handle_unknown='ignore')
        covariates = self.cov_matrix.columns[1:]
        for covariate in covariates:
            cov_col = self.cov_matrix[covariate]
            #print(covariate, cov_col.dtype.name)

            if np.issubdtype(cov_col.dtype.name, np.number):
                pass

            if np.issubdtype(cov_col.dtypes.name, np.object0):

                unique_category = cov_col.unique()
                #This will convert binary data to 0/1 variable
                if(len(unique_category)==2):
                    bin_op = cov_col[0]
                    tmp_var = np.where(cov_col == f'{bin_op}', 1, 0)
                    self.cov_matrix.drop([covariate], axis=1, inplace=True)
                    self.cov_matrix[covariate] = tmp_var

                else:
                    cov_col = cov_col.values.reshape(-1, 1)
                    enc.fit(cov_col)
                    onehotlabels_names = enc.categories_
                    onehotlabels_names = list(map(( lambda x: f'{covariate}' + x), onehotlabels_names))[0]

                    onehotlabels = enc.transform(cov_col).toarray()

                    self.cov_matrix.drop([covariate], axis=1, inplace=True)
                    #self.cov_matrix[[onehotlabels_names]] = onehotlabels
                    self.cov_matrix = self.cov_matrix.join(pd.DataFrame(onehotlabels, index=self.cov_matrix.index,
                                                          columns=onehotlabels_names))


    def run(self):
        """
        main class of metabolite association analysis, this will load in metabolite matrix, phenotype matrix, and
        covariates of each simulations required.

        :return: outputs a summary stats file for the output created. will be outputted based on what the user enters
        in the -o flag.
        """
        #first we will load and check the input files for the metabolite analysis
        self.load_files()
        self.check_matrix()
        self.check_categorical_data()


        metabolite_id = self.m_matrix.columns[1:]
        summ_stats_table_output = []

        for metabolite in metabolite_id:
            cols = [self.m_matrix.columns[0], metabolite]
            x = pd.merge(self.m_matrix[cols], self.cov_matrix, on=self.m_matrix.columns[0])
            x = x.iloc[: , 1:]
            y = self.p_arr.iloc[:,1]

            if(len(np.unique(x[metabolite])) < 10):
                print(f"not enough unqiue metabolite values {metabolite}")
                row = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
                summ_stats_table_output.append(row)
                continue

            #  Fit logistic regression model for single metabolite
            log_reg = sm.Logit(y,x).fit()

            #Get beta, se, p_val, z_score, ci_start, and ci_end for each estimate.
            beta = log_reg.params[0]

            #Standard error
            cov = log_reg.cov_params()
            std_err = np.sqrt(np.diag(cov))

            p_val = log_reg.pvalues[0]

            #Calculate Z score from the standard error of the beta coefficient estimate.
            z_values = log_reg.params / std_err
            z = z_values[0]

            #Caluclate 95% confident interval
            ci_est = log_reg.conf_int(alpha=0.05, cols=None)
            ci_start = ci_est[0][0]
            ci_end = ci_est[1][0]

            #Finally we will calulcate the fold change between the two groups.
            # Fold Change = mean ( group with phenotype status) / mean (group not with phenotype status)
            x['phenotype'] = y
            group1_mean = np.mean(x[x['phenotype'] == 1][metabolite])
            group0_mean = np.mean(x[x['phenotype'] == 0][metabolite])
            fold_change = group1_mean / group0_mean

            row = [metabolite, fold_change, beta, std_err[0], z, p_val, ci_start, ci_end]
            summ_stats_table_output.append(row)


        #Finally after the for loop, we convert the output list to a df and output based on user input output_prefix.
        summ_stats_col_names = ['metabolite_id','fold_change','beta_coef','standard_error',
                                'z_score','p_value','ci_start','ci_end']
        summ_stats_table_output = pd.DataFrame(summ_stats_table_output, columns=summ_stats_col_names)

        txt_filepath = f'{self.output_prefix_filepath}_summ_stats.txt'
        summ_stats_table_output.to_csv(txt_filepath, sep='\t', index=False)

def load_args():
    """
    This function is used to read in and initalize the user parameters
    entered from the command line.
    This wont check the parameters in but just reads them in the filepaths.
    metabolite_assocation_analysis class will check the files and their dimensions.

    USER MUST INPUT THESE PARAMETERS to make any of the simulations run:
        -m metabolite_matrix
        -p phenotype_array
        -o output prefix

    Covariance matrix does not need to be included, but should!

    This function is used outside of the metabolite_association_analysis class. not apart of the class's internal
    function.

    :return: returns the users command line inputs for each of the parameters defined
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--metabolite_matrix')
    parser.add_argument('-p', '--phenotype_array')
    parser.add_argument('-c', '--covariance_matrix')
    parser.add_argument('-o', "--output_prefix")

    #This will check if the user inputted these three filename, will not check the contents,
    user_args = parser.parse_args()
    if user_args.metabolite_matrix is None:
        sys.exit("InputError: metabolite_matrix not inputted")
    if user_args.phenotype_array is None:
        sys.exit("InputError: phenotype array no inputted")
    if user_args.output_prefix is None:
        sys.exit("InputError: Output prefix not inputted, please try again")

    #Now we will check if these files exist
    if not os.path.isfile(f'{user_args.metabolite_matrix}'):
        sys.exit("FileError: metabolite matrix file not found")
    if not os.path.isfile(f'{user_args.phenotype_array}'):
        sys.exit("FileError: phenotype array not found")

    m_matrix = os.path.abspath(f"{user_args.metabolite_matrix}")
    p_arr = os.path.abspath(f"{user_args.phenotype_array}")
    output_prefix = f'{os.getcwd()}/{user_args.output_prefix}'

    if os.path.isfile(f'{user_args.covariance_matrix}'):
        cov_matrix = os.path.abspath(f'{user_args.covariance_matrix}')
    else:
        cov_matrix = None

    return m_matrix,p_arr,cov_matrix,output_prefix


if __name__ == '__main__':
    """
    Quick start to run this code once the conda enviroment is downloaded. 
    
    python run_association_analysis.py 
    -m ../tpn_metabolite_data/tolsurf_post_tpn_m_30missing.txt 
    -p ../tpn_metabolite_data/tolsurf_post_tpn_p.txt 
    -c ../tpn_metabolite_data/tolsurf_post_tpn_cov.txt 
    -o tpn_tolsurf
    """

    m_matrix, p_arr, cov_matrix, output_prefix = load_args()
    metabolite_association_analysis(m_matrix=m_matrix,
                                    p_arr=p_arr,
                                    cov_matrix=cov_matrix,
                                    output_prefix=output_prefix)
