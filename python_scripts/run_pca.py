import pandas as pd
import argparse
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import umap

"""
This script will is a command line python script which will preform PCA and uMAP analysis from the 
resulting input matrix. This script will use scikit learn to preform the pca analysis. 

"""


def find_pca_dim_umap(pc_matrix, eigen_cumsum):
    """
    Qucik filter to determine the number of components to use for pca plot initialization. Either I will extract the top
    10 dimensions, or will extract the number of dimensions that explains the variance up to 90%.
    :param pc_matrix: pca component matrix
    :param eigen_cumsum: array of eigenvalue cumulative sum for each dimension of the pc_matrix.
    :return: pc_matrix
    """
    if(eigen_cumsum[10] > 0.90):
        for i in range(0,10):
            if(eigen_cumsum[i] > 0.90):
                return pc_matrix.iloc[:,0:i+1]
    else: # Now we will just extract the first ten dimensions of the pca_matrix
        return pc_matrix.iloc[:,0:11]



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--metabolite_matrix')
    parser.add_argument('-o', "--output_prefix")
    user_args = parser.parse_args()


    raw_metabolite_data = pd.read_table(f'{user_args.metabolite_matrix}')
    sample_id = raw_metabolite_data.iloc[:,0]
    raw_metabolite_data.drop(columns=raw_metabolite_data.columns[0], axis=1, inplace=True)
    raw_metabolite_data = StandardScaler().fit_transform(raw_metabolite_data)

    #Returns the smallest dimension of the metabolite matrix, used as the number of components that are generated
    num_comp = min(raw_metabolite_data.shape)

    #Set the number of components to be the number of samples from the metabolite matrix
    pca = PCA(n_components= num_comp)
    principal_components = pca.fit_transform(raw_metabolite_data)
    eigen_cumsum = pca.explained_variance_ratio_.cumsum()
    pc_cols = list(range(1,num_comp+1))
    pc_cols = ["PC" + str(pc) for pc in pc_cols]

    pca_df = pd.DataFrame(data=principal_components, columns=pc_cols)
    pca_df.insert(loc=0, column='SampleID', value=sample_id)
    pca_df.to_csv(f"{user_args.output_prefix}_pca.txt", sep='\t', index=False)

    #Now with the PCA components, we will additionally preform uMAP analysis of the data.
    # pca_df = find_pca_dim_umap(pca_df, eigen_cumsum)
    # reducer = umap.UMAP()
    # # First we will scale the data, which is done from the manuel
    # embedding = reducer.fit_transform(pca_df)
    # umap_df = pd.DataFrame(list(zip(sample_id, embedding[:,0],embedding[:,1])), columns=['SampleID','UMAP1','UMAP2'])
    # umap_df.to_csv(f"{user_args.output_prefix}_umap.txt", sep='\t', index=False)

