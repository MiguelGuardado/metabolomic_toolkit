import numpy as np
import pandas as pd
import argparse
from sklearn.model_selection import KFold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_auc_score

def subset_testing_metabolites():
    pass


def check_arguments(user_args):
    print(user_args)
    pass


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--metabolite_matrix')
    parser.add_argument('-c', '--covariates')
    parser.add_argument('-l', '--testing_metabolites')
    parser.add_argument('-t', '--true_status')
    parser.add_argument('-o', '--output_prefix')

    user_args = parser.parse_args()
    metab_matrix = pd.read_csv(str(user_args.metabolite_matrix), sep='\t')
    true_status = pd.read_csv(str(user_args.true_status), sep='\t')
    metab_to_test = np.loadtxt(str(user_args.testing_metabolites), dtype='str')
    joint_cov = pd.read_csv(str(user_args.covariates), sep='\t')

    #To make my life easier, I convert Maternal Race/Ethnicity to a numerical catergorical data.
    joint_cov['Maternal_Race_Ethnicity'].replace({'NHW': 0, 'AA': 1, 'HL': 2}, inplace=True)
    metab_to_test_matrix = metab_matrix[metab_to_test]

    X = pd.concat([metab_to_test_matrix.reset_index(drop=True),joint_cov.reset_index(drop=True)], axis=1)
    y = true_status['TPN_status'].values



    k = 10
    kf = KFold(n_splits=k, random_state=None)
    model = LogisticRegression(solver='liblinear')

    acc_score = []
    auc_score = []

    for train_index, test_index in kf.split(X):
        X_train, X_test = X.iloc[train_index, :], X.iloc[test_index, :]
        y_train, y_test = y[train_index], y[test_index]
        print("-------------------------------------")

        clf = model.fit(X_train, y_train)
        y_score = clf.decision_function(X_test)
        pred_values = clf.predict(X_test)
        roc_auc = roc_auc_score(y_test, y_score)
        print(roc_auc)
        auc_score.append(roc_auc)
        print("-------------------------------------")

        acc = accuracy_score(pred_values, y_test)
        acc_score.append(acc)

    avg_acc_score = sum(acc_score) / k
    avg_roc_score = sum(auc_score) / k

    print(f'Avg accuracy : {avg_acc_score}')
    print(f'ROC area under the curve accuracy : {avg_roc_score}')