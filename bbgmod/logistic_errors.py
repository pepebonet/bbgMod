import os
import click
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
from sklearn.linear_model import LogisticRegression


def get_data(freq_treat, freq_untreat):
    df_treat = pd.read_csv(freq_treat)
    df_treat['treatment'] = 1
    df_treat = df_treat.sample(n=2000000) #TODO remove and run for whole set

    df_untreat = pd.read_csv(freq_untreat)
    df_untreat['treatment'] = 0
    df_untreat = df_untreat.sample(n=2000000) #TODO remove and run for whole set

    data = pd.concat([df_treat, df_untreat]).reset_index(drop=True)
    
    X = data[data.columns[4:-1]].values
    Y = data[data.columns[-1]].values
    return data, X, Y 


def train_logistic_regression(X, y):
    return LogisticRegression(max_iter=1000).fit(X, y)


def test_logreg(logreg, x):
    predict = logreg.predict_proba(x)
    return predict[:, 0], predict[:, 1]


#TODO improve to exactly select what you want to test
def get_testset(data, cpg='cpg'):
    data_cpg = get_kmer_CpGs(data)
    
    if cpg == 'cpg':
        import pdb;pdb.set_trace()
        test = data_cpg[data_cpg['CpG_status'] == 1]
        Y_test = test['treatment'].values
    elif cpg == 'no_cpg':
        test = data_cpg[data_cpg['treatment'] == 1]
        Y_test = test['treatment'].values
    else:
        test = data_cpg[data_cpg['treatment'] == 1]
        Y_test = test[test.columns['CpG_status']].values

    return test[test.columns[4:-8]].values, Y_test


def get_test_sig_positions(data, df_sig):
    data = data[data['treatment'] == 1]
    df_sig['id'] = df_sig['CHROM'] + '_' + df_sig['pos'].astype(str)
    
    positions = data['Window'].str.split(':', expand=True)
    data['id'] = 'chr' + data['Ref'] + '_' + positions[2].astype(str)
    sig_test = pd.merge(df_sig, data, how='inner', on='id')

    Y_test = sig_test['treatment'].values
    X_test = sig_test[sig_test.columns[12:-1]]

    return X_test, Y_test


def split_cpg_data(df2):
    cpgs = df2[(df2['base3'] == 'C') & (df2['base4'] == 'G')]
    cpgs['CpG_status'] = 1
    non_cpgs = df2[(df2['base3'] != 'C') | (df2['base4'] != 'G')]
    non_cpgs['CpG_status'] = 0

    return pd.concat([cpgs, non_cpgs])


def get_kmer_CpGs(df):
    df['kmer'] = df['#Kmer'].apply(lambda x: list(x))
    df1 = pd.DataFrame(df['kmer'].values.tolist(), 
        columns=['base1', 'base2', 'base3', 'base4', 'base5'])
    df2 = df.join(df1) 

    return split_cpg_data(df2)


def plot_ROC (y_test, probas, fig_out, kn='Linear'):
    fpr, tpr, thresholds = roc_curve(y_test, probas)

    roc_auc = auc(fpr,tpr)
    plt.plot (fpr, tpr, 
        label = 'Logistic Regression (Area under ROC = %0.2f)'% roc_auc)
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Multiple features (5mC Ecoli)')
    plt.legend(loc="lower right")
    plt.savefig(fig_out)


#TODO select only kmers that interest the user. Focus for instance in Gs 
def main(treat_errors, untreat_errors, signif_pos, output):
    data, X, Y = get_data(treat_errors, untreat_errors)
    if signif_pos:
        df_sig = pd.read_csv(signif_pos,sep='\t')
        test_X, test_Y = get_test_sig_positions(data, df_sig)
    else:
        test_X, test_Y = get_testset(data)

    logreg = train_logistic_regression(X, Y)
    
    unm, mod = test_logreg(logreg, test_X)

    if signif_pos:
        import numpy as np
        print(np.where(mod > 0.8)[0].shape)
        import pdb;pdb.set_trace()
    else:
        out_fig = os.path.join(output, 'AUC_logreg_mod.png')
        plot_ROC(test_Y, mod, out_fig)


if __name__ == "__main__":
    main()