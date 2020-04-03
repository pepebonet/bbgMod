#!/usr/bin/env python3
import os
import click
import pandas as pd

import plots as pl


def get_data(freq_treat, freq_untreat, cpg_ecoli):
    df_treat = pd.read_csv(freq_treat)

    if cpg_ecoli:
        df_treat['treatment'] = 't'
    else:
        treatment = freq_treat.rsplit('/')[-2]
        df_treat['treatment'] = '{}'.format(treatment)

    df_untreat = pd.read_csv(freq_untreat)
    df_untreat['treatment'] = 'unt'

    return pd.concat([df_treat, df_untreat]).reset_index(drop=True)


def split_cpg_data(df2):
    cpgs = df2[(df2['base3'] == 'C') & (df2['base4'] == 'G')]
    cpgs['ID'] = cpgs['treatment'] + '_CpG'
    non_cpgs = df2[(df2['base3'] != 'C') | (df2['base4'] != 'G')]
    non_cpgs['ID'] = non_cpgs['treatment'] + '_no_CpG'

    return pd.concat([cpgs, non_cpgs])


def get_kmer_CpGs(df):
    df['kmer'] = df['#Kmer'].apply(lambda x: list(x))
    df1 = pd.DataFrame(df['kmer'].values.tolist(), 
        columns=['base1', 'base2', 'base3', 'base4', 'base5'])
    df2 = df.join(df1) 

    return split_cpg_data(df2)


def get_single_CpGs(df):
    cpgs = df[df['base'] == 'C']
    cpgs['ID'] = cpgs['treatment'] + '_C'
    non_cpgs = df[df['base'] != 'C']
    non_cpgs['ID'] = non_cpgs['treatment'] + '_no_C'

    return pd.concat([cpgs, non_cpgs])

    
def merge_sig_positions(df, df_sig, kmer_study, output):
    df_sig['id'] = df_sig['CHROM'] + '_' + df_sig['pos'].astype(str)
    
    if kmer_study:
        positions = df['Window'].str.split(':', expand=True)
        df['id'] = 'chr' + df['Ref'] + '_' + positions[2].astype(str)
    else:
        df['id'] = 'chr' + df['#Ref'] + '_' + df['pos'].astype(str)

    df = pd.merge(df_sig, df, how='outer', on='id')
    no_sig = df[(df['strand'] != '+') & (df['strand'] != '-')]
    no_sig['treatment'] = no_sig['treatment'] + '_no_sig'
    sig = df[(df['strand'] == '+') | (df['strand'] == '-')]
    sig['treatment'] = sig['treatment'] + '_is_sig'
    df = pd.concat([no_sig, sig])
    df = df.rename(columns={'treatment':'ID'})

    if kmer_study:
        pl.kmer_advanced_plots(df, output, label='sig', stat=True)
    else:
        pl.advanced_plots(df, output, label='sig', stat=True)

# ------------------------------------------------------------------------------
# ARGPARSER
# ------------------------------------------------------------------------------
@click.command(short_help='Obtain list of reads that passed q-score')
@click.option(
    '-sp', '--signif-pos', default='', 
    help='Significant positions or motif to do restrict the error comparison'
)
@click.option(
    '-ft', '--freq_treat', default='', 
    help='Frequency file of the treated sample'
)
@click.option(
    '-fu', '--freq_untreat', default='', 
    help='Frequency file of the untreated sample'
)
@click.option(
    '-o', '--output', default='', 
    help='output directory to save plots and files'
)
@click.option(
    '-ks', '--kmer_study', default=False, 
    help='whether to do analysis on kmer or single nucleotide'
)
@click.option(
    '-cpg', '--cpg-ecoli', default=False, 
    help='Whether to compare profiles of only CpG islands'
)

def main(signif_pos, freq_treat, freq_untreat, output, kmer_study, cpg_ecoli):
    df = get_data(freq_treat, freq_untreat, cpg_ecoli)
    if kmer_study: 
        if cpg_ecoli: 
            df = get_kmer_CpGs(df)
            pl.kmer_advanced_plots(df, output)
        elif signif_pos:
            df_sig = pd.read_csv(signif_pos,sep='\t')
            merge_sig_positions(df, df_sig, kmer_study, output)
        else:
            pl.kmer_general_plots(df, output)
    else:
        if cpg_ecoli: 
            df = get_single_CpGs(df)
            pl.advanced_plots(df, output)
        elif signif_pos:
            df_sig = pd.read_csv(signif_pos,sep='\t')
            merge_sig_positions(df, df_sig, kmer_study, output)
        else:
            pl.general_plots(df, output)
  

if __name__ == '__main__':
    main()