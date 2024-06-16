import pandas as pd
import numpy as np
import os
import sys
from tcrdist.repertoire import TCRrep


#user fold
datset_directory='data/small_enosh/enosh_small.csv'
light_chain='empty'# fill empty string for each cdr in light chain 

def process_df(df,light_chain):
    """
    preprocessing df, return clena df for tcrrep
    """
    df["v_a_gene"]=np.nan
    df['v_b_gene']=np.nan
    df['pmhc_a_aa']=''
    df['pmhc_b_aa']=''
    if light_chain=='empty':
        df['cdr1_a_aa']=''
        df['cdr2_a_aa']=''
        df['cdr3_a_aa']=''

    df['subject']=df.index
    return df


# initialize tcrep object with tcr sequences (only for initialization)
df_tcr = pd.read_csv("data/dash/dash.csv")

#define which chains we want to calcualte (beta or alpha+beta)
chains=['beta']
if light_chain!="empty":
    chains = ['alpha','beta']

bcrrep = TCRrep( cell_df=df_tcr,
            organism = 'mouse',
            chains = chains,
            db_file = 'alphabeta_gammadelta_db.tsv',imgt_aligned=True)
bcrrep.kargs_b['cdr3_b_aa']['fixed_gappos']=False
bcrrep.kargs_b['cdr2_b_aa']['fixed_gappos']=False
bcrrep.kargs_b['cdr1_b_aa']['fixed_gappos']=False
bcrrep.kargs_b['cdr2_b_aa']['gap_penalty']=4



#import dataset
df_bcr=pd.read_csv(datset_directory)
df_bcr=process_df(df_bcr,light_chain)

#compute bcr distance
bcrrep.clone_df=df_bcr
bcrrep.compute_distances()

df_matrix=pd.DataFrame(bcrrep.pw_beta)
df_matrix.to_csv(f"{os.path.dirname(datset_directory)}/distnance_matrix.csv")