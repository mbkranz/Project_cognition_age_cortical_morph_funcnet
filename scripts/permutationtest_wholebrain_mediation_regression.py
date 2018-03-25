
# coding: utf-8

# __run permutation tests from `regression_mediation_pervertex_wholebrain` script__

# In[1]:

#general packages used for data import etc
from nibabel import freesurfer as fs
import pandas as pd
import numpy as np
import os,glob
import re
import cPickle as pkl #using csv now for file formats (except cv indices currently)
from sklearn.feature_selection import variance_threshold
from sklearn.utils import resample
from joblib import Parallel,delayed
idx=pd.IndexSlice


# In[2]:

measure_list=['area','thickness']
np_name_list=['Memory','ExFunction']
hemilist=['rh','lh']
bootlist=['bootscore','bootmean']
meas_key={'thickness':np.int(0),'area':np.int(1)}
hemi_key={'lh':np.int(0),'rh':np.int(1)}
n=1000


# In[3]:

def getannot(annotname):
    #initiate DataFrame
    #may want to make concatenation/join (instead of append) 
    #so can have one column per annotation/set of labels
    annot_df=[]
    for hemi in hemilist:
        annot_data=fs.read_annot(
            '/Applications/freesurfer/'
            'subjects/fsaverage/'
            'label/'+
            hemi + 
            '.' + 
            annotname + 
             '.annot'
        )
        annot_hemi=pd.DataFrame(
            {"annot_label" : annot_data[0],
             "annot_name": annotname, 
             "vertex_index" : range(
                 len(annot_data[0])
             ), 
             "hemi": hemi_key[hemi]})
        annot_df.append(annot_hemi)
    annots=pd.concat(annot_df).set_index(
        ['hemi','vertex_index']
    )
    return annots


# In[4]:

def test_permutation(pklfile,
                     num_perm=100000,
                     threshList=[1.96,2.58,3.3,3.9]):
    def get_network_summary(df,iteration='computed'):
        def compute_summary(group):
            total_sigvertices=[
                (group>thresh).sum()+
                (group<(-thresh)).sum()
                for thresh in threshList
            ]
            return np.array(total_sigvertices)
        df_copy=df.copy()
        #for permutation test: shuffle annotation labels
        if iteration!='computed': 
            df_copy['annot_label']=(
                df_copy['annot_label']
                .sample(frac=1)
                .values
            )
        select_summary=(
            df_copy
            .set_index(['annot_label'])
            ['bootscore']
            .groupby(level=['annot_label'])
            .apply(compute_summary))
        return select_summary
    bootdf=pd.read_pickle(pklfile)
    multi_i=(
        pd.MultiIndex.from_tuples(
            bootdf['hemi_vertex'],
            names=[0,'hemi','vertex_index']
        )
        .droplevel([0])
    )
    del bootdf['hemi_vertex']
    bootdf=pd.DataFrame(
        data=bootdf.values,
        index=multi_i,
        columns=bootdf.columns
    )
    annots=getannot('Yeo2011_7Networks_N1000')
    bootdf_wannot=(
        bootdf
        .join(annots)
        .query('annot_label!=0')
    )
    del bootdf_wannot['annot_name']
    empirical_vals=np.vstack(
        get_network_summary(
            bootdf_wannot
        )
    )
    perm_greaterthan=(
        np.zeros(empirical_vals.shape)
    )
    perm_totalsigvertices=(
        np.zeros(empirical_vals.shape)
    )
    #get permutations 
    #(shuffling network labels to make inferences)
    for i in xrange(num_perm):
        perm_vals=np.vstack(
            get_network_summary(
                bootdf_wannot,
                iteration=i
            )
        )
        perm_greaterthan+=perm_vals>=empirical_vals
        perm_totalsigvertices+=perm_vals
    perm_meansigvertices=np.divide(
        perm_totalsigvertices,
        num_perm
    )
    perm_pvals=np.divide(
        perm_greaterthan,num_perm
    )
    multiindex=pd.MultiIndex.from_tuples(
        tuples=zip(np.arange(1,8),
                   [pklfile]*len(empirical_vals)),
        names=['network','pklfile'])
    return_vars=[
        'perm_meansigvertices',
        'perm_totalsigvertices',
        'perm_greaterthan',
        'perm_pvals',
        'empirical_vals'
    ]
    return_df=pd.concat(
        [(pd.DataFrame(
                data=eval(x),
                index=multiindex,
                columns=threshList)
          .assign(var=x)
          .set_index(['var'],append=True))
         for x in return_vars])
    return return_df


# In[ ]:

bootstrap_pklfile = get_ipython().getoutput(u'ls ../data/wholebrain_bootstrap/*.pkl')
summary_perms_all=pd.concat(
    Parallel(5)(
        delayed(test_permutation)
        (pkl) 
        for pkl in bootstrap_pklfile
    )
)


# In[ ]:

summary_perms_all.head()


# In[ ]:

summary_perms_unstack=(
    summary_perms_all
    .stack()
    .unstack(['var'])
)
summary_perms_unstack.head()


# In[ ]:

(summary_perms_unstack
 .to_csv(
     '../data/wholebrain_bootstrap/'
     'wholebrain_bootstrap_
     'bootscore_permtests.csv'
 )
)

