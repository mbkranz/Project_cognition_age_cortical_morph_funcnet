
# coding: utf-8

# ## Modules and subject directory

# __TO DO: test with code improvements__

# In[ ]:

#general packages used for data import etc
from nibabel import freesurfer as fs
import nibabel as nib
import pandas as pd
import numpy as np
import os,glob
import re
from joblib import Parallel,delayed
import cPickle as pkl
idx=pd.IndexSlice
datapath=os.getwd().replace(
    'scripts/','data/'
)
#subject directory and subject list
SUBJECTS_DIR=(
    '/Volumes/Users/mbkranz/'
    'projects/ACT_Freesurfer_NewProc/'
)
os.chdir(SUBJECTS_DIR)
sublist = get_ipython().getoutput(u'ls -d ACT*_1')
#npdata from make_npdata_from_pca.R 
npdata_all=pd.read_csv(
    datapath+
    'np_all.csv',
    na_values='NA'
)
#filter out all subs without behav data
sublist_filt=[
    x 
    for x in sublist 
    if int(re.sub("ACT0|ACT|_1","",x)) 
    in npdata_all['subs'].values 
    and x!='ACT6142_1' #bad sub (see notes)
]
os.chdir(gitpath)


# In[ ]:

meas_key={'thickness':np.int(0),
          'area':np.int(1)}
hemi_key={'lh':np.int(0),
          'rh':np.int(1)}
net_key={'network'+str(net):net 
         for net in [1,2,3,4,5,6,7]}
net_key['wholebrain']=[1,2,3,4,5,6,7]


# In[ ]:

def getmorph(sub,measure,fwhm,hemi):
    path_data=(
        SUBJECTS_DIR+
        sub+'/surf/'+
        hemi+'.'+measure+
        '.'+'fwhm'+fwhm+
        '.fsaverage.mgh'
    )
    morph_data=(nib.load(path_data)
                .get_data()
                .flatten()
               )
    morph_df=pd.DataFrame({
        'subs':np.int(re.sub("ACT0|ACT|_1",
                             "",
                             sub)),
        'measure':meas_key[measure],
        'hemi':hemi_key[hemi],
        'value':morph_data
    })
    return (morph_df
            .rename_axis('vertex_index')
            .set_index(['hemi'],append=True)


# In[ ]:

def getannot(annotname):
    hemilist=['lh','rh']
    #initiate DataFrame
    #may want to make concatenation/join (instead of append) 
    #so can have one column per annotation/set of labels
    annot_df=[]
    for hemi in hemilist:
        annot_data=fs.read_annot(
            '/Applications/freesurfer/'
            'subjects/fsaverage/label/' +
            hemi + '.' + annotname + '.annot'
        )
        annot_hemi=pd.DataFrame({
            "annot_label" : annot_data[0],
            "hemi": hemi_key[hemi]
        })
        annot_df.append(annot_hemi)
    annots=pd.concat(annot_df)
    return (annots
            .rename_axis('vertex_index')
            .set_index(['hemi'],append=True)
            .reorder_levels(['vertex_index','hemi'])
           )


# In[ ]:

morph_df_all=pd.concat(
    Parallel(n_jobs=4)(
        delayed(getmorph)
        (sub,meas,'10',hemi) 
        for sub in sublist_filt 
        for hemi in ['lh','rh']
        for meas in ['thickness','area']
    )
)


# In[ ]:

annots=getannot('Yeo2011_7Networks_N1000')
morph_df_annot=(
    morph_df_all
    .join(annots)
    .set_index(
        ['subs','measure','annot_label'],
        append=True
    )
    .reorder_levels(
        ['annot_label',
         'measure',
         'hemi',
         'subs',
         'vertex_index']
    )
    .sort_index()
)


# In[ ]:

def write_morph_data(meas,imeas,netname,inet):
    imorph_data=(
        morph_df_annot
        .loc[idx[inet,imeas,:,:,:],:]
    )
    #drop as unstacking only hemi and vertex
    imorph_data.index=(
        imorph_data
        .index
        .droplevel('annot_label')
    )
    hemi_vertex=['hemi','vertex_index']
    imorph_data=(imorph_data
                 .sort_index(level=hemi_vertex)
                 .unstack(hemi_vertex)
                )
    netpath=(
        SUBJECTS_DIR+
        '/networks_ML/'+
        meas+'_fwhm10_'+
        netname+
        'fsaverage_df.pkl' 
    )
    imorph_data.to_pickle(netpath)
    print('Done with '+netname+' for '+meas)


# In[ ]:

#save individual network and wholebrain files
for meas,imeas in meas_key.iteritems():
    for netname,inet in net_key.iteritems():
        write_morph_data(meas,imeas,netname,inet)


# #### add wholebrain average values to npdata

# In[ ]:

npdata_all.index.names


# In[ ]:

npdata_all.query(
    'subs in @morph_df_annot'
    '.index.levels[3].values'
)


# In[ ]:

wb_area=pd.read_pickle(
    '{}networks_ML/'
    '{}_fwhm{}_{}fsaverage_df.pkl'
    .format(SUBJECTS_DIR,
            'area',
            str(10),
            'wholebrain')
)
wb_thick=pd.read_pickle(
    '{}networks_ML/'
    '{}_fwhm{}_{}fsaverage_df.pkl'
    .format(SUBJECTS_DIR,
            'thickness',
            str(10),
            'wholebrain')
)
wb_thick_avg=(wb_thick
              .mean(axis=1)
              .reset_index(
                  ['measure'],
                  drop=True)
              .rename('wholebrain_thickness')
             )
wb_area_avg=(wb_area
             .sum(axis=1)
             .reset_index(['measure'],
                          drop=True)
             .rename('wholebrain_area')
            )
wb_both=pd.concat(
    [wb_thick_avg,wb_area_avg],
    axis=1
)
npdata=(npdata_all
        .set_index(['subs'])
        .join(wb_both,how='right')
       )
npdata['Gender']=np.where(
    npdata['Gender']=='Female',
    0,
    1
)


# In[ ]:

npdata.to_csv(datapath+
              'np_filter_wb_gendernum.csv')


# In[ ]:



