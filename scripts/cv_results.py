
# coding: utf-8

# __TO DO__
# - make cv_run fxn into a class with functions (rather than function with nested functions)
# - make separate scripts for running cv pipeline and saving as df

# ___

# In[ ]:




# In[2]:

#general packages used for data import etc
from nibabel import freesurfer as fs
import pandas as pd
import numpy as np
import os,glob
import re
import cPickle as pkl #using csv now for file formats (except cv indices currently)
from joblib import Parallel,delayed
#machine learning tools used
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import FunctionTransformer, StandardScaler
from sklearn.feature_selection import (f_regression, 
                                       SelectFpr, 
                                       SelectPercentile)
from sklearn.pipeline import Pipeline, FeatureUnion
from sklearn.model_selection import KFold, LeaveOneOut
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.feature_selection.univariate_selection import (check_X_y,
                                                            safe_sparse_dot,
                                                            issparse,
                                                            row_norms,stats)
from sklearn.base import BaseEstimator, TransformerMixin 
from pandas.core.algorithms import rank
import itertools
import time
#slice indexing
idx=pd.IndexSlice


# In[3]:

#matplotlib to see distribution
get_ipython().magic(u'matplotlib inline')
import matplotlib.pyplot as plt


# In[20]:

###################directories###################
SUBJECTS_DIR=(
    '/Volumes/Users/mbkranz/'
    'projects/ACT_Freesurfer_NewProc/'
datapath='../data/'
pkldirList=glob.glob(
    SUBJECTS_DIR+
    'networks_ML/*network*.pkl'
)

##model fitting parameters/variables##
feature_select_type='fpr'
pval_fprList=[.05,.01,.001]
perc_List=np.hstack(
    [np.linspace(.0001,.9,50),
     np.linspace(1,100,50)]
)
npList=['Memory','ExFunction']
allcovarcombosList=[
    ['Gender'],
    ['Gender','wholebrain']
]
##data frame names##
index_names=[
    'cviter',
    'cvsplit',
    'np_measure_name',
    'pkldir',
    'covariates',
    'threshold'
]
predict_names=[
    'test','subs','prediction'
]
results_names=[
    'test_rsq',
    'train_rsq',
    'test_mse',
    'train_mse',
    'num_features_selected',
    'time',
    'coef',
    'intercept'
]


# In[5]:

def _split_data(X,indices):
    dummyvals=(
        X[:,[i for i, x 
             in enumerate(indices) 
             if x]]
    )
    notdummyvals=(
        X[:,[i for i, x 
             in enumerate(indices) 
             if not x]]
    )
    return notdummyvals,dummyvals
class StandardScaler_notdummy(BaseEstimator, 
                              TransformerMixin):
    """sklearn preprocessing
    - standardizes all variables except dummy variables 
    - i.e., (vars that are 0 and 1)
    - e.g., Gender """
    def _reset(self):
        if hasattr(self, 'dummy_indices_'):
            del self.dummy_indices_
            del self.notdummy_mean_
            del self.notdummy_sd_
    def fit(self, X, y=None):
        self._reset()
        self.dummy_indices_=[np.array_equal(np.unique(x),[0,1]) 
                   for x in X.T]
        notdumm,dumm=_split_data(X,self.dummy_indices_)
        self.notdummy_mean_=np.mean(notdumm,axis=0)
        self.notdummy_sd_=np.std(notdumm,axis=0)
        return self
    def transform(self,X,y=None):
        notdumm,dumm=_split_data(X,self.dummy_indices_)
        if notdumm.shape[1]>0:
            notdumm -= self.notdummy_mean_
            notdumm /= self.notdummy_sd_
            return np.hstack((notdumm,dumm))
        else: #only dummy variables
            return dumm


# In[6]:

def f_corr(X,y,spearman=True,center=True):
    '''
    this is taken from sklearn f_regression except option to rank
    variables and outputs correlations instead of f values
    '''
    if spearman==True:
        X=rank(X)
        y=rank(y)

    X, y = check_X_y(X, y, ['csr', 'csc', 'coo'], dtype=np.float64)
    n_samples = X.shape[0]

    # compute centered values
    # note that E[(x - mean(x))*(y - mean(y))] = E[x*(y - mean(y))], so we
    # need not center X
    if center:
        y = y - np.mean(y)
        if issparse(X):
            X_means = X.mean(axis=0).getA1()
        else:
            X_means = X.mean(axis=0)
        # compute the scaled standard deviations via moments
        X_norms = np.sqrt(row_norms(X.T, squared=True) -
                          n_samples * X_means ** 2)
    else:
        X_norms = row_norms(X.T)

    # compute the correlation
    corr = safe_sparse_dot(y, X)
    corr /= X_norms
    corr /= np.linalg.norm(y)

    # convert to p-value
    degrees_of_freedom = y.size - (2 if center else 1)
    F = corr ** 2 / (1 - corr ** 2) * degrees_of_freedom
    pv = stats.f.sf(F, 1, degrees_of_freedom)
    return corr, pv


# In[7]:

def cv_run(pkldir,cvtype,
           feature_select_type=None,
           threshList=[1],
           covarListList=[''],
           npdata_features=False,
           npdata_features_names=None,
           cviter=1,nfolds=None,
           permtest=False, 
           extra_savefile_str='',
           run_features=False,fit_inter=True):

    ########cv split method#######
    if cvtype=='loo':
        cv=[{cog:LeaveOneOut() for cog in npList}]
    else:
        
        cv=list(cogcvs)
    ##########load MRI data########
    if npdata_features==False:
        print(' loading '+pkldir)
        i_network_df=pd.read_pickle(pkldir)
        i_networkvals=i_network_df.values
        print('  starting '+pkldir)
    else:
        print('running cv for {}'.format(','.join(npdata_features_names)))
        i_networkvals=npdata[npdata_features_names].values
    ##########run and return cv_features if selected########
    def cv_features(np_measure_name,train,test):
        y_target=npdata[np_measure_name].values
        X_train=i_networkvals[train]
        y_train=y_target[train]
        X_test=i_networkvals[test]
        y_test=y_target[test]
        corrs,pvals=f_corr(X_train,y_train)
        return tuple([(np_measure_name,pkldir),(corrs,pvals)])
    if run_features==True:
        cv_features=[
            (i,splitnum,
             cv_features(cog,splits[0],splits[1]))
            for cog in npList 
            for i in xrange(cviter) 
            for splitnum,splits 
            in enumerate(cv[i][cog].split(npdata))
        ]
        return cv_features
 
    #######define and run cv pipeline##########
    def cv_evaluation(np_measure_name,train,test,
                      thresh,covarList,scaley=True):
        ####make some of the feature selection,transform, estimate objects####
        starttime=time.time()
        
        y_target=npdata[np_measure_name].values
        if permtest==True:
            np.random.shuffle(y_target)
        ######create train and test datasets######
        X_train=i_networkvals[train]
        y_train=y_target[train]
        X_test=i_networkvals[test]
        y_test=y_target[test]
        if scaley==True:
            y_mean=np.mean(y_train)
            y_sd=np.std(y_train)
            y_train -= y_mean
            y_train /= y_sd
            y_test -= y_mean
            y_test /= y_sd
        #concatenate covariate data if covariates
        if covarList!='': 
            measure=re.search(
                'area|thickness',
                pkldir
            ).group(0)
            covarList=[x+'_'+measure 
                       if x=='wholebrain' 
                       else x 
                       for x in covarList]
            X_train=np.hstack(
                (npdata[covarList].values[train],
                 X_train)
            )
            X_test=np.hstack(
                (npdata[covarList].values[test],
                 X_test)
            )
        ######construct pipeline##############
        #define linear regression estimator with scaling
        
        linear_est=[
            ('scale2',StandardScaler_notdummy()),
            ('linear',LinearRegression(
                fit_intercept=fit_inter
            ))
        ]
        #if running npdata as features, then only need estimator in pipeline
        if npdata_features==True:
            cv_pipeline=Pipeline(linear_est)
            num_features_selected=(-1)
        else:
            #define transform fxns
            def meansum_transform(x,key):
                if any([s in key for s in ['thickness']]):
                    trans=np.mean(x,axis=1).reshape(-1,1)
                elif 'area' in key:
                    trans=np.sum(x,axis=1).reshape(-1,1)
                return trans
            def selectfeat(x,key,ncov):
                if any([s in key for s in ['area','thickness']]):
                    trans=x[:,ncov:]
                elif key=='covars':
                    trans=x[:,:ncov-1]
                return trans
            
            numcovar=len(covarList) #used in selectfeat
            selectvertices=FunctionTransformer(
                selectfeat,kw_args={'key':pkldir,'ncov':numcovar})
            selectcovars=FunctionTransformer(
                selectfeat,kw_args={'key':'covars','ncov':numcovar})
            if feature_select_type=='fpr':
                anova_select=SelectFpr(f_corr,alpha=thresh)
            elif feature_select_type=='percentile':
                anova_select=SelectPercentile(f_corr,percentile=thresh)
                
            combine_morph=FunctionTransformer(meansum_transform,
                                              kw_args={'key':pkldir})
            vertex_pipe=[('selector_vert',selectvertices),
                        # ('scale1',RobustScaler()),
                         ('anovaselect',anova_select),
                         ('combine', combine_morph)]
            if numcovar>0:
                covars_pipe=[('selector_cov',selectcovars)]
                vertex_covars_tuple=[('morph',Pipeline(vertex_pipe)),
                                     ('covars',Pipeline(covars_pipe))]
                vertex_covars_union=[('union',FeatureUnion(
                    transformer_list=vertex_covars_tuple))]
            
            #if covariates,
            #then include the FeatureUnion object
            if numcovar>0:
                cv_pipeline=Pipeline(vertex_covars_union+linear_est)
             #if not covar list 
            #then just do:
                ##feature selection, 
                ##avg/sum vertices, 
                ##run estimator
            else:
                cv_pipeline=Pipeline(vertex_pipe+linear_est) 
            only_vertices_train=selectvertices.fit_transform(
                X_train
            )
            num_features_selected=anova_select.fit_transform(
                X=only_vertices_train,
                y=y_train
            ).shape[1]
        ####run pipeline if greater than 0 features or npdata as features####
        ##if there are features select -->
        ##there are no features selected (can't fit)
        #no features selected for npdata as features so num_feat=NA
        if num_features_selected>0 or npdata_features==True:
            #fit
            cv_pipeline.fit(X=X_train,y=y_train)
            #predictions
            X_train_pred=cv_pipeline.predict(X_train)
            X_test_pred=cv_pipeline.predict(X_test)
            #metrics
            train_rsq=r2_score(y_train,X_train_pred)
            train_mse=mean_squared_error(y_train,X_train_pred)
            test_rsq=r2_score(y_test,X_test_pred)
            test_mse=mean_squared_error(y_test,X_test_pred)
            cv_results=[
                test_rsq,
                train_rsq, 
                test_mse, 
                train_mse,
                num_features_selected
            ]
            cv_predicts=X_test_pred.reshape(-1,1)
            
            coefs=list(
                cv_pipeline
                .named_steps
                ['linear']
                .coef_
            )
            intercept=[
                cv_pipeline
                .named_steps
                ['linear']
                .intercept_
            ]
        else:
            cv_results=[np.nan]*4+[0]
            cv_predicts=np.array([[np.nan]]*len(test))
            coefs=[np.nan]
            intercept=[np.nan]
        ##return a list containing 
        ###1. list of indices 
        ###2. train/test performance metrics 
        ###3. prediction values with ids (e.g., subs)
        ids_test=ids[test].reshape(-1,1)
        index_results=[np_measure_name,
                       pkldir,'_'.join(covarList),
                       str(thresh)]
        return [
            index_results,
            cv_results+
            [time.time()-starttime]+
            [coefs]+
            intercept,
            np.hstack((test.reshape(-1,1),
                       ids_test,cv_predicts))
        ]
    
    ##runs cv train test (cv_evaluation) ##
    '''
    - for each target (i.e., cognitive var)
    - feature selection threshold
    - covariate
    - cv split
    - train/test cv performance results for:
       - target variables, 
       - features selection thresholds
       - covariates
       - for eaceh cv split and iteration
    - returns a list of tuples with: 
        - iteration, split, 
        - cv_results/predictions/index vars
    '''
    cv_results=[
        (i,splitnum,
         [cv_evaluation(
             cog,
             splits[0],splits[1],
             thresh,covarList
         ) 
          for thresh in threshList
          for covarList in covarListList])
        for cog in npList 
        for i in xrange(cviter) 
        for splitnum,splits 
        in enumerate(cv[i][cog].split(npdata))]
    pkldir_save=(
        datapath+
        'cv_results/tmp/'
        '_tmp_cv_results_'+
        extra_savefile_str+
        '_'.join(
            [y 
             for x in covarListList 
             for y in x]
        )
    )
    if npdata_features==True:
        cv_results_pkldir=(pkldir_save+
                           pkldir+
                           '.pkl')
    else:
        cv_results_pkldir=pkldir.replace(
            SUBJECTS_DIR+
            'networks_ML/',pkldir_save
        )
    print('saving '+ cv_results_pkldir)
    pkl.dump(cv_results,open(cv_results_pkldir,'wb'))
    ###################################
    return pkldir


# In[8]:

#save dataframes of cv results
def save_cv_dfs(cv_list,
                featselectstr,
                endstr='_df.csv',
                return_data=False,
                save_data=True):

    cv_list_flat=[
        (i,s,y) 
        for r in cv_list 
        for i,s,x in r 
        for y in x
    ]
    predict_index=list(
        itertools
        .chain
        .from_iterable(
            itertools.repeat(
                [i,s]+x[0], len(x[2])) 
            for i,s,x in cv_list_flat
        )
    )
    results_index=[
        [i,s]+x[0] 
        for i,s,x in cv_list_flat
    ]
    def make_indices(index):
        df=pd.DataFrame(index,columns=index_names)
        return df.set_index(index_names).index
    predict_df=pd.DataFrame(
        np.vstack([x[2] for i,s,x in cv_list_flat]),
        index=make_indices(predict_index),
        columns=predict_names
    )
    results_df=pd.DataFrame(
        [x[1] for i,s,x in cv_list_flat],
        index=make_indices(results_index),
        columns=results_names
    )
    if save_data==True:
        predict_df.to_csv(
            datapath+
            'cv_results/'+
            'cv_predict_'+
            featselectstr+
            endstr,
            na_rep='NA'
        )
        results_df.to_csv(
            datapath+
            'cv_results/'+
            'cv_results_'+
            featselectstr+
            endstr,
            na_rep='NA'
        )
    if return_data==True:
        return (results_df,predict_df)


# In[9]:

#import behavioral data (target)
#np_all filtered to include good MRI data in MakeDataFrames script
npdata=pd.read_csv(
    datapath+
    'np_filter_wb_gendernum.csv',
    na_values='NA',
    index_col="subs"
)
ids=npdata.index.values


# In[10]:

#create or read in sklearn validation object list
#if create makes a unique file with time.time
nfolds=5
read_folds=True
create_folds=False
write_folds=False
if create_folds==True:
    cogcvs=[
        {cog:KFold(n_splits=nfolds,
                   shuffle=True) 
         for cog in npList} 
        for i in xrange(100)
    ]
if write_folds==True:
    pkl.dump(
        cogcvs,open(
        '../data/cv_results/'
        'cvfold_indices_{}fold_{}.pkl'
        .format(str(nfolds),str(time.time()))
        ),
        'wb')
if read_folds==True:
    cogcvs=pkl.load(open(
        '../data/cv_results/'
        'cvfold_indices_{}fold.pkl'
        .format(str(nfolds)),
        'rb'
    ))


# In[ ]:

#run cross validation with covariates (gender and wholebrain)
Parallel(5)(
    delayed(cv_run)(
        pkldir=pkldir,
        feature_select_type='fpr',
        threshList=pval_fprList,
        covarListList=allcovarcombosList,
        cvtype='kfold',
        cviter=100
    )
    for pkldir in pkldirList
)


# In[ ]:

#run cross validation with NPDATA covariates (gender and wholebrain)
Parallel(5)(
    delayed(cv_run)(
        pkldir='npdata_'+'_'.join(names),
        npdata_features=True,
        npdata_features_names=names,
        cvtype='kfold',
        cviter=100
    )
    for names in [
        ['Gender'],
        ['Gender','wholebrain_thickness'],
        ['Gender','wholebrain_area']
    ]
)


# In[ ]:

#run cross validation
Parallel(5)(delayed(cv_run)
               (pkldir=pkldir,
                feature_select_type='fpr',
                threshList=pval_fprList,
                covarListList=[''],
                cvtype='kfold',
                cviter=100)
               for pkldir in pkldirList)


# In[ ]:

#run cross validation more conservative p values
Parallel(5)(delayed(cv_run)
               (pkldir=pkldir,
                extra_savefile_str='extremepvals',
                feature_select_type='fpr',
                threshList=[.0005,.0001,.00001],
                covarListList=[''],
                cvtype='kfold',
                cviter=100)
               for pkldir in pkldirList)


# In[ ]:

#running network 2 thickness as Kramerlab is running 3-7 (and done with area)
pkldirList_thicknet2=[x for x in pkldirList if re.search('thickness_.*network2',x)]
#run cross validation with percentile selection
Parallel(5)(delayed(cv_run)
               (pkldir=pkldir,
                extra_savefile_str='percentile',
                feature_select_type='percentile',
                threshList=perc_List,
                covarListList=[''],
                cvtype='kfold',
                cviter=100)
               for pkldir in pkldirList_thicknet2)


# In[21]:

#run cross validation with percentile selection
Parallel(5)(delayed(cv_run)
               (pkldir=pkldir,
                extra_savefile_str='percentile',
                feature_select_type='percentile',
                threshList=perc_List,
                covarListList=[''],
                cvtype='kfold',
                cviter=100)
               for pkldir in pkldirList)


# In[ ]:

#run cross validation permutations
for permiter in range(10,21):
    Parallel(5)(delayed(cv_run)
                (pkldir=pkldir,
                 feature_select_type='fpr',
                 threshList=pval_fprList,
                 covarListList=[''],
                 cvtype='kfold',
                 cviter=100,
                 extra_savefile_str='permtest'+str(permiter),
                 permtest=True)
                for pkldir in pkldirList)


# In[14]:

tmpfiledir=(
    datapath+
    'cv_results/tmp/'
    '_tmp_cv_results_{}_'
    'fwhm10_network*fsaverage_df.pkl'
)


# In[15]:

#save cv results as dataframe
cv_fpr_tmpList=(
    glob.glob(
        tmpfiledir.format('extremepvalsarea')
    )+
    glob.glob(
        tmpfiledir.format('extremepvalsthickness')
    )
)
cv_fpr=[pkl.load(open(tmp,'rb')) 
        for tmp in cv_fpr_tmpList]
save_cv_dfs(
    cv_fpr,
    'fpr',
    endstr='_df_extremepvals_5fold.csv'
)


# In[ ]:

#save cv results as dataframe
cv_fpr_tmpList=(
    glob.glob(
        tmpfiledir.format('area')
    )+
    glob.glob(
        tmpfiledir.format('thickness')
    )
)
cv_fpr=[pkl.load(open(tmp,'rb')) 
        for tmp in cv_fpr_tmpList]
save_cv_dfs(
    cv_fpr,
    'fpr',
    endstr='_df_5fold.csv'
)


# In[14]:

#save NPDATA cv results with covariates as dataframe
cv_fpr_npdata_tmpList=glob.glob(
    datapath+
    'cv_results/tmp/*npdata*.pkl'
)
cv_fpr_npdata=[pkl.load(open(tmp,'rb')) 
               for tmp in cv_fpr_npdata_tmpList]
save_cv_dfs(
    cv_fpr_npdata,
    'fpr',
    endstr='_df_5fold_npdata.csv'
)


# In[ ]:

#save cv results with covariates as dataframe
cv_fpr_covars_tmpList=glob.glob(
    tmpfiledir.format('Gender*')
)
cv_fpr_covars=[pkl.load(open(tmp,'rb')) 
               for tmp in cv_fpr_covars_tmpList]
save_cv_dfs(
    cv_fpr_covars,
    'fpr',
    endstr='_df_5fold_covars.csv'
)


# In[29]:

#save perm cv_results as dataframe
for permiter in range(10):
    cv_perm_tmpList=glob.glob(
        tmpfiledir.format('permtest'+
                          str(permiter)+'*')
    )
    cv_fpr_perm=[pkl.load(open(tmp,'rb')) 
                 for tmp in cv_perm_tmpList]
    save_cv_dfs(
        cv_fpr_perm,
        'fpr',
        endstr=('_df_5fold_permtest{}.csv'
                .format(str(permiter)))
    )


# ###### make cv features for better interpretation

# In[ ]:

#run for whole brain thickness and area
cv_features=[cv_run(pkldir=pkldir,cvtype='kfold',run_features=True)


# In[ ]:

def make_cv_feature_dfs(cv_list):
    i_list=[z[0] 
            for df in cv_list 
            for x,y,z in df]
    i_index=pd.MultiIndex.from_tuples(
        i_list,
        names=['np_measure_name','pkldir']
    )
    corrs=[z[1][0] 
           for df in cv_list 
           for x,y,z in df]
    pvals=[z[1][1] 
           for df in cv_list 
           for x,y,z in df]
    return (pd.DataFrame(corrs,i_index),
            pd.DataFrame(pvals,i_index))


# In[ ]:

corrs_df,pvals_df=make_cv_feature_dfs(cv_features)


# In[ ]:

corrs_df.to_pickle(
    datapath+
    'cv_results/'
    'corrs_df_features_5fold.pkl'
)
pvals_df.to_pickle(
    datapath+
    'cv_results/'
    'pvals_df_features_5fold.pkl'
)

