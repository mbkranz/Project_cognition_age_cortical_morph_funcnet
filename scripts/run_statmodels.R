library(tidyverse)
library(parallel)
library(multidplyr)
library(parallel)
library(magrittr)
library(kableExtra)
library(knitr)
library(formattable)

###Keys for plot colors and misc defined variables####
setwd('../figures')
NetworkColors<- c(
  rgb(120,18,134,maxColorValue = 255),
  rgb(70,130,180,maxColorValue = 255),
  rgb(0,118,14,maxColorValue = 255),
  rgb(196,58,250,maxColorValue = 255),
  rgb(128,128,128,maxColorValue = 255),
  rgb(230,148,34,maxColorValue = 255),
  rgb(205,62,78,maxColorValue = 255)
  )
TypeAlphaColors<- c(
  "average"='black',
  "alpha: 0.01" = 'white',
  "alpha: 0.05" = 'white',
  "alpha: 0.001" = 'white'
  )
corr_method <- 'spearman'
#using parallel clusters with dplyr
#so set inner bootstrap loops to serial
parallel_boot='no'
parallel_boot_cores=1
cols_group_by<-c(
  'np_measure_name',
  'network',
  'measure',
  'type_alpha',
  'threshold',
  'covariates',
  'index',
  'pkldir'
  )
######Demographics####
#get full dataset (with all participants) and demog vars
#get final dataset with list of included subjects
data_raw<- 
  readxl::read_excel("../data/ACT_NP.xlsx",
                     sheet='All Subjects')
data_filt<-
  data_raw %>% 
  filter(MRIPre==1)
subs_final<-
  read_csv(
    '../data/np_filter_wb_gendernum.csv'
  )[['subs']]
demogs_final<-
  data_raw %>% 
  filter(Subject_ID %in% subs_final) %>%
  select(Gender,matches('Age|Years.*Education'))

comp_baseline<-table(data_raw$EligiblePre)['1']
education<-psych::describe(demogs_final$Years_of_Education)
age<-psych::describe(demogs_final$Age_at_start_of_study)
gender<-table(demogs_final$Gender)
#####Principal Component Analysis####
#This code is taken from the make_npdata_frompca.R script
windsorizefxn<-function(data,
                        variablelist,
                        outlieroutput=FALSE,
                        replaceval=TRUE,
                        std=3){
  'outliers (outputs a list: 
  [1] windsorized value data set and 
  [2]makes a list of subject outliers per task)'
  #if data.table from dplyr can't pull numeric values
  ddata<-as.data.frame(data)
  wdata<-as.data.frame(data) 
  outliers<-data.frame(subs=data[,1])
  for(variable in variablelist){
    #probablity values for outliers 
    #(for mean=0,sd=1-->i.e., zscore probablity)
    uwind<-pnorm(std)
    lwind<-pnorm(-std)
    #gets values of variable for probabilities
    cut_point_top <-qnorm(
      p=uwind,
      mean=mean(ddata[,variable],na.rm=TRUE),
      sd=sd(ddata[,variable],na.rm=TRUE))
    cut_point_bottom <-qnorm(
      p=lwind,
      mean=mean(ddata[,variable],na.rm=TRUE),
      sd=sd(ddata[,variable],na.rm=TRUE))
    #replaces outlier variable values within new wdata dataframe
    if(replaceval==TRUE){
      wdata[,variable]<-ifelse(
        ddata[,variable]>cut_point_top,cut_point_top,
        ifelse(ddata[,variable]<cut_point_bottom,
               cut_point_bottom,
               ddata[,variable]))
    }else{
      wdata[,variable]<-ifelse(
        ddata[,variable]>cut_point_top,NA,
        ifelse(ddata[,variable]<cut_point_bottom,
               NA,
               ddata[,variable]))}
    #makes dataset with outlier values replaced with windsorized values 
    #(calling from inputted data dataframe)
    outliers[,variable]<-ifelse(
      (ddata[,variable]>cut_point_top | 
         ddata[,variable]<cut_point_bottom),
      ddata[,variable],NA)}
  if(outlieroutput==TRUE){outputdata<-outliers
  }else{outputdata<-wdata}  
  return(outputdata)}

pr.model<-
  data_filt %>%
  #make RT/error measures negative so positive means good for all measures
  mutate(Flanker_Pre_Incongruent_RT=-Flanker_Pre_Incongruent_RT,
         Time_Trial_B_Pre=-Time_Trial_B_Pre,
         DotComparison_RT_Pre=-DotComparison_RT_Pre,
         Imove_Swaps_Pre=-Imove_Swaps_Pre) %>%
  transmute(
    `Flanker Incongruent`=
      Flanker_Pre_Incongruent_RT,
    `Trail Making Part B`=
      Time_Trial_B_Pre,
    `Spatial Working Memory`=
      .[] %>% 
      select(matches("SWM_ACC")) %>% 
      rowMeans(na.rm=TRUE),
    `Letter Nback (2 Back)`=
      ACC_2back_Pre,
    `Category Fluency (Animals + Vegetables)`=
      CategoryFluency_Fruits_Veg_Pre+
      CategoryFluency_Animals_Pre,
    `Face-Scene (D')`=
      Dprime_Pre, 
    `Spatial Reconstruction (Swaps)`=
      Imove_Swaps_Pre,
    `Dot Comparison`=
      DotComparison_RT_Pre,
    `Digit Symbol Coding`=
      Digit_Symbol_Correct_Pre, 
    `CVLT Free Recall (Delayed)`=
      Delayed_Free_Recall_Pre, 
    `Story Free Recall (Delayed)`=
      Story_Unit_Recall_Delayed_Pre
    ) %>% 
  windsorizefxn(.,names(.),replaceval=TRUE) %>%
  as_data_frame()

pr.model.comp2<-psych::principal(
  pr.model, 
  nfactors= 2 , 
  rotate= "varimax" , covar= FALSE , scores=TRUE)

pca_tbl_vals<-
  rbind(
    pr.model.comp2$loadings[,1:2],
    pr.model.comp2$Vaccounted['Proportion Var',]
    )
pca_tbl_names<-
  pr.model.comp2$loadings %>% 
  row.names() %>% 
  c(.,c('% Variance explained'))

cell_bold<-function(x){ 
  cell_spec(x,bold = T,format = 'html')
}
pca_table<-
  pca_tbl_vals %>%
  as_data_frame() %>%
  transmute(
    Task=pca_tbl_names,
    `Component 1: Executive Function`=round(RC1,2) %>%
      as.character %>%
      ifelse(RC1>RC2,cell_bold(.),.),
    `Component 2: Declarative Memory`=round(RC2,2) %>%
      as.character() %>%
      ifelse(RC2>RC1,cell_bold(.),.)
    )

pca_table %>%
  kable(format = 'html',
        escape = F) %>%
  kable_styling(
    bootstrap_options = "striped", 
    # full_width = F, 
    position = "float_right"
  )
######Training models summary stats####
str_meansd<- function(df,x,num) {
  mean_str=str_c(x,'_mean')
  sd_str=str_c(x,'_sd')
  str_c(round(df[[mean_str]],num),
        ' (',round(df[[sd_str]],num),')'
  )}

cv_results_all<- 
  make_cv_results_df(
    cv_results_filestr = 
      '../data/cv_results/cv_results_fpr_df_5fold.csv'
  ) %>%
  transmute(
    np_measure_name,
    measure,
    network,
    threshold=ifelse(
      round(threshold,2)==.02,
      'average',
      round(threshold,3) %>% as.character()
    ),
    coef=ifelse(
      coef_mean>0,
      str_meansd(.,'coef',2) %>% 
        cell_spec(.,color='red'),
      str_meansd(.,'coef',2) %>%
        cell_spec(.,color='blue')
    ),
    num_features=str_meansd(.,'num_features_selected'),
    nofeatures=str_c(
      round(nofeatures_selected_mean,2)*100,' %'
    )) %>%
  gather(stat,value,
         -np_measure_name,-measure,-network,-threshold) %>%
  mutate(network=network %>% 
           str_remove('network'),
         np_measure_name=np_measure_name %>%
           str_replace('ExFunction','Executive Function'))

#grouped cols:np_measure_name,measure,threshold,network,stat
#unite: grouped cols you want across columns
#group_by/spread: grouped cols you want down rows
#drop_vars<-c('stat','threshold')
drop_vars<-c('threshold')
row_tbl_vars<-c('np_measure_name','network')
col_tbl_vars<-c('stat','measure')
#col_tbl_remove_regex='(area|thickness)'
col_tbl_remove_regex='_.*$'
threshold_vars<-c('average')

row_tbl_vars_num<-row_tbl_vars %>%
  length()
col_vars_num<-cv_results_all[[col_tbl_vars[2]]] %>% 
  unique() %>% 
  length()

cv_results_table <- 
  cv_results_all %>%
  filter(threshold %in% threshold_vars) %>%
  unite(colnames,one_of(col_tbl_vars)) %>%
  select(-one_of(drop_vars)) %>%
  group_by_at(row_tbl_vars) %>%
  spread(colnames,value) 

col_var_names<-
  cv_results_table %>% 
  names(.) %>% 
  str_extract('area|thickness') %>% 
  str_replace_all(c(
    'area'='Surface Area',
    'thickness'='Thickness'
  )) %>%
  str_replace_na('')

cv_results_table %>%
  kable(format = 'html',
        escape = F,
        col.names = col_var_names) %>%
  kable_styling(
    bootstrap_options = "striped", 
    # full_width = F, 
    position = "float_right"
  ) %>%
  column_spec(1:row_tbl_vars_num, bold = T) %>%
  collapse_rows(1:row_tbl_vars_num) %>%
  add_header_above(c(
    " " = row_tbl_vars_num,
    "Coefficient M (SD)" = col_vars_num,
    "% of Folds with No Features Selected" = col_vars_num,
    "Number of Features Selected M (SD)" = col_vars_num
  )) %>%
  footnote(str_c(
    'Networks are labeled by index in Yeo et al., 2011: ',
    '1 = Visual, ',
    '2 = Somatomotor, ',
    '3 = Dorsal Attention, ',
    '4 = Salience, ',
    '5 = Limbic, ',
    '6 = Control, ',
    '7 = Default; ',
    'Linear regression training coefficient weights are color coded ',
    'by positive (red) or negative (blue) ',
    'for easier reading .'
  )
  )
####Define Model Functions######
pcorrfxn_noboot<-function(
  d,str_lmcovariates,simnum){
  'partial correlation function used in permutation tests
  i.e., controls for a set of variables and determines 
  correlation of predicted versus residuals of observed'
  lm_eval<-eval(
    parse(
      text=str_c(
        "lm(formula=",
        str_lmcovariates,
        ",data=d)"
      )
    )
  )
  corr=cor(lm_eval$residuals,
           scale(d[['predict_avg']]),
           method=corr_method)
  return(corr)
}
i_pcorrfxn_boot<-function(
  str_lmcovariates,data,i){
  'partial correlation user defined bootstrap function.
  i.e., controls for a set of variables and determines 
  correlation of predicted versus residuals of observed'
  d<-data[i,]
  bootlmobj<-
    eval(parse(text=str_c(
      "lm(formula=",
      str_lmcovariates,
      ",data=d)"
    )))
  corr=cor(bootlmobj$residuals,
           scale(d[['predict_avg']]),
           method=corr_method)
  return(corr)
}
pcorrfxn_boot<-function(
  data,str_lmcovariates,simnum){
  'get bootstrapped stats using the
  corr_residualsfxn_boot function'
  corr_boot<-boot::boot(
    data=data,
    statistic=i_pcorrfxn_boot,
    R=simnum,
    str_lmcovariates=str_lmcovariates,
    parallel = parallel_boot,
    ncpus=parallel_boot_cores
  )
  confboot<-boot::boot.ci(corr_boot,type="bca")$bca
  #bootstrap stats
  ##lower confidence interval (2.5% is default)
  bootlower<-confboot[4]
  ##upper confidence interval (97.5% is default)
  bootupper<-confboot[5]  
  ##observed value of correlation (entire dataset)
  corr_obs<- corr_boot$t0
  return(data_frame(
    corr=corr_obs[1],
    lower=bootlower,
    upper=bootupper
  ))
}
corrfxn_onlycorr_noboot<-function(data,varlist,simnum){
  'bivariate correlation function 
  used in permutation tests'
  d<-data[,varlist]
  corr=cor(d[1],d[2],
           method=corr_method,
           use = 'pairwise.complete.obs')
  return(corr)
}
i_corrfxn_onlycorr_boot<-function(varlist,data,i){
  'bivariate correlation for bootstrap function'
  'varlist=a character vector of 2 variables in data'
  d<-data[i,varlist]
  corr=cor(d[1],d[2],
           method=corr_method,
           use = 'pairwise.complete.obs')
  return(corr)
}
corrfxn_onlycorr_boot<-function(data,varlist,simnum){
  'gets bootstrapped stats for bivariate correlation
  of 2 variables specified'
  'used for Age and predicted'
  corr_boot<-boot::boot(
    data=data,
    statistic=i_corrfxn_onlycorr_boot,
    R=simnum,
    parallel = parallel_boot,
    ncpus=parallel_boot_cores,
    varlist=varlist
  )
  confboot<-boot::boot.ci(corr_boot,type="bca")$bca
  #bootstrap stats
  ##lower confidence interval (2.5% is default)
  bootlower<-confboot[4] 
  ##upper confidence interval (97.5% is default)
  bootupper<-confboot[5]
  ##observed value of correlation (entire dataset)
  corr_obs<- corr_boot$t0 
  return(data_frame(
    corr=corr_obs[1],
    lower=bootlower,
    upper=bootupper
  )
  )
}
med_fxn<- function(
  x,model.y_str,model.m_str,simnum){
  'get mediation stats for indirect effect 
  using following variables:
  Y=observ
  X=age
  C=Gender,wholebrain prediction
  M=predicted
  - (see mediation fxn doc for more info)
  - for permutation testing, 
  use this function but enter simnum=1'
  med_treat_str<- gsub('\\+.*','',model.m_str)
  treat_str <- gsub('^.*~','',med_treat_str)
  med_str <- gsub('~.*$','',model.m_str)
  model.y<-eval(parse(text=str_c('lm(',model.y_str,',x)')))
  model.m<-eval(parse(text=str_c('lm(',model.m_str,',x)')))
  results <- 
    mediation::mediate(
      model.m = model.m, 
      model.y = model.y, 
      treat=treat_str, 
      mediator=med_str, 
      boot=TRUE, 
      boot.ci.type='bca',
      sims=simnum
    )
  med_model <- summary(results)
  med_coefs<- data_frame(
    prop_ex=med_model[['d0']],
    prop_ex_p=med_model[['d0.p']],
    lower=med_model[['d0.ci']][[1]],
    upper=med_model[['d0.ci']][[2]]
  )
  return(med_coefs)
}
tryCatch_medfxn <- function(x,model.y_str,model.m_str){
  'to catch errors for mediation models with some NAs....'
  tryCatch(
    med_fxn(x,model.y_str,model.m_str),
    error=function(e){
      data_frame(prop_ex=NA,
                 prop_ex_p=NA,
                 lower=NA,
                 upper=NA)
    }
  )
}
map_corrstats<-function(
  df,
  permtest=FALSE,
  runboot=TRUE){
  if(runboot==TRUE){
    pcorrfxn_fxn<-pcorrfxn_boot
    corrfxn_onlycorr_fxn<-corrfxn_onlycorr_boot
    simnum=5000
  }else{
    pcorrfxn_fxn<-pcorrfxn_noboot
    corrfxn_onlycorr_fxn<-corrfxn_onlycorr_noboot
    simnum=1
  }
  #if(permtest==TRUE){print('starting calc')}
  df %>%
    mutate(wholebrain_obs_scale=scale(wholebrain_obs)) %>%
    summarise(
      corr_onlycorr=corrfxn_onlycorr_fxn(
        .[],
        varlist=c('predict_avg','observed_perf'),
        simnum=simnum
      ) %>% 
        list(),
      corr_cntrlgender=pcorrfxn_fxn(
        simnum=simnum,
        .[],
        str_lmcovariates='observed_perf~Gender'
      ) %>% 
        list(),
      corr_cntrlgender_wb=pcorrfxn_fxn(
        simnum=simnum,
        .[],
        str_lmcovariates='observed_perf~Gender+wholebrain_obs_scale'
      ) %>%
        list(),
      corr_onlycorr_age=corrfxn_onlycorr_fxn(
        simnum=simnum,
        .[],
        varlist=c('predict_avg','negAge')
      ) %>% 
        list(),
      corr_age_cntrlgender=pcorrfxn_fxn(
        simnum=simnum,
        .[],
        str_lmcovariates='negAge~Gender'
      ) %>% 
        list(),
      corr_age_cntrlgender_wb=pcorrfxn_fxn(
        simnum=simnum,
        .[],
        str_lmcovariates='negAge~Gender+wholebrain_obs_scale'
      ) %>%
        list(),
      med_cntrlgender=med_fxn(
        .[],
        simnum = simnum,
        model.y_str=str_c('observed_perf~',
                          'predict_avg+',
                          'negAge+',
                          'Gender'),
        model.m_str='predict_avg~negAge'
      ) %>%
        list(),
      med_cntrlgender_wb=med_fxn(
        .[],
        simnum = simnum,
        model.y_str=str_c('observed_perf~',
                          'negAge+',
                          'predict_avg+',
                          'Gender+wholebrain_obs_scale'),
        model.m_str='predict_avg~negAge'
      ) %>%
        list()
    )
}

####Make Parallel dataframe######
#loads make_cv_predict_df function
source('../scripts/make_cv_dfs.R')
make_parallel_df <- function(df,pkgs,fxns,cl){
  'preps df for parallel analyses with a list of
  packages and functions to import into each cluster core
  useful for bootstrap analyses and permutation testing'
  #cl=detectCores()-10
  clusters=create_cluster(cores = cl)
  groups <- nrow(df) %>% 
    rep(1:cl,length.out=.)
  str1 <- str_c(
    'df %>% ',
    'mutate(groups=groups) %>% ',
    'partition(groups,cluster=clusters) %>% ', 
    'cluster_library(pkgs) %>% '
  )
  str2 <- paste(sapply(fxns,FUN = function(x) 
    paste('cluster_assign_value(',"'",x,"',",x,')',' %>% ',
          sep='')),
    collapse='')
  str_all <- paste(str1,str2,sep = '') %>% 
    gsub('%>% $','',.)
  df_parallel<-eval(parse(text = str_all))
  return(df_parallel)
}
gather_npdata <- function(var){
  d<-npdata
  d_g<- d %>% 
    gather(key,value,-one_of(c(var))) %>%
    mutate(key2=var)
  names(d_g)[grepl(var,names(d_g))]<- 'value2'
  return(d_g)
}
cv_predict_avg_all_nested <- 
  make_cv_predict_df(
    cv_predict_filestr = '../data/cv_results/cv_predict_fpr_df_5fold.csv',
    npdata_filestr = '../data/np_filter_wb_gendernum.csv',
    permtest = FALSE) %>%
  group_by_at(cols_group_by) %>%
  nest()

cv_predict_avg_all_nested_parallel<- 
  make_parallel_df(
    df = cv_predict_avg_all_nested,
    pkgs = c('tidyverse'),
    fxns = ls()[str_detect(
      ls(),'corr|fxn|boot|make|cols_group_by|vars|map|scale'
    )],
    cl
    = 10
  )
npdata<-
  read_csv(
    '../data/np_filter_wb_gendernum.csv'
  ) %>% 
  select(matches('Age|ExF|Mem|wholebrain|Gender')) %>%
  mutate_at(vars(matches('wholebrain')),
            function(x) as.numeric(scale(x))) 
#combos of vars to run bootstraps for each combo in npdata
npdata_gather_nest<- bind_rows(list(
  gather_npdata('Age'),
  gather_npdata('wholebrain_area') %>% 
    filter(key!='Age'),
  gather_npdata('wholebrain_thickness') %>% 
    filter(!(key %in% c('Age','wholebrain_area'))),
  gather_npdata('Gender') %>% 
    filter(grepl('Mem|ExF',key))
)) %>%
  group_by(key,key2) %>%
  nest()

npdata_gather_nest_parallel<- 
  make_parallel_df(
    df = npdata_gather_nest,
    pkgs = c('tidyverse'),
    fxns = ls()[str_detect(
      ls(),'corr|fxn|boot|make|cols_group_by|vars|map|scale'
    )],
    cl
    = 10
  )

#####Run Models#####
boostrap_npdata_corrmodels<-
  npdata_gather_nest_parallel %>%
  mutate(corr_stats=data %>% map(function(x) {
    corrfxn_onlycorr_boot(
      data=x,
      varlist = c('value','value2'),
      simnum = 2000)
    }))

time1<-Sys.time()
bootstrap_corrmodels<-
  cv_predict_avg_all_nested_parallel %>%
  mutate(corr_stats=data %>%
           map(function(x){map_corrstats(x)})) %>%
  collect()
time2<-Sys.time()-time1
print(time2)
save(bootstrap_corrmodels,
     file = 
       '../data/stats/bootstrap_corrmodels_fromcv_predict.RData')
#####Permutation testing######
map_permstats<- function(
  x,
  npdata_filestr,
  pkgs=c('tidyverse'),
  fxns
){
  #print(fxns)
  df <- 
    make_cv_predict_df(
      cv_predict_filestr=x,
      npdata_filestr,
      permtest
    ) %>%
    group_by_at(
      vars(cviter,one_of(cols_group_by))
    ) %>%
    nest() %>%
    make_parallel_df(
      df=.,
      pkgs=c('tidyverse'),
      fxns = fxns,
      cl = 6
    )
  df_return <- df %>%
    mutate(
      corr_stats=data %>% 
        map(function(x){
          map_corrstats(df=x,permtest=TRUE,runboot=FALSE)})) %>%
    select(-data)
  return(df_return %>% collect())
}

permtest=TRUE
runboot=FALSE
env_vars<-ls()
scalevec<- function(x) as.numeric(scale(x))
permfile_list<-list.files(
  '../data/cv_results','cv_predict_fpr_df_5fold_perm',
  full.names = TRUE)

permfile_df<- data_frame(
  filestr=permfile_list
  ) %>%
  mutate(
    permfilenum=permfile_list %>% 
      str_remove_all('(^.*permtest|\\.csv$)') %>%
      as.integer()
  ) %>%
  arrange(desc(filestr))
  
time1<-Sys.time()
#get individual permutation stats#
perm_summaries<-
  permfile_df %>%
  group_by(permfilenum) %>%
  mutate(
    perm_stats=filestr %>% 
      map_permstats(
        x=.,
        npdata_filestr = 
          '../data/np_filter_wb_gendernum.csv',
        pkgs=c('tidyverse'),
        fxns = env_vars
      ) %>% 
      list()
    )
time2<-Sys.time()-time1
#take around 1 and half hours for 1000 permutations
print(time2)
save(perm_summaries,
     file = 
       '../data/perm_summaries_fromcv_predict.RData'
     )
#get and unnest empirical stats
empirical_summaries_unnest <-
  cv_predict_avg_all_nested_parallel %>%
  mutate(corr_stats=data %>%
           map(function(x){
             map_corrstats(
               x,runboot = FALSE)
             }
             )) %>%
  collect() %>%
  #filter(type_alpha=='average') %>%
  select(-data) %>%
  unnest(corr_stats) %>% 
  ungroup() %>% 
  unnest(.sep = '_') %>%
  select(-matches(
    '_p$|upper|lower|groups|pkldir|index|covar'
    )) %>%
  gather(analysis,value,
         -matches(
           paste(cols_group_by,collapse = '|')
           )) %>%
  mutate(analysis=gsub('_corr$','',analysis))
#compute permutation numbers 
    #each file consisted of 100 permutations
#join with empirical stats
#compute if perm stat >= empirical
perm_summaries_unnest<-
  perm_summaries %>% 
  unnest(perm_stats) %>% 
  ungroup() %>% 
  mutate(permnum=cviter*permfilenum) %>%
  unnest(corr_stats) %>%
  ungroup() %>% 
  unnest(.sep = '_') %>%
  select(-matches(
    '_p$|upper|lower|groups|pkldir|index|filestr|covar'
    )) %>%
  gather(analysis,value,
         -matches(paste(cols_group_by,collapse = '|')),
         -matches('cviter|perm'))
perm_summaries_w_empir<-
  perm_summaries_unnest %>%
  left_join(empirical_summaries_unnest,
             by=c(cols_group_by[1:5],'analysis'),
             suffix=c('_perm','_empirical')) %>%
  mutate(perm_greaterthanequal_emp=value_perm>=value_empirical)
save(perm_summaries_w_empir,
     file = 
       '../data/stats/perm_summaries_unnest_withemp_fromcv_predict.RData')
perm_pvals<-
  perm_summaries_w_empir %>% 
  group_by_at(c(cols_group_by[1:5],'analysis')) %>% 
  summarise(
    perm_greaterthanequal_total=perm_greaterthanequal_emp %>% 
      sum(),
    perm_total=n(),
    pval=perm_greaterthanequal_total/perm_total,
    perm_mean=value_perm %>% mean(na.rm=TRUE),
    perm_sd=value_perm %>% sd(na.rm=TRUE))

save(perm_pvals,
     file = 
       '../data/stats/perm_pvals_fromcv_predict.RData')

perm_pvals %>% 
  filter(type_alpha=='average') %>% 
  select(-perm_greaterthanequal_total,-perm_total) %>%
  ggplot(aes(x=network,y=perm_mean))+
  facet_grid(analysis~np_measure_name+measure,scale='free')+   
  geom_bar(
    stat='identity',
    aes(fill=network),
    position = position_dodge(width = 1)
    )+
  geom_errorbar(
    aes(ymax = perm_mean+perm_sd, 
        ymin = perm_mean-perm_sd,
        group=type_alpha),
    width=.2,color='black',
    position = position_dodge(width = 1)
  )
 # geom_label(aes(y=.75,label=value),size=3)


###Figures#####
make_ggplot <- function(bootstrap_df,modeltype,sizenum=1){
  'ggplot code for all correlation stats'
  graph_df<-bootstrap_df %>%
    ungroup() %>%
    unnest(corr_stats) %>%
    select(np_measure_name:threshold,one_of(modeltype)) %>%
    unnest() %>%
    mutate(np_measure_name=np_measure_name %>% 
             str_replace('ExFunction','Executive Function'),
           type_alpha=type_alpha %>%
             str_replace_all(c('alpha: '='','0\\.'='\\.')),
           measure=measure %>% 
             str_replace_all(c('area'='Surface Area',
                               'thickness'='Thickness')),
           network=network %>% 
             str_replace('network',''))
  names(graph_df) <- 
    str_replace(names(graph_df),'prop_ex','corr')
  graph_df %>%
    ggplot(aes(x=network,y=corr))+
    geom_bar(
      stat='identity',
      aes(group=type_alpha,
          fill=network,
          color=type_alpha),
      size=sizenum,
      position = position_dodge(width = 1)
    )+
    geom_errorbar(
      aes(ymax = upper, 
          ymin = lower,
          group=type_alpha),
      width=.2,color='black',
      position = position_dodge(width = 1)
    ) +
    facet_grid(np_measure_name~measure) +
    guides(size=FALSE) +
    scale_fill_manual(values=NetworkColors,guide=FALSE)+
    scale_color_manual(values = TypeAlphaColors,guide=FALSE) +
    theme_bw()+
    theme(
      text=element_text(family="Helvetica"),
      panel.grid=element_blank(),
      strip.text= element_text(size=10),
      strip.background = element_rect(colour="black",fill="white"),
      axis.text=element_text(size=10),
      axis.text.x = element_text(),
      legend.position=c(.8,.85)
    )
}
#####Yeo Functional Network Parcellation######
#####Cross validation workflow#######
#####Bivariate Correlations####
corr_ggplots<- 
  cv_predict_avg_all_nested %>%
  unnest() %>%
  filter(type_alpha=='average') %>%
  ggplot(aes(x=predict_avg,y=observed_perf)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  facet_grid(measure+np_measure_name~network)

bootstrap_corrmodels %>%
  make_ggplot(modeltype = 'corr_onlycorr') +
  labs(x='Functional Networks (Cortical Morphometry)',
       y='Bivariate Correlation Coefficient')+
  scale_y_continuous(limits=c(-0.3,.52),
                     breaks=round(seq(-.3,.5,.1),3))

ggsave("Correlations.png",width=6,height=6)

#####Correlations (Cntrl Gender)####
bootstrap_corrmodels %>% 
  make_ggplot(modeltype = 'corr_cntrlgender') +
  labs(x='Functional Networks (Cortical Morphometry)',
       y='Partial Correlation Coefficient')+
  scale_y_continuous(limits=c(-0.32,.52),
                     breaks=round(seq(-.3,.5,.1),3))
ggsave("Correlations_CntrlGender.png",width=6,height=6)
#####Correlations (Cntrl Whole Brain and Gender)####
bootstrap_corrmodels %>% 
  make_ggplot(modeltype = 'corr_cntrlgender_wb') +
  labs(x='Functional Networks (Cortical Morphometry)',
       y='Partial Correlation Coefficient')+
  scale_y_continuous(limits=c(-0.32,.52),
                     breaks=round(seq(-.3,.5,.1),3))

ggsave("Correlations_GenderWholebrain.png",width=2.4,height=6)

####Bivariate Correlations with Age#######
bootstrap_corrmodels %>% 
  make_ggplot(modeltype = 'corr_onlycorr_age') +
  labs(x='Functional Networks (Cortical Morphometry)',
       y='Partial Correlation Coefficient')+
  scale_y_continuous(limits=c(-0.32,.6),
                     breaks=round(seq(-.3,.5,.1),3))
ggsave("CorrelationsAge.png",width=6,height=6)
####Correlations with Age (Cntrl Gender)#####
bootstrap_corrmodels %>% 
  make_ggplot(modeltype = 'corr_age_cntrlgender') +
  labs(x='Functional Networks (Cortical Morphometry)',
       y='Partial Correlation Coefficient')+
  scale_y_continuous(limits=c(-0.32,.52),
                     breaks=round(seq(-.3,.5,.1),3))

ggsave("CorrelationsAge_GenderWholebrain.png",width=2.4,height=6)
####Correlations with Age (Cntrl Gender + WholeBrain)#####
bootstrap_corrmodels %>% 
  make_ggplot(modeltype = 'corr_age_cntrlgender_wb') +
  labs(x='Functional Networks (Cortical Morphometry)',
       y='Partial Correlation Coefficient')+
  scale_y_continuous(limits=c(-0.32,.52),
                     breaks=round(seq(-.3,.5,.1),3))

ggsave("CorrelationsAge_GenderWholebrain.png",width=2.4,height=6)

####Mediation of Age~cognition#####
#http://data.library.virginia.edu/introduction-to-mediation-analysis/
#controlling for gender
bootstrap_corrmodels %>% 
  ungroup() %>%
  make_ggplot(sizenum=0,modeltype='med_cntrlgender') +
  labs(x='Functional Networks (Cortical Morphometry)',
       y='Mediation Effect (controlling for gender)')

ggsave("Mediation_Gender.png",width=6,height=6)
#controlling for gender+wholebrain
bootstrap_corrmodels %>% 
  ungroup() %>%
  make_ggplot(sizenum=0,modeltype='med_cntrlgender_wb') +
  labs(x='Functional Networks (Cortical Morphometry)',
       y='Mediation Effect')
ggsave("Mediation_Gender_wholebrain.png",width=6,height=6)
#####Number of vertices per Network in Whole Brain Analysis###
######Whole Brain Analysis Figure#####
######*********Percentile Selection*******########
###Make Parallel DataFrame#####
cv_predict_avg_all_nested_extreme <- 
  make_cv_predict_df(
    cv_predict_filestr = '../data/cv_results/cv_predict_fpr_df_extremepvals_5fold.csv',
    npdata_filestr = '../data/np_filter_wb_gendernum.csv',
    permtest = FALSE) %>%
  group_by_at(cols_group_by) %>%
  nest()
cv_predict_avg_all_nested_extreme_parallel<- 
  make_parallel_df(
    df = cv_predict_avg_all_nested_extreme,
    pkgs = c('tidyverse'),
    fxns = ls()[str_detect(
      ls(),'corr|fxn|boot|make|cols_group_by|vars|map|scale'
    )],
    cl = 10
  )
#####Run models######
time1<-Sys.time()
bootstrap_corrmodels_extreme<-
  cv_predict_avg_all_nested_extreme_parallel %>%
  mutate(corr_stats=data %>%
           map(function(x){map_corrstats(x)})) %>%
  collect()
time2<-Sys.time()-time1
print(time2)
####Bivariate Correlations####
bootstrap_corrmodels_extreme %>%
  make_ggplot(modeltype = 'corr_onlycorr') +
  labs(x='Functional Networks (Cortical Morphometry)',
       y='Bivariate Correlation Coefficient')+
  scale_y_continuous(limits=c(-0.3,.52),
                     breaks=round(seq(-.3,.5,.1),3))

ggsave("Correlations.png",width=6,height=6)

####Correlations (Cntrl Gender)######
####Correlations (Cntrl Gender+Wholebrain)#######
bootstrap_corrmodels_extreme %>% 
  make_ggplot(modeltype = 'corr_cntrlgender_wb') +
  labs(x='Functional Networks (Cortical Morphometry)',
       y='Partial Correlation Coefficient')+
  scale_y_continuous(limits=c(-0.32,.52),
                     breaks=round(seq(-.3,.5,.1),3))
####Mediation of Age~Cognition####
#controlling for gender+wholebrain
bootstrap_corrmodels_extreme %>% 
  ungroup() %>%
  make_ggplot(sizenum=0,modeltype='med_cntrlgender_wb') +
  labs(x='Functional Networks (Cortical Morphometry)',
       y='Mediation Effect')
