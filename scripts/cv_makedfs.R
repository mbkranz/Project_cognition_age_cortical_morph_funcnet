library(tidyverse)
library(magrittr)
setwd('../data')

cv_fpr_predict_file <- 'cv_results/cv_predict_fpr_df_5fold.csv'
cv_fpr_results_file <- 'cv_results/cv_results_fpr_df_5fold.csv'

cv_npdata_predict_file <- 'cv_results/cv_predict_fpr_df_5fold_npdata.csv'
cv_npdata_results_file <- 'cv_results/cv_results_fpr_df_5fold_npdata.csv'

npdata_file <- 'np_filter_wb_gendernum.csv'

cv_fpr_covars_results_file <- 'cv_results/cv_results_fpr_df_covarsTESTING.csv'

####import and format data####
#cv results with covariates
cv_results_fpr_covars <- read_csv(cv_fpr_covars_results_file) %>%
  bind_rows(read_csv(cv_npdata_results_file))

cv_results_fpr_covars %<>%
  mutate(measure=pkldir %>%
           str_extract('area|thickness'),
         network=pkldir %>% 
           str_extract('network[:digit:]'),
         covariates=ifelse(grepl('npdata',pkldir),
                           str_replace(pkldir,'npdata_',''),
                           covariates)) %>%
  filter(grepl('Gender',covariates))

cv_results_covars <- cv_results_fpr_covars

#NP data and add 
npdata <- read_csv(npdata_file) %>%
  select(-matches('CVLT')) %>%
  gather(np_measure_name,observed_perf,
         -subs,-Gender,-Age,-matches('wholebrain')) %>%
  gather(measure,wholebrain_obs,
         -subs,-Gender,-Age,-np_measure_name,-observed_perf) %>%
  mutate(measure= measure %>%
           str_extract('thickness|area'))

#False discovery rate
cv_predict_fpr<-read_csv(cv_fpr_predict_file,
                         col_types = cols(
                           threshold=col_double()))
cv_results_fpr<-read_csv(cv_fpr_results_file,
                         col_types = cols(
                           threshold=col_double()))

#some splits did not select features (.001 and .01 cut offs) 
#so imputed mean cog perf (i.e., 0)
cv_predict_avg_fpr <- cv_predict_fpr %>%
  group_by(np_measure_name,
           pkldir,
           covariates,
           threshold) %>%
  mutate(prediction_wzero=ifelse(
    is.na(prediction),
    0,prediction)) %>% 
  group_by(np_measure_name,
           pkldir,
           covariates,
           threshold,
           subs) %>%
  summarise(predict_avg=mean(prediction_wzero),
            netpredict_nanum=prediction %>% 
              is.na() %>% 
              sum(),
            type='alpha')

#some cv splits did not have any vertices selected. 
#made a cut off of 100 
to_include <- cv_predict_avg_fpr %>% 
  group_by(np_measure_name,pkldir,covariates,threshold) %>% 
  summarise(nanum=sum(netpredict_nanum)) %>% 
  filter(nanum<100) %>%
  unite(index,
        np_measure_name,pkldir,covariates,threshold,
        sep='_',remove=FALSE)

cv_predict_avg_fpr %>% 
  group_by(np_measure_name,pkldir,covariates,threshold) %>% 
  summarise(nanum=sum(netpredict_nanum)) %>% 
  write_csv(
    'cv_results/number_of_subs_without_featuresselected.csv'
    )

cv_predict_avg <- cv_predict_avg_fpr %>%
  mutate(type_alpha=ifelse(type=='alpha',
                           paste('alpha:',threshold),
                           'percentile'),
         measure=pkldir %>%
           str_extract('area|thickness'),
         network=pkldir %>% 
           str_extract('network[:digit:]')) %>%
  left_join(npdata,by=c('subs','np_measure_name','measure'))

scalevec<- function(x) as.numeric(scale(x))
cv_predict_avg_alpha_scaled <- cv_predict_avg %>% 
  unite(index,
        np_measure_name,pkldir,covariates,threshold,
        sep='_',remove=FALSE) %>%
  filter(subs %in% (npdata[['subs']] %>% unique())) %>%
  filter(!is.na(predict_avg) & !is.na(observed_perf)) %>%
  filter(sd(predict_avg)>0) %>%
  mutate(Gender=as.factor(Gender)) %>%
  dplyr::select(-netpredict_nanum) %>%
  group_by(index) %>%
  #mutate_at(vars(observed_perf,predict_avg),funs(scalevec)) %>%
  mutate(wholebrain_pred= NA)

cv_predict_avg_allthresh <- cv_predict_avg_alpha_scaled %>%
  filter(threshold<1) %>%
  group_by(np_measure_name,
           Gender,
           Age,
           pkldir,
           covariates,
           measure,
           network,
           subs) %>%
  summarise(observed_perf= observed_perf %>% mean(),
            wholebrain_pred= NA, 
            wholebrain_obs= wholebrain_obs %>% mean(),
            predict_avg= predict_avg %>% mean(),
            threshold=NA,
            type_alpha='average') %>% 
  ungroup()

str_ext_del <- function(str,ex_str,del_str)
{
  str %>% 
    str_extract(ex_str) %>% 
    str_replace(del_str,'')
  }

cv_predict_avg_all_nested <- cv_predict_avg_alpha_scaled %>%
  bind_rows(cv_predict_avg_allthresh) %>%
  ungroup() %>%
  mutate(negAge=-Age,
         Gender=ifelse(Gender=='1',1,0)) %>%
  group_by_at(vars(np_measure_name,
                   pkldir,
                   network,
                   measure,
                   type_alpha,
                   threshold,
                   covariates)) %>% 
  nest() %>%
  ungroup()
