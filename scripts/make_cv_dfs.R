make_cv_predict_df <- function(cv_predict_filestr,
                       npdata_filestr,
                       permtest=FALSE){
  require(tidyverse)
  require(magrittr)
  #as each iteration is random and a separate perm
  if(permtest==TRUE){
    include_cviter=TRUE
    id_c<-c('subs','cviter')
  }else{
    include_cviter=FALSE
    id_c<-c('subs')
  }
  #NP data and add
  npdata <- read_csv(npdata_filestr) %>%
    select(-matches('CVLT')) %>%
    gather(np_measure_name,observed_perf,
           -subs,-Gender,-Age,-matches('wholebrain')) %>%
    gather(measure,wholebrain_obs,
           -subs,-Gender,-Age,-np_measure_name,-observed_perf) %>%
    mutate(measure= measure %>%
             str_extract('thickness|area'))
  #False discovery rate
  cv_predict_fpr<-read_csv(
    cv_predict_filestr,
    col_types = cols(threshold=col_double(),
                     prediction=col_double())
    )
  #some splits did not select features (.001 and .01 cut offs) 
  #so imputed mean cog perf (i.e., 0)
  cv_predict_impute0_fpr <- 
    cv_predict_fpr %>%
    mutate(prediction_wzero=ifelse(
      is.na(prediction),
      0,prediction)) 

  if (include_cviter==TRUE){
    cv_predict_avg_fpr <-
    cv_predict_impute0_fpr %>%
      transmute(
        np_measure_name,
        pkldir,
        covariates,
        threshold,
        subs,
        cviter,
        predict_avg=prediction_wzero,
        netpredict_nanum=prediction %>%
          is.na(),
        type='alpha')
    }else{
      cv_predict_avg_fpr <- 
      cv_predict_impute0_fpr %>%
        group_by_at(
          .vars=vars(
            np_measure_name,
            pkldir,
            covariates,
            threshold,
            subs)) %>%
        summarise(
          predict_avg=mean(prediction_wzero),
          netpredict_nanum=prediction %>% 
            is.na() %>% 
            sum(),
          type='alpha')
        }

  # #some cv splits did not have any vertices selected. 
  # #made a cut off of 100 
  # to_include <- cv_predict_avg_fpr %>% 
  #   group_by(np_measure_name,pkldir,covariates,threshold) %>% 
  #   summarise(nanum=sum(netpredict_nanum)) %>% 
  #   filter(nanum<100) %>%
  #   unite(index,
  #         np_measure_name,pkldir,covariates,threshold,
  #         sep='_',remove=FALSE)
  # cv_predict_avg_fpr %>% 
  #   group_by(np_measure_name,pkldir,covariates,threshold) %>% 
  #   summarise(nanum=sum(netpredict_nanum)) %>% 
  #   write_csv(
  #     'cv_results/number_of_subs_without_featuresselected.csv'
  #     )
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
  cv_predict_avg_alpha_scaled <- 
    cv_predict_avg %>% 
    unite(index,
          np_measure_name,pkldir,covariates,threshold,
          sep='_',remove=FALSE) %>%
    filter(subs %in% (npdata[['subs']] %>% 
                        unique())) %>%
    filter(!is.na(predict_avg) & !is.na(observed_perf)) %>%
    filter(sd(predict_avg)>0) %>%
    mutate(Gender=as.factor(Gender)) %>%
    dplyr::select(-netpredict_nanum) %>%
    #group_by(index) %>%
    #mutate_at(.vars=vars(observed_perf,predict_avg),funs(scalevec)) %>%
    mutate(wholebrain_pred= NA)
  
  if(permtest==TRUE){print(str_c(
    'starting avg across thresh calc:',cv_predict_filestr
    ))}
  cv_predict_avg_allthresh <- 
    cv_predict_avg_alpha_scaled %>%
    filter(threshold<1) %>%
    ungroup() %>%
    group_by_at(
      .vars=vars(
        np_measure_name,
        Gender,
        Age,
        pkldir,
        covariates,
        measure,
        network,
        one_of(id_c)
        )) %>%
    summarise(observed_perf= observed_perf %>% 
                mean(),
              wholebrain_pred= NA, 
              wholebrain_obs= wholebrain_obs %>% 
                mean(),
              predict_avg= predict_avg %>%
                mean(),
              threshold=NA,
              type_alpha='average') %>% 
    ungroup()
  if(permtest==TRUE){print(str_c(
    'finished avg across thresh calc:',cv_predict_filestr
  ))}
  str_ext_del <- function(str,ex_str,del_str)
  {
    str %>% 
      str_extract(ex_str) %>% 
      str_replace(del_str,'')
  }
  #not nested (perm testing doesn't require bootstrap)
  cv_predict_avg_all <- 
    cv_predict_avg_alpha_scaled %>%
    bind_rows(cv_predict_avg_allthresh) %>%
    ungroup() %>%
    mutate(negAge=-Age,
           Gender=ifelse(Gender=='1',1,0))
  if(permtest==TRUE){print(str_c(
    'COMPLETED makecvs_dfs:',cv_predict_filestr
  ))}
  return(cv_predict_avg_all)
  }
           
make_cv_results_df <- function(cv_results_filestr,
                               permtest=FALSE){
  require(tidyverse)
  require(magrittr)
  #as each iteration is random and a separate perm
  if(permtest==TRUE){
    include_cviter=TRUE
    id_c<-c('cviter')
  }else{
    include_cviter=FALSE
    id_c<-c('thisvarisjustplaceholder')
  }
  #False discovery rate
  cv_results_fpr<-read_csv(
    cv_results_filestr,
    col_types = cols(threshold=col_double())
  )
  #some splits did not select features (.001 and .01 cut offs) 
  #so imputed mean cog perf (i.e., 0)
  cv_results_impute0_fpr <- 
    cv_results_fpr %>%
    mutate(coef=ifelse(
      is.na(coef),
      0,coef)) 
  
  if (include_cviter==TRUE){
    cv_results_avg_fpr <-
      cv_results_impute0_fpr %>%
      select(
        np_measure_name,
        pkldir,
        covariates,
        threshold,
        cviter,
        coef,
        num_features_selected,
        type='alpha'
        )
  }else{
    cv_results_avg_fpr <- 
      cv_results_impute0_fpr %>%
      mutate(
        nofeatures_selected=(num_features_selected==0)
        ) %>%
      group_by_at(
        .vars=vars(
          np_measure_name,
          pkldir,
          covariates,
          threshold
          )) %>%
      summarise_all(funs(mean=mean,sd=sd),na.rm=TRUE)
  }
  
  # #some cv splits did not have any vertices selected. 
  # #made a cut off of 100 
  # to_include <- cv_results_avg_fpr %>% 
  #   group_by(np_measure_name,pkldir,covariates,threshold) %>% 
  #   summarise(nanum=sum(netpredict_nanum)) %>% 
  #   filter(nanum<100) %>%
  #   unite(index,
  #         np_measure_name,pkldir,covariates,threshold,
  #         sep='_',remove=FALSE)
  # cv_results_avg_fpr %>% 
  #   group_by(np_measure_name,pkldir,covariates,threshold) %>% 
  #   summarise(nanum=sum(netpredict_nanum)) %>% 
  #   write_csv(
  #     'cv_results/number_of_subs_without_featuresselected.csv'
  #     )
  if(permtest==TRUE){print(str_c(
    'starting avg across thresh calc:',cv_results_filestr
  ))}
  cv_results_avg_allthresh <- 
    cv_results_avg_fpr %>%
    filter(threshold<1) %>%
    ungroup() %>%
    group_by_at(
      .vars=vars(
        np_measure_name,
        pkldir,
        covariates,
        matches(id_c)
      )) %>%
    summarise_all(mean)
  
  if(permtest==TRUE){print(str_c(
    'finished avg across thresh calc:',cv_results_filestr
  ))}
  str_ext_del <- function(str,ex_str,del_str)
  {
    str %>% 
      str_extract(ex_str) %>% 
      str_replace(del_str,'')
  }
  #not nested (perm testing doesn't require bootstrap)
  cv_results_avg_all <- 
    cv_results_avg_fpr %>%
    bind_rows(cv_results_avg_allthresh) %>%
    ungroup() %>%
    mutate(
      type='alpha',
      type_alpha=ifelse(
        type=='alpha',
        paste('alpha:',threshold),
        'percentile'
        ),
      measure=pkldir %>%
        str_extract('area|thickness'),
      network=pkldir %>% 
        str_extract('network[:digit:]')) %>%
    select(matches('np_meas|network|measure|thresh'),
           matches('coef|test_rsq|features')) 

  if(permtest==TRUE){print(str_c(
    'COMPLETED make_cv_results_df:',cv_results_filestr
  ))}
  return(cv_results_avg_all)
}


