library(tidyverse)
'this script takes a measure from each task and explores PCA structure. 
Used the 2 component structure as it is:
1.) highly interpretable (Declarative Memory vs. Executive Function)
2.) inflection point on scree plot

- note: excluded MEPS as it is not widely used in cognitive research 
and
task switching to give a good balance of # of EF vs. Memory tasks
 '
'on bottom is PCA exploring structure of just PCA variables
which shows two variables that emerge from EF:
1.) Executive Function
2.) Processing Speed
- interpretation emerges from Trail Making task loading
higher on EF when Part B used and higher on PSpeed when 
Part A used
'
#####import behavioral data######
taskdir=str_c(
  '/Volumes/Users/mbkranz/projects/',
  'ACT_Freesurfer_NewProc/scripts/Rproject_paper'
  )
setwd(taskdir)

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

####import behavioral data#####
dat.act.pre<- 
  readxl::read_excel("ACT_NP.xlsx",
                  sheet='All Subjects')
ACT_PreInterventionNP.sub<-
  dat.act.pre %>% 
  filter(MRIPre==1)
####Principal Component Analysis of Tasks in ACT Battery#####
pr.model<-
  ACT_PreInterventionNP.sub %>% 
  transmute(
    Flanker_Pre_Incongruent_RT,
    #Time_BminusA_Pre=Time_Trial_B_Pre-Time_Trial_A_Pre,
    Time_Trial_B_Pre,
    SPWMall=.[] %>% select(matches("SWM_ACC")) %>% 
      rowMeans(na.rm=TRUE),
    ACC_2back_Pre,
    CategoryFluency=CategoryFluency_Fruits_Veg_Pre+
      CategoryFluency_Animals_Pre,
    Dprime_Pre, 
    Imove_Swaps_Pre,
    DotComparison_RT_Pre,
    #TS_Switch_RT_Pre,
    #LocalSwitchCost_RT_Pre,
    Digit_Symbol_Correct_Pre, 
    Delayed_Free_Recall_Pre, 
    #MEPS_Total_Pre, 
    Story_Unit_Recall_Delayed_Pre) %>% 
  #make RT/error measures negative so positive means good for all measures
  mutate(Flanker_Pre_Incongruent_RT=-Flanker_Pre_Incongruent_RT,
         Time_Trial_B_Pre=-Time_Trial_B_Pre,
         DotComparison_RT_Pre=-DotComparison_RT_Pre,
         Imove_Swaps_Pre=-Imove_Swaps_Pre) %>%
  windsorizefxn(.,names(.),replaceval=TRUE) %>%
  as_data_frame()

###histograms of vars####
pr.model %>%
  gather(task,value) %>%
  group_by(task) %>%
  # mutate(value_scale=scale(value)) %>%
   #ggplot(aes(value_scale)) +
  ggplot(aes(value)) +
  geom_histogram() +
  facet_wrap(~task,ncol = 3,scales = 'free')

psych::scree(pr.model)

pr.model.comp2<-psych::principal(
  pr.model, 
  nfactors= 2 , 
  rotate= "varimax" , covar= FALSE , scores=TRUE)
pr.model.comp3<-psych::principal(
  pr.model, 
  nfactors= 3 , 
  rotate= "varimax" , covar= FALSE , scores=TRUE )
pr.model.comp4<-psych::principal(
  pr.model, 
  nfactors= 4 , 
  rotate= "varimax" , covar= FALSE , scores=TRUE )

np.all<- 
  ACT_PreInterventionNP.sub %>% 
  as_data_frame() %>%
  transmute(subs=as.character(Subject_ID),
            Age=Age_at_start_of_study,
            Gender=ifelse(Gender=="Female","Female","Male"),
            ExFunction=pr.model %>% 
              select(matches(
                "Flanker_Pre_Incongruent_RT|Time_Trial_B|SPWM|ACC_2back|DotComp|Digit_Symbol"
                )) %>%
              scale() %>% 
              rowMeans(na.rm=TRUE),
            Memory=pr.model %>% select(matches(
              "CategoryFluency|Dprime|Imove|Recall"
              )) %>% 
              scale() %>% 
              rowMeans(na.rm=TRUE))
names(np.all)<-gsub("_Pre","",names(np.all))

###histograms of construct scores####
np.all %>%
  select(-subs,-Gender) %>%
  gather(task,value) %>%
  group_by(task) %>%
  # mutate(value_scale=scale(value)) %>%
  #ggplot(aes(value_scale)) +
  ggplot(aes(value)) +
  geom_histogram() +
  facet_wrap(~task,ncol = 3,scales = 'free')

###Executive Function PCA#####
efonly_pr.model<-
  dat.act.pre %>% 
  filter(MRIPre==1) %>%
  transmute(Flanker_Pre_Incongruent_RT,
            Time_Trial_A_Pre,
            Time_Trial_B_Pre,
            SPWMall=.[] %>% select(matches("SWM_ACC")) %>% 
              rowMeans(na.rm=TRUE),
            ACC_2back_Pre,
            Digit_Symbol_Correct_Pre,
            DotComparison_RT_Pre) %>% 
  #make RT/error measures negative so + = good for all measures
  mutate(Flanker_Pre_Incongruent_RT=-Flanker_Pre_Incongruent_RT,
         Time_Trial_B_Pre=-Time_Trial_B_Pre,
         Time_Trial_A_Pre=-Time_Trial_A_Pre,
         DotComparison_RT_Pre=-DotComparison_RT_Pre) %>%
  as_data_frame()
windsorizefxn(.,names(.),replaceval=TRUE)

efonly_pr.model_TrailsA<-psych::principal(
  efonly_pr.model %>% select(-matches('Trial_B|^subs|Age|Gender')), 
  nfactors= 2 , 
  rotate= "varimax" , 
  covar= FALSE , 
  scores=TRUE)
efonly_pr.model_TrailsB<-psych::principal(
  efonly_pr.model %>% select(-matches('Trial_A|^subs|Age|Gender')), 
  nfactors= 2 , 
  rotate= "varimax" , 
  covar= FALSE , 
  scores=TRUE)

np.ef<- 
  ACT_PreInterventionNP.sub %>% 
  as_data_frame() %>%
  transmute(subs=as.character(Subject_ID),
            Age=Age_at_start_of_study,
            Gender=ifelse(Gender=="Female","Female","Male"),
            EF_ExFunction=efonly_pr.model %>% 
              select(matches("Time_Trial_B|SPWM|ACC_2back")) %>%
              scale() %>% rowMeans(na.rm=TRUE),
            EF_Speed=efonly_pr.model %>% 
              select(matches(
              "Flanker_Pre_Incongruent_RT|Time_Trial_A|DotComp|Digit_Symbol"
              )) %>% 
              scale() %>% rowMeans(na.rm=TRUE))
names(np.ef)<-gsub("_Pre","",names(np.ef))
write_csv(np.ef,'np_ef.csv')


