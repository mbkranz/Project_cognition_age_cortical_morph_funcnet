library(tidyverse)
####This script makes a figure displaying the number of vertices selected at the two thresholds (3.3 and 4.25)
####To be displayed with surface image of 3.3 threshold? Or potentially masks overlayed on annotation?
#To put CI bars on plot: extract lower and upper CIs from BCa string and make numeric
NetworkColors<- c(rgb(120,18,134,maxColorValue = 255),
                  rgb(70,130,180,maxColorValue = 255),
                  rgb(0,118,14,maxColorValue = 255),
                  rgb(196,58,250,maxColorValue = 255),
                  rgb(217, 201, 144,maxColorValue = 255),
                  rgb(230,148,34,maxColorValue = 255),
                  rgb(205,62,78,maxColorValue = 255))
NetworkLabels=c('1'='Visual',
                '2'='Somato-\nmotor',
                '3'='Dorsal\nAttention',
                '4'='Ventral\nAttention',
                '5'='Limbic',
                '6'='Fronto-\nparietal',
                '7'='Default')
cv_perm<- read_csv(str_c('../data/',
                   'wholebrain_bootstrap/',
                   'wholebrain_bootstrap_bootscore_permtests.csv'))
cv_perm_clean <- cv_perm %>%
  rename(threshold=X3) %>%
  mutate(anal= pklfile %>%
           str_extract('regression|mediation'),
         measure= pklfile %>%
           str_extract('thickness|area') %>%
           str_replace_all(c('thickness'='Thickness','area'='Surface Area')),
         np_measure_name= pklfile %>%
           str_extract('ExFunction|Memory') %>%
           str_replace('ExFunction','Executive Function'),
         annot_label= network %>%
           as.character()) %>%
  select(-pklfile)
graph_fxn <- function(analtype)
  {cv_perm_clean %>%
    filter(anal==analtype) %>%
    ggplot(aes(x=annot_label,
           fill=annot_label,
           color=annot_label,
           group=interaction(threshold,annot_label))) +
    geom_bar(stat='identity',
             position='dodge',
             aes(y=empirical_vals))+
    geom_bar(stat='identity',
             position = 'dodge',
             alpha=.2,fill='black',
             color='grey',
             aes(y=perm_meansigvertices)) +
    facet_grid(np_measure_name~measure) +
    scale_fill_manual(values=NetworkColors,guide=FALSE)+
    scale_color_manual(values=NetworkColors,guide=FALSE)+
    labs(y='Number of Vertices',x='Network') +
    theme_bw()+
    theme(text=element_text(family="Helvetica"),
          panel.grid=element_blank(),
          strip.text= element_text(size=12,face = 'bold'),
          strip.background = element_rect(colour="black",fill="white"),
          axis.text=element_text(size=10),
          axis.title=element_text(face='bold'),
          axis.text.x = element_text(angle = 0, hjust = .5),
          plot.title = element_text(hjust = 0.5))
  }

graph_fxn('mediation')+
  ggtitle('Mediation of the Age ~ Cognition Relationship')
ggsave('../figures/numvertices_mediation_bootstrap.png',dpi = 300,width = 8,height = 6)
graph_fxn('regression') +
  ggtitle('Relationship with Cognition')
ggsave('../figures/numvertices_regression_bootstrap.png',dpi = 300,width = 8,height = 6)
