## ----setup, include = FALSE---------------------------------------------------
  library(dplyr)  
  library(tibble)
  library(rlang)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(echoice2)
  knitr::opts_chunk$set(fig.align  = "center",
                        fig.height = 3.5,
                        warning    = FALSE,
                        error      = FALSE,
                        message    = FALSE)

## -----------------------------------------------------------------------------
  data(icecream)
  icecream %>% head

## -----------------------------------------------------------------------------
icecream %>% ec_summarize_attrlvls

## -----------------------------------------------------------------------------
icecream %>% 
  arrange(id,task,alt) %>%
  group_by(id,task) %>%  
  summarise(units=sum(x)) %>%
  group_by(id) %>% summarise(mean_units=mean(units)) %>%
  ggplot(aes(x=mean_units)) + geom_histogram()

## -----------------------------------------------------------------------------
#select holdout tasks
  set.seed(8438453)
      ho_tasks=
      icecream %>%
        distinct(id,task) %>%
        mutate(id=as.integer(id))%>%
        group_by(id) %>%
        summarise(task=sample(task,1), .groups = 'drop')
  set.seed(NULL)

#split data
  ice_cal= icecream %>% mutate(id=as.integer(id)) %>%
    anti_join(ho_tasks, by=c('id','task'))
  
  ice_ho= icecream %>% mutate(id=as.integer(id)) %>%
    semi_join(ho_tasks, by=c('id','task'))

## -----------------------------------------------------------------------------
 icecream_est = ice_cal %>% vd_est_vdm(R=20000)

## -----------------------------------------------------------------------------
icecream_est %>% ec_estimates_MU 

## -----------------------------------------------------------------------------
icecream_est %>% ec_trace_MU()

## -----------------------------------------------------------------------------
icecream_est %>% ec_lmd_NR

## -----------------------------------------------------------------------------
#generate predictions
ho_demand=
    ice_ho %>%   
      prep_newprediction(ice_cal) %>% 
        vd_dem_vdm(icecream_est %>% vd_thin_draw) 

## -----------------------------------------------------------------------------
ho_demand

## -----------------------------------------------------------------------------
ho_demand %>%
  ec_dem_eval

## -----------------------------------------------------------------------------
my_scenario = tibble(
  id=1L,task=1L,alt=1:10, 
  Brand=  'HaagenDa',
  Flavor=c("Chocolate", "ChocChip", "ChocDough", "CookieCream", "Neapolitan", "Oreo", "RockyRoad", "Vanilla", "VanillaBean", "VanillaFudge"),
  Size=8,
  p=4)

## -----------------------------------------------------------------------------
my_scenario <- my_scenario %>% prep_newprediction(ice_cal)

## -----------------------------------------------------------------------------
N = n_distinct(icecream$id)

my_market=
  tibble(
      id = rep(seq_len(N),each=nrow(my_scenario)),
      task = 1,
      alt = rep(1:nrow(my_scenario),N)
    ) %>% 
  bind_cols(
      my_scenario[rep(1:nrow(my_scenario),N),-(1:3)]
    )

## -----------------------------------------------------------------------------
my_market_demand <-  my_market %>%
                     vd_dem_vdm(est = icecream_est)

## -----------------------------------------------------------------------------
my_market_demand %>% head

## -----------------------------------------------------------------------------
my_market_demand %>%
  ec_dem_aggregate(c('Brand','Flavor','Size'))

## -----------------------------------------------------------------------------
my_market_demand %>%
  ec_dem_aggregate(c('Brand','Flavor','Size')) %>%
  ec_dem_summarise()

## -----------------------------------------------------------------------------
my_market_demand %>%
  ec_dem_aggregate(c('Brand')) %>%
  ec_dem_summarise()

## ----echo=FALSE, message=FALSE, warning=FALSE---------------------------------
eps_not<-my_market %>% ec_gen_err_ev1(icecream_est, 99887766)

my_demcurve =
  my_market %>%
    ec_demcurve((my_market %>% transmute(focal=Flavor=='Oreo') %>% pull(focal)),
                c(.8,.9,1,1.1,1.2),
                vd_dem_vdm,
                icecream_est,
                eps_not) 

## -----------------------------------------------------------------------------
my_demcurve %>% do.call('rbind',.) %>%
  ggplot(aes(x=scenario,y=`E(demand)`,color=Flavor)) + geom_line()

