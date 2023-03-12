# echoice2


# Notes -------------------------------------------------------------------


# ec_ functions work for both discrete demand (dd) and volumetric demand (vd)
# vd_ functions are specific to volumetric demand models
# dd_ functions are specific to discrete choice models






# Namespace ---------------------------------------------------------------

#' @importFrom magrittr %>%
NULL

#' @importFrom dplyr group_by
NULL

#' @importFrom dplyr select
NULL

#' @importFrom dplyr select_if
NULL

#' @importFrom dplyr filter
NULL

#' @importFrom dplyr mutate
NULL

#' @importFrom dplyr mutate_if
NULL

#' @importFrom dplyr mutate_all
NULL

#' @importFrom dplyr transmute
NULL


#' @importFrom dplyr n_distinct
NULL

#' @importFrom dplyr left_join
NULL

#' @importFrom dplyr arrange
NULL

#' @importFrom dplyr summarise
NULL

#' @importFrom dplyr summarise_all
NULL

#' @importFrom dplyr summarise_if
NULL

#' @importFrom dplyr across
NULL


#' @importFrom dplyr pull
NULL

#' @importFrom dplyr one_of
NULL

#' @importFrom dplyr arrange
NULL

#' @importFrom dplyr bind_rows
NULL

#' @importFrom dplyr bind_cols
NULL

#' @importFrom dplyr relocate
NULL

#' @importFrom dplyr rename_all
NULL

#' @importFrom dplyr group_split
NULL

#' @importFrom dplyr group_keys
NULL

#' @importFrom dplyr ntile
NULL

#' @importFrom dplyr rename
NULL




#' @importFrom tidyr pivot_longer
NULL

#' @importFrom tidyr pivot_wider
NULL



#' @importFrom stringr str_subset
NULL

#' @importFrom stringr str_remove
NULL

#' @importFrom stringr str_extract
NULL

#' @importFrom forcats fct_relabel
NULL

#' @importFrom forcats fct_recode
NULL


#' @importFrom purrr map
NULL

#' @importFrom purrr map_dfr
NULL

#' @importFrom purrr map_df
NULL

#' @importFrom purrr imap_dfr
NULL

#' @importFrom purrr map2
NULL

#' @importFrom purrr map_dbl
NULL

#' @importFrom purrr reduce
NULL


#' @importFrom tidyselect any_of
NULL

#' @importFrom tidyselect all_of
NULL

#' @importFrom tidyselect contains
NULL

#' @importFrom tidyselect last_col
NULL

#' @importFrom tidyselect everything
NULL




#' @importFrom tibble tibble
NULL

#' @importFrom tibble as_tibble
NULL

#' @importFrom tibble add_row
NULL

#' @importFrom tibble add_column
NULL

#' @importFrom tibble rowid_to_column
NULL

#' @importFrom tibble enframe
NULL

#' @importFrom stats cov2cor
NULL

#' @importFrom stats model.matrix.lm
NULL

#' @importFrom stats quantile
NULL

#' @importFrom stats rnorm
NULL

#' @importFrom stats sd
NULL

#' @importFrom graphics par
NULL



#' @importFrom ggplot2 ggplot
NULL

#' @importFrom ggplot2 geom_boxplot
NULL

#' @importFrom ggplot2 geom_line
NULL

#' @importFrom ggplot2 coord_flip
NULL

#' @importFrom ggplot2 facet_wrap
NULL

#' @importFrom ggplot2 aes
NULL

#' @importFrom ggplot2 guides
NULL

#' @importFrom utils combn
NULL




utils::globalVariables(c(".", ".MAE", ".MSE", ".bias", ".demdraws", ".hp", ".hpall",
                         "lvl","id","reference", "reference_lvl",
                         "rowid", "x",
                         ".isfocal", ".prodid", ".s", "alt", "attr_level",
                         "attribute", "draw", "n", "task", "value", "xp", "part",
                         "p", "MAE", ".screendraws"))



# Utilities ---------------------------------------------------------------


#' Get the attribute of an object
#'
#' @param obj The object to get the attribute from.
#' @param attrname The name of the attribute to get.
#'
#' @return The attribute of the object.
#'
#' @examples
#' obj <- list(a = 1, b = 2)
#' attributes(obj)$test="hello"
#' `%.%`(obj, "test")
#'
#' @export
`%.%` <- function(obj,attrname) (attributes(obj))[[attrname]]


#' Log Marginal Density (Newton-Raftery)
#' 
#' This function uses the quick-and-dirty Newton-Raftery approximation for log-marginal-density.
#'
#' Approximation of LMD based on Newton-Raftery.
#' It is not the most accurate, but a very fast method.
#'
#' @param ll A vector of log-likelihood values (i.e., draws)
#' @return A single numeric value representing the log marginal density
#'
#'
#' @examples
#' logll_values <- c(-4000, -4001, -4002)
#' logMargDenNRu(logll_values)
#' @export
logMargDenNRu=function(ll) 
{
  med = stats::median(ll)
  return(med - log(mean(exp(-ll + med))))
}

#' Obtain Log Marginal Density from draw objects
#' 
#' This is a helper function to quickly obtain log marginal density from a draw object
#' 
#' Draws are split in 4 equal parts from start to finish, and LMD
#' is computed for each part. This helps to double-check convergence.
#'
#'
#' @param est 'echoice2' draw object
#' @return tibble with LMDs (first 25% of draws, next 25% of draws, ...) 
#' @examples
#' data(icecream)
#' #run MCMC sampler (use way more than 50 draws for actual use)
#' icecream_est <- icecream %>% dplyr::filter(id<100) %>% vd_est_vdm(R=50)
#' #obtain LMD by quartile of draws
#' ec_lmd_NR(icecream_est)
#'
#' @export
ec_lmd_NR=function(est){
  return(
    drop(est$loglike) %>% 
      enframe(name = 'draw') %>% 
      mutate(part=ntile(draw,4)) %>% 
      group_by(part) %>% summarise(lmd=logMargDenNRu(value), .groups='drop') %>% 
      mutate(part=part/4)
  )
}

# Data manipulation -------------------------------------------------------


#clean data
vd_janitor=function(vd, 
                    maxquant=999){
  `%!in%` <- Negate(`%in%`)
  
  #remove too high volumes
  fid_toomuch   = vd %>% dplyr::filter(x>maxquant) %>% dplyr::pull(id) %>% unique
  
  #remove all0s
  fid_toolittle = vd %>% group_by(id) %>% summarise(.s=sum(x)) %>% dplyr::filter(.s==0) %>% pull(id)
  
  #combine filter
  fid_all = base::intersect(fid_toolittle, fid_toomuch)
  
  filter_summary=
    tibble(
      `overall`=fid_all%>%n_distinct,
      `large quantities`=fid_toolittle%>%n_distinct,
      `all 0` = fid_toolittle %>% n_distinct())
  
  #filter and return
  vd = vd %>% dplyr::filter(id %!in% fid_all) 
  attributes(vd)$filter=filter_summary
  return(vd)
}


#' Dummy-code a categorical variable
#'
#'
#' @param data one column of categorical data to be dummy-coded
#' @return tibble with dummy variables
#' @examples
#' mytest=data.frame(attribute=factor(c('a','a','b','c','c')))
#' dummyvar(mytest)
#' @export
dummyvar<-function(data){
  #note: retains missing values as NA
  out = model.matrix.lm(~.-1,
                        data=data, 
                        na.action = "na.pass") 
  attributes(out)[c("assign", "contrasts")]=NULL
  return(as_tibble(out))
}


#' Create dummy variables within a tibble
#'
#'
#' @param dat A \code{tibble} with the data.
#' @param sel A character vector with the name(s) of the variables to be dummied.
#' @return tibble with dummy variables
#' @examples
#' mytest=data.frame(A=factor(c('a','a','b','c','c')), B=1:5)
#' dummify(mytest,"A")
#'
#' @export
dummify=function(dat, sel){
  dat=as_tibble(dat)
  for(i in seq_along(sel)){
    selv=sel[i]
    dummidata <- dat %>% select(all_of(selv)) %>% dummyvar()
    dat<- dat %>% add_column(dummidata, .before = selv) %>% 
      select(-all_of(selv)) %>% 
      rename_all(gsub, 
                 pattern = paste0("^(",selv,")"), 
                 replacement = paste0(selv,":"))
  }
  return(dat)
}


#' Obtain attributes and levels from tidy choice data with dummies
#'
#'
#' @param tdc A tibble with choice data
#' @return tibble
#' @examples
#' mytest=data.frame(A=factor(c('a','a','b','c','c')), B=1:5)
#' dummied_data = dummify(mytest,"A")
#' get_attr_lvl(dummied_data)
#'
#' @export
get_attr_lvl=function(tdc){
  tdc %>%
    select(-any_of(c('id','task','alt','p','x')))%>% 
    names %>% tibble::enframe() %>% 
      mutate(attribute=stringr::str_extract(.$value,"^.*?(?=\\:)")) %>%
      mutate(lvl=stringr::str_remove(.$value, .$attribute)) %>%
      mutate(lvl=stringr::str_remove(.$lvl,"^(:)")) %>%
      group_by(across("attribute")) %>%
      mutate(reference_lvl=dplyr::first(lvl)) %>%
      mutate(reference=ifelse(lvl==reference_lvl,1,0))%>%
      mutate(lvl_abbrv=abbreviate(lvl))%>%
      rename(attr_level=value)
}



#' Generate tidy choice data with dummies from long-format choice data
#'
#'
#' @param longdata tibble
#' @return tibble
#' @examples 
#' data(icecream)
#' vd_long_tidy(icecream)
#'
#' @export
vd_long_tidy<-function(longdata){

  #find categorical variables
  catvars <-
    longdata %>% select(tidyselect::where(is.factor)) %>% names
  
  #dummify categorical variables
  dummified <-
    longdata %>% 
      dummify(catvars) 
  
  #get list of attribute levels
  attrs <-
    dummified %>% get_attr_lvl
  
  #find reference categories in dummy coding
  refcats <-
    attrs %>% dplyr::filter(reference==1) %>% pull(attr_level)
  
  #generate output, then add attributes
  out <- 
    dummified %>% select(-any_of(refcats)) %>% add_column(int=1,.after='p')
  
  attributes(out)$ec_data = list(choice_type='volumetric',
                                 data_type='vd_tidy_choice',
                                 attributes=attrs)
  
  attributes(out)$Af = dummified %>% select_if(!(colnames(.)%in%c('id','task','alt','p','x')))
  
  return(out)
}







#' Prepare choice data for analysis
#'
#' This utility function prepares tidy choice data for fast MCMC samplers.
#'
#' Note: This function is only exported because it makes it easier to tinker with this package.
#' This function re-arranges choice data for fast access in highly-optimized MCMC samplers.
#' It Pre-computes task-wise total expenditures `sumpsx` and generates indices `xfr`,`xto`,`lfr`,`lto` for fast data access.
#'
#' @param dt tidy choice data (columns: id, task, alt, x, p, attributes)
#' @param Af (optional) contains a full design matrix (for attribute-based screening), or, more generally, a design matrix used for attribute-based screening
#'
#' @return list containing information for estimation functions
#'
#' @examples
#' #minimal data example
#' dt <- structure(list(id = c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 
#'                             2L, 2L), 
#'                      task = c(1L, 1L, 1L, 2L, 2L, 2L, 1L, 1L, 1L, 2L, 2L, 2L), 
#'                      alt = c(1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L), 
#'                      x = c(1, 0, 2, 1, 0, 1, 2, 3, 1, 1, 0, 1), 
#'                      p = c(0, 1, 1, 1, 2, 0, 2, 2, 1, 2, 1, 1), 
#'                      attr2 = c(1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0), 
#'                      attr1 = c(0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1)), 
#'                  class = c("tbl_df", "tbl", "data.frame"), row.names = c(NA,-12L))
#' #run prep function
#' test <- dt %>% vd_prepare 
#' 
#' @export
vd_prepare <- function(dt, Af=NULL){
  
  #arrange
  dt <- dt %>% 
    mutate(id=as.integer(id)) %>%
    arrange(as.numeric(id),task,alt) %>% mutate(xp=x*p)
  
  #sumpxs+nalts
  sumpxs_nalts <- dt %>% group_by(id,task) %>% 
    summarise(sumpxs=sum(xp),nalts=n(), .groups="drop")
  
  #ntasks
  ntasks<-dt %>% group_by(id) %>% 
    summarise(ntasks=n_distinct(task), .groups="drop") %>% pull(ntasks)
  
  #xfr-xto
  xfr_xto<-dt %>% group_by(id) %>% rowid_to_column() %>%
    summarise(xfr=min(rowid),xto=max(rowid), .groups="drop") 
  
  #lfr-lto
  lfr_lto<-sumpxs_nalts %>% group_by(id) %>%
    rowid_to_column() %>% 
    summarise(lfr=min(rowid),lto=max(rowid), .groups="drop")
  
  #tlens
  tlens<-dt %>% group_by(id) %>% summarise(tlens=n_distinct(task), .groups='drop') %>% pull(tlens)
  
  out=list(XX=dt$x,
           PP=dt$p, 
           AA=dt %>% select(-id,-task,-alt,-x,-p,-xp)%>% as.matrix(),
           nalts=sumpxs_nalts$nalts,
           sumpxs=sumpxs_nalts$sumpxs,
           ntasks=ntasks,
           xfr=xfr_xto$xfr,xto=xfr_xto$xto,
           lfr=lfr_lto$lfr,lto=lfr_lto$lto,
           tlens=tlens,
           idx=select(dt,id,task,alt))
  
  
  #meta-information
  ec_data = attributes(dt)$ec_data
  
  ec_data$data_type="vdm_choice"
  
  
  #full attribute data, screening lower limit
  if(!is.null(Af)){
    
    if(all(c("id","task","alt") %in% colnames(Af) )){
      out$AAf=Af%>%select(-id,-task,-alt)%>% as.matrix()
      foo <- dt %>% select(id,task,alt,x) %>% left_join(Af,by=c("id","task","alt")) %>% arrange(as.numeric(id),task,alt) 
      foo_check<-foo %>% select(-id,-task,-alt,-x) %>% is.na %>% sum
      
      if(foo_check==0){
        tauconst=1-(foo %>% dplyr::filter(x>0) %>%group_by(id) %>% summarise_all(max) %>% select(-id,-task,-alt,-x))
      }else{
        stop("Could not match full attribute tibble with choice data")
      }
      
      #bind_cols(Af)%>% filter(x>0) %>%group_by(id) %>% summarise_all(max) %>% select(-id,-task,-alt,-x)
    }else{
      message("Af does not contain id, task, alt columns. Assuming that attribute columns are properly sorted...")
      out$AAf=Af%>% as.matrix()
      tauconst=1-(dt %>% select(id,task,alt,x) %>%bind_cols(Af)%>% dplyr::filter(x>0) %>%group_by(id) %>% summarise_all(max) %>% select(-id,-task,-alt,-x))
    }
    
    out$tauconst=tauconst
    
    #meta-information
    ec_data$screening='attribute-based'
  }
  
  attributes(out)$ec_data = ec_data
  attributes(out)$Af      = attributes(dt)$Af
  
  return(out)
}



#' Prepare choice data for analysis (without x being present)
#'
#' This utility function prepares tidy choice data (without x) for fast data access.
#'
#' Note: This function is only exported because it makes it easier to tinker with this package.
#' This function re-arranges choice data for fast access, mainly for demand prediction.
#'
#' @param dt tidy choice data (columns: id, task, alt, p, attributes)
#' @param Af (optional) contains a full design matrix (for attribute-based screening), or, more generally, a design matrix used for attribute-based screening
#'
#' @return list containing information for prediction functions
#'
#' @examples
#' #Minimal example:
#' #One attribute with 3 levels, 2 subjects, 3 alternatives, 2 tasks
#' dt <- structure(list(id = c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 
#'                             2L, 2L), 
#'                      task = c(1L, 1L, 1L, 2L, 2L, 2L, 1L, 1L, 1L, 2L, 2L, 2L), 
#'                      alt = c(1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L), 
#'                      x = c(1, 0, 2, 1, 0, 1, 2, 3, 1, 1, 0, 1), 
#'                      p = c(0, 1, 1, 1, 2, 0, 2, 2, 1, 2, 1, 1), 
#'                      attr2 = c(1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0), 
#'                      attr1 = c(0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1)), 
#'                  class = c("tbl_df", "tbl", "data.frame"), row.names = c(NA,-12L))
#' test <- dt %>% dplyr::select(-all_of("x")) %>% vd_prepare_nox() 
#' 
#' @export
vd_prepare_nox <- function(dt, Af=NULL){
  
  `%!in%` = Negate(`%in%`)
  #arrange
  dt <- dt %>% 
    mutate(id=as.integer(id)) %>%
    arrange(as.numeric(id),task,alt)
  
  #ntasks
  ntasks<-dt %>% group_by(id) %>% 
    summarise(ntasks=n_distinct(task), .groups="drop") %>% pull(ntasks)
  
  #xfr-xto
  xfr_xto<-dt %>% group_by(id) %>% rowid_to_column() %>%
    summarise(xfr=min(rowid),xto=max(rowid), .groups="drop") 
  
  #sumpxs+nalts
  sump_nalts <- dt %>% group_by(id,task) %>% 
    summarise(sump=sum(p),nalts=n(), .groups="drop")
  
  #lfr-lto
  lfr_lto<-sump_nalts %>% group_by(id) %>%
    rowid_to_column() %>% 
    summarise(lfr=min(rowid),lto=max(rowid), .groups="drop")
  
  
  #tlens
  tlens<-dt %>% group_by(id) %>% summarise(tlens=n_distinct(task), .groups='drop') %>% pull(tlens)
  
  
  #select_if(colnames(.) %in% c("a","c","d","e"))
  
  out=list(PP=dt$p, 
           AA=dt %>% select_if((colnames(.) %!in% c('id','task','alt','p','x')))%>% as.matrix(),
           nalts=sump_nalts$nalts,
           ntasks=ntasks,
           xfr=xfr_xto$xfr,xto=xfr_xto$xto,
           lfr=lfr_lto$lfr,lto=lfr_lto$lto,
           tlens=tlens,
           idx=select(dt,id,task,alt))
  

  
  #meta-information
  ec_data = attributes(dt)$ec_data
  
  #full attribute data, screening lower limit
  if(!is.null(Af)){
    
    if(all(c("id","task","alt") %in% colnames(Af) )){
      out$AAf=Af%>%select(-id,-task,-alt)%>% as.matrix()
      foo <- dt %>% select(id,task,alt,x) %>% left_join(Af,by=c("id","task","alt")) %>% arrange(as.numeric(id),task,alt) 
      foo_check<-foo %>% select(-id,-task,-alt,-x) %>% is.na %>% sum
      
      if(foo_check==0){
        tauconst=1-(foo %>% dplyr::filter(x>0) %>%group_by(id) %>% summarise_all(max) %>% select(-id,-task,-alt,-x))
      }else{
        stop("Could not match full attribute tibble with choice data")
      }
      
      #bind_cols(Af)%>% filter(x>0) %>%group_by(id) %>% summarise_all(max) %>% select(-id,-task,-alt,-x)
    }else{
      message("Af does not contain id, task, alt columns. Assuming that attribute columns are properly sorted...")
      out$AAf=Af%>% as.matrix()
      tauconst=1-(dt %>% select(id,task,alt,x) %>%bind_cols(Af)%>% dplyr::filter(x>0) %>%group_by(id) %>% summarise_all(max) %>% select(-id,-task,-alt,-x))
    }
    
    out$tauconst=tauconst
    
    #meta-information
    ec_data$screening='attribute-based'
  }
  
  
  # if('attributes_levels' %in% names(attributes(dt))){
  #   out$attributes_levels=attributes(dt)$attributes_levels
  # }
  # ec_data$attributes=
  
  attributes(out)$ec_data = ec_data
  attributes(out)$Af =  attributes(dt)$Af
  
  return(out)
}








#utility function
#input checker
vd_check_long=function(dat){
  check1 <- all(c("id","task",'alt',"x","p") %in% colnames(dat))
  return(check1)
}


dd_check_long=function(dat){
  check1 <- all(c("id","task",'alt',"x","p") %in% colnames(dat))
  return(check1)
}




#' Summarize attributes and levels
#' 
#' Summarize attributes and levels in tidy choice data containing categorical attributes (before dummy-coding)
#' 
#' This functions looks for categorical attributes and summaries their levels
#' This is helpful when evaluating a new choice data file.
#'
#'
#' @param data_in A tibble, containing long-format choice data
#' @return A tibble with one row per attribute, and a list of the levels
#' @examples
#' data(icecream)
#' ec_summarize_attrlvls(icecream)
#'
#' @export
ec_summarize_attrlvls<-function(data_in){
  return(
    data_in %>% select(-any_of(c('id','task','alt','p','x'))) %>% 
      map(table) %>% 
      map(names) %>% 
      map(paste,collapse=', ') %>% 
      as_tibble() %>% 
      pivot_longer(everything()) %>% rlang::set_names(c('attribute','levels')) )
}
#' @rdname ec_summarize_attrlvls
#' @export
ec_summarise_attrlvls <- ec_summarize_attrlvls




# Working with estimates and draw objects ---------------------------------


#' Obtain upper level model estimates
#'
#'
#' @param est is an 'echoice2' draw object (list)
#' @param quantiles quantile for CI
#' @return tibble with MU (upper level) summaries
#' @examples
#' data(icecream)
#' #run MCMC sampler (use way more than 20 draws for actual use)
#' icecream_est <- icecream %>% dplyr::filter(id<20) %>% vd_est_vdm(R=20, cores=2)
#' #Upper-level summary
#' icecream_est %>% ec_estimates_MU
#' @export
ec_estimates_MU=function(est, quantiles=c(.05,.95)){
  
  quantiles_name=paste0("CI-",quantiles*100,"%")
  parnames=est$parnames
  
  estimates=  
    est$MUDraw %>% 
    as_tibble(.name_repair = ~make.names(seq_along(.), unique=TRUE)) %>% 
    rlang::set_names(parnames) %>%
    pivot_longer(everything(),names_to='par') %>%
    mutate(par=factor(par,levels=parnames)) %>%
    group_by(par) %>%
    summarise(mean=mean(value), 
              sd=sd(value), 
              !!(quantiles_name[1]):=quantile(value,probs=quantiles[1]),
              !!(quantiles_name[2]):=quantile(value,probs=quantiles[2]),
              sig=(prod(quantile(value,probs=quantiles))>0),
              .groups='drop')
  
  estimates$model=est$ec_type_short
  estimates$error=est$error_specification
  

  #add attribute groups
  chk<-attributes(est)$ec_data$attributes %>% pull(attribute) %>% n_distinct()
  if(chk>0){
    estimates<-
      estimates %>% 
      left_join(attributes(est)$ec_data$attributes %>%
                  select(attr_level,attribute,lvl,reference_lvl), 
                by=c('par'='attr_level')) %>%
      relocate(attribute,.before=par) %>%
      relocate(lvl,.before=par) %>% mutate(parameter=ifelse(!is.na(lvl),lvl,par))
  }
  #attr
  
  return(estimates)
}



#' Obtain posterior mean estimates of upper level correlations
#'
#'
#' @param est is an 'echoice2' draw object (list)
#' @return estimates of upper level correlations
#' @examples
#' data(icecream)
#' #run MCMC sampler (use way more than 20 draws for actual use)
#' icecream_est <- icecream %>% dplyr::filter(id<20) %>% vd_est_vdm(R=20, cores=2)
#' icecream_est %>% ec_estimates_SIGMA_corr %>% round(2)
#' @export
ec_estimates_SIGMA_corr=function(est){
  parnames=est$parnames
  rownames(est$SIGMADraw)=parnames
  colnames(est$SIGMADraw)=parnames
  return(est$SIGMADraw %>% apply(1:2,mean) %>% cov2cor)
}



#' Obtain posterior mean estimates of upper level covariance
#'
#'
#' @param est is an 'echoice2' draw object (list)
#' @return estimates of upper level covariance
#' @examples
#' data(icecream)
#' #run MCMC sampler (use way more than 20 draws for actual use)
#' icecream_est <- icecream %>% dplyr::filter(id<50) %>% vd_est_vdm(R=20, cores=2)
#' icecream_est %>% ec_estimates_SIGMA %>% round(2)
#' @export
ec_estimates_SIGMA=function(est){
  parnames=est$parnames
  rownames(est$SIGMADraw)=parnames
  colnames(est$SIGMADraw)=parnames
  return(est$SIGMADraw %>% apply(1:2,mean))
}


#' Summarize attribute-based screening parameters
#'
#' Summarize attribute-based screening parameters from an attribute-based screening model in 'echoice2'
#' 
#'
#' @param est is an 'echoice2' draw object (list) from a model with attribute-based screening
#' @param quantiles quantile for CI
#' @return tibble with screening summaries
#' @importFrom rlang :=
#' @examples
#' #run MCMC sampler (use way more than 20 draws for actual use)
#' data(icecream)
#' est_scr_icecream <- vd_est_vdm_screen(icecream%>%dplyr::filter(id<30), R=20, cores=2)
#' #summarise draws of screening probabilities
#' ec_estimates_screen(est_scr_icecream)
#' #Note: There is no variance in this illustrative example - more draws are needed
#' 
#' @export
ec_estimates_screen=function(est,quantiles=c(.05,.95)){
  
  quantiles_name=paste0("CI-",quantiles*100,"%")
  
  out<-
    est$deltaDraw %>% 
    as_tibble(.name_repair = ~make.names(seq_along(.), unique=TRUE)) %>%
    rlang::set_names(colnames(attributes(est)$Af)) %>% 
    pivot_longer(cols = everything(), names_to = 'par') %>%
    group_by(par) %>% 
    summarise(mean=mean(value), 
              sd=sd(value), 
              !!(quantiles_name[1]):=quantile(value,probs=quantiles[1],na.rm=TRUE),
              !!(quantiles_name[2]):=quantile(value,probs=quantiles[2],na.rm=TRUE),
              .groups='drop') 
  
  #add limits (maximum possible screening probability)
  if(!is.null(est$dat)){
    out<- out %>% 
      left_join(est$dat$tauconst %>% 
                  colMeans() %>% enframe(name = 'par', value = 'limit'),
                by = "par")
  }
  
  #add attribute groups
  chk<-attributes(est)$ec_data$attributes %>% pull(attribute) %>% n_distinct()
  if(chk>0){
    out<-
      out %>% 
      left_join(attributes(est)$ec_data$attributes %>%
                  select(attr_level,attribute,lvl), 
                by=c('par'='attr_level')) %>%
      relocate(attribute,.before=par) %>%
      relocate(lvl,.before=par)
  }
  
  return(out)
} 





# The Models --------------------------------------------------------------

# vd_ functions are specific to volumetric demand models
# dd_ functions are specific to discrete choice models
# _est_ are estimation functions, i.e. they run the corresponding MCMC sampler
# _LL_ obtain the log lilelihood
# _dem_ generates demand prediction

# Key model variants
# - vdm: standard volumetric demand model
# - vdm_screen: attribute-based screening
# - vdm_screenpr: including price
# - vdm_ss: set size variation



# Volumetric Demand Estimation --------------------------------------------



#' Estimate volumetric demand model
#'
#'
#' @param vd A tibble, containing volumetric demand data (long format)
#' @param tidy A logical, whether to apply 'echoice2' tidier function (default: TRUE)
#' @param R A numeric, no of draws
#' @param keep A numeric, thinning factor
#' @param cores An integer, no of CPU cores to use (default: auto-detect)
#' @param error_dist A string defining the error term distribution, 'EV1' or 'Normal'
#' @param control A list containing additional settings
#' 
#' @return An 'echoice2' draw object, in the form of a list
#' 
#' @seealso [vd_dem_vdm()] to generate demand predictions based on this model
#' @seealso [vd_est_vdm_screen()] to estimate a volumetric demand model with screening
#' 
#' @examples
#' data(icecream)
#' #run MCMC sampler (use way more than 10 draws for actual use)
#' icecream_est <- icecream %>% dplyr::filter(id<50) %>% vd_est_vdm(R=10, cores=2)
#' @export
vd_est_vdm=
  function(vd,
           tidy=TRUE,
           R=100000, 
           keep=10,
           cores=NULL,
           error_dist="EV1",
           control=list(include_data=TRUE)){
    
    #check input data
    if(!vd_check_long(vd)) stop("Check data")
    
    #error dist: either Normal or EV1
    if(!(error_dist=="Normal")){
      error_dist="EV1"
    }
    
    #integer variables
    vd<-vd %>% mutate(task=as.integer(.$task),
                      alt =as.integer(.$alt))
    
    #Multicore settings
    if(is.null(cores)){
      cores=parallel::detectCores(logical=FALSE)
    }
    message(paste0("Using ",cores," cores"))
    
    #re-arrange data
    if(tidy){
      dat <- 
        vd %>% 
        vd_long_tidy %>% vd_prepare
    }else{
      dat <- 
        vd %>% vd_prepare    
    }
    
    #Prior
    Bbar=matrix(rep(0,ncol(dat$AA)+3), ncol=ncol(dat$AA)+3)
    A=0.01*diag(1)
    nu=ncol(dat$AA)+9
    V=(ncol(dat$AA)+9)*diag(ncol(dat$AA)+3)
    Prior=list(Bbar=Bbar,A=A,nu=nu,V=V)
    
    #Run model
    if(error_dist=="EV1"){
    out=
      loop_vd2_RWMH( dat$XX, 
                     dat$PP,
                     dat$AA,
                     dat$nalts,
                     dat$sumpxs,  
                     dat$ntasks,  
                     dat$xfr-1,  
                     dat$xto-1,  
                     dat$lfr-1,  
                     dat$lto-1,
                     p=ncol(dat$AA)+3, 
                     N=length(dat$xfr),
                     R=R, 
                     keep=keep, 
                     Bbar=Bbar, 
                     A=A, 
                     nu=nu,  
                     V=V, 
                     tuneinterval = 30, steptunestart=.15, tunelength=10000, tunestart=500, 
                     progressinterval=100, cores=cores)
    }else{
      out=
        loop_vdn_RWMH( dat$XX, 
                       dat$PP,
                       dat$AA,
                       dat$nalts,
                       dat$sumpxs,  
                       dat$ntasks,  
                       dat$xfr-1,  
                       dat$xto-1,  
                       dat$lfr-1,  
                       dat$lto-1,
                       p=ncol(dat$AA)+3, 
                       N=length(dat$xfr),
                       R=R, 
                       keep=keep, 
                       Bbar=Bbar, 
                       A=A, 
                       nu=nu,  
                       V=V, 
                       tuneinterval = 30, steptunestart=.15, tunelength=10000, tunestart=500, 
                       progressinterval=100, cores=cores)
    }
    
    #Add data information
    out$A_names<-colnames(dat$AA)
    out$parnames=c(colnames(dat$AA),'sigma','gamma','E')
    out$model_parnames=c('sigma','gamma','E')
    
    #Add model information
    out$ec_type="volumetric-compensatory"
    out$error_specification=error_dist
    out$ec_type_short="VD-comp"
    out$Prior=Prior
    
    attributes(out)$Af<-attributes(dat)$Af
    attributes(out)$ec_data<-attributes(dat)$ec_data

    #Add training data
    if(control$include_data){
      out$dat=dat
    }
    
    return(out)
  }



# note on message suppression
# invisible(capture.output(icecream_est <- icecream %>% dplyr::filter(id<100) %>% vd_est_vdm(R=50, error_dist="Normal"), type = "message"))


#' Estimate volumetric demand model with attribute-based conjunctive screening
#'
#'
#'
#' @param vd volumetric demand data (long format)
#' @param R draws
#' @param keep thinning
#' @param cores no of CPU cores to use (default: auto-detect)
#' @param error_dist A string defining the error term distribution, 'EV1' or 'Normal' (default: 'EV1')
#' @param price_screen A logical, indicating whether price tag screening should be estimated (default: TRUE)
#' @param control list containing additional settings
#' 
#' 
#' @return est ec-draw object (List)
#' 
#' @examples
#' data(icecream)
#' icecream_est <- icecream %>% vd_est_vdm_screen(R=10, cores=2)
#' @export
vd_est_vdm_screen = function(vd,
                      R=100000, 
                      keep=10,
                      cores=NULL,
                      error_dist="EV1",
                      price_screen=TRUE,
                      control=list(include_data=TRUE)){
  
  #check input data
  if(!vd_check_long(vd)) stop("Check data")
  
  #integer variables
  vd<-vd %>% mutate(task=as.integer(.$task),
                    alt=as.integer(.$alt))
  
  #Multicore settings
  if(is.null(cores)){
    cores=parallel::detectCores(logical=FALSE)
  }
  message(paste0("Using ",cores," cores"))
  
  #re-arrange data
  dat <- 
    vd %>% 
    vd_long_tidy %>% vd_prepare
  
  #screening-relevant data
  dat$Af <- vd %>% vd_long_tidy %>%attributes() %>% `[[`('Af') %>% as.matrix()
  
  dat$tauconst=
    1-(
      vd %>% 
        mutate(x=sign(x)) %>%
        select(all_of(c("id","x")))%>%
        bind_cols(dat$Af%>%as_tibble(.name_repair = ~make.names(seq_along(.), unique=TRUE)))  %>%
        mutate(across(-any_of(c('id','x')),
                      ~.*x)) %>%
        select(-any_of("x")) %>%
        group_by(across("id")) %>%
        summarise(across(everything(), max)) %>%
        select(-any_of("id")) %>%
        as.matrix())
  
  
  #Prior
  Bbar=matrix(rep(0,ncol(dat$AA)+3), ncol=ncol(dat$AA)+3)
  A=0.01*diag(1)
  nu=ncol(dat$AA)+9
  V=(ncol(dat$AA)+9)*diag(ncol(dat$AA)+3)
  Prior=list(Bbar=Bbar,A=A,nu=nu,V=V)
  
  #Run model
  
  if(error_dist!="Normal"){
    if(!price_screen){
      out=
        loop_vdsr_e_RWMH( dat$XX, 
                          dat$PP,
                          dat$AA,
                          dat$Af,
                          t(dat$tauconst),
                          dat$nalts,
                          dat$sumpxs,  
                          dat$ntasks,  
                          dat$xfr-1,  
                          dat$xto-1,  
                          dat$lfr-1,  
                          dat$lto-1,
                          p=ncol(dat$AA)+3, N=length(dat$xfr),
                          R=R, keep=keep, 
                          Bbar=matrix(rep(0,ncol(dat$AA)+3),ncol=ncol(dat$AA)+3), A=0.01*diag(1), 
                          nu=ncol(dat$AA)+9,  V=(ncol(dat$AA)+9)*diag(ncol(dat$AA)+3), 
                          tuneinterval = 30, steptunestart=.15, tunelength=10000, tunestart=500, 
                          progressinterval=100, cores=cores)
    }else{
      out=
        loop_vdsrpr_e_RWMH( dat$XX, 
                            dat$PP,
                            dat$AA,
                            dat$Af,
                            t(dat$tauconst),
                            dat$nalts,
                            dat$sumpxs,  
                            dat$ntasks,  
                            dat$xfr-1,  
                            dat$xto-1,  
                            dat$lfr-1,  
                            dat$lto-1,
                            p=ncol(dat$AA)+3, N=length(dat$xfr),
                            R=R, keep=keep, 
                            Bbar=matrix(rep(0,ncol(dat$AA)+3),ncol=ncol(dat$AA)+3), A=0.01*diag(1), 
                            nu=ncol(dat$AA)+9,  V=(ncol(dat$AA)+9)*diag(ncol(dat$AA)+3), 
                            tuneinterval = 30, steptunestart=.15, tunelength=10000, tunestart=500, 
                            progressinterval=100, cores=cores)
    }
  } else{
  
    if(!price_screen){
    out=
      loop_vdsr_n_RWMH( dat$XX, 
                       dat$PP,
                       dat$AA,
                       dat$Af,
                       t(dat$tauconst),
                       dat$nalts,
                       dat$sumpxs,  
                       dat$ntasks,  
                       dat$xfr-1,  
                       dat$xto-1,  
                       dat$lfr-1,  
                       dat$lto-1,
                       p=ncol(dat$AA)+3, N=length(dat$xfr),
                       R=R, keep=keep, 
                       Bbar=matrix(rep(0,ncol(dat$AA)+3),ncol=ncol(dat$AA)+3), A=0.01*diag(1), 
                       nu=ncol(dat$AA)+9,  V=(ncol(dat$AA)+9)*diag(ncol(dat$AA)+3), 
                       tuneinterval = 30, steptunestart=.15, tunelength=10000, tunestart=500, 
                       progressinterval=100, cores=cores)
    }else{
    out=
      loop_vdsrpr_n_RWMH( dat$XX, 
                        dat$PP,
                        dat$AA,
                        dat$Af,
                        t(dat$tauconst),
                        dat$nalts,
                        dat$sumpxs,  
                        dat$ntasks,  
                        dat$xfr-1,  
                        dat$xto-1,  
                        dat$lfr-1,  
                        dat$lto-1,
                        p=ncol(dat$AA)+3, N=length(dat$xfr),
                        R=R, keep=keep, 
                        Bbar=matrix(rep(0,ncol(dat$AA)+3),ncol=ncol(dat$AA)+3), A=0.01*diag(1), 
                        nu=ncol(dat$AA)+9,  V=(ncol(dat$AA)+9)*diag(ncol(dat$AA)+3), 
                        tuneinterval = 30, steptunestart=.15, tunelength=10000, tunestart=500, 
                        progressinterval=100, cores=cores)
    }
  
  }

  
  #Add data information
  out$A_names<-colnames(dat$AA)
  out$parnames=c(colnames(dat$AA),'sigma','gamma','E')
  out$model_parnames=c('sigma','gamma','E')
  
  #Add model information
  if(!price_screen){
    out$ec_type="volumetric-conjunctive"
    out$ec_type_short="VD-conj"
  }else{
    out$ec_type="volumetric-conjunctive-pr"
    out$ec_type_short="VD-conj-pr"
  }
  if(error_dist!="Normal"){
    out$error_specification="EV1"
  }else{
    out$error_specification="Normal"
  }
  out$Prior=Prior
  
  attributes(out)$Af<-attributes(dat)$Af
  attributes(out)$ec_data<-attributes(dat)$ec_data

  #Add training data
  if(control$include_data){
    out$dat=dat
  }
  
  return(out)
}





#' Estimate volumetric demand model accounting for set size variation (1st order)
#' 
#' This model REQUIRES variation in choice-set size
#' 
#'
#' @param vd volumetric demand data (long format) with set size variation
#' @param order integer, either 1 or 2 (for now), indicating linear or quadratic set-size effect
#' @param R draws
#' @param keep thinning
#' @param cores no of CPU cores to use (default: auto-detect)
#' @param control list containing additional settings
#' @examples
#' data(icecream)
#' #note that for this example dataset, the model is not identified
#' #because the data lacks variation in set size
#' icecream_est <- icecream %>% vd_est_vdm_ss(R=10, cores=2)
#' @return est ec-draw object (List)
#' 
#' 
#' 
#' @export
vd_est_vdm_ss = function(vd,
                         order=1,
                         R=100000, 
                         keep=10,
                         cores=NULL,
                         control=list(include_data=TRUE)){
  
  #only supporting 1st and 2nd order
  if(!order==2){
    order=1
  }
  
  #check input data
  if(!vd_check_long(vd)) stop("Check data")
  
  #integer variables
  vd<-vd %>% mutate(task=as.integer(.$task),alt=as.integer(.$alt))
  
  #Multicore settings
  if(is.null(cores)){
    cores=parallel::detectCores(logical=FALSE)
  }
  message(paste0("Using ",cores," cores"))
  
  #re-arrange data
  dat <- 
    vd %>% 
    vd_long_tidy %>% vd_prepare
  
  
  #Prior
  Bbar=matrix(rep(0,ncol(dat$AA)+3+order), ncol=ncol(dat$AA)+3+order)
  Bbar[,ncol(dat$AA)+1]=-3
  
  A=0.01*diag(1)
  nu=ncol(dat$AA)+9
  V=(ncol(dat$AA)+9)*diag(ncol(dat$AA)+3+order)
  Prior=list(Bbar=Bbar,A=A,nu=nu,V=V)
  
  #Run model
  if(order==1){
  out=
    loop_vd_ss_RWMH(dat$XX, 
                   dat$PP,
                   dat$AA,
                   dat$nalts,
                   dat$sumpxs,  
                   dat$ntasks,  
                   dat$xfr-1,  
                   dat$xto-1,  
                   dat$lfr-1,  
                   dat$lto-1,
                   p=ncol(dat$AA)+4, 
                   N=length(dat$xfr),
                   R=R, 
                   keep=keep, 
                   Bbar=Bbar, 
                   A=A, 
                   nu=nu,  
                   V=V, 
                   tuneinterval = 30, steptunestart=.15, tunelength=10000, tunestart=500, 
                   progressinterval=100, cores=cores)
  }else{
    out=
      loop_vd_ssQ_RWMH(dat$XX, 
                       dat$PP,
                       dat$AA,
                       dat$nalts,
                       dat$sumpxs,  
                       dat$ntasks,  
                       dat$xfr-1,  
                       dat$xto-1,  
                       dat$lfr-1,  
                       dat$lto-1,
                       p=ncol(dat$AA)+5, 
                       N=length(dat$xfr),
                       R=R, 
                       keep=keep, 
                       Bbar=Bbar, 
                       A=A, 
                       nu=nu,  
                       V=V, 
                       tuneinterval = 30, steptunestart=.15, tunelength=10000, tunestart=500, 
                       progressinterval=100, cores=cores)    
  }
  
  #Add data information
  out$A_names<-colnames(dat$AA)
  if(order==1){
    out$parnames=c(colnames(dat$AA),'xi','sigma','gamma','E')
    out$model_parnames=c('xi','sigma','gamma','E')
  }
  
  if(order==2){
    out$parnames=c(colnames(dat$AA),'tau','xi','sigma','gamma','E')
    out$model_parnames=c('tau','xi','sigma','gamma','E')
  }
  

  #Add model information
  out$error_specification="EV1"
  if(order==1) out$ec_type="volumetric-compensatory-setsize_linear"
  if(order==2) out$ec_type="volumetric-compensatory-setsize_quadratic"
  if(order==1) out$ec_type_short="VD-comp-ssl"
  if(order==2) out$ec_type_short="VD-comp-ssq"
  
  out$Prior=Prior
  
  attributes(out)$Af<-attributes(dat)$Af
  attributes(out)$ec_data<-attributes(dat)$ec_data

  #Add training data
  if(control$include_data){
    out$dat=dat
  }
  
  return(out)
}




# Log-Likelihoods for entire datasets -------------------------------------


#' Log-Likelihood for compensatory volumetric demand model
#' 
#' 
#' @param draw A list, 'echoice2' draws object
#' @param vd A tibble, tidy choice data (before dummy-coding)
#' @param fromdraw An integer, from which draw onwards to compute LL (i.e., excl. burnin)
#' @examples
#' data(icecream)
#' #fit model
#' icecream_est <- icecream %>% vd_est_vdm(R=10, keep=1, cores=2)
#' #compute likelihood for each subject in each draw
#' loglls<-vd_LL_vdm(icecream_est, icecream, fromdraw = 2)
#' dim(loglls)
#' @return N x Draws Matrix of log-Likelihood values
#' @export
vd_LL_vdm <- function(draw, vd, fromdraw=1){
  
  #Number of draws total
  R=dim(draw$thetaDraw)[3]
  
  #starting point
  fromdraw=as.integer(fromdraw)
  
  #error dist
  error_specification = "EV1"
  if(!is.null(draw$error_specification)){
    error_specification=draw$error_specification
  }
  
  #prepare data
  dat <- 
    vd %>% 
    vd_long_tidy %>% vd_prepare
  
  #get LLs
  if(error_specification=="EV1"){
  out<- 
    vd2LLs(draw$thetaDraw[,,seq(fromdraw,R)],
             dat$XX,
             dat$PP,
             dat$AA, 
             dat$nalts, 
             dat$sumpxs, 
             dat$ntasks, 
             dat$xfr-1, dat$xto-1,
             dat$lfr-1, dat$lto-1, 
             ncol(draw$MUDraw), 
             ncol(draw$thetaDraw))
  }
  
  if(error_specification=="Normal"){
    out<- 
      vdnLLs(draw$thetaDraw[,,seq(fromdraw,R)],
             dat$XX,
             dat$PP,
             dat$AA, 
             dat$nalts, 
             dat$sumpxs, 
             dat$ntasks, 
             dat$xfr-1, dat$xto-1,
             dat$lfr-1, dat$lto-1, 
             ncol(draw$MUDraw), 
             ncol(draw$thetaDraw))
  }
  return(out) 
}


#' Log-Likelihood for conjunctive-screening volumetric demand model
#' 
#' 
#' @param draw A list, 'echoice2' draws object
#' @param vd A tibble, tidy choice data (before dummy-coding)
#' @param fromdraw An integer, from which draw onwards to compute LL (i.e., excl. burnin)
#' @return N x Draws Matrix of log-Likelihood values
#' @examples
#' data(icecream)
#' #fit model
#' icecream_est <- icecream %>% filter(id<20) %>% vd_est_vdm_screen(R=10, keep=1, cores=2)
#' #compute likelihood for each subject in each draw
#' loglls<-vd_LL_vdm_screen(icecream_est, icecream%>% filter(id<20), fromdraw = 2)
#' dim(loglls)
#' @export
vd_LL_vdm_screen <- function(draw, vd, fromdraw=1){
  
  #Number of draws total
  R=dim(draw$thetaDraw)[3]
  
  #starting point
  fromdraw=as.integer(fromdraw)
  
  #error dist
  error_specification = "EV1"
  if(!is.null(draw$error_specification)){
    error_specification=draw$error_specification
  }

  #screening model type "olumetric-conjunctive" or "volumetric-conjunctive-pr"
  screening_model_type <- draw$ec_type 

  #prepare data
  dat <- 
    vd %>% 
    vd_long_tidy %>% vd_prepare
  
  #screening-relevant data
  dat$Af <- vd %>% vd_long_tidy %>% 
    attributes() %>% `[[`('Af') %>% as.matrix()

  out<-NULL
  
  if(error_specification == "Normal"){
    if(screening_model_type=="volumetric-conjunctive"){
      out<- 
        vdsr2LLs(draw$thetaDraw[,,seq(fromdraw,R)],
                 draw$tauDraw[,,seq(fromdraw,R)],
                 dat$XX,
                 dat$PP,
                 dat$AA, 
                 attributes(dat)$Af%>%as.matrix(), 
                 dat$nalts, 
                 dat$sumpxs, 
                 dat$ntasks, 
                 dat$xfr-1, dat$xto-1,
                 dat$lfr-1, dat$lto-1, 
                 ncol(draw$MUDraw), 
                 ncol(draw$thetaDraw))
    }
    if(screening_model_type=="volumetric-conjunctive-pr"){
      out<- 
        vdsrprLLs(draw$thetaDraw[,,seq(fromdraw,R)],
                  draw$tauDraw[,,seq(fromdraw,R)],
                  draw$tau_pr_draw[,seq(fromdraw,R)], 
                  dat$XX,
                  dat$PP,
                  dat$AA, 
                  attributes(dat)$Af%>%as.matrix(), 
                  dat$nalts, 
                  dat$sumpxs, 
                  dat$ntasks, 
                  dat$xfr-1, dat$xto-1,
                  dat$lfr-1, dat$lto-1, 
                  ncol(draw$MUDraw), 
                  ncol(draw$thetaDraw))
    }
    
  }else{
    if(screening_model_type=="volumetric-conjunctive"){
      out<- 
        vdsreLLs(draw$thetaDraw[,,seq(fromdraw,R)],
                 draw$tauDraw[,,seq(fromdraw,R)],
                 dat$XX,
                 dat$PP,
                 dat$AA, 
                 attributes(dat)$Af%>%as.matrix(), 
                 dat$nalts, 
                 dat$sumpxs, 
                 dat$ntasks, 
                 dat$xfr-1, dat$xto-1,
                 dat$lfr-1, dat$lto-1, 
                 ncol(draw$MUDraw), 
                 ncol(draw$thetaDraw))
    }
    if(screening_model_type=="volumetric-conjunctive-pr"){
      out<-
        vdsrpreLLs(draw$thetaDraw[,,seq(fromdraw,R)],
                   draw$tauDraw[,,seq(fromdraw,R)],
                   draw$tau_pr_draw[,seq(fromdraw,R)], 
                   dat$XX,
                   dat$PP,
                   dat$AA, 
                   attributes(dat)$Af%>%as.matrix(), 
                   dat$nalts, 
                   dat$sumpxs, 
                   dat$ntasks, 
                   dat$xfr-1, dat$xto-1,
                   dat$lfr-1, dat$lto-1, 
                   ncol(draw$MUDraw), 
                   ncol(draw$thetaDraw))
    }

  }
  return(out) 
}



#' Log-Likelihood for volumetric demand model with set-size variation
#' 
#' 
#' @param draw A list, 'echoice2' draws object
#' @param vd A tibble, tidy choice data (before dummy-coding)
#' @param fromdraw An integer, from which draw onwards to compute LL (i.e., excl. burnin)
#' @return N x Draws Matrix of log-Likelihood values
#' @examples
#' data(icecream)
#' #fit model
#' #note: this is just for demo purposes
#' #on this demo dataset, the model is not identified
#' #due to a lack of set size variation
#' icecream_est <- icecream %>% vd_est_vdm_ss(R=10, keep=1, cores=2)
#' #compute likelihood for each subject in each draw
#' loglls<-vd_LL_vdmss(icecream_est, icecream, fromdraw = 2)
#' #300 respondents, 10 draws
#' dim(loglls)
#' @export
vd_LL_vdmss <- function(draw, vd, fromdraw=1){
  
  #Number of draws total
  R=dim(draw$thetaDraw)[3]
  
  #starting point
  fromdraw=as.integer(fromdraw)
  
  dat <- 
    vd %>% 
    vd_long_tidy %>% vd_prepare
  
  out<- 
    vdss_LLs(draw$thetaDraw[,,seq(fromdraw,R)],
           dat$XX,
           dat$PP,
           dat$AA, 
           dat$nalts, 
           dat$sumpxs, 
           dat$ntasks, 
           dat$xfr-1, dat$xto-1,
           dat$lfr-1, dat$lto-1, 
           ncol(draw$MUDraw), 
           ncol(draw$thetaDraw))
  
  return(out) 
}







# Volumetric Predictions --------------------------------------------------


#' Match factor levels between two datasets
#' 
#' Makes sure the factor levels in `data_new` are aligned with `data_old`
#' This is helpful for demand simulations.
#'
#'
#' @param data_new New long-format choice data 
#' @param data_old Old long-format choice data
#' @examples 
#' data(icecream)
#' prep_newprediction(icecream, icecream)
#' 
#' @return long-format choice data 
#' 
#' @export
prep_newprediction <- function(data_new,data_old){
  names_old <-
    data_old %>% 
    select(-any_of(c('id','task','alt','x','p','respnum'))) %>% names
  
  names_new <-
    data_new %>% select(-any_of(c('id','task','alt','x','p','respnum'))) %>% names
  
  check_allvars=all(names_old %in% names_new)
  if(!check_allvars) stop("mismatching product attributes")
  
  #relocate
  data_new <- 
    data_new %>% relocate(!!(names_old), .after='p') %>% select(-any_of('respnum'))
  
  #adjust factor levels
  for(k in seq_along(names_new)){
    data_new[[names_old[k]]] = 
      factor(  data_new[[names_old[k]]], 
               levels=levels(data_old[[names_old[k]]]))
  }
  
  return(data_new)
}




#' Demand Prediction (Volumetric Demand Model)
#'
#' Generating demand predictions for volumetric demand model. 
#' Reminder: there is no closed-form solution for demand, thus we need to integrate not only over the posterior distribution of parameters and the error distribution.
#' The function outputs a tibble containing id, task, alt, p, attributes, draws from the posterior of demand.
#' Error realizations can be pre-supplied to the `epsilon_not`. This helps create smooth demand curves or conduct optimization.
#'
#'
#' @param vd data
#' @param est ec-model draws 
#' @param epsilon_not (optional) error realizations
#' @param error_dist (optional) A string defining the error term distribution (default: 'EV1')
#' @param tidy (optional) apply 'echoice2' tidier (default: TRUE)
#' @param cores (optional) cores (default: auto-detect)
#' @return Draws of expected demand
#' @examples
#' data(icecream)
#' #run MCMC sampler (use way more than 10 draws for actual use)
#' icecream_est <- icecream %>% dplyr::filter(id<20) %>% vd_est_vdm(R=10, keep=1, cores=2)
#' #Generate demand predictions
#' icecream_predicted_demand=
#'  icecream %>% dplyr::filter(id<20) %>%   
#'    vd_dem_vdm(icecream_est, cores=2)
#' #column .demdraws contains draws from posterior of predicted demand
#' 
#' @seealso [prep_newprediction()] to match `vd`'s factor levels,
#'   [ec_gen_err_ev1()] for pre-generating error realizations and
#'   [vd_est_vdm()] for estimating the corresponding model
#' 
#' @export
vd_dem_vdm=function(vd,
                    est,
                    epsilon_not=NULL,
                    error_dist=NULL,
                    tidy=TRUE,
                    cores=NULL){
  
  #read error dist from draws if not supplied
  if(is.null(error_dist)){
    error_dist = est$error_specification
  }
  
  #error dist: either Normal or EV1
  if(!(error_dist=="Normal")){
    error_dist="EV1"
  }
  
  #cores  
  if(is.null(cores)){
    cores=parallel::detectCores(logical=FALSE)
  }
  message(paste0("Using ",cores," cores"))
  
  #re-arrange data
  if(tidy){
    dat <- 
      vd %>% 
      vd_long_tidy %>% vd_prepare_nox()
  }else{
    dat <- 
      vd  %>% vd_prepare_nox()
  }
  
  #demand sim
  if(is.null(epsilon_not)){
    
    if(error_dist=="Normal"){
    out=
      des_dem_vdmn(dat$PP,
                   dat$AA,
                   dat$nalts,
                   dat$ntasks,  
                   dat$xfr-1,
                   dat$xto-1,  
                   dat$lfr-1,
                   dat$lto-1,
                   dat$tlens,
                   est$thetaDraw,
                   cores=cores)      
    }else{
      #EV1
    out=
      des_dem_vdm( dat$PP,
                   dat$AA,
                   dat$nalts,
                   dat$ntasks,  
                   dat$xfr-1,
                   dat$xto-1,  
                   dat$lfr-1,
                   dat$lto-1,
                   dat$tlens,
                   est$thetaDraw,
                   cores=cores)
    }
    
    
  }else{
    #pre-supplied
    out=
      der_dem_vdm( dat$PP,
                   dat$AA,
                   dat$nalts,
                   dat$ntasks,  
                   dat$xfr-1,
                   dat$xto-1,  
                   dat$lfr-1,
                   dat$lto-1,
                   dat$tlens,
                   est$thetaDraw,
                   epsilon_not,
                   cores=cores)
  }
  
  attributes(out)=NULL
  #add draws to data tibble
  vd=as_tibble(vd)
  vd$.demdraws<-map(out,drop)  
  
  #add attributes
  attributes(vd)$attr_names <- vd %>% colnames %>% setdiff(c("id","task","alt","x","p" )) %>% str_subset('^\\.', negate = TRUE)
  # attributes(vd)$ec_model   <- attributes(est)$ec_model
  attributes(vd)$model <- list(ec_type=est$ec_type,
                               ec_type_short=est$ec_type_short,
                               error_specification=est$error_specification)
    
  return(vd)
}





#' Demand Prediction (Volumetric demand, attribute-based screening)
#'
#' Generating demand predictions for volumetric demand model with attribute-based screening. 
#' Reminder: there is no closed-form solution for demand, thus we need to integrate not only over the posterior distribution of parameters and the error distribution.
#' The function outputs a tibble containing id, task, alt, p, attributes, draws from the posterior of demand.
#' Eerror realisations can be pre-supplied to the `epsilon_not`. This helps create smooth demand curves or conduct optimization.
#'
#'
#' @param vd data
#' @param est ec-model draws 
#' @param epsilon_not (optional) error realizations
#' @param error_dist (optional) A string defining the error term distribution (default: 'EV1')
#' @param cores (optional) cores
#' @examples
#' data(icecream)
#' #run MCMC sampler (use way more than 20 draws for actual use)
#' icecream_est <- icecream %>% dplyr::filter(id<20) %>% vd_est_vdm_screen(R=20, keep=1, cores=2)
#' #Generate demand predictions
#' icecream_predicted_demand=
#'  icecream %>% dplyr::filter(id<20) %>%   
#'    vd_dem_vdm_screen(icecream_est, cores=2)
#' #column .demdraws contains draws from posterior of predicted demand
#' @return Draws of expected demand
#' 
#' 
#' @seealso [prep_newprediction()] to match `vd`'s factor levels,
#'   [ec_gen_err_normal()] for pre-generating error realizations and
#'   [vd_est_vdm_screen()] for estimating the corresponding model
#' 
#' @export
vd_dem_vdm_screen=function(vd,
                           est,
                           epsilon_not=NULL,
                           error_dist=NULL,
                           cores=NULL){
  
  #read error dist from draws if not supplied
  if(is.null(error_dist)){
    error_dist = est$error_specification
  }
  
  #error dist: either Normal or EV1
  if(!(error_dist=="Normal")){
    error_dist="EV1"
  }
  
  #screening model type "olumetric-conjunctive" or "volumetric-conjunctive-pr"
  screening_model_type <- est$ec_type
  
  #cores  
  if(is.null(cores)){
    cores=parallel::detectCores(logical=FALSE)
  }
  message(paste0("Using ",cores," cores"))
  
  
  #re-arrange data
  dat <- 
    vd %>% 
    vd_long_tidy %>% vd_prepare_nox()
  
  #screening-relevant data
  dat$Af <- vd %>% vd_long_tidy %>%attributes() %>% `[[`('Af') %>% as.matrix()
  
  #demand sim
  if(is.null(epsilon_not)){
    
    if(error_dist=="Normal"){
    if(screening_model_type=="volumetric-conjunctive"){
    out=
      des_dem_vdm_screen(dat$PP,
               dat$AA,
               dat$Af,
               dat$nalts,
               dat$tlens,  
               dat$ntasks,  
               dat$xfr-1,
               dat$xto-1,  
               dat$lfr-1,  
               dat$lto-1,
               est$thetaDraw,
               est$tauDraw, 
               cores=cores)
    }
    
    if(screening_model_type=="volumetric-conjunctive-pr"){
      out=
        des_dem_vdm_screenpr(dat$PP,
                             dat$AA,
                             dat$Af,
                             dat$nalts,
                             dat$tlens,  
                             dat$ntasks,  
                             dat$xfr-1,
                             dat$xto-1,  
                             dat$lfr-1,  
                             dat$lto-1,
                             est$thetaDraw,
                             est$tauDraw, 
                             est$tau_pr_draw,
                             cores=cores)
    }
    }
  
    if(error_dist=="EV1"){
      if(screening_model_type=="volumetric-conjunctive"){
        out=
          des_ev_dem_vdm_screen(dat$PP,
                                dat$AA,
                                dat$Af,
                                dat$nalts,
                                dat$tlens,  
                                dat$ntasks,  
                                dat$xfr-1,
                                dat$xto-1,  
                                dat$lfr-1,  
                                dat$lto-1,
                                est$thetaDraw,
                                est$tauDraw, 
                                cores=cores)
      }
      
      if(screening_model_type=="volumetric-conjunctive-pr"){
        out=
          des_ev_dem_vdm_screenpr(dat$PP,
                                  dat$AA,
                                  dat$Af,
                                  dat$nalts,
                                  dat$tlens,  
                                  dat$ntasks,  
                                  dat$xfr-1,
                                  dat$xto-1,  
                                  dat$lfr-1,  
                                  dat$lto-1,
                                  est$thetaDraw,
                                  est$tauDraw, 
                                  est$tau_pr_draw,
                                  cores=cores)
      }
    }  
    
    
  }else{
    
    if(screening_model_type=="volumetric-conjunctive"){
    out=
      der_dem_vdm_screen(dat$PP,
             dat$AA,
             dat$Af,
             dat$nalts,
             dat$tlens,  
             dat$ntasks,  
             dat$xfr-1,
             dat$xto-1,  
             dat$lfr-1,  
             dat$lto-1,
             est$thetaDraw,
             est$tauDraw,
             epsilon_not,
             cores=cores)
    }
    
    if(screening_model_type=="volumetric-conjunctive-pr"){
      out=
        der_dem_vdm_screenpr(dat$PP,
                             dat$AA,
                             dat$Af,
                             dat$nalts,
                             dat$tlens,  
                             dat$ntasks,  
                             dat$xfr-1,
                             dat$xto-1,  
                             dat$lfr-1,  
                             dat$lto-1,
                             est$thetaDraw,
                             est$tauDraw, 
                             est$tau_pr_draw,
                             epsilon_not,
                             cores=cores)
    }
    
    
  }
  
  attributes(out)=NULL
  #add draws to data tibble
  vd=as_tibble(vd)
  vd$.demdraws<-map(out,drop)   
  
  #add attributes
  attributes(vd)$attr_names <- 
    vd %>% colnames %>% setdiff(c("id","task","alt","x","p")) %>% 
      str_subset('^\\.', negate = TRUE)

  attributes(vd)$model <- list(ec_type=est$ec_type,
                               ec_type_short=est$ec_type_short,
                               error_specification=est$error_specification)
  return(vd)
}




#' Demand Prediction (Volumetric demand, accounting for set-size variation, EV1 errors)
#'
#' Generating demand predictions for volumetric demand model with set-size adjustment.
#' Reminder: there is no closed-form solution for demand, thus we need to integrate not only over the posterior distribution of parameters and the error distribution.
#' The function outputs a tibble containing id, task, alt, p, attributes, draws from the posterior of demand.
#' Eerror realizations can be pre-supplied to the `epsilon_not`. This helps create smooth demand curves or conduct optimization.
#'
#'
#' @param vd data
#' @param est ec-model draws 
#' @param epsilon_not (optional) error realizations
#' @param cores (optional) cores
#' @return Draws of expected demand
#' 
#' @examples
#' data(icecream)
#' #run MCMC sampler (use way more than 10 draws for actual use)
#' icecream_est <- icecream %>% dplyr::filter(id<10) %>% vd_est_vdm_ss(R=10, keep=1, cores=2)
#' #Generate demand predictions
#' icecream_predicted_demand=
#'  icecream %>% dplyr::filter(id<10) %>%   
#'    vd_dem_vdm_ss(icecream_est, cores=2)
#' #column .demdraws contains draws from posterior of predicted demand
#' 
#' @seealso [prep_newprediction()] to match `vd`'s factor levels,
#'   [ec_gen_err_ev1()] for pre-generating error realizations and
#'   [vd_est_vdm_ss()] for estimating the corresponding model
#' 
#' @export
vd_dem_vdm_ss=function(vd,
                       est,
                       epsilon_not=NULL,
                       cores=NULL){
  
  #model type
  model_type = est$ec_type 
  # "volumetric-compensatory-setsize_linear"
  # "volumetric-compensatory-setsize_quadratic"
  
  #cores  
  if(is.null(cores)){
    cores=parallel::detectCores(logical=FALSE)
  }
  message(paste0("Using ",cores," cores"))
  
  #re-arrange data
  dat <- 
    vd %>% 
    vd_long_tidy %>% vd_prepare_nox()
  
  #demand sim
  if(is.null(epsilon_not)){
    if(model_type=="volumetric-compensatory-setsize_linear"){
    out=
      des_dem_vdm_ss( dat$PP,
                   dat$AA,
                   dat$nalts,
                   dat$ntasks,  
                   dat$xfr-1,
                   dat$xto-1,  
                   dat$lfr-1,
                   dat$lto-1,
                   dat$tlens,
                   est$thetaDraw,
                   cores=cores)
    }
    
    if(model_type=="volumetric-compensatory-setsize_quadratic"){
      des_dem_vdm_ssq( dat$PP,
                       dat$AA,
                       dat$nalts,
                       dat$sumpxs,  
                       dat$ntasks,  
                       dat$xfr-1,
                       dat$xto-1,  
                       dat$lfr-1,
                       dat$lto-1,
                       dat$tlens,
                       est$thetaDraw,
                       cores=cores)
    }
    
  }else{
    if(model_type=="volumetric-compensatory-setsize_linear"){
    out=
      der_dem_vdm_ss( dat$PP,
                      dat$AA,
                      dat$nalts,
                      dat$ntasks,  
                      dat$xfr-1,
                      dat$xto-1,  
                      dat$lfr-1,
                      dat$lto-1,
                      dat$tlens,
                      est$thetaDraw,
                      epsilon_not,
                      cores=cores)
    }
    if(model_type=="volumetric-compensatory-setsize_quadratic"){
    out=
      der_dem_vdm_ssq(dat$PP,
                      dat$AA,
                      dat$nalts,
                      dat$ntasks,  
                      dat$xfr-1,
                      dat$xto-1,  
                      dat$lfr-1,
                      dat$lto-1,
                      dat$tlens,
                      est$thetaDraw,
                      epsilon_not,
                      cores=cores)
    }
  }
  

  attributes(out)=NULL
  #add draws to data tibble
  vd=as_tibble(vd)
  vd$.demdraws<-map(out,drop)  
  
  #add attributes
  attributes(vd)$attr_names <- vd %>% colnames %>% setdiff(c("id","task","alt","x","p" )) %>% str_subset('^\\.', negate = TRUE)
  attributes(vd)$model <- list(ec_type=est$ec_type,
                               ec_type_short=est$ec_type_short,
                               error_specification=est$error_specification)
  return(vd)
}









# Discrete Demand ---------------------------------------------------------


#' Estimate discrete choice model (HMNL)
#'
#'
#' @param dd discrete choice data (long format)
#' @param R draws
#' @param keep thinning
#' @param cores no of CPU cores to use (default: auto-detect)
#' @param control list containing additional settings
#' @examples 
#' data(icecream_discrete)
#' icecream_est <- icecream_discrete %>% dd_est_hmnl(R=20, cores=2)
#' 
#' @return est ec-draw object (List)
#' 
#' @seealso [dd_dem()] to generate demand predictions based on this model
#' 
#' @export
dd_est_hmnl = function(dd,
                       R=100000, 
                       keep=10,
                       cores=NULL,
                       control=list(include_data=TRUE)){
  
  #check input data
  if(!dd_check_long(dd)) stop("Check data")
  
  #integer variables
  dd<-dd %>% mutate(dplyr::across(all_of(c("task","alt")),as.integer))
  
  #Multicore settings
  if(is.null(cores)){
    cores=parallel::detectCores(logical=FALSE)
  }
  message(paste0("Using ",cores," cores"))
  
  #outside good present or not?
  outside_good_check <- 
    dd  %>% 
      mutate(x=sign(x)) %>%
      group_by(dplyr::across(all_of(c("id","task")))) %>% 
      summarise(dplyr::across(all_of(c("x")), sum)) %>% pull(all_of("x")) %>% mean
  
  #re-arrange data
  if(outside_good_check==1){
    #no outside
    dat <- 
      dd %>% 
      vd_long_tidy %>% select(-all_of("int")) %>% vd_prepare
  }else{
    #outside good
    dat <- 
      dd %>% 
      vd_long_tidy %>% vd_prepare
  }
  
  #Prior
  Bbar=matrix(rep(0,ncol(dat$AA)+1), ncol=ncol(dat$AA)+1)
  A=0.01*diag(1)
  nu=ncol(dat$AA)+9
  V=(ncol(dat$AA)+9)*diag(ncol(dat$AA)+1)
  Prior=list(Bbar=Bbar,A=A,nu=nu,V=V)
  
  #Run model
  out=
    loop_dd_RWMH( dat$XX, 
                   dat$PP,
                   dat$AA,
                   dat$nalts,
                   dat$ntasks,  
                   dat$xfr-1,  
                   dat$xto-1,  
                   dat$lfr-1,  
                   dat$lto-1,
                   p=ncol(dat$AA)+1, 
                   N=length(dat$xfr),
                   R=R, 
                   keep=keep, 
                   Bbar=Bbar, 
                   A=A, 
                   nu=nu,  
                   V=V, 
                   tuneinterval = 30, steptunestart=.15, tunelength=10000, tunestart=500, 
                   progressinterval=100, cores=cores)
  
  #Add data information
  out$A_names<-colnames(dat$AA)
  out$parnames<-c(colnames(dat$AA),'ln_beta_p')
  colnames(out$MUDraw)=c(colnames(dat$AA),'ln_beta_p')
  
  #Add model information
  out$ec_type="discrete-compensatory"
  out$ec_type_short="hmnl-comp"
  out$error_specification="EV1"
  out$model_parnames=c('beta_p')
  out$Prior=Prior
  
  attributes(out)$Af<-attributes(dat)$Af
  attributes(out)$ec_data<-attributes(dat)$ec_data
  
  #Add training data
  if(control$include_data){
    out$dat=dat
  }
  
  return(out)
}





#' Estimate discrete choice model (HMNL, attribute-based screening (not including price))
#'
#'
#' @param dd discrete choice data (long format)
#' @param price_screen A logical, indicating whether price tag screening should be estimated
#' @param R draws
#' @param keep thinning
#' @param cores no of CPU cores to use (default: auto-detect)
#' @param control list containing additional settings
#' @examples 
#' data(icecream_discrete)
#' icecream_est <- icecream_discrete %>% dplyr::filter(id<20) %>% 
#'   dd_est_hmnl_screen(R=20, cores=2)
#' @return est ec-draw object (List)
#' 
#' @seealso [dd_dem_sr()] to generate demand predictions based on this model
#' 
#' @export
dd_est_hmnl_screen = function(dd,
                       price_screen=TRUE,
                       R=100000,
                       keep=10,
                       cores=NULL,
                       control=list(include_data=TRUE)){
    
    #check input data
    if(!dd_check_long(dd)) stop("Check data")
    
    #integer variables, double variables
    dd<-dd %>% 
      mutate(dplyr::across(all_of(c("task","alt")) ,as.integer)) %>%
      mutate(dplyr::across(all_of(c("x")),as.double))
  
    #sorting
    dd<-dd %>% 
      mutate(dplyr::across(all_of(c("id")), as.numeric)) %>%
      arrange(dplyr::across("id"), dplyr::across("task"), dplyr::across("alt"))
    
    #Multicore settings
    if(is.null(cores)){
      cores=parallel::detectCores(logical=FALSE)
    }
    message(paste0("Using ",cores," cores"))
    
    #outside good present or not?
    outside_good_check <- 
      dd %>% 
      mutate(x=sign(x)) %>%
      group_by(dplyr::across(all_of(c("id","task")))) %>% 
      summarise(dplyr::across(all_of(c("x")), sum)) %>% pull(all_of("x")) %>% mean
    
    #re-arrange data
    if(outside_good_check==1){
      #no outside
      dat <- 
        dd %>% 
        vd_long_tidy %>% 
          #remove intercept
          select(-all_of("int")) %>% vd_prepare
    }else{
      #outside good
      dat <- 
        dd %>% 
        vd_long_tidy %>% vd_prepare
    }
    
    dat$Af <- dd %>% vd_long_tidy %>% 
      attributes() %>% `[[`('Af') %>% as.matrix()
  
    #boundaries for screning
    dat$tauconst=
      1-(
        dd %>% 
          mutate(x=sign(x)) %>%
          select(all_of(c("id","x")))%>%
          bind_cols(dat$Af%>%as_tibble(.name_repair = ~make.names(seq_along(.), unique=TRUE)))  %>%
          mutate(across(-any_of(c('id','x')),
                        ~.*x)) %>%
          select(-any_of("x")) %>%
          group_by(across("id")) %>%
          summarise(across(everything(), max)) %>%
          select(-any_of("id")) %>%
          as.matrix())
        
    
    #Prior
    Bbar=matrix(rep(0,ncol(dat$AA)+1), ncol=ncol(dat$AA)+1)
    A=0.01*diag(1)
    nu=ncol(dat$AA)+9
    V=(ncol(dat$AA)+9)*diag(ncol(dat$AA)+1)
    Prior=list(Bbar=Bbar,A=A,nu=nu,V=V)
    
    
    #Run model
    if(!price_screen){
    out=
      loop_ddrs_RWMH(dat$XX, 
                    dat$PP,
                    dat$AA,
                    dat$Af,
                    t(dat$tauconst),
                    dat$nalts,
                    dat$ntasks,  
                    dat$xfr-1,  
                    dat$xto-1,  
                    dat$lfr-1,  
                    dat$lto-1,
                    p=ncol(dat$AA)+1, 
                    N=length(dat$xfr),
                    R=R, 
                    keep=keep, 
                    Bbar=Bbar, 
                    A=A, 
                    nu=nu,  
                    V=V, 
                    tuneinterval = 30, steptunestart=.15, tunelength=10000, tunestart=500, 
                    progressinterval=100, cores=cores)
    }else{
      out=
        loop_ddrspr_RWMH(dat$XX, 
                         dat$PP,
                         dat$AA,
                         dat$Af,
                         t(dat$tauconst),
                         dat$nalts,
                         dat$ntasks,  
                         dat$xfr-1,  
                         dat$xto-1,  
                         dat$lfr-1,  
                         dat$lto-1,
                         p=ncol(dat$AA)+1, 
                         N=length(dat$xfr),
                         R=R, 
                         keep=keep, 
                         Bbar=Bbar, 
                         A=A, 
                         nu=nu,  
                         V=V, 
                         tuneinterval = 30, steptunestart=.15, tunelength=10000, tunestart=500, 
                         progressinterval=100, cores=cores) 
    }
    
  
    #Add data information
    out$A_names<-colnames(dat$AA)
    out$Af_names<-colnames(dat$AAf)
    out$parnames=c(colnames(dat$AA),'ln_beta_p')
    colnames(out$MUDraw)=c(colnames(dat$AA),'ln_beta_p')
    
    #Add model information
    if(!price_screen){
      out$ec_type="discrete-conjunctive"
      out$ec_type_short="hmnl-conj"
    }else{
      out$ec_type="discrete-conjunctive-pr"
      out$ec_type_short="hmnl-conj-pr"
    }
    
    out$error_specification="EV1"
    out$model_parnames=c('beta_p')
    out$Prior=Prior
    
    attributes(out)$Af<-attributes(dat)$Af
    attributes(out)$ec_data<-attributes(dat)$ec_data
    
    #Add training data
    if(control$include_data){
      out$dat=dat
    }
    
    return(out)
  }







# Log-Likelihoods for entire datasets -------------------------------------



#' Log-Likelihood for compensatory hmnl model
#' 
#' 
#' @param draw A list, 'echoice2' draws object
#' @param dd A tibble, tidy choice data (before dummy-coding)
#' @param fromdraw An integer, from which draw onwards to compute LL (i.e., excl. burnin)
#' @examples
#' data(icecream_discrete)
#' #fit model
#' icecream_est <- icecream_discrete %>% dd_est_hmnl(R=10, keep=1, cores=2)
#' #compute likelihood for each subject in each draw
#' loglls<-dd_LL(icecream_est, icecream_discrete, fromdraw = 2)
#' @return N x Draws Matrix of log-Likelihood values
#' @export
dd_LL <- function(draw, dd, fromdraw=1){
  
  R=dim(draw$thetaDraw)[3]
  
  dat <- 
    dd %>% 
      vd_long_tidy %>% 
        vd_prepare
  
  out<- 
    ddLLs(draw$thetaDraw[,,seq(fromdraw,R)],
           dat$XX,
           dat$PP,
           dat$AA, 
           dat$nalts, 
           dat$ntasks, 
           dat$xfr-1, dat$xto-1,
           dat$lfr-1, dat$lto-1, 
           ncol(draw$MUDraw), 
           ncol(draw$thetaDraw))
  
  return(out) 
}


#' Log-Likelihood for screening hmnl model
#' 
#' 
#' @param draw A list, 'echoice2' draws object
#' @param dd A tibble, tidy choice data (before dummy-coding)
#' @param fromdraw An integer, from which draw onwards to compute LL (i.e., excl. burnin)
#' @examples
#' data(icecream_discrete)
#' #fit model
#' icecream_est <- icecream_discrete %>% dd_est_hmnl_screen(R=10, keep=1, cores=2)
#' #compute likelihood for each subject in each draw
#' loglls<-dd_LL_sr(icecream_est, icecream_discrete, fromdraw = 2)
#' @return N x Draws Matrix of log-Likelihood values
#' @export
dd_LL_sr <- function(draw, dd, fromdraw=1){
  
  model_type = draw$ec_type
  #"discrete-conjunctive" "discrete-conjunctive-pr"
  
  R=dim(draw$thetaDraw)[3]
  
  out=NULL
  
  dat <- 
    dd %>% 
    vd_long_tidy %>% 
    vd_prepare
  
  if(model_type=="discrete-conjunctive"){
  out<- 
    ddsrLLs(draw$thetaDraw[,,seq(fromdraw,R)],
            draw$tauDraw[,,seq(fromdraw,R)],
            dat$XX,
            dat$PP,
            dat$AA, 
            attributes(dat)$Af%>%as.matrix(), 
            dat$nalts, 
            dat$ntasks, 
            dat$xfr-1, dat$xto-1,
            dat$lfr-1, dat$lto-1, 
            ncol(draw$MUDraw), 
            ncol(draw$thetaDraw))
  }
  
  
  if(model_type=="discrete-conjunctive-pr"){
    out<- 
      ddsrprLLs(draw$thetaDraw[,,seq(fromdraw,R)],
              draw$tauDraw[,,seq(fromdraw,R)],
              draw$tau_pr_draw[,seq(fromdraw,R)], 
              dat$XX,
              dat$PP,
              dat$AA, 
              attributes(dat)$Af%>%as.matrix(), 
              dat$nalts, 
              dat$ntasks, 
              dat$xfr-1, dat$xto-1,
              dat$lfr-1, dat$lto-1, 
              ncol(draw$MUDraw), 
              ncol(draw$thetaDraw))
  }
  return(out) 
}






# Discrete Demand Prediction ----------------------------------------------


#' Discrete Choice Predictions (HMNL)
#'
#'
#' @param dd tibble with long-format choice data
#' @param est estimation object
#' @param prob logical, report probabilities instead of demand
#' @param cores cores
#' @return Draws of expected choice
#' 
#' @examples 
#' data(icecream_discrete)
#' icecream_est <- icecream_discrete %>% filter(id<20) %>% dd_est_hmnl(R=10)
#' #demand prediction
#' icecream_dempred <- icecream_discrete %>% filter(id<20) %>% 
#'   dd_dem(icecream_est, cores=2)
#' 
#' @seealso [dd_est_hmnl()] to generate demand predictions based on this model
#' 
#' @export
dd_dem=function(dd,
                est,
                prob=FALSE,
                cores=NULL){
  
  #cores  
  if(is.null(cores)){
    cores=parallel::detectCores(logical=FALSE)
  }
  message(paste0("Using ",cores," cores"))
  
  
  #re-arrange data
  dat <- 
    dd %>% 
    vd_long_tidy %>% vd_prepare_nox()
  
  
  #demand sim
  if(prob){
  out=
    ddprob( dat$PP,
            dat$AA,
            dat$nalts,
            dat$tlens,
            dat$ntasks,  
            dat$xfr-1,
            dat$xto-1,  
            dat$lfr-1,
            dat$lto-1,
            est$thetaDraw,
            cores=cores)    
  }else{
  out=
    dddem( dat$PP,
           dat$AA,
           dat$nalts,
           dat$tlens,
           dat$ntasks,  
           dat$xfr-1,
           dat$xto-1,  
           dat$lfr-1,
           dat$lto-1,
           est$thetaDraw,
           cores=cores)
  }
  
  attributes(out)=NULL
  #add draws to data tibble
  dd=as_tibble(dd)
  dd$.demdraws<-map(out,drop)  
  
  #add attributes
  attributes(dd)$attr_names <- dd %>% colnames %>% setdiff(c("id","task","alt","x","p" )) %>% str_subset('^\\.', negate = TRUE)
  # attributes(dd)$ec_model   <- attributes(est)$ec_model
  attributes(dd)$model <- list(ec_type=est$ec_type,
                               ec_type_short=est$ec_type_short,
                               error_specification=est$error_specification)
  
  return(dd)
}






#' Discrete Choice Predictions (HMNL with attribute-based screening)
#'
#'
#' @param dd data
#' @param est est
#' @param prob logical, report probabilities instead of demand
#' @param cores cores
#' @return Draws of expected choice
#' 
#' @examples 
#' data(icecream_discrete)
#' icecream_est <- icecream_discrete %>% filter(id<20) %>% dd_est_hmnl_screen(R=10)
#' #demand prediction
#' icecream_dempred <- icecream_discrete %>% filter(id<20) %>% 
#'  dd_dem_sr(icecream_est, cores=2)
#' 
#' @seealso [dd_est_hmnl_screen()] to generate demand predictions based on this model
#' 
#' @export
dd_dem_sr=function(dd,
                   est,
                   prob=FALSE,
                   cores=NULL){
  
  #screening model type
  screening_model_type <- est$ec_type 
  
  #cores  
  if(is.null(cores)){
    cores=parallel::detectCores(logical=FALSE)
  }
  message(paste0("Using ",cores," cores"))
  
  
  #re-arrange data
  dat <- 
    dd %>% 
    vd_long_tidy %>% vd_prepare_nox()
  
  #screening-relevant data
  dat$Af <- dd %>% vd_long_tidy %>%attributes() %>% `[[`('Af') %>% as.matrix()
  
  if(screening_model_type=="discrete-conjunctive"){
    if(prob){
      out=
        ddsrprob(dat$PP,
                 dat$AA,
                 dat$Af,
                 dat$nalts,
                 dat$tlens,  
                 dat$ntasks,  
                 dat$xfr-1,
                 dat$xto-1,  
                 dat$lfr-1,  
                 dat$lto-1,
                 est$thetaDraw,
                 est$tauDraw, 
                 cores=cores)      
    }else{
      out=
        ddsrdem(dat$PP,
                dat$AA,
                dat$Af,
                dat$nalts,
                dat$tlens,  
                dat$ntasks,  
                dat$xfr-1,
                dat$xto-1,  
                dat$lfr-1,  
                dat$lto-1,
                est$thetaDraw,
                est$tauDraw, 
                cores=cores)
    }
    
  }
  
  
  if(screening_model_type=="discrete-conjunctive-pr"){
    if(prob){
      out=
        ddsrprprob(dat$PP,
                   dat$AA,
                   dat$Af,
                   dat$nalts,
                   dat$tlens,  
                   dat$ntasks,  
                   dat$xfr-1,
                   dat$xto-1,  
                   dat$lfr-1,  
                   dat$lto-1,
                   est$thetaDraw,
                   est$tauDraw, 
                   est$tau_pr_draw,
                   cores=cores)    
    }else{
      out=
        ddsrprdem(dat$PP,
                  dat$AA,
                  dat$Af,
                  dat$nalts,
                  dat$tlens,  
                  dat$ntasks,  
                  dat$xfr-1,
                  dat$xto-1,  
                  dat$lfr-1,  
                  dat$lto-1,
                  est$thetaDraw,
                  est$tauDraw, 
                  est$tau_pr_draw,
                  cores=cores)
    }
  }
  
  
  attributes(out)=NULL
  #add draws to data tibble
  dd=as_tibble(dd)
  dd$.demdraws<-map(out,drop)  
  
  #add attributes
  attributes(dd)$attr_names <- dd %>% colnames %>% setdiff(c("id","task","alt","x","p" )) %>% str_subset('^\\.', negate = TRUE)
  # attributes(dd)$ec_model   <- attributes(est)$ec_model
  attributes(dd)$model <- list(ec_type=est$ec_type,
                               ec_type_short=est$ec_type_short,
                               error_specification=est$error_specification)
  return(dd)
}




# working with demand -----------------------------------------------------






#' Add product id to demand draws
#' 
#' This adds a unique product identifier to  demand draw objects.
#'
#'
#' @param de demand draws
#' @return est
#' 
#' @examples
#' data(icecream)
#' #run MCMC sampler (use way more than 10 draws for actual use)
#' icecream_est <- icecream %>% dplyr::filter(id<50) %>% vd_est_vdm(R=10, keep=1)
#' #Generate demand predictions
#' icecream_predicted_demand=
#'  icecream %>% dplyr::filter(id<50) %>%   
#'    vd_dem_vdm(icecream_est)
#' #add prodid
#' icecream_predicted_demand_w_id<-icecream_predicted_demand %>% vd_add_prodid
#' 
#' @export
vd_add_prodid<-function(de){
  
  prodids<-
    de %>% 
    select(one_of(attributes(de)$attr_names)) %>%
    unique %>% 
    mutate(.prodid=1:n())
  
  attributes(de)$prodids=prodids
  
  de %>% 
    left_join(prodids, by=attributes(de)$attr_names) %>%
    relocate(.prodid,.before=attributes(de)$attr_names[1]) %>%
    return()
}


##utility function
ec_named_group_split <- function(.tbl, ...) {
  grouped <- group_by(.tbl, ...)
  names <- rlang::eval_bare(rlang::expr(paste(!!!group_keys(grouped), sep = ":/:")))
  
  grouped %>% 
    group_split() %>% 
    rlang::set_names(names)
}



#' Aggregate posterior draws of demand
#'
#' Aggregate demand draws, e.g. from individual-choice occasion-alternative level to individual level.
#' (using the new demand draw format)
#'
#' @usage ec_dem_aggregate(de,groupby)
#'
#' @param de demand draws
#' @param groupby groupby grouping variables (as (vector of) string(s))
#' @return Aggregated demand predictions
#' @examples
#' data(icecream)
#' #run MCMC sampler (use way more than 50 draws for actual use)
#' icecream_est <- icecream %>% dplyr::filter(id<20) %>% vd_est_vdm(R=10, keep=1)
#' #Generate demand predictions
#' icecream_predicted_demand=
#'  icecream %>% dplyr::filter(id<20) %>%   
#'    vd_dem_vdm(icecream_est)
#' #aggregate
#' brand_lvl_pred_demand <-
#'  icecream_predicted_demand %>% ec_dem_aggregate("Brand")
#' 
#' @export
ec_dem_aggregate = function(de, groupby){
  return(
  de %>% 
    group_by(!!!rlang::syms(groupby)) %>% 
    summarise(.demdraws=list(purrr::reduce(.demdraws,`+`))) )
}



#' Summarize posterior draws of demand
#'
#' Adds summaries of posterior draws of demand to tibble.
#' (using the new demand draw format)
#'
#' @usage ec_dem_summarise(de,quantiles)
#'
#' @param de demand draws
#' @param quantiles Quantiles for Credibility Intervals (default: 90% interval)
#' @return Summary of demand predictions
#' @examples
#' \donttest{
#' data(icecream)
#' #run MCMC sampler (use way more than 10 draws for actual use)
#' icecream_est <- icecream %>% dplyr::filter(id<10) %>% vd_est_vdm(R=10, keep=1)
#' #Generate demand predictions
#' icecream_predicted_demand=
#'  icecream %>% dplyr::filter(id<10) %>%   
#'    vd_dem_vdm(icecream_est)
#' #aggregate
#' brand_lvl_pred_demand <-
#'  icecream_predicted_demand %>% ec_dem_aggregate("Brand")
#' #summarise
#' brand_lvl_pred_demand %>% ec_dem_summarise()
#' }
#' 
#' @importFrom rlang :=
#' @export
ec_dem_summarise = function(de, quantiles=c(.05,.95)){
  quantiles_name=paste0("CI-", quantiles*100, "%")
  
  return(
  de %>% 
    mutate(
      `E(demand)`=map_dbl(.demdraws, mean),
      `S(demand)`=map_dbl(.demdraws, sd),
      !!(quantiles_name[1]):=map_dbl(.demdraws, quantile, probs=quantiles[1]),
      !!(quantiles_name[2]):=map_dbl(.demdraws, quantile, probs=quantiles[2])
    ) )
}
#' @rdname ec_dem_summarise
#' @export
ec_dem_summarize <- ec_dem_summarise




#' Summarize posterior draws of screening
#'
#' Adds summaries of posterior draws of demand to tibble.
#' (using the new demand draw format)
#'
#'
#' @param sc tibble containing screening draws in .screendraws  
#' @param quantiles Quantiles for Credibility Intervals (default: 90% interval)
#' @return Summary of screening probabilities
#' @examples
#' \donttest{
#' data(icecream)
#' icecream_est <- icecream %>% vd_est_vdm_screen(R=20,  price_screen=TRUE)
#' #consideration set by respondent
#' cons_ss <- 
#' ec_screenprob_sr(icecream, icecream_est) %>%
#' group_by(id, task)  %>%
#'   summarise(.screendraws=list(purrr::reduce(.screendraws ,`+`))) %>%
#'   ec_screen_summarise() %>%
#'   group_by(id) %>%
#'   summarise(n_screen=mean(`E(screening)`))
#'   }
#' @importFrom rlang :=
#' @export
ec_screen_summarise = function(sc, quantiles=c(.05,.95)){
  quantiles_name=paste0("screening-CI-", quantiles*100, "%")
  
  return(
  sc %>% 
    mutate(
      `E(screening)`=map_dbl(.screendraws, mean),
      `S(screening)`=map_dbl(.screendraws, sd),
      !!(quantiles_name[1]):=map_dbl(.screendraws, quantile, probs=quantiles[1]),
      !!(quantiles_name[2]):=map_dbl(.screendraws, quantile, probs=quantiles[2])
    ) )
}
#' @rdname ec_screen_summarise
#' @export
ec_screen_summarize <- ec_screen_summarise



#' Summarize posterior draws of demand (volumetric models only)
#'
#' Adds summaries of posterior draws of demand to tibble.
#' (using the new demand draw format)
#'
#'
#' @param de demand draws
#' @param quantiles Quantiles for Credibility Intervals (default: 90% interval)
#' @return Summary of demand predictions
#' @examples
#' \donttest{
#' data(icecream)
#' #run MCMC sampler (use way more than 10 draws for actual use)
#' icecream_est <- icecream %>% dplyr::filter(id<10) %>% vd_est_vdm(R=10, keep=1)
#' #Generate demand predictions
#' icecream_predicted_demand=
#'  icecream %>% dplyr::filter(id<10) %>%   
#'    vd_dem_vdm(icecream_est)
#' #aggregate
#' brand_lvl_pred_demand <-
#'  icecream_predicted_demand %>% ec_dem_aggregate("Brand")
#' #summarise
#' brand_lvl_pred_demand %>% vd_dem_summarise()
#' }
#' @importFrom rlang :=
#' @export
vd_dem_summarise = function(de, quantiles=c(.05,.95)){
  quantiles_name=paste0("CI-", quantiles*100, "%")
  
  return(
  de %>% 
    mutate(
      `E(demand)`=map_dbl(.demdraws, mean),
      `S(demand)`=map_dbl(.demdraws, sd),
      !!(quantiles_name[1]):=map_dbl(.demdraws, quantile, probs=quantiles[1]),
      !!(quantiles_name[2]):=map_dbl(.demdraws, quantile, probs=quantiles[2]),
      `E(interior)`= map_dbl(map(.demdraws, sign),mean),
      `S(interior)`= map_dbl(map(.demdraws, sign),sd)
    ) )
}
#' @rdname vd_dem_summarise
#' @export
vd_dem_summarize <- vd_dem_summarise






#' Evaluate (hold-out) demand predictions
#'
#' This function obtains proper posterior fit statistics. 
#' It computes the difference between true demand and each draw from the demand posterior. 
#' Then, fit statistics are obtained.
#'
#'
#' @param de demand draws (output from vd_dem_x function)
#' @return Predictive fit statistics (MAE, MSE, RAE, bias, hit-probability)
#' 
#' data(icecream)
#' #run MCMC sampler (use way more than 50 draws for actual use)
#' icecream_est <- icecream %>% dplyr::filter(id<100) %>% vd_est_vdm(R=100, keep=1)
#' #Generate demand predictions
#' icecream_predicted_demand=
#'  icecream %>% dplyr::filter(id<100) %>%   
#'    vd_dem_vdm(icecream_est)
#' #evaluate in-sample fit (note: too few draws for good results)
#' ec_dem_eval(icecream_predicted_demand)
#' 
#' @export
ec_dem_eval = function(de){
  `%!in%` = Negate(`%in%`)
  
  is_this_discrete=attributes(de)$model$ec_type %>% stringr::str_detect('discrete')
  
  if(is_this_discrete){
    out <- de %>%
      mutate(.MSE=map_dbl(purrr::map2(.demdraws,x,function(draws,x)(draws-x)^2 ),mean),
             .MAE=map_dbl(purrr::map2(.demdraws,x,function(draws,x)abs(draws-x)),mean),
             .bias=map_dbl(purrr::map2(.demdraws,x,function(draws,x)(draws-x)),mean),
             .R=abs(x-mean(x)),
             .hp=map_dbl(purrr::map2(.demdraws,x,function(draws,x)(draws==x) ),mean)) %>%
      summarise(MSE=mean(.$.MSE),
                MAE=mean(.$.MAE),
                bias=mean(.$.bias),
                RAE=MAE/mean(.$.R),
                hitrate=mean(.$.hp))
  }else{
    out <- de %>%
      mutate(.MSE=purrr::map_dbl(purrr::map2(.demdraws,x,function(draws,x)(draws-x)^2 ),mean),
             .MAE=purrr::map_dbl(purrr::map2(.demdraws,x,function(draws,x)abs(draws-x)),mean),
             .bias=purrr::map_dbl(purrr::map2(.demdraws,x,function(draws,x)(draws-x)),mean),
             .R=abs(x-mean(x))) %>%
      summarise(MSE=mean(.$.MSE),
                MAE=mean(.$.MAE),
                bias=mean(.$.bias),
                RAE=MAE/mean(.$.R))
  }
  
  return(out)
}





#' Create demand curves
#' 
#' This helper function creates demand curves
#'
#'
#' @param ec_long choice scenario (discrete or volumetric)
#' @param focal_product Logical vector picking the focal product for which to create a demand curve
#' @param rel_pricerange Price range, relative to base case price; this is used to create demand curve
#' @param dem_fun demand function (e.g., `dd_prob` for HMNL or `vd_dem_vdm` for volumetric demand). For discrete choice, use choice probabilities instead of choice predictions.
#' @param draws ec-draws object (e.g., output from `dd_est_hmnl` or `vd_est_vd`)
#' @param epsilon_not (optional) error realisatins (this helps make curves look smother for voumetric models)
#' @return List containing aggregate demand quantities for each scenario defined by `rel_pricerange`
#' 
#' @examples
#' \donttest{
#' data(icecream)
#' #run MCMC sampler (use way more than 50 draws for actual use)
#' icecream_est <- icecream %>% dplyr::filter(id<100) %>% 
#' vd_est_vdm(R=20, keep=1)
#' #demand at different price points
#' dem_scenarios<-
#' ec_demcurve(icecream%>% dplyr::filter(id<100),
#'  icecream%>% dplyr::filter(id<100) %>% pull('Brand')=="Store",
#'  c(.75,1,1.25),vd_dem_vdm,icecream_est)
#' #optional plot
#' # dem_scenarios %>% 
#' #   do.call('rbind',.) %>%
#' #   ggplot(aes(x=scenario,y=`E(demand)`,color=Flavor)) + geom_line()
#' }
#' 
#' @seealso [ec_gen_err_normal()] to generate error realization from Normal distribution,
#' [ec_gen_err_ev1()] to generate error realization from EV1 distribution
#' 
#' @export
ec_demcurve=function(ec_long,
                     focal_product,
                     rel_pricerange,
                     dem_fun,
                     draws,
                     epsilon_not=NULL){
  
  #select correct demand function
  if(is.function(dem_fun)){
    demfun={dem_fun}
  }else{
    demfun=eval(parse(text=dem_fun))
  }
  
  #get key attributes
  attr_names=ec_long %>% 
    select_if(!(colnames(.)%in%c('id','task','alt','p'))) %>% colnames
  
  #output
  out=list()
  
  if(is.null(epsilon_not)){
    for(kk in seq_along(rel_pricerange)){
      testmarket_temp=ec_long
      testmarket_temp$p[focal_product]=testmarket_temp$p[focal_product]*rel_pricerange[kk]
      
      out[[kk]]=
        testmarket_temp %>% 
        demfun(est=draws) %>% 
        ec_dem_aggregate(attr_names) %>% ec_dem_summarise %>% select(-.demdraws) %>%
        bind_cols( scenario=rel_pricerange[kk])
      
    }
  }else{
    for(kk in seq_along(rel_pricerange)){
      testmarket_temp=ec_long
      testmarket_temp$p[focal_product]=testmarket_temp$p[focal_product]*rel_pricerange[kk]
      
      out[[kk]]=
        testmarket_temp %>% 
        demfun(est=draws) %>% 
        ec_dem_aggregate(attr_names )%>% ec_dem_summarise  %>% select(-.demdraws) %>%
        bind_cols( scenario=rel_pricerange[kk])
    }  
  }
  
  return(out)
}



#' Create demand-incidence curves
#' 
#' This helper function creates demand curves
#'
#'
#' @param ec_long choice scenario (discrete or volumetric)
#' @param focal_product Logical vector picking the focal product for which to create a demand curve
#' @param rel_pricerange Price range, relative to base case price; this is used to create demand curve
#' @param dem_fun demand function (e.g., `dd_prob` for HMNL or `vd_dem_vdm` for volumetric demand). For discrete choice, use choice probabilities instead of choice predictions.
#' @param draws ec-draws object (e.g., output from `dd_est_hmnl` or `vd_est_vd`)
#' @param epsilon_not (optional) error realisatins (this helps make curves look smother for voumetric models)
#' @return List containing aggregate demand quantities for each scenario defined by `rel_pricerange`
#' @examples
#' \donttest{
#' data(icecream)
#' #run MCMC sampler (use way more than 50 draws for actual use)
#' icecream_est <- icecream %>% dplyr::filter(id<50) %>% 
#' vd_est_vdm(R=20, keep=1)
#' #demand at different price points
#' inci_scenarios<-
#' ec_demcurve_inci(icecream%>% dplyr::filter(id<50),
#'  icecream%>% dplyr::filter(id<50) %>% pull('Brand')=="Store",
#'  c(.75,1,1.25),vd_dem_vdm,icecream_est)
#' }
#' @seealso [ec_gen_err_normal()] to generate error realization from Normal distribution,
#' [ec_gen_err_ev1()] to generate error realization from EV1 distribution
#' 
#' @export
ec_demcurve_inci=function(ec_long,
                          focal_product,
                          rel_pricerange,
                          dem_fun,
                          draws,
                          epsilon_not=NULL){
  
  #select correct demand function
  if(is.function(dem_fun)){
    demfun={dem_fun}
  }else{
    demfun=eval(parse(text=dem_fun))
  }
  
  #get key attributes
  attr_names=ec_long %>% 
    select_if(!(colnames(.)%in%c('id','task','alt','p'))) %>% colnames
  
  #output
  out=list()
  
  if(is.null(epsilon_not)){
    for(kk in seq_along(rel_pricerange)){
      testmarket_temp=ec_long
      testmarket_temp$p[focal_product]=testmarket_temp$p[focal_product]*rel_pricerange[kk]
      
      out[[kk]]=
        testmarket_temp %>% 
        demfun(est=draws) %>% mutate(.demdraws=map(.demdraws, sign))%>%
        ec_dem_aggregate(attr_names) %>% ec_dem_summarise %>% select(-.demdraws) %>%
        bind_cols( scenario=rel_pricerange[kk])
      
    }
  }else{
    for(kk in seq_along(rel_pricerange)){
      testmarket_temp=ec_long
      testmarket_temp$p[focal_product]=testmarket_temp$p[focal_product]*rel_pricerange[kk]
      
      out[[kk]]=
        testmarket_temp %>% 
        demfun(est=draws) %>% mutate(.demdraws=map(.demdraws, sign))%>%
        ec_dem_aggregate(attr_names )%>% ec_dem_summarise  %>% select(-.demdraws) %>%
        bind_cols( scenario=rel_pricerange[kk])
    }  
  }
  
  return(out)
}



#' Create demand-incidence curves
#' 
#' This helper function creates demand curves
#'
#'
#' @param ec_long choice scenario (discrete or volumetric)
#' @param focal_product Logical vector picking the focal product for which to create a demand curve
#' @param rel_pricerange Price range, relative to base case price; this is used to create demand curve
#' @param dem_fun demand function (e.g., `dd_prob` for HMNL or `vd_dem_vdm` for volumetric demand). For discrete choice, use choice probabilities instead of choice predictions.
#' @param draws ec-draws object (e.g., output from `dd_est_hmnl` or `vd_est_vd`)
#' @param epsilon_not (optional) error realisatins (this helps make curves look smother for voumetric models)
#' @return List containing aggregate demand quantities for each scenario defined by `rel_pricerange`
#' @examples
#' \donttest{
#' data(icecream)
#' #run MCMC sampler (use way more draws for actual use)
#' icecream_est <- icecream %>% dplyr::filter(id<20) %>% 
#' vd_est_vdm(R=2, keep=1)
#' #demand at different price points
#' conddem_scenarios<-
#' ec_demcurve_cond_dem(icecream%>% dplyr::filter(id<20),
#'  icecream%>% dplyr::filter(id<20) %>% pull('Brand')=="Store",
#'  c(.75,1),vd_dem_vdm,icecream_est)
#' }
#' 
#' @seealso [ec_gen_err_normal()] to generate error realization from Normal distribution,
#' [ec_gen_err_ev1()] to generate error realization from EV1 distribution
#' 
#' @export
ec_demcurve_cond_dem=function(ec_long,
                              focal_product,
                              rel_pricerange,
                              dem_fun,
                              draws,
                              epsilon_not=NULL){
  
  #select correct demand function
  if(is.function(dem_fun)){
    demfun={dem_fun}
  }else{
    demfun=eval(parse(text=dem_fun))
  }
  
  #get key attributes
  attr_names=ec_long %>% 
    select_if(!(colnames(.)%in%c('id','task','alt','p'))) %>% colnames
  
  #output
  out=list()
  
  if(is.null(epsilon_not)){
    for(kk in seq_along(rel_pricerange)){
      
      testmarket_temp=ec_long
      testmarket_temp$p[focal_product]=testmarket_temp$p[focal_product]*rel_pricerange[kk]
      
      # out[[kk]]=
      #   testmarket_temp %>% 
      #   demfun(draws) %>% 
      #   mutate(.demdraws=map(.demdraws, sign))%>%
      #   ec_dem_aggregate(attr_names) %>% ec_dem_summarise %>% select(-.demdraws) %>%
      #   bind_cols( scenario=rel_pricerange[kk])
      # 
      
      dem_temp=
      testmarket_temp %>%
        demfun(est=draws) %>%
        tibble::add_column(.isfocal=focal_product) %>%
        dplyr::filter(.isfocal)%>% select(-.isfocal) 
      
      demm_split=
        dem_temp%>%
        split(.[c(attr_names)]) 
      
      demm_split=demm_split[(demm_split %>% map_dbl(nrow))>0]
      #%>% map_dbl(length) %>% drop

      res=demm_split %>%
        map_dfr(.  %>% pull(.demdraws)  %>% 
                  do.call('rbind',.) %>% apply(2,function(x)mean(x[x>0])))%>%
        summarise_all(mean)
        
      out[[kk]]=bind_cols(demm_split %>% map(.%>%select(attr_names)%>%unique),
                res %>% rlang::set_names('condem')) %>%
        bind_cols( scenario=rel_pricerange[kk]) 

        
       
    }
  }else{
    for(kk in seq_along(rel_pricerange)){
      testmarket_temp=ec_long
      testmarket_temp$p[focal_product]=testmarket_temp$p[focal_product]*rel_pricerange[kk]
      
      # out[[kk]]=
      #   testmarket_temp %>% 
      #   demfun(draws) %>% mutate(.demdraws=map(.demdraws, sign))%>%
      #   ec_dem_aggregate(attr_names )%>% ec_dem_summarise  %>% select(-.demdraws) %>%
      #   bind_cols( scenario=rel_pricerange[kk])
      # 
      
      dem_temp=
        testmarket_temp %>%
        demfun(est=draws) %>%
        tibble::add_column(.isfocal=focal_product) %>%
        dplyr::filter(.isfocal)%>% select(-.isfocal) 
      
      demm_split=
        dem_temp%>%
        split(.[c(attr_names)]) 
      
      demm_split=demm_split[(demm_split %>% map_dbl(nrow))>0]
      #%>% map_dbl(length) %>% drop
      
      res=demm_split %>%
        map_dfr(.  %>% pull(.demdraws)  %>% 
                  do.call('rbind',.) %>% apply(2,function(x)mean(x[x>0])))%>%
        summarise_all(mean)
      
      out[[kk]]=bind_cols(demm_split %>% map(.%>%select(attr_names)%>%unique),
                    res %>% rlang::set_names('condem')) %>%
        bind_cols( scenario=rel_pricerange[kk]) 
      
      
      
    }  
  }
  
  return(out)
}






#' Simulate error realization from Normal distribution
#'
#'
#'
#' @param ec_dem discrete or volumetric choice data, with or without x
#' @param draws draws from volumetric demand model
#' @param seed seed for reproducible error realisations; seet is automatically reset of running this function
#' @examples
#' data(icecream)
#' #run MCMC sampler (use way more than 50 draws for actual use)
#' icecream_est <- icecream %>% dplyr::filter(id<100) %>% 
#' vd_est_vdm(R=100, keep=1, error_dist = "Normal")
#' #generate error realizations
#' errs<- ec_gen_err_normal(icecream %>% dplyr::filter(id<100), icecream_est, seed=123)
#' 
#' @return error realizations
#' 
#' @export
ec_gen_err_normal = function(ec_dem, draws, seed=NULL){
  R=dim(draws$thetaDraw)[3]
  set.seed(seed)
  out=matrix(rnorm(nrow(ec_dem)*R),nrow(ec_dem),R)
  set.seed(NULL)
  return(out)
}


#' Simulate error realization from EV1 distribution
#'
#'
#'
#' @param ec_dem discrete or volumetric choice data, with or without x
#' @param draws draws from volumetric demand model
#' @param seed seed for reproducible error realisations; seet is automatically reset of running this function
#' @examples
#' data(icecream)
#' #run MCMC sampler (use way more than 50 draws for actual use)
#' icecream_est <- icecream %>% dplyr::filter(id<100) %>% 
#' vd_est_vdm(R=100, keep=1)
#' #generate error realizations
#' errs<- ec_gen_err_ev1(icecream %>% dplyr::filter(id<100), icecream_est, seed=123)
#' @return error realizations
#' 
#' @export
ec_gen_err_ev1 = function(ec_dem, draws, seed=NULL){
  R=dim(draws$thetaDraw)[3]
  set.seed(seed)
  out=matrix(revd0(nrow(ec_dem)*R,1),nrow(ec_dem),R)
  set.seed(NULL)
  return(out)
}



#' Screening probabilities of choice alternatives
#'
#' Obtain draws of screening probabilities of choiec alternatives
#'
#' @usage ec_screenprob_sr(xd, est, cores=NULL)
#'
#' @param xd data
#' @param est ec-model draws 
#' @param cores (optional) cores
#' 
#' @return Draws of screening probabilities of choice alternatives
#' @examples
#' data(icecream)
#' icecream_est <- icecream %>% filter(id<10) %>% vd_est_vdm_screen(R=10,  price_screen=TRUE)
#' ec_screenprob_sr(icecream %>% filter(id<10), icecream_est) 
#' 
#' @export
ec_screenprob_sr=function(xd,
                          est,
                          cores=NULL){
  
  #screening model type "olumetric-conjunctive" or "volumetric-conjunctive-pr"
  screening_model_type <- est$ec_type
  with_price_screening <- stringr::str_detect(screening_model_type, "-pr")
  
  #cores  
  if(is.null(cores)){
    cores=parallel::detectCores(logical=FALSE)
  }
  message(paste0("Using ",cores," cores"))
  
  
  #re-arrange data
  dat <- 
    xd %>% 
    vd_long_tidy %>% vd_prepare_nox()
  
  #screening-relevant data
  dat$Af <- xd %>% vd_long_tidy %>%attributes() %>% `[[`('Af') %>% as.matrix()
  
  #demand sim
  if(with_price_screening){
  out=
    ec_screenpr_prob_cpp(dat$PP,
                         dat$AA,
                         dat$Af,
                         dat$nalts,
                         dat$tlens,  
                         dat$ntasks,  
                         dat$xfr-1,
                         dat$xto-1,  
                         dat$lfr-1,  
                         dat$lto-1,
                         est$thetaDraw,
                         est$tauDraw,
                         est$tau_pr_draw,
                         cores=cores)  
  }else{
  out=
    ec_screen_prob_cpp(dat$PP,
                   dat$AA,
                   dat$Af,
                   dat$nalts,
                   dat$tlens,  
                   dat$ntasks,  
                   dat$xfr-1,
                   dat$xto-1,  
                   dat$lfr-1,  
                   dat$lto-1,
                   est$thetaDraw,
                   est$tauDraw, 
                   cores=cores)
  }
  
  attributes(out)=NULL
  
  #add draws to data tibble
  xd=as_tibble(xd)
  xd$.screendraws<-map(out,drop)
  
  #add attributes
  attributes(xd)$attr_names <- xd %>% colnames %>% setdiff(c("id","task","alt","x","p" )) %>% str_subset('^\\.', negate = TRUE)
  attributes(xd)$ec_model   <- attributes(est)$ec_model
  
  return(xd)
}




# Inspect Draws --------------------------------------------------------------

#' Obtain MU_theta draws
#'
#'
#' @param draws A list, 'echoice2' draws object
#' @return  A tibble, long format, draws of MU
#' @examples
#' data(icecream)
#' #run MCMC sampler (use way more than 20 draws for actual use
#' icecream_est <- icecream %>% dplyr::filter(id<50) %>% vd_est_vdm(R=20)
#' ec_draws_MU(icecream_est)
#' 
#' @seealso [ec_draws_screen()] to obtain screening parameter draws (where applicable),
#' [ec_trace_MU()] to generate a traceplot of MU_theta draws
#' 
#' @export
ec_draws_MU <- function(draws){
  
  draws$MUDraw %>% as_tibble %>%
    rlang::set_names(draws$parnames) %>% 
    rowid_to_column(var = 'draw') %>%
    pivot_longer(cols = -any_of('draw'), 
                 names_to = 'attribute_level') 
  
}


#' Obtain Screening probability draws
#'
#'
#' @param draws A list, 'echoice2' draws object
#' @return  A tibble, long format, draws of MU
#' @examples
#' data(icecream)
#' #run MCMC sampler (use way more than 20 draws for actual use
#' icecream_scr_est <- icecream %>% dplyr::filter(id<50) %>% vd_est_vdm_screen(R=20)
#' ec_draws_screen(icecream_scr_est)
#' 
#' @seealso [ec_draws_MU()] to obtain MU_theta draws,
#' [ec_trace_screen()] to generate a traceplot of screening draws
#' 
#' @export
ec_draws_screen <- function(draws){
  
  draws$deltaDraw %>% as_tibble %>%
    rlang::set_names(colnames(attributes(draws)$Af)) %>% 
    rowid_to_column(var = 'draw') %>%
    pivot_longer(cols = -any_of('draw'), 
                 names_to = 'attribute_level') 
  
}


#' Generate MU_theta traceplot
#'
#'
#' @param draws A list, 'echoice2' draws object
#' @param burnin burn-in to remove
#' @return A ggplot2 plot containing traceplots of draws
#' @examples
#' data(icecream)
#' #run MCMC sampler (use way more than 20 draws for actual use
#' icecream_est <- icecream %>% dplyr::filter(id<100) %>% vd_est_vdm(R=20)
#' ec_trace_MU(icecream_est)
#' 
#' @seealso [ec_boxplot_MU()] to obtain boxplot
#' 
#' @export
ec_trace_MU <- function(draws, burnin=100){
  
  myplot<-
  draws$MUDraw %>% data.frame %>% as_tibble %>%
    rlang::set_names(draws$parnames) %>% 
    rowid_to_column(var = 'draw') %>%
    pivot_longer(cols = -any_of('draw'), 
                 names_to = 'par') %>%
    left_join(attributes(draws)$ec_data$attributes %>%
                select(attr_level,attribute,lvl,reference_lvl), 
              by=c('par'='attr_level')) %>%
    ggplot(aes(x=draw, y=value)) +
    geom_line() +guides(color='none')+facet_wrap(~par,scales='free_y')
  myplot
}



#' Generate Screening probability traceplots
#'
#'
#' @param draws A list, 'echoice2' draws object, from a model with attribute-based screening
#' @param burnin burn-in to remove
#' @return A ggplot2 plot containing traceplots of draws
#' @examples
#' data(icecream)
#' #run MCMC sampler (use way more than 20 draws for actual use
#' icecream_scr_est <- icecream %>% dplyr::filter(id<20) %>% vd_est_vdm_screen(R=20)
#' ec_trace_screen(icecream_scr_est, burnin=1)
#' 
#' 
#' @seealso [ec_draws_MU()] to obtain MU_theta draws,
#' [ec_boxplot_screen()] to generate boxplot
#' 
#' @export
ec_trace_screen <- function(draws, burnin=100){
  
  myplot<-
  draws$deltaDraw %>% 
    as_tibble(.name_repair = ~make.names(seq_along(.), unique=TRUE)) %>%
    rlang::set_names(colnames(attributes(draws)$Af)) %>% 
    rowid_to_column(var = 'draw') %>%
    pivot_longer(cols = -any_of('draw'), 
                 names_to = 'attribute_level') %>%
    left_join(attributes(draws)$ec_data$attributes %>%
                select(attr_level,attribute,lvl,reference_lvl), 
              by=c('attribute_level'='attr_level')) %>%
    dplyr::filter(draw>burnin) %>%
    ggplot(aes(x=draw, y=value)) +
    geom_line() +guides(color='none')+facet_wrap(~attribute_level,scales='free_y')
  
  myplot
  
}



#' Generate MU_theta boxplot
#'
#'
#' @param draws A list, 'echoice2' draws object
#' @param burnin burn-in to remove
#' @return A ggplot2 plot containing traceplots of draws
#' @examples
#' data(icecream)
#' #run MCMC sampler (use way more than 50 draws for actual use
#' icecream_est <- icecream %>% dplyr::filter(id<100) %>% vd_est_vdm(R=50)
#' ec_boxplot_MU(icecream_est, burnin=1)
#' 
#' @seealso [ec_trace_MU()] to obtain traceplot
#' 
#' @export
ec_boxplot_MU <- function(draws, burnin=100){
  
  myplot<-
  draws$MUDraw %>% 
    as_tibble(.name_repair = ~make.names(seq_along(.), unique=TRUE)) %>%
    rlang::set_names(draws$parnames) %>% 
    rowid_to_column(var = 'draw') %>%
    pivot_longer(cols = -any_of('draw'), 
                 names_to = 'par') %>%
    left_join(x=.,
              y = attributes(draws)$ec_data$attributes %>%
                dplyr::select(all_of(c('attr_level','attribute','lvl','reference_lvl'))), 
              by=c("par"="attr_level")) %>%
    mutate(par=stringr::str_replace_all(par,paste0(.$attribute,":"),"")) %>%
    select(all_of(c('draw','par','value','attribute'))) %>%
    dplyr::filter(.$draw>burnin) %>%
    ggplot2::ggplot(aes(x=par,y=value)) + geom_boxplot() + coord_flip() +
    facet_wrap(~attribute, scales='free')
  myplot
  
}


#' Generate Screening probability boxplot
#'
#'
#' @param draws A list, 'echoice2' draws object, from a model with attribute-based screening
#' @param burnin burn-in to remove
#' @return A ggplot2 plot containing traceplots of draws
#' @examples
#' data(icecream)
#' #run MCMC sampler (use way more than 20 draws for actual use
#' icecream_scr_est <- icecream %>% dplyr::filter(id<100) %>% vd_est_vdm_screen(R=20)
#' ec_boxplot_screen(icecream_scr_est, burnin = 1)
#' 
#' @seealso [ec_draws_MU()] to obtain MU_theta draws,
#' [ec_trace_screen()] to generate traceplot
#' 
#' @export
ec_boxplot_screen <- function(draws, burnin=100){
  
  myplot<-
  draws$deltaDraw %>% 
    as_tibble(.name_repair = ~make.names(seq_along(.), unique=TRUE)) %>%
    rlang::set_names(colnames(attributes(draws)$Af)) %>% 
    rowid_to_column(var = 'draw') %>%
    pivot_longer(cols = -any_of('draw'), 
                 names_to = 'attribute_level') %>%
    left_join(attributes(draws)$ec_data$attributes %>%
                select(all_of(c('attr_level','attribute','lvl','reference_lvl'))), 
              by=c('attribute_level'='attr_level')) %>%
    dplyr::filter(.$draw>burnin) %>%
    ggplot(aes(x=lvl,y=value)) + geom_boxplot() + coord_flip() +
    facet_wrap(~attribute, scales='free_y')
  
  myplot
}




#' Thin 'echoice2'-vd draw objects
#'
#'
#' @param est is an 'echoice2' draw object (list)
#' @param burnin_perc how much burn-in to remove
#' @param total_draws how many draws to keep after thinning
#' @examples
#' data(icecream)
#' #run MCMC sampler (use way more than 50 draws for actual use)
#' icecream_est <- icecream %>% dplyr::filter(id<100) %>% vd_est_vdm(R=50, keep = 1)
#' #without thinning, yields R=50 draWs
#' dim(icecream_est$MUDraw)
#' icecream_est_thinned <- vd_thin_draw(icecream_est,.5)
#' #26 draws left after thinning about half
#' dim(icecream_est_thinned$MUDraw)
#' 
#' @return thinned 'echoice2' draw object (list)
#'
#' @export
vd_thin_draw=function(est, 
                      burnin_perc=.5, 
                      total_draws=NULL){
  
  R = dim(est$thetaDraw)[3]
  
  #burnin
  first_draw=floor(burnin_perc*R)
  usableR   =length(first_draw:R)
  message("Draws post burnin: ", usableR)
  
  check_replicate = F
  if(!is.null(total_draws)){
    check_replicate <- usableR<total_draws
  }
  
  
  if(check_replicate){
    message("More desired draws than original draws")
    total_draws=usableR
    keepdraws=sample(first_draw:R, size = total_draws, replace = TRUE)
  }else{
    all_draws=first_draw:R
    
    #try to thin as much as possible
    chosen_draws = unique(round(seq.int(first_draw, R, length.out=total_draws),0))
    
    check_truncate =F
    if(!is.null(total_draws)){
      check_truncate <- length(chosen_draws)<total_draws
    }
    
    #if missing some draws, randomly grab them from remaining draws
    if(check_truncate){
      message("draws in quence: ",length(chosen_draws), " adding more")
      
      added_draws=sample(all_draws[all_draws%in%chosen_draws],total_draws-length(chosen_draws))
      keepdraws=c(chosen_draws,added_draws)
    } else{
      keepdraws=chosen_draws
    }
  
  }
  
  if(is.null(est$tauDraw)){
    est$thetaDraw=  est$thetaDraw[,,keepdraws]
    est$MUDraw=     est$MUDraw[keepdraws,]
    est$SIGMADraw=  est$SIGMADraw[,,keepdraws]
    est$loglike=    est$loglike[keepdraws,,drop=FALSE]
    est$logpost=    est$logpost[keepdraws,,drop=FALSE]
    est$reject=     est$reject[keepdraws,,drop=FALSE]
  }else{
    est$thetaDraw=  est$thetaDraw[,,keepdraws]
    est$MUDraw=     est$MUDraw[keepdraws,]
    est$SIGMADraw=  est$SIGMADraw[,,keepdraws]
    est$tauDraw=    est$tauDraw[,,keepdraws]
    est$deltaDraw=  est$deltaDraw[keepdraws,]
    est$loglike=    est$loglike[keepdraws,,drop=FALSE]
    est$logpost=    est$logpost[keepdraws,,drop=FALSE]
    est$reject=     est$reject[keepdraws,,drop=FALSE]
  }
  
  if(is.null(est$tau_pr_draw)){
    est$tau_pr_draw       = est$tau_pr_draw[,keepdraws,drop=FALSE]
    est$prscreenMuSigDraw = est$prscreenMuSigDraw[keepdraws,,drop=FALSE]
  }
          
  return(est)
}




# Special helper functions---------------------------------------------------------

#' Convert a vector of choices to long format
#'
#' Converts a vector of choices into a long format data frame, where each row represents 
#' a single choice and contains the choice status for each alternative.
#'
#' @param myvec A vector of choices, where each element represents the index of the chosen alternative.
#' @param all_index A vector of all the possible alternative indices.
#' 
#' @return A tibble with columns 'x', 'task', and 'alt', where 'x' is a binary indicator of whether
#' the alternative was chosen or not, 'task' is the task index, and 'alt' is the alternative index.
#' 
#' @examples
#' #There are 3 alternatives in this task. 
#' #Since there are 3 observations in myvec, there are 3 tasks total.
#' ec_util_choice_to_long(c(1, 2, 1), c(1, 2, 3))
#'
#' @export
ec_util_choice_to_long <- function(myvec, all_index) {
  # Create a matrix of zeros with dimensions n x p, 
  # where n is the length of myvec and p is the length of all_index
  dummy_mat <- matrix(0, nrow = length(myvec), ncol = length(all_index))
  
  # Loop through each index in myvec and set the corresponding element in dummy_mat to 1
  for (i in 1:length(myvec)) {
    dummy_mat[i, myvec[i]] <- 1
  }
  long_choices= as.vector(t(dummy_mat))
  long_task   = rep(seq_along(myvec), each=length(all_index))
  long_alt    = rep(seq_along(all_index),length(myvec))
  
  return(tibble(x    = long_choices,
                task = long_task,
                alt  = long_alt))
}


#' Convert "list of lists" format to long "tidy" format
#'
#' @param data_lol A list of data frames containing design matrices and response vectors
#' @param X The column name of the design matrix, default: "X"
#' @param y The column name of the response vector, default: "y"
#'
#' @return A tidy data frame with columns for each design matrix column, the response vector,
#' and an id column indicating which data frame the row came from
#' 
#' @examples
#' loldata<-list()
#' loldata[[1]]=list()
#' loldata[[1]]$y = c(1,2)
#' loldata[[1]]$X= data.frame(brand1=c(1,0, 1,0),brand2=c(0,1, 0,1),price=c(1,2))
#' loldata[[2]]=list()
#' loldata[[2]]$y = c(1,1)
#' loldata[[2]]$X= data.frame(brand1=c(1,0, 1,0),brand2=c(0,1, 0,1),price=c(1,2))
#' ec_lol_tidy1(loldata)
#' 
#' @export
ec_lol_tidy1 <- function(data_lol, X="X", y="y"){
  
  all_possible_responses <-
    map(data_lol,`[[`,y) %>% do.call("c",.) %>% unique()
  
  combined_design_matrices <- 
    map(data_lol,`[[`, X) %>% 
    map_dfr(~data.frame(.),.id = "id")
  
  combined_responses <-  
    map(data_lol,`[[`,y) %>%
    map_dfr(~ec_util_choice_to_long(.,all_possible_responses))
  
  combined_design_matrices %>%
    add_column(combined_responses) %>% 
    relocate(c(task,alt),.after=id) %>%
    as.data.frame() %>% as_tibble() %>% return()
}



#' Find mutually exclusive columns
#'
#' This function finds pairs of columns in a data frame that are mutually exclusive, i.e., that never have positive values at the same time.
#'
#' @param data_in A data frame containing the data.
#' @param filtered A logical value indicating whether to return only the mutually exclusive pairs (TRUE) or all pairs (FALSE). Default is TRUE.
#' 
#' @return A tibble containing all pairs of mutually exclusive columns in the data frame.
#' 
#' @examples
#' minidata=structure(list(id = c("1", "1", "1", "1", "2", "2", "2", "2"), 
#' task = c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L), 
#' alt = c(1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L), 
#' brand1 = c(1, 0, 1, 0, 1, 0, 1, 0), 
#' brand2 = c(0, 1, 0, 1, 0, 1, 0, 1), 
#' price = c(1, 2, 1, 2, 1, 2, 1, 2), 
#' x = c(1, 0, 0, 1, 1, 0, 1, 0)), 
#' class = c("tbl_df", "tbl", "data.frame"), row.names = c(NA, -8L))
#' ec_util_dummy_mutualeclusive(minidata)
#'
#' 
#' @export
ec_util_dummy_mutualeclusive= function(data_in, filtered=TRUE){
  
  #helper  
  ec_util_mut_exc = function(a,b) max(sign(a)+sign(b))<=1
  
  #select relevant columns
  dat<-data_in %>% select(-any_of(c("id","task","alt","x")))
  
  combs       <- t(utils::combn(seq_len(ncol(dat)),2))
  combs_name  <- t(utils::combn(colnames(dat),2)) %>% as_tibble()
  
  m_ex=rep(NA,nrow(combs))
  for(kk in seq_len(nrow(combs))){
    m_ex[kk]<-ec_util_mut_exc(dat[combs[kk,1]], dat[combs[kk,2]])
  }
  if(filtered){
    combs_name %>% add_column(mut_ex=m_ex) %>% filter(mut_ex=TRUE) %>% return()
  }else{
    combs_name %>% add_column(mut_ex=m_ex) %>% return()
  }
} 



#' Converts a set of dummy variables into a single categorical variable
#'
#' Given a set of dummy variables, this function converts them into a single
#' categorical variable. The categorical variable is created by determining
#' which variables are active (i.e. have a value of 1) for each observation and
#' assigning a category based on the set of active variables. If necessary, a
#' reference level can be specified to ensure that all possible categories are
#' represented. Often, all brands of a brand attribute are added as brand 
#' intercepts, while other categorical attributes are coded with respect to a
#' reference level.
#'
#' @param data_in a data frame containing the dummy variables
#' @param set_members a character vector of the names of the dummy variables
#' @param attribute_name a character string representing the name of the new
#'   categorical variable
#' @param ref_level a character string representing the name of the reference
#'   level. If specified, a new dummy variable will be created for this level,
#'   and it will be used as the reference category for the categorical variable.
#'   Defaults to NULL.
#'
#' @return a data frame with the same columns as \code{data_in}, except for the
#'   dummy variables in \code{set_members}, which are replaced with the new
#'   categorical variable \code{attribute_name}
#'
#' @examples
#' minidata=structure(list(id = c("1", "1", "1", "1", "2", "2", "2", "2"), 
#' task = c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L), 
#' alt = c(1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L), 
#' brand1 = c(1, 0, 1, 0, 1, 0, 1, 0), 
#' brand2 = c(0, 1, 0, 1, 0, 1, 0, 1), 
#' price = c(1, 2, 1, 2, 1, 2, 1, 2), 
#' x = c(1, 0, 0, 1, 1, 0, 1, 0)), 
#' class = c("tbl_df", "tbl", "data.frame"), row.names = c(NA, -8L))
#' 
#' minidata %>% ec_undummy(c('brand1','brand2'),"brand")
#'
#' 
#' @export
ec_undummy=function(data_in, 
                    set_members, 
                    attribute_name, 
                    ref_level=NULL){
  
  datax=data_in %>% select(all_of(set_members)) 
  row_sum_check<-datax %>% rowSums() %>% min
  
  if(row_sum_check!=1){
    if(is.null(ref_level)) stop("You need to define the name of the reference level.")
    datax<-datax %>% mutate(!!(ref_level):= 1-rowSums(select(.,everything())))
    set_members=c(set_members,ref_level)
  }else{
    if(!is.null(ref_level)) warning("Reference level is not needed.")
  }
  
  data_in %>%
    add_column(
      !!(attribute_name):=
      (sign(data.matrix(datax))>0)  %>% which(arr.ind = TRUE) %>% 
        as.data.frame%>%arrange(across("row")) %>% pull("col") %>% {set_members[.]} %>% factor(),
      .before = set_members[1]
    ) %>% 
    select(-any_of(set_members)) %>%
    return()
  
}



#' Convert dummy-coded variables to low/high factor
#' 
#' @param vec_in A vector of dummy-coded variables (0/1) 
#' @return A factor vector with levels "low" and "high"
#' 
#' @examples
#' ec_undummy_lowhigh(c(0,1,0,1,1))
#'
#'
#' @export
ec_undummy_lowhigh=function(vec_in){
  vec_in %>% factor() %>% fct_recode("low"="0","high"="1") %>% return()
}


#' Convert dummy-coded variables to yes/no factor
#' 
#' @param vec_in A vector of dummy-coded variables (0/1) 
#' @return A factor vector with levels "no" and "yes"
#' 
#' @examples
#' ec_undummy_yesno(c(0,1,0,1,1))
#'
#'
#' @export
ec_undummy_yesno=function(vec_in){
  vec_in %>% factor() %>% fct_recode("no"="0","yes"="1") %>% return()
}


#' Convert dummy-coded variables to low/medium/high factor
#' 
#' @param vec_in A vector of dummy-coded variables (0/1/2) 
#' @return A factor vector with levels "low", "medium" and "high"
#' 
#' @examples
#' ec_undummy_lowmediumhigh(c(0,1,2,1,0,2))
#'
#'
#' @export
ec_undummy_lowmediumhigh=function(vec_in){
  vec_in %>% factor() %>% fct_recode("low"="0","medium"="1","high"="2") %>% return()
}
