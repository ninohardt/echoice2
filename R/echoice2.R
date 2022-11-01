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




#' @importFrom tidyr pivot_longer
NULL




#' @importFrom stringr str_subset
NULL

#' @importFrom stringr str_remove
NULL

#' @importFrom stringr str_extract
NULL




#' @importFrom purrr map
NULL

#' @importFrom purrr map_dfr
NULL


#' @importFrom tidyselect any_of
NULL

#' @importFrom tidyselect all_of
NULL

#' @importFrom tidyselect contains
NULL

#' @importFrom tidyselect last_col
NULL





#' @importFrom tibble tibble
NULL

#' @importFrom tibble as_tibble
NULL

#' @importFrom tibble add_row
NULL

#' @importFrom tibble rowid_to_column
NULL


# Utilities ---------------------------------------------------------------


#' Generate huxtable for model comparison
#'
#' This requires huxtable package
#' Pass a list of dd or vd models in a named list to this function, and it will generate a huxtable.
#' This will include posterior mean and standard deviations of the upper-level parameters
#' The table can easily by converted to latex using huxtable functions.
#'
#' @param model_list model_list
#'
#' @usage ht_modelMuCompare(model_list)
#
#' @export
ht_modelMuCompare=function(model_list){
  
  n_models=length(model_list)
  
  paras=
    model_list %>%
    map_dfr(ec_estimates_MU, .id = 'model') %>% select(attribute,parameter,mean,sd,sig,model) %>% 
    pivot_wider(names_from = model,values_from = all_of(c('mean','sd','sig')))
  
  #arrangement
  paras_order=names(paras) %>% enframe %>%  mutate(modelname=str_remove(str_remove(str_remove(value,'mean_'),'sd_'),'sig_')) %>% filter(value!='attribute') %>% filter(value!='parameter') %>% mutate(mean=str_detect(value,'mean'),sd=str_detect(value,'sd'),sig=str_detect(value,'sig')) %>% filter(sig==F) %>% arrange(modelname) %>%pull(name)
  paras_order=c(1,2,paras_order)
  
  #significant paras
  paras_order_sig=names(paras) %>% enframe %>%  mutate(modelname=str_remove(str_remove(str_remove(value,'mean_'),'sd_'),'sig_')) %>% filter(value!='attribute') %>% filter(value!='parameter') %>% mutate(mean=str_detect(value,'mean'),sd=str_detect(value,'sd'),sig=str_detect(value,'sig')) %>% filter(sig==T) %>% arrange(modelname) %>%pull(name)
  paras2=paras[,paras_order]
  paras2_sig=paras[,paras_order_sig]
  
  #a few steps to build
  
  #add row identifier
  prep_1 = paras2 %>% rowid_to_column()
  
  #don't repeat the grouping attribute name
  prep_2=
    prep_1  %>% split(.$attribute) %>% 
    map((function(x){x$attribute=NA;x})) %>% map(add_row,.before=1) %>%   
    imap_dfr((function(x,y){x$attribute[1]=y;x}))
  
  #clean up
  mpar=prep_2[NA,][1,]
  mpar[1,2]='model'
  
  #put together
  prep_3=
    bind_rows(
      prep_1[1,],
      prep_2,
      mpar,
      prep_1[base::setdiff(prep_1$rowid,prep_2$rowid),][-1,]) %>% select(-rowid)
  
  #some arrangement
  prep_4=prep_3[,seq(1,2+n_models*2)]
  
  #generate huxtable
  ht=prep_4 %>% 
    huxtable::huxtable() %>%
    set_all_padding(0) %>%
    set_row_height(everywhere,0)%>% set_align(everywhere,2:ncol(ht),'right')
  
  
  #bold-face significant paras
  pickrows= (ht[,1]=="") %>% replace_na(T)
  huxtable::bold(ht)[pickrows,2+seq(1,n_models*2,by=2)] <- (paras2_sig %>% as.matrix %>% replace_na(F))
  
  #number format
  huxtable::number_format(ht)[-1,-(1:2)]= "%.2f"
  
  #merging attribute-name cells row-wise for compact look
  mergerows=which(!is.na(ht$attribute))[-1]
  for(i in seq_along(mergerows)){
    ht=merge_across(ht,mergerows[i],seq_len(ncol(ht)))
  }
  ht[1,1]=""
  
  #attr column in bold
  huxtable::bold(ht)[,1]=TRUE
  
  #modelname headers
  modelnames=sort(names(model_list))
  ht=ht %>% insert_row(c("","",rep(modelnames,each=2)))
  headers=seq(1,n_models*2,by=2)+2
  
  ht[2,-(1:2)]=rep(c('',''),n_models)
  
  for(i in seq_len(n_models)) {
    ht=ht %>% merge_cells(c(1,1),c(headers[i],headers[i]+1)) %>% set_align(1,headers[i],'centre')
  }
  
  #standard devs in parentheses
  huxtable::number_format(ht)[-(1:2),headers+1]= "(%.2f)"
  
  #
  ht[2,2]=""
  ht=ht[-2,]
  
  #padding  
  bottom_padding(ht)=1
  top_padding(ht)=1
  right_padding(ht)[,2+seq(1,n_models*2,by=2)]=1
  right_padding(ht)[,2+seq(2,n_models*2,by=2)]=10
  left_padding(ht)[,2+seq(2,n_models*2,by=2)]=1
  
  #borders
  for(i in seq_len(n_models)) {
    bottom_border(ht)[1,seq(headers[i],headers[i]+1)] = brdr(0.4, "solid", "black")
  }
  top_border(ht)[1,] = brdr(1, "solid", "black")
  bottom_border(ht)[nrow(ht),] = brdr(1, "solid", "black")
  
  return(ht)
}


#' Generate huxtable for model comparison - screening vs no screening
#'
#' This requires huxtable package
#' Pass a list of dd or vd models in a named list to this function, and it will generate a huxtable.
#' This will include posterior mean and standard deviations of the upper-level parameters
#' The table can easily by converted to latex using huxtable functions.
#'
#' @param vdm vdm object
#' @param vdmsrpr vdmsrpr object
#
#' @export
ht_screenCompare=
  function(vdm,vdmsrpr){
    
    paras = vdmsrpr %>% ec_estimates_MU() %>%transmute(attribute,lvl,par,mean,sd,sig,para='theta',model='vdsr') %>%
      bind_rows(vdmsrpr %>% ec_estimates_screen()%>%transmute(attribute,lvl,par,mean,sd,sig=NA,para='delta',model='vdsr')) %>%
      bind_rows(vdm %>% ec_estimates_MU() %>%transmute(attribute,lvl,par,mean,sd,sig,para='theta',model='vd'))%>%
      pivot_wider(names_from = c(model,para), id_cols = 1:3, values_from = c(mean,sd, sig)) %>% arrange(par!='int',(attribute))%>% relocate(contains('vdsr'),.after=last_col()) %>% relocate(contains('sig_'),.after = last_col()) %>% `[`(
        c("attribute", "lvl", "par", 
          "mean_vd_theta", "sd_vd_theta", 
          "mean_vdsr_theta", "sd_vdsr_theta", "mean_vdsr_delta", "sd_vdsr_delta", 
          "sig_vd_theta", "sig_vdsr_theta", "sig_vdsr_delta"))
    
    paras=
      rbind(
        paras, 
        c(NA, NA, 'MU_delta_p', rep(NA,4), mean(vdmsrpr$prscreenMuSigDraw[,1]), sd(vdmsrpr$prscreenMuSigDraw[,1]), NA, NA, NA),
        c(NA, NA, 'SIG_delta_p', rep(NA,4), mean(vdmsrpr$prscreenMuSigDraw[,2]), sd(vdmsrpr$prscreenMuSigDraw[,2]), NA, NA, NA))
    
    
    #arrangement
    paras2=paras %>% select(-contains('sig_'))
    paras2_sig= paras %>% select(contains('sig_')) %>% mutate_all(as.logical)
    
    #add row identifier
    prep_1 = paras2 %>% rowid_to_column()
    prep_1$par[!is.na(prep_1$attribute)]=prep_1$lvl[!is.na(prep_1$attribute)]
    
    
    #don't repeat the grouping attribute name
    prep_2=
      prep_1  %>% split(.$attribute) %>%
      map((function(x){x$attribute=NA;x})) %>% map(add_row,.before=1) %>%
      imap_dfr((function(x,y){x$attribute[1]=y;x}))
    
    
    prep_2$par[is.na(prep_2$rowid)]= prep_2$attribute %>% unique() %>% na.omit()
    
    
    #clean up
    mpar=prep_2[NA,][1,]
    mpar[1,2]='model'
    
    
    #put together
    prep_3=
      bind_rows(
        prep_1[1,],
        prep_2,
        prep_1[base::setdiff(prep_1$rowid,prep_2$rowid),][-1,]) %>% select(-rowid)
    
    
    prep_3$par[!is.na(prep_3$lvl)] = paste0("~",prep_3$par[!is.na(prep_3$lvl)])
    
    #some arrangement
    prep_4=prep_3[,-(1:2)]
    
    
    #generate huxtable
    ht=prep_4 %>% huxtable::huxtable() %>%  
      set_all_padding(0) %>%
      set_row_height(everywhere,0)
    
    
    #bold-face significant paras
    #pickrows= (ht[,1]=="") %>% replace_na(T)
    
    pickrows=c(FALSE,drop(!is.na(prep_4[-c(nrow(prep_4)-1,nrow(prep_4)),2])))
    
    
    huxtable::bold(ht[-c(nrow(ht)-1,nrow(ht)),])[pickrows,  c(2,4) ] <- (paras2_sig[,1:2] %>% as.matrix %>%na.omit() %>% set_attrs("na.action"=NULL) )
    
    #number format
    huxtable::number_format(ht)[-1,-(1)]= "%.2f"
    
    #attributes bold
    #huxtable::bold(ht)[drop(is.na(ht[,2])&is.na(ht[,6])),1]=T
    huxtable::italic(ht)[drop(is.na(ht[,2])&is.na(ht[,6])),1]=T
    
    
    
    ht[1,1]=""
    
    
    ht = ht %>% set_align(everywhere,2:ncol(ht),'right')
    
    ht=ht %>% 
      insert_row(c("",rep('theta',each=2), rep('theta',each=2), rep('delta',each=2)   )) %>%
      insert_row(c("",rep('vd',each=2), rep('vd-sr',each=4)   ))
    
    
    ht[3,-(1)]=rep(c('mean','sd'),3)
    ht=ht %>% merge_cells(c(1,1),c(2,3)) %>% set_align(1,2:3,'centre')
    ht=ht %>% merge_cells(c(1,1),c(4,7)) %>% set_align(1,4:7,'centre')
    
    ht=ht %>% merge_cells(c(2,2),c(2,3)) %>% set_align(2,2:3,'centre')
    ht=ht %>% merge_cells(c(2,2),c(4,5)) %>% set_align(2,4:5,'centre')
    ht=ht %>% merge_cells(c(2,2),c(6,7)) %>% set_align(2,6:7,'centre')
    
    ht=ht[-3,]
    
    
    #standard devs in parentheses
    huxtable::number_format(ht)[-(1:2),c(3,5,7)]= "(%.2f)"
    
    
    ht=
      ht %>% 
      set_right_padding(everywhere,3,3) %>% set_right_padding(everywhere,3,2) %>%  set_right_padding(everywhere,5,2) %>%
      set_bottom_border(1,2:3,1) %>% set_bottom_border(1,4:7,1) %>%
      set_bottom_border(2,2:3,1) %>% set_bottom_border(2,4:5,1) %>% set_bottom_border(2,6:7,1) 
    
    
    
    return(ht)
  }




#' Attribute pipe
#'
#' Read attributes of an object, pipe-style
#'
#' @param obj object
#' @param attrname Attribute to pull
#
#' @export
`%.%` <- function(obj,attrname) (attributes(obj))[[attrname]]



# data manipulation -------------------------------------------------------


#clean data
vd_janitor=function(vd, maxquant=999){
  `%!in%` <- Negate(`%in%`)
  
  #remove too high volumes
  fid_toomuch   = vd %>% filter(x>maxquant) %>% pull(id) %>% unique
  
  #remove all0s
  fid_toolittle = vd %>% group_by(id) %>% summarise(.s=sum(x)) %>% filter(.s==0) %>% pull(id)
  
  #combine filter
  fid_all = base::intersect(fid_toolittle, fid_toomuch)
  
  filter_summary=
    tibble(
      `overall`=fid_all%>%n_distinct,
      `large quantities`=fid_toolittle%>%n_distinct,
      `all 0` = fid_toolittle %>% n_distinct())
  
  
  #filter and return
  vd = vd %>% filter(id %!in% fid_all) 
  
  print("Filter summary:")
  print(filter_summary)
  
  print("Data:")
  attributes(vd)$filter=filter_summary
  
  return(vd)
}


#utility function
#tidying a single element of a list-of-lists choice dataset
vd_lol_tidyelement=function(da){
  return(
    bind_cols(
      tibble(
      task=rep(1:nrow(da$X),each=ncol(da$X)),
      alt=rep(1:ncol(da$X),nrow(da$X)),
      x=as.vector(t(da$X)),
      p=as.vector(t(da$P))),
    as_tibble(da$A))
  )
}


#' Convert list-of-lists choice format to stacked tidy data
#'
#'
#' @param datalist is a list of lists (one list per unit) with elements A, X, P
#'
#' @return tibble/data.frame in stacked format
#' Columns of output data.frame:
#' id; task: ask/choice occasion number; alt: alternative/product number; x: quantity; p: price; attributes
#'
#' @export
vd_lol_tidy=function(datalist){
  out = mutate(
    map_dfr(datalist,vd_lol_tidyelement, .id = 'id'),
    id=as.integer(id))
  
  attributes(out)$ec_data = list(choice_type='volumetric',
                                 screening='none',
                                 data_type='vd_tidy_choice')
  return(out)
}




#' Dummy-code a categorical variable
#'
#'
#' @param data one column of categorical data to be dummy-coded
#'
#' @return tibble with dummy variables
#'
#' @export
dummyvar<-function(data){
  #retains missing values as NA
  out = model.matrix.lm(~.-1,
                        data=data, 
                        na.action = "na.pass") 
  attributes(out)[c("assign", "contrasts")]=NULL
  return(as_tibble(out))
}


#' Dummy-code select variables in a tibble
#'
#'
#' @param dat tibble
#' @param sel vector of strings, variables to be replaced by dummy coding
#'
#' @return tibble with dummy variables
#'
#' @export
dummify=function(dat, sel){
  dat=as_tibble(dat)
  for(i in seq_along(sel)){
    selv=sel[i]
    dummidata <- dat %>% select(all_of(selv)) %>% dummyvar()
    dat<- dat %>% add_column(dummidata, .before = selv) %>% 
      select(-all_of(selv)) %>% 
      rename_all(gsub, pattern = paste0("^(",selv,")"), replacement = paste0(selv,":"))
  }
  
  return(dat)
}



#' Obtain attributes and levels from tidy choice data
#'
#'
#' @param tdc tibble
#'
#' @return tibble
#'
#' @export
get_attr_lvl=function(tdc){
  tdc %>%
    select(-any_of(c('id','task','alt','p','x')))%>% 
    names %>% enframe %>% 
    mutate(attribute=str_extract(value,"^.*?(?=\\:)")) %>%
    mutate(lvl=str_remove(value,attribute)) %>%
    mutate(lvl=str_remove(lvl,"^(:)")) %>%
    group_by(attribute) %>%
    mutate(reference_lvl=first(lvl)) %>%
    mutate(reference=ifelse(lvl==reference_lvl,1,0))%>%
    mutate(lvl_abbrv=abbreviate(lvl))%>%
    rename(attr_level=value)
}




#' vd_tidy_choice data from long-format choice data
#'
#'
#' @param longdata tibble
#'
#' @return tibble
#' 
#'
#' @export
vd_long_tidy<-function(longdata){

  catvars = longdata %>% select_if(is.factor) %>% names
  dummified =
    longdata %>% 
    dummify(catvars) 
  attrs=dummified %>%get_attr_lvl
  refcats=attrs %>% filter(reference==1) %>% pull(attr_level)
  
  out<-dummified %>% select(-any_of(refcats)) %>% add_column(int=1,.after='p')
  
  attributes(out)$ec_data = list(choice_type='volumetric',
                                 data_type='vd_tidy_choice',
                                 attributes=attrs)
  
  attributes(out)$Af = dummified %>% select_if(!(colnames(.)%in%c('id','task','alt','p','x')))
  
  return(out)
}







#' Prepare choice data for analysis
#'
#' Pre-computing task-wise total expenditures `sumpsx` and indices `xfr`,`xto`,`lfr`,`lto` speeds up computation
#'
#'
#' @param dt is tidy choice data (columns: id, task, alt, x, p, attributes)
#' @param Af (optional) contains a full design matrix (for attribute-based screening), or, more generally, a design matrix used for attribute-based screening
#'
#' @return list containing information for estimation functions
#'
#' @examples
#' \dontrun{
#' #Minimal example:
#' #One attribute with 3 levels, 2 subjects, 3 alternatives, 2 tasks
#'Af=
#'  rbinom(12,2,.5) %>% factor %>% fct_relabel(~paste0('attr',.x)) %>% 
#'  enframe %>% add_column(val=1) %>% 
#'  pivot_wider(names_from = value, values_from=val) %>% map_df(replace_na,0) %>%
#'  select(-name) %>%
#'  add_column(  id=c(rep(1L,6),rep(2L,6)),
#'               task=c(1L,1L,1L, 2L,2L,2L, 1L,1L,1L, 2L,2L,2L),
#'               alt=c(1L,2L,3L, 1L,2L,3L, 1L,2L,3L, 1L,2L,3L), .before = 1)
#'
#Choice data
#'dt=tibble( id=c(rep(1L,6),rep(2L,6)),
#'           task=c(1L,1L,1L, 2L,2L,2L, 1L,1L,1L, 2L,2L,2L),
#'           alt=c(1L,2L,3L, 1L,2L,3L, 1L,2L,3L, 1L,2L,3L),
#'           x=c(1,0,2, 1,0,1, 2,3,1, 1,0,1),
#'           p=c(0, 1, 1, 1, 2, 0, 2, 2, 1, 2, 1, 1),
#'           Af%>%select(attr2:attr1))
#' 
#' dt %>% vd_prepare 
#' dt %>% vd_prepare(Af) #with data for attribute-based screening
#' }
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
        tauconst=1-(foo %>% filter(x>0) %>%group_by(id) %>% summarise_all(max) %>% select(-id,-task,-alt,-x))
      }else{
        stop("Could not match full attribute tibble with choice data")
      }
      
      #bind_cols(Af)%>% filter(x>0) %>%group_by(id) %>% summarise_all(max) %>% select(-id,-task,-alt,-x)
    }else{
      message("Af does not contain id, task, alt columns. Assuming that attribute columns are properly sorted...")
      out$AAf=Af%>% as.matrix()
      tauconst=1-(dt %>% select(id,task,alt,x) %>%bind_cols(Af)%>% filter(x>0) %>%group_by(id) %>% summarise_all(max) %>% select(-id,-task,-alt,-x))
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
#' This is mostly for internal use
#'
#'
#' @param dt is tidy choice data (columns: id, task, alt, x, p, attributes)
#' @param Af Af
#'
#' @return list containing information for estimation functions
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
        tauconst=1-(foo %>% filter(x>0) %>%group_by(id) %>% summarise_all(max) %>% select(-id,-task,-alt,-x))
      }else{
        stop("Could not match full attribute tibble with choice data")
      }
      
      #bind_cols(Af)%>% filter(x>0) %>%group_by(id) %>% summarise_all(max) %>% select(-id,-task,-alt,-x)
    }else{
      message("Af does not contain id, task, alt columns. Assuming that attribute columns are properly sorted...")
      out$AAf=Af%>% as.matrix()
      tauconst=1-(dt %>% select(id,task,alt,x) %>%bind_cols(Af)%>% filter(x>0) %>%group_by(id) %>% summarise_all(max) %>% select(-id,-task,-alt,-x))
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


#' Log Marginal Density (Newton-Raftery)
#'
#'
#' @param ll vector of likelihood draws
#'
#' @return LMD (double)
#'
#'
#' @examples
#' \dontrun{
#' logMargDenNR(ll)dd_dem_sr
#' }
#' @export
logMargDenNR=function(ll) 
{
  med = median(ll)
  return(med - log(mean(exp(-ll + med))))
}


#' Log Marginal Density (Newton-Raftery)
#' 
#' Very rough approximation of LMD based on Newton-Raftery.
#' It is not the most accurate, but a very fast method.
#' Draws are split in 4 equal parts from start to finish, and LMD
#' is computed for each part. This helps to double-check convergence.
#'
#'
#' @param M echoice draw object
#'
#' @return tibble with LMDs (first 25% of draws, next 25% of draws, ...) 
#'
#' @export
ec_lmd_NR=function(M){
  return(
    drop(M$loglike) %>% 
      enframe(name = 'draw') %>% 
      mutate(part=ntile(draw,4)) %>% 
      group_by(part) %>% summarise(lmd=logMargDenNR(value), .groups='drop') %>% 
      mutate(part=part/4)
  )
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



#' Summarise attributes and levels
#' 
#' This functions looks for categorical attributes and summaries their levels
#' This is helpful when evaluating a new choice data file.
#'
#'
#' @param data_in long-format choice data
#'
#' @return tibble with summary
#'
#' @export
ec_summarize_attrlvls<-function(data_in){
  return(
  data_in %>% select(-any_of(c('id','task','alt','p','x'))) %>% map(table) %>% 
    map(names) %>% map(paste,collapse=', ') %>% as_tibble() %>% 
    pivot_longer(everything()) %>% set_names(c('attribute','levels')) )
}


# Estimates ---------------------------------------------------------------


#' Obtain upper level model estimates
#'
#'
#' @param est is an echoice draw object (list)
#' @param quantiles quantile for CI
#' @return tibble with MU (upper level) summaries
#' @examples
#' \dontrun{
#' est %>% ec_estimates_MU
#' est %>% ec_estimates_MU %>% 
#' ggplot(aes(x=par,y=mean))+geom_bar(stat='identity')+
#' facet_wrap(~attribute,scales = 'free') + coord_flip()
#' }
#' @export
ec_estimates_MU=function(est, quantiles=c(.05,.95)){
  
  quantiles_name=paste0("CI-",quantiles*100,"%")
  parnames=est$parnames
  
  estimates=  
    est$MUDraw %>% 
    as_tibble() %>% 
    set_names(parnames) %>%
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
#' @param est is an echoice draw object (list)
#' @return estimates of upper level correlations
#' @examples
#' \dontrun{
#' #obtain correlations and visualize using corrr package
#' est %>% ec_estimates_SIGMA_corr %>% 
#' corrr::as_cordf() %>% 
#' corrr::rplot(print_cor = TRUE, legend = FALSE)+ 
#'   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#' }
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
#' @param est is an echoice draw object (list)
#' @return estimates of upper level covariance
#' @export
ec_estimates_SIGMA=function(est){
  parnames=est$parnames
  rownames(est$SIGMADraw)=parnames
  colnames(est$SIGMADraw)=parnames
  return(est$SIGMADraw %>% apply(1:2,mean))
}


#' Summarize attribute-based screening parameters
#'
#'
#' @param est is an echoice draw object (list)
#' @param quantiles quantile for CI
#' @return tibble with screening summaries
#' @export
ec_estimates_screen=function(est,quantiles=c(.05,.95)){
  
  quantiles_name=paste0("CI-",quantiles*100,"%")
  
  out<-
    est$deltaDraw %>% as_tibble %>%
    set_names(colnames(attributes(est)$Af)) %>% 
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
      left_join(est$dat$tauconst %>% colMeans() %>% enframe(name = 'par', value = 'limit'))
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





# Volumetric Demand Estimation --------------------------------------------



#' Estimate volumetric demand model, EV1 errors
#'
#'
#' @param vd volumetric demand data (long format)
#' @param tidy apply echoice tidier
#' @param R draws
#' @param keep thinning
#' @param cores no of CPU cores to use (default: auto-detect)
#' @param control list containing additional settings
#' 
#' @return est ec-draw object (List)
#' 
#' @seealso [vd_dem_vdm()] to generate demand predictions based on this model
#' 
#' @examples
#' \dontrun{
#' data(icecream)
#' icecream_est <- icecream %>% vd_est_vdm(R=50000)
#' }
#' @export
vd_est_vdm=
  function(vd,
           tidy=T,
           R=100000, 
           keep=10,
           cores=NULL,
           control=list(include_data=TRUE)){
    
    #check input data
    if(!vd_check_long(vd)) stop("Check data")
    
    #integer variables
    vd<-vd %>% mutate(task=as.integer(task),
                      alt =as.integer(alt))
    
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
    
    #Add data information
    out$A_names<-colnames(dat$AA)
    out$parnames=c(colnames(dat$AA),'sigma','gamma','E')
    if(!is.null(dat$attributes_levels)){
      out$attributes_levels=dat$attributes_levels
    }
    
    #Add model information
    out$ec_type="volumetric-compensatory"
    out$error_specification="EV1"
    out$ec_type_short="VD-comp"
    out$Prior=Prior
    
    attributes(out)$Af<-attributes(dat)$Af
    attributes(out)$ec_data<-attributes(dat)$ec_data
    
    
    ec_model = list(model_name_full="volumetric-compensatory-EV1",
                    model_name_short="VD-comp",
                    error_specification="EV1",
                    model_parnames=c('sigma','gamma','E'),
                    Prior=Prior)
    
    attributes(out)$ec_model=ec_model
    
    
    #Add training data
    if(control$include_data){
      out$dat=dat
    }
    
    return(out)
  }




#' Estimate volumetric demand model (Normal Error)
#'
#'
#' @param vd volumetric demand data (long format)
#' @param R draws
#' @param keep thinning
#' @param cores no of CPU cores to use (default: auto-detect)
#' @param control list containing additional settings
#' 
#' @return est ec-draw object (List)
#' 
#' @examples
#' \dontrun{
#' icecream_est <- icecream %>% vd_est_vdmn(R=50000)
#' }
#' @export
vd_est_vdmn = function(vd,
                       R=100000, 
                       keep=10,
                       cores=NULL,
                       control=list(include_data=TRUE)){
  
  #check input data
  if(!vd_check_long(vd)) stop("Check data")
  
  #integer variables
  vd<-vd %>% mutate(task=as.integer(task),
                    alt =as.integer(alt))
  
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
  Bbar=matrix(rep(0,ncol(dat$AA)+3), ncol=ncol(dat$AA)+3)
  A=0.01*diag(1)
  nu=ncol(dat$AA)+9
  V=(ncol(dat$AA)+9)*diag(ncol(dat$AA)+3)
  Prior=list(Bbar=Bbar,A=A,nu=nu,V=V)
  
  #Run model
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
  
  #Add data information
  out$A_names<-colnames(dat$AA)
  out$parnames=c(colnames(dat$AA),'sigma','gamma','E')
  if(!is.null(dat$attributes_levels)){
    out$attributes_levels=dat$attributes_levels
  }
  
  #Add model information
  out$ec_type="volumetric-compensatory"
  out$error_specification="Normal"
  out$ec_type_short="VD-comp-norm"
  out$Prior=Prior
  
  attributes(out)$Af<-attributes(dat)$Af
  attributes(out)$ec_data<-attributes(dat)$ec_data
  
  
  ec_model = list(model_name_full="volumetric-compensatory-Normal",
                  model_name_short="VD-comp-norm",
                  error_specification="Normal",
                  model_parnames=c('sigma','gamma','E'),
                  Prior=Prior)
  
  attributes(out)$ec_model=ec_model
  
  
  #Add training data
  if(control$include_data){
    out$dat=dat
  }
  
  return(out)
}




#' Estimate volumetric demand model with attribute-based conjunctive screening
#'
#' See https://dx.doi.org/10.2139/ssrn.2770025 for more details
#'
#'
#' @param vd volumetric demand data (long format)
#' @param R draws
#' @param keep thinning
#' @param cores no of CPU cores to use (default: auto-detect)
#' @param control list containing additional settings
#' 
#' 
#' @return est ec-draw object (List)
#' 
#' @examples
#' \dontrun{
#' icecream_est <- icecream %>% vd_est_vdm_screen(R=50000)
#' }
#' @export
vd_est_vdm_screen = function(vd,
                      R=100000, 
                      keep=10,
                      cores=NULL,
                      control=list(include_data=TRUE)){
  
  #check input data
  if(!vd_check_long(vd)) stop("Check data")
  
  #integer variables
  vd<-vd %>% mutate(task=as.integer(task),alt=as.integer(alt))
  
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
      1-(vd %>% 
        select(id,task,alt,x) %>% 
        bind_cols(dat$Af%>%as_tibble)  %>%
        mutate_if(is.double,function(col){vd$x*col})%>%
        group_by(id) %>% summarise_if(is.double,max) %>%
        arrange(as.numeric(id)) %>%
        select(-any_of(c('id','x'))) %>% as.matrix %>% sign)
  
  
  #Prior
  Bbar=matrix(rep(0,ncol(dat$AA)+3), ncol=ncol(dat$AA)+3)
  A=0.01*diag(1)
  nu=ncol(dat$AA)+9
  V=(ncol(dat$AA)+9)*diag(ncol(dat$AA)+3)
  Prior=list(Bbar=Bbar,A=A,nu=nu,V=V)
  
  #Run model
  out=
    loop_vdrs2_RWMH( dat$XX, 
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
  
  
  

  
  #Add data information
  out$A_names<-colnames(dat$AA)
  out$parnames=c(colnames(dat$AA),'sigma','gamma','E')
  if(!is.null(dat$attributes_levels)){
    out$attributes_levels=dat$attributes_levels
  }
  
  #Add model information
  out$ec_type="volumetric-conjunctive"
  out$error_specification="Normal"
  out$ec_type_short="VD-conj"
  out$Prior=Prior
  
  attributes(out)$Af<-attributes(dat)$Af
  attributes(out)$ec_data<-attributes(dat)$ec_data
  
  
  ec_model = list(model_name_full="volumetric-conjunctive-Normal",
                  model_name_short="VD-conj",
                  error_specification="Normal",
                  model_parnames=c('sigma','gamma','E'),
                  Prior=Prior)
  
  attributes(out)$ec_model=ec_model
  
  
  #Add training data
  if(control$include_data){
    out$dat=dat
  }
  
  return(out)
}


#' Estimate volumetric demand model with attribute-based conjunctive screening (w/ price)
#'
#' See https://dx.doi.org/10.2139/ssrn.2770025 for more details
#'
#'
#' @param vd volumetric demand data (long format)
#' @param R draws
#' @param keep thinning
#' @param cores no of CPU cores to use (default: auto-detect)
#' @param control list containing additional settings
#' 
#' 
#' @return est ec-draw object (List)
#' 
#' @export
vd_est_vdm_screenpr = function(vd,
                             R=100000, 
                             keep=10,
                             cores=NULL,
                             control=list(include_data=TRUE)){
  
  #check input data
  if(!vd_check_long(vd)) stop("Check data")
  
  #integer variables
  vd<-vd %>% mutate(task=as.integer(task),alt=as.integer(alt))
  
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
    1-(vd %>% 
         select(id,task,alt,x) %>% 
         bind_cols(dat$Af%>%as_tibble)  %>%
         mutate_if(is.double,function(col){vd$x*col})%>%
         group_by(id) %>% summarise_if(is.double,max) %>%
         arrange(as.numeric(id)) %>%
         select(-any_of(c('id','x'))) %>% as.matrix %>% sign)
  
  
  
  #Prior
  Bbar=matrix(rep(0,ncol(dat$AA)+3), ncol=ncol(dat$AA)+3)
  A=0.01*diag(1)
  nu=ncol(dat$AA)+9
  V=(ncol(dat$AA)+9)*diag(ncol(dat$AA)+3)
  Prior=list(Bbar=Bbar,A=A,nu=nu,V=V)
  
  #Run model
  out=
    loop_vdrspr_RWMH( dat$XX, 
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
  
  
  
  
  
  #Add data information
  out$A_names<-colnames(dat$AA)
  out$parnames=c(colnames(dat$AA),'sigma','gamma','E')
  if(!is.null(dat$attributes_levels)){
    out$attributes_levels=dat$attributes_levels
  }
  
  #Add model information
  out$ec_type="volumetric-conjunctive"
  out$error_specification="Normal"
  out$ec_type_short="VD-conj"
  out$Prior=Prior
  
  attributes(out)$Af<-attributes(dat)$Af
  attributes(out)$ec_data<-attributes(dat)$ec_data
  
  
  ec_model = list(model_name_full="volumetric-conjunctive-Normal",
                  model_name_short="VD-conj",
                  error_specification="Normal",
                  model_parnames=c('sigma','gamma','E'),
                  Prior=Prior)
  
  attributes(out)$ec_model=ec_model
  
  
  #Add training data
  if(control$include_data){
    out$dat=dat
  }
  
  return(out)
}







#' Estimate volumetric demand model accounting for set size variation (1st order)
#' 
#' For more details on the model: https://dx.doi.org/10.2139/ssrn.3418383
#' This model REQUIRES variation in choice-set size
#' 
#'
#' @param vd volumetric demand data (long format) with set size variation
#' @param R draws
#' @param keep thinning
#' @param cores no of CPU cores to use (default: auto-detect)
#' @param control list containing additional settings
#' 
#' @return est ec-draw object (List)
#' 
#' 
#' 
#' @export
vd_est_vdm_ss = function(vd,
                          R=100000, 
                          keep=10,
                          cores=NULL,
                          control=list(include_data=TRUE)){
  
  #check input data
  if(!vd_check_long(vd)) stop("Check data")
  
  #integer variables
  vd<-vd %>% mutate(task=as.integer(task),alt=as.integer(alt))
  
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
  Bbar=matrix(rep(0,ncol(dat$AA)+4), ncol=ncol(dat$AA)+4)
  Bbar[,ncol(dat$AA)+1]=-3
  
  A=0.01*diag(1)
  nu=ncol(dat$AA)+9
  V=(ncol(dat$AA)+9)*diag(ncol(dat$AA)+4)
  Prior=list(Bbar=Bbar,A=A,nu=nu,V=V)
  
  #Run model
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
  
  #Add data information
  out$A_names<-colnames(dat$AA)
  out$parnames=c(colnames(dat$AA),'xi','sigma','gamma','E')
  if(!is.null(dat$attributes_levels)){
    out$attributes_levels=dat$attributes_levels
  }
  
  #Add model information
  out$ec_type="volumetric-compensatory-setsize_linear"
  out$error_specification="EV1"
  out$ec_type_short="VD-comp-ssl"
  out$Prior=Prior
  
  attributes(out)$Af<-attributes(dat)$Af
  attributes(out)$ec_data<-attributes(dat)$ec_data
  
  
  ec_model = list(model_name_full="volumetric-compensatory-EV1-setsize_linear",
                  model_name_short="VD-comp-ssl",
                  error_specification="EV1",
                  model_parnames=c('xi','sigma','gamma','E'),
                  Prior=Prior)
  
  attributes(out)$ec_model=ec_model
  
  
  #Add training data
  if(control$include_data){
    out$dat=dat
  }
  
  return(out)
}





#' Estimate volumetric demand model accounting for set size variation (2nd order)
#' 
#' For more details on the model: https://dx.doi.org/10.2139/ssrn.3418383
#' This model REQUIRES variation in choice-set size
#' 
#'
#' @param vd volumetric demand data (long format)
#' @param R draws
#' @param keep thinning
#' @param cores no of CPU cores to use (default: auto-detect)
#' @param control list containing additional settings
#' 
#' @return est ec-draw object (List)
#' 
#' @export
vd_est_vdm_ssq = function(vd,
                         R=100000, 
                         keep=10,
                         cores=NULL,
                         control=list(include_data=TRUE)){
  
  #check input data
  if(!vd_check_long(vd)) stop("Check data")
  
  #integer variables
  vd<-vd %>% mutate(task=as.integer(task),alt=as.integer(alt))
  
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
  Bbar=matrix(rep(0,ncol(dat$AA)+5), ncol=ncol(dat$AA)+5)
  A=0.01*diag(1)
  nu=ncol(dat$AA)+9
  V=(ncol(dat$AA)+9)*diag(ncol(dat$AA)+5)
  Prior=list(Bbar=Bbar,A=A,nu=nu,V=V)
  
  #Run model
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
  
  #Add data information
  out$A_names<-colnames(dat$AA)
  out$parnames=c(colnames(dat$AA),'tau','xi','sigma','gamma','E')
  if(!is.null(dat$attributes_levels)){
    out$attributes_levels=dat$attributes_levels
  }
  
  #Add model information
  out$ec_type="volumetric-compensatory-setsize_linear"
  out$error_specification="EV1"
  out$ec_type_short="VD-comp-ssq"
  out$Prior=Prior
  
  attributes(out)$Af<-attributes(dat)$Af
  attributes(out)$ec_data<-attributes(dat)$ec_data
  
  
  ec_model = list(model_name_full="volumetric-compensatory-EV1-setsize_linear",
                  model_name_short="VD-comp-ssq",
                  error_specification="EV1",
                  model_parnames=c('tau','xi','sigma','gamma','E'),
                  Prior=Prior)
  
  attributes(out)$ec_model=ec_model
  
  
  #Add training data
  if(control$include_data){
    out$dat=dat
  }
  
  return(out)
}





# Log-Likelihoods for entire datasets -------------------------------------


#logll


#' Log-Likelihood for compensatory volumetric demand model
#' 
#' @usage vd_LL_vdm(draw, vd, fromdraw=1)
#' 
#' @param draw draws
#' @param vd vd object
#' @param fromdraw from which draw (excl. burnin)
#' 
#' @return Log-Likelihood
#' 
#' @export
vd_LL_vdm <- function(draw, vd, fromdraw=1){
  
  R=dim(draw$thetaDraw)[3]
  
  dat <- 
    vd %>% 
    vd_long_tidy %>% vd_prepare
  
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
  
  return(out) 
}

#' Log-Likelihood for volumetric demand model with set-size variation
#' 
#' 
#' @param draw draws
#' @param vd vd object
#' @param fromdraw from which draw (excl. burnin)
#' 
#' @return Log-Likelihood
#' 
#' @export
vd_LL_vdmss <- function(draw, vd, fromdraw=1){
  
  R=dim(draw$thetaDraw)[3]
  
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

#' Log-Likelihood for conjunctive-screening volumetric demand model
#' 
#' 
#' @param draw draws
#' @param vd vd object
#' @param fromdraw from which draw (excl. burnin)
#' 
#' @return Log-Likelihood
#' 
#' @export
vd_LL_vdm_screen <- function(draw, vd, fromdraw=1){
  
  R=dim(draw$thetaDraw)[3]
  
  dat <- 
    vd %>% 
    vd_long_tidy %>% vd_prepare
  
  dat$Af <- vd %>% vd_long_tidy %>% 
    attributes() %>% `[[`('Af') %>% as.matrix()
  
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
  
  return(out) 
}

#' Log-Likelihood for conjunctive-screening (with price) volumetric demand model
#' 
#' 
#' @param draw draws
#' @param vd vd object
#' @param fromdraw from which draw (excl. burnin)
#' 
#' @return Log-Likelihood
#' 
#' @export
vd_LL_vdm_screenpr <- function(draw, vd, fromdraw=1){
  
  R=dim(draw$thetaDraw)[3]
  
  dat <- 
    vd %>% 
    vd_long_tidy %>% vd_prepare
  
  dat$Af <- vd %>% vd_long_tidy %>% 
    attributes() %>% `[[`('Af') %>% as.matrix()
  
  out<- 
    vdsrprLLs(draw$thetaDraw[,,seq(fromdraw,R)],
             draw$tauDraw[,,seq(fromdraw,R)],
             draw$tau_pr_draw[,,seq(fromdraw,R)], 
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
  
  return(out) 
}





# Volumetric Predictions --------------------------------------------------


#' Match factor levels between two datasets
#' 
#' Makes sure the factor levels in `data_new` are aligned with `data_old`
#' This is helpful for demand simulations.
#'
#' @usage prep_newprediction(data_new,data_old)
#'
#' @param data_new New long-format choice data 
#' @param data_old Old long-format choice data
#' 
#' @return data_new (adjusted)
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




#' Demand Prediction (Volumetric demand, EV1 errors)
#'
#' Generating demand predictions for volumetric demand model. 
#' Reminder: there is no closed-form solution for demand, thus we need to integrate not only over the posterior distribution of parameters and the error distribution.
#' The function outputs a tibble containing id, task, alt, p, attributes, draws from the posterior of demand.
#' Eerror realisations can be pre-supplied to the `epsilon_not`. This helps create smooth demand curves or conduct optimization.
#'
#'
#' @param vd data
#' @param tidy apply echoice tidier
#' @param est ec-model draws 
#' @param epsilon_not (optional) error realizations
#' @param cores (optional) cores
#' 
#' @return Draws of expected demand
#' 
#' 
#' @seealso [prep_newprediction()] to match `vd`'s factor levels,
#'   [ec_gen_err_ev1()] for pre-generating error realizations and
#'   [vd_est_vdm()] for estimating the corresponding model
#' 
#' @export
vd_dem_vdm=function(vd,
                    tidy=T,
                    est,
                    epsilon_not=NULL,
                    cores=NULL){
  
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
  }else{
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
  attributes(vd)$ec_model   <- attributes(est)$ec_model
  
  return(vd)
}






#' Demand Prediction (Volumetric demand, Normal errors)
#'
#' Generating demand predictions for volumetric demand model. 
#' Reminder: there is no closed-form solution for demand, thus we need to integrate not only over the posterior distribution of parameters and the error distribution.
#' The function outputs a tibble containing id, task, alt, p, attributes, draws from the posterior of demand.
#' Eerror realisations can be pre-supplied to the `epsilon_not`. This helps create smooth demand curves or conduct optimization.
#'
#'
#' @param vd data
#' @param est ec-model draws 
#' @param epsilon_not (optional) error realizations
#' @param cores (optional) cores
#' 
#' @return Draws of expected demand
#' 
#' 
#' @seealso [prep_newprediction()] to match `vd`'s factor levels,
#'   [ec_gen_err_normal()] for pre-generating error realizations and
#'   [vd_est_vdmn()] for estimating the corresponding model
#' 
#' @export
vd_dem_vdmn=function(vd,
                    est,
                    epsilon_not=NULL,
                    cores=NULL){
  
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
    out=
      des_dem_vdmn( dat$PP,
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
  attributes(vd)$ec_model   <- attributes(est)$ec_model
  
  return(vd)
}




#' Demand Prediction (Volumetric demand, attribute-based screening, Normal)
#'
#' Generating demand predictions for volumetric demand model with attribute-based screening. 
#' Reminder: there is no closed-form solution for demand, thus we need to integrate not only over the posterior distribution of parameters and the error distribution.
#' The function outputs a tibble containing id, task, alt, p, attributes, draws from the posterior of demand.
#' Eerror realisations can be pre-supplied to the `epsilon_not`. This helps create smooth demand curves or conduct optimization.
#'
#' @usage vd_dem_vdmsr(vd, est, epsilon_not=NULL, cores=NULL)
#'
#' @param vd data
#' @param est ec-model draws 
#' @param epsilon_not (optional) error realizations
#' @param cores (optional) cores
#' 
#' @return Draws of expected demand
#' 
#' 
#' @seealso [prep_newprediction()] to match `vd`'s factor levels,
#'   [ec_gen_err_normal()] for pre-generating error realizations and
#'   [vd_est_vdm_screen()] for estimating the corresponding model
#' 
#' @export
vd_dem_vdmsr=function(vd,
                    est,
                    epsilon_not=NULL,
                    cores=NULL){
  
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
  }else{
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
  
  attributes(out)=NULL
  #add draws to data tibble
  vd=as_tibble(vd)
  vd$.demdraws<-map(out,drop)   
  
  #add attributes
  attributes(vd)$attr_names <- vd %>% colnames %>% setdiff(c("id","task","alt","x","p" )) %>% str_subset('^\\.', negate = TRUE)
  attributes(vd)$ec_model   <- attributes(est)$ec_model
  
  return(vd)
}



#' Demand Prediction (Volumetric demand, attribute-based screening including price, Normal)
#'
#' Generating demand predictions for volumetric demand model with attribute-based screening including price. 
#' Reminder: there is no closed-form solution for demand, thus we need to integrate not only over the posterior distribution of parameters and the error distribution.
#' The function outputs a tibble containing id, task, alt, p, attributes, draws from the posterior of demand.
#' Eerror realisations can be pre-supplied to the `epsilon_not`. This helps create smooth demand curves or conduct optimization.
#'
#' @usage vd_dem_vdmsrpr(vd, est, epsilon_not=NULL, cores=NULL)
#'
#' @param vd data
#' @param est ec-model draws 
#' @param epsilon_not (optional) error realizations
#' @param cores (optional) cores
#' 
#' @return Draws of expected demand
#' 
#' 
#' @seealso [prep_newprediction()] to match `vd`'s factor levels,
#'   [ec_gen_err_normal()] for pre-generating error realizations and
#'   [vd_est_vdm_screenpr()] for estimating the corresponding model
#' 
#' @export
vd_dem_vdmsrpr=function(vd,
                      est,
                      epsilon_not=NULL,
                      cores=NULL){
  
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
  }else{
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
  
  attributes(out)=NULL
  #add draws to data tibble
  vd=as_tibble(vd)
  vd$.demdraws<-map(out,drop)  
  
  #add attributes
  attributes(vd)$attr_names <- vd %>% colnames %>% setdiff(c("id","task","alt","x","p" )) %>% str_subset('^\\.', negate = TRUE)
  attributes(vd)$ec_model   <- attributes(est)$ec_model
  
  return(vd)
}



#' Demand Prediction (Volumetric demand, accounting for set-size variation (1st order), EV1 errors)
#'
#' Generating demand predictions for volumetric demand model with set-size adjustment.
#' Reminder: there is no closed-form solution for demand, thus we need to integrate not only over the posterior distribution of parameters and the error distribution.
#' The function outputs a tibble containing id, task, alt, p, attributes, draws from the posterior of demand.
#' Eerror realizations can be pre-supplied to the `epsilon_not`. This helps create smooth demand curves or conduct optimization.
#'
#' @usage vd_dem_vdmss(vd, est, epsilon_not=NULL, cores=NULL)
#'
#' @param vd data
#' @param est ec-model draws 
#' @param epsilon_not (optional) error realizations
#' @param cores (optional) cores
#' 
#' @return Draws of expected demand
#' 
#' 
#' @seealso [prep_newprediction()] to match `vd`'s factor levels,
#'   [ec_gen_err_ev1()] for pre-generating error realizations and
#'   [vd_est_vdm_ss()] for estimating the corresponding model
#' 
#' @export
vd_dem_vdmss=function(vd,
                    est,
                    epsilon_not=NULL,
                    cores=NULL){
  
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
  }else{
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
  
  attributes(out)=NULL
  #add draws to data tibble
  vd=as_tibble(vd)
  vd$.demdraws<-map(out,drop)  
  
  #add attributes
  attributes(vd)$attr_names <- vd %>% colnames %>% setdiff(c("id","task","alt","x","p" )) %>% str_subset('^\\.', negate = TRUE)
  attributes(vd)$ec_model   <- attributes(est)$ec_model
  
  return(vd)
}







#' Demand Prediction (Volumetric demand, accounting for set-size variation (2nd order), EV1 errors)
#'
#' Generating demand predictions for volumetric demand model with set-size adjustment.
#' Reminder: there is no closed-form solution for demand, thus we need to integrate not only over the posterior distribution of parameters and the error distribution.
#' The function outputs a tibble containing id, task, alt, p, attributes, draws from the posterior of demand.
#' Eerror realizations can be pre-supplied to the `epsilon_not`. This helps create smooth demand curves or conduct optimization.
#'
#'
#' @param vd data
#' @param est ec-model draws 
#' @param cores (optional) cores
#' 
#' @return Draws of expected demand
#' 
#' 
#' @seealso [prep_newprediction()] to match `vd`'s factor levels,
#'   [ec_gen_err_ev1()] for pre-generating error realizations and
#'   [vd_est_vdm_ssq()] for estimating the corresponding model
#' 
#' @export
vd_dem_vdmssq=function(vd,
                      est,
                      cores=NULL){
  
  #cores  
  if(is.null(cores)){
    cores=parallel::detectCores(logical=FALSE)
  }
  message(paste0("Using ",cores," cores"))
  
  
  #re-arrange data
  dat <- 
    vd %>% 
    vd_long_tidy %>% vd_prepare
  
  
  #demand sim
  out=
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
  
  attributes(out)=NULL
  #add draws to data tibble
  vd=as_tibble(vd)
  vd$.demdraws<-map(out,drop)  
  
  #add attributes
  attributes(vd)$attr_names <- vd %>% colnames %>% setdiff(c("id","task","alt","x","p" )) %>% str_subset('^\\.', negate = TRUE)
  attributes(vd)$ec_model   <- attributes(est)$ec_model
  
  return(vd)
}









# Discrete Demand ---------------------------------------------------------


#' Estimate discrete choice model (HMNL)
#'
#' @usage dd_est_hmnl(dd, R=100000, keep=10, cores=NULL, control=list(include_data=TRUE))
#'
#' @param dd discrete choice data (long format)
#' @param R draws
#' @param keep thinning
#' @param cores no of CPU cores to use (default: auto-detect)
#' @param control list containing additional settings
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
  dd<-dd %>% mutate(task=as.integer(task),alt=as.integer(alt))
  
  #Multicore settings
  if(is.null(cores)){
    cores=parallel::detectCores(logical=FALSE)
  }
  message(paste0("Using ",cores," cores"))
  
  #outside good present or not?
  outside_good_check <- dd %>% group_by(id,task) %>% summarise(x=sum(x)) %>% pull(x) %>% mean
  
  #re-arrange data
  if(outside_good_check==1){
    #no outside
    dat <- 
      dd %>% 
      vd_long_tidy %>% select(-int) %>% vd_prepare
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
  out$parnames=c(colnames(dat$AA),'ln_beta_p')
  if(!is.null(dat$attributes_levels)){
    out$attributes_levels=dat$attributes_levels
  }
  
  colnames(out$MUDraw)=c(colnames(dat$AA),'ln_beta_p')
  
  #Add model information
  out$ec_type="discrete-compensatory"
  out$error_specification="EV1"
  out$ec_type_short="hmnl-comp"
  out$Prior=Prior
  
  attributes(out)$Af<-attributes(dat)$Af
  attributes(out)$ec_data<-attributes(dat)$ec_data
  
  
  ec_model = list(model_name_full="discrete-compensatory",
                  model_name_short="hmnl-comp",
                  error_specification="EV1",
                  model_parnames=c('beta_p'),
                  Prior=Prior)
  
  attributes(out)$ec_model=ec_model
  
  
  #Add training data
  if(control$include_data){
    out$dat=dat
  }
  
  return(out)
}





#' Estimate discrete choice model (HMNL, attribute-based screening (not including price))
#'
#' @usage dd_est_hmnl_screen(dd, R=100000, keep=10, cores=NULL, control=list(include_data=TRUE))
#'
#' @param dd discrete choice data (long format)
#' @param R draws
#' @param keep thinning
#' @param cores no of CPU cores to use (default: auto-detect)
#' @param control list containing additional settings
#' 
#' @return est ec-draw object (List)
#' 
#' @seealso [dd_dem_sr()] to generate demand predictions based on this model
#' 
#' @export
dd_est_hmnl_screen = function(dd,
                     R=100000,
                     keep=10,
                     cores=NULL,
                     control=list(include_data=TRUE)){
  
  #check input data
  if(!dd_check_long(dd)) stop("Check data")
  
  #integer variables, double variables
  dd<-dd %>% 
    mutate(task=as.integer(task),
           alt=as.integer(alt),
           x=as.double(x))
  
  #sorting
  dd<-dd%>%arrange(as.numeric(id),task,alt)
  
  #Multicore settings
  if(is.null(cores)){
    cores=parallel::detectCores(logical=FALSE)
  }
  message(paste0("Using ",cores," cores"))
  
  #outside good present or not?
  outside_good_check <- dd %>% group_by(id,task) %>% summarise(x=sum(x)) %>% pull(x) %>% mean
  
  #re-arrange data
  if(outside_good_check==1){
    #no outside
    dat <- 
      dd %>% 
      vd_long_tidy %>% select(-int) %>% vd_prepare
  }else{
    #outside good
    dat <- 
      dd %>% 
      vd_long_tidy %>% vd_prepare
  }
  
  dat$Af <- dd %>% vd_long_tidy %>% 
    attributes() %>% `[[`('Af') %>% as.matrix()


  dat$tauconst=
    1-(
  dd %>% 
    select(id,task,alt,x) %>% 
    bind_cols(dat$Af%>%as_tibble)  %>%
    mutate_if(is.double,function(col){dd$x*col})%>%
    group_by(id) %>% summarise_if(is.double,max) %>%
    arrange(as.numeric(id)) %>%
    select(-any_of(c('id','x'))) %>% as.matrix)
  

  #Prior
  Bbar=matrix(rep(0,ncol(dat$AA)+1), ncol=ncol(dat$AA)+1)
  A=0.01*diag(1)
  nu=ncol(dat$AA)+9
  V=(ncol(dat$AA)+9)*diag(ncol(dat$AA)+1)
  Prior=list(Bbar=Bbar,A=A,nu=nu,V=V)
  
  
  #Run model
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
  
  

  #Add data information
  out$A_names<-colnames(dat$AA)
  out$Af_names<-colnames(dat$AAf)
  
  out$parnames=c(colnames(dat$AA),'ln_beta_p')
  
  if(!is.null(dat$attributes_levels)){
    out$attributes_levels=dat$attributes_levels
  }
  
  colnames(out$MUDraw)=c(colnames(dat$AA),'ln_beta_p')
  
  #Add model information
  out$ec_type="discrete-conjunctive"
  out$error_specification="EV1"
  out$ec_type_short="hmnl-conj"
  out$Prior=Prior
  
  attributes(out)$Af<-attributes(dat)$Af
  attributes(out)$ec_data<-attributes(dat)$ec_data
  
  
  ec_model = list(model_name_full="discrete-conjunctive",
                  model_name_short="hmnl-conj",
                  error_specification="EV1",
                  model_parnames=c('ln_beta_p'),
                  Prior=Prior)
  
  attributes(out)$ec_model=ec_model
  
  
  #Add training data
  if(control$include_data){
    out$dat=dat
  }
  
  return(out)
}




#' Estimate discrete choice model (HMNL, attribute-based screening (including price))
#'
#'
#' @param dd discrete choice data (long format)
#' @param R draws
#' @param keep thinning
#' @param cores no of CPU cores to use (default: auto-detect)
#' @param control list containing additional settings
#' 
#' @return est ec-draw object (List)
#' 
#' @seealso [dd_dem_sr()] to generate demand predictions based on this model
#' 
#' @export
dd_est_hmnl_screenpr = function(dd,
                                R=100000,
                                keep=10,
                                cores=NULL,
                                control=list(include_data=TRUE)){
  
  #check input data
  if(!dd_check_long(dd)) stop("Check data")
  
  #integer variables, double variables
  dd<-dd %>% 
    mutate(task=as.integer(task),
           alt=as.integer(alt),
           x=as.double(x))
  
  #sorting
  dd<-dd%>%arrange(as.numeric(id),task,alt)
  
  #Multicore settings
  if(is.null(cores)){
    cores=parallel::detectCores(logical=FALSE)
  }
  message(paste0("Using ",cores," cores"))
  
  #outside good present or not?
  outside_good_check <- dd %>% group_by(id,task) %>% summarise(x=sum(x)) %>% pull(x) %>% mean
  
  #re-arrange data
  if(outside_good_check==1){
    #no outside
    dat <- 
      dd %>% 
      vd_long_tidy %>% select(-int) %>% vd_prepare
  }else{
    #outside good
    dat <- 
      dd %>% 
      vd_long_tidy %>% vd_prepare
  }
  
  dat$Af <- dd %>% vd_long_tidy %>% 
    attributes() %>% `[[`('Af') %>% as.matrix()
  
  
  dat$tauconst=
    1-(
      dd %>% 
        select(id,task,alt,x) %>% 
        bind_cols(dat$Af%>%as_tibble)  %>%
        mutate_if(is.double,function(col){dd$x*col})%>%
        group_by(id) %>% summarise_if(is.double,max) %>%
        arrange(as.numeric(id)) %>%
        select(-any_of(c('id','x'))) %>% as.matrix)
  
  
  #Prior
  Bbar=matrix(rep(0,ncol(dat$AA)+1), ncol=ncol(dat$AA)+1)
  A=0.01*diag(1)
  nu=ncol(dat$AA)+9
  V=(ncol(dat$AA)+9)*diag(ncol(dat$AA)+1)
  Prior=list(Bbar=Bbar,A=A,nu=nu,V=V)
  
  
  #Run model
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
  
  
  
  #Add data information
  out$A_names<-colnames(dat$AA)
  out$Af_names<-colnames(dat$AAf)
  
  out$parnames=c(colnames(dat$AA),'ln_beta_p')
  
  if(!is.null(dat$attributes_levels)){
    out$attributes_levels=dat$attributes_levels
  }
  
  colnames(out$MUDraw)=c(colnames(dat$AA),'ln_beta_p')
  
  #Add model information
  out$ec_type="discrete-conjunctive-price"
  out$error_specification="EV1"
  out$ec_type_short="hmnl-conjpr"
  out$Prior=Prior
  
  attributes(out)$Af<-attributes(dat)$Af
  attributes(out)$ec_data<-attributes(dat)$ec_data
  
  
  ec_model = list(model_name_full="discrete-conjunctive-price",
                  model_name_short="hmnl-conjpr",
                  error_specification="EV1",
                  model_parnames=c('ln_beta_p'),
                  Prior=Prior)
  
  attributes(out)$ec_model=ec_model
  
  
  #Add training data
  if(control$include_data){
    out$dat=dat
  }
  
  return(out)
}



#logll
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


#logll
dd_LL_sr <- function(draw, dd, fromdraw=1){
  
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






# Discrete Demand Prediction ----------------------------------------------


#' Discrete Choice Predictions (HMNL)
#'
#' @usage dd_dem(dd, est, cores)
#'
#' @param dd data
#' @param est est
#' @param cores cores
#' 
#' @return Draws of expected choice
#' 
#' @seealso [dd_est_hmnl()] to generate demand predictions based on this model
#' 
#' @export
dd_dem=function(dd,
                est,
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
  
  
  attributes(out)=NULL
  #add draws to data tibble
  dd=as_tibble(dd)
  dd$.demdraws<-map(out,drop)  
  
  #add attributes
  attributes(dd)$attr_names <- dd %>% colnames %>% setdiff(c("id","task","alt","x","p" )) %>% str_subset('^\\.', negate = TRUE)
  attributes(dd)$ec_model   <- attributes(est)$ec_model
  
  return(dd)
}





#' Discrete Choice Probabilities (HMNL)
#'
#'
#' @param dd data
#' @param est est
#' @param cores cores
#' 
#' @return Draws of choice probabilities
#' 
#' @seealso [dd_est_hmnl()] to generate demand predictions based on this model
#' 
#' @export
dd_prob=function(dd,
                est,
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
  
  attributes(out)=NULL
  #add draws to data tibble
  dd=as_tibble(dd)
  dd$.demdraws<-map(out,drop)  
  
  #add attributes
  attributes(dd)$attr_names <- dd %>% colnames %>% setdiff(c("id","task","alt","x","p" )) %>% str_subset('^\\.', negate = TRUE)
  attributes(dd)$ec_model   <- attributes(est)$ec_model
  
  return(dd)
}




#' Discrete Choice Predictions (HMNL with attribute-based screening)
#'
#' @usage dd_dem_sr(dd, est, cores)
#'
#' @param dd data
#' @param est est
#' @param cores cores
#' 
#' @return Draws of expected choice
#' 
#' @seealso [dd_est_hmnl_screen()] to generate demand predictions based on this model
#' 
#' @export
dd_dem_sr=function(dd,
                   est,
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
  
  #screening-relevant data
  dat$Af <- dd %>% vd_long_tidy %>%attributes() %>% `[[`('Af') %>% as.matrix()
  
  #demand sim
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
  
  
  attributes(out)=NULL
  #add draws to data tibble
  dd=as_tibble(dd)
  dd$.demdraws<-map(out,drop)  
  
  #add attributes
  attributes(dd)$attr_names <- dd %>% colnames %>% setdiff(c("id","task","alt","x","p" )) %>% str_subset('^\\.', negate = TRUE)
  attributes(dd)$ec_model   <- attributes(est)$ec_model
  
  return(dd)
}


#' Discrete Choice Probabilities (HMNL with attribute-based screening)
#'
#' @usage dd_prob_sr(dd, est, cores)
#'
#' @param dd data
#' @param est est
#' @param cores cores
#' 
#' @return Draws of expected choice
#' 
#' @seealso [dd_est_hmnl_screen()] to generate demand predictions based on this model
#' 
#' @export
dd_prob_sr=function(dd,
                   est,
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
  
  #screening-relevant data
  dat$Af <- dd %>% vd_long_tidy %>%attributes() %>% `[[`('Af') %>% as.matrix()
  
  #demand sim
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
  
  
  attributes(out)=NULL
  #add draws to data tibble
  dd=as_tibble(dd)
  dd$.demdraws<-map(out,drop)  
  
  #add attributes
  attributes(dd)$attr_names <- dd %>% colnames %>% setdiff(c("id","task","alt","x","p" )) %>% str_subset('^\\.', negate = TRUE)
  attributes(dd)$ec_model   <- attributes(est)$ec_model
  
  return(dd)
}




#' Discrete Choice Predictions (HMNL with attribute-based screening w/price)
#'
#' @usage dd_dem_srpr(dd, est, cores)
#'
#' @param dd data
#' @param est est
#' @param cores cores
#' 
#' @return Draws of expected choice
#' 
#' @seealso [dd_est_hmnl_screen()] to generate demand predictions based on this model
#' 
#' @export
dd_dem_srpr=function(dd,
                   est,
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
  
  #screening-relevant data
  dat$Af <- dd %>% vd_long_tidy %>%attributes() %>% `[[`('Af') %>% as.matrix()
  
  #demand sim
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
  
  
  attributes(out)=NULL
  #add draws to data tibble
  dd=as_tibble(dd)
  dd$.demdraws<-map(out,drop)   
  
  #add attributes
  attributes(dd)$attr_names <- dd %>% colnames %>% setdiff(c("id","task","alt","x","p" )) %>% str_subset('^\\.', negate = TRUE)
  attributes(dd)$ec_model   <- attributes(est)$ec_model
  
  return(dd)
}


#' Discrete Choice Probabilities (HMNL with attribute-based screening w/ price)
#'
#' @usage dd_prob_srpr(dd, est, cores)
#'
#' @param dd data
#' @param est est
#' @param cores cores
#' 
#' @return Draws of expected choice
#' 
#' @seealso [dd_est_hmnl_screen()] to generate demand predictions based on this model
#' 
#' @export
dd_prob_srpr=function(dd,
                    est,
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
  
  #screening-relevant data
  dat$Af <- dd %>% vd_long_tidy %>%attributes() %>% `[[`('Af') %>% as.matrix()
  
  #demand sim
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
  
  
  attributes(out)=NULL
  #add draws to data tibble
  dd=as_tibble(dd)
  dd$.demdraws<-map(out,drop)  
  
  #add attributes
  attributes(dd)$attr_names <- vd %>% colnames %>% setdiff(c("id","task","alt","x","p" )) %>% str_subset('^\\.', negate = TRUE)
  attributes(dd)$ec_model   <- attributes(est)$ec_model
  
  return(dd)
}



# working with demand -----------------------------------------------------






#' Add product id to demand draws
#' 
#' This adds a unique product identifier to  demand draw objects.
#'
#' @usage vd_add_prodid(de)
#'
#' @param de demand draws
#' 
#' @return est
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
#' 
#' @return Aggregated demand predictions
#' 
#' @export
ec_dem_aggregate = function(de, groupby){
  return(
  de %>% 
    group_by(!!!syms(groupby)) %>% 
    summarise(.demdraws=list(reduce(.demdraws,`+`))) )
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
#' 
#' @return Summary of demand predictions
#' 
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




#' Summarize posterior draws of screening
#'
#' Adds summaries of posterior draws of demand to tibble.
#' (using the new demand draw format)
#'
#'
#' @param sc tibble containing screening draws in .screendraws  
#' @param quantiles Quantiles for Credibility Intervals (default: 90% interval)
#' 
#' @return Summary of screening probabilities
#' 
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


#' Summarize posterior draws of demand (volumetric models only)
#'
#' Adds summaries of posterior draws of demand to tibble.
#' (using the new demand draw format)
#'
#'
#' @param de demand draws
#' @param quantiles Quantiles for Credibility Intervals (default: 90% interval)
#' 
#' @return Summary of demand predictions
#' 
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










#' Evaluate (hold-out) demand predictions
#'
#' This function obtains proper posterior fit statistics. 
#' It computes the difference between true demand and each draw from the demand posterior. 
#' Then, fit statistics are obtained.
#'
#'
#' @param de demand draws (output from vd_dem_x function)
#' 
#' @return Predictive fit statistics (MAE, MSE, RAE, bias, hit-probability)
#' 
#' @export
ec_dem_eval = function(de){
  `%!in%` = Negate(`%in%`)
  
  is_this_discrete=attributes(de)$ec_model$model_name_full %>% str_detect('discrete')
  
  if(is_this_discrete){
    out <- de %>%
      mutate(.MSE=map_dbl(map2(.demdraws,x,function(draws,x)(draws-x)^2 ),mean),
             .MAE=map_dbl(map2(.demdraws,x,function(draws,x)abs(draws-x)),mean),
             .bias=map_dbl(map2(.demdraws,x,function(draws,x)(draws-x)),mean),
             .R=abs(x-mean(x)),
             .hp=map_dbl(map2(.demdraws,x,function(draws,x)(draws==x) ),mean)) %>%
      summarise(MSE=mean(.MSE),
                MAE=mean(.MAE),
                bias=mean(.bias),
                RAE=MAE/mean(.R),
                hitrate=mean(.hp))
  }else{
    out <- de %>%
      mutate(.MSE=map_dbl(map2(.demdraws,x,function(draws,x)(draws-x)^2 ),mean),
             .MAE=map_dbl(map2(.demdraws,x,function(draws,x)abs(draws-x)),mean),
             .bias=map_dbl(map2(.demdraws,x,function(draws,x)(draws-x)),mean),
             .R=abs(x-mean(x))) %>%
      summarise(MSE=mean(.MSE),
                MAE=mean(.MAE),
                bias=mean(.bias),
                RAE=MAE/mean(.R))
  }
  
  return(out)
}


#' Evaluate (hold-out) demand predictions based on posterior mean predictions
#'
#' This function obtains proper posterior fit statistics. 
#' 
#'
#' @param de demand draws (output from vd_dem_x function)
#' 
#' @return Predictive fit statistics (MAE, MSE, RAE, bias, hit-probability)
#' 
#' @export
ec_dem_eval_pm=
  function(de){
    `%!in%` = Negate(`%in%`)
    
    is_this_discrete=attributes(de)$ec_model$model_name_full %>% str_detect('discrete')
    
    if(is_this_discrete){
      out <- 
        de %>% ec_dem_summarise() %>%
        mutate(.MSE=(`E(demand)`-x)^2,
               .MAE=abs(`E(demand)`-x),
               .bias=(`E(demand)`-x),
               .R=abs(x-mean(x)),
               .hr=as.numeric(`E(demand)`==x),
               .hp=ifelse(x==1,`E(demand)`,NA)) %>%
        group_by(id,task) %>% 
        mutate(.outside=1-sum(x),.outsideprob=1-sum(`E(demand)`)) %>%
        ungroup() %>%
        mutate(.hpoutside=ifelse(.outside==1,.outsideprob,NA)) %>%
        mutate(.hpall=ifelse(!is.na(.hp),.hp,.hpoutside))%>%
        summarise(MSE=mean(.MSE),
                  MAE=mean(.MAE),
                  bias=mean(.bias),
                  RAE=MAE/mean(.R),
                  hitrate=mean(.hr),
                  hitprob=mean(.hpall,na.rm=T))
    }else{
      out <- 
        de %>% ec_dem_summarise() %>%
        mutate(.MSE=(`E(demand)`-x)^2,
               .MAE=abs(`E(demand)`-x),
               .bias=(`E(demand)`-x),
               .R=abs(x-mean(x)))%>%
        summarise(MSE=mean(.MSE),
                  MAE=mean(.MAE),
                  bias=mean(.bias),
                  RAE=MAE/mean(.R))
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
#' 
#' @return List containing aggregate demand quantities for each scenario defined by `rel_pricerange`
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
#' 
#' @return List containing aggregate demand quantities for each scenario defined by `rel_pricerange`
#' 
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
#' 
#' @return List containing aggregate demand quantities for each scenario defined by `rel_pricerange`
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
        add_column(.isfocal=focal_product) %>%
        filter(.isfocal)%>% select(-.isfocal) 
      
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
                res %>% set_names('condem')) %>%
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
        add_column(.isfocal=focal_product) %>%
        filter(.isfocal)%>% select(-.isfocal) 
      
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
                    res %>% set_names('condem')) %>%
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
#' 
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
#' 
#' 
#' @export
ec_screenprob_sr=function(xd,
                          est,
                          cores=NULL){
  
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
  
  attributes(out)=NULL
  
  #add draws to data tibble
  xd=as_tibble(xd)
  xd$.screendraws<-map(out,drop)
  
  #add attributes
  attributes(xd)$attr_names <- xd %>% colnames %>% setdiff(c("id","task","alt","x","p" )) %>% str_subset('^\\.', negate = TRUE)
  attributes(xd)$ec_model   <- attributes(est)$ec_model
  
  return(xd)
}



#' Screening probabilities of choice alternatives (w/ price)
#'
#' Obtain draws of screening probabilities of choiec alternatives
#'
#'
#' @param xd data
#' @param est ec-model draws 
#' @param cores (optional) cores
#' 
#' @return Draws of screening probabilities of choice alternatives
#' 
#' 
#' @export
ec_screenprprob_sr=function(xd,
                            est,
                            cores=NULL){
  
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
#' @usage ec_draws_MU(draws)
#'
#' @param draws ec-draws, output from any ec-model
#' 
#' @return  draws
#' 
#' @seealso [ec_draws_screen()] to obtain screening parameter draws (where applicable),
#' [ec_trace_MU()] to generate a traceplot of MU_theta draws
#' 
#' @export
ec_draws_MU <- function(draws){
  
  draws$MUDraw %>% as_tibble %>%
    set_names(draws$parnames) %>% 
    rowid_to_column(var = 'draw') %>%
    pivot_longer(cols = -any_of('draw'), 
                 names_to = 'attribute_level') 
  
}

#' Obtain Screening probability draws
#'
#'
#' @usage ec_draws_screen(draws)
#'
#' @param draws ec-draws, output from any ec-model with screening
#' 
#' @return  draws
#' 
#' @seealso [ec_draws_MU()] to obtain MU_theta draws,
#' [ec_trace_screen()] to generate a traceplot of screening draws
#' 
#' @export
ec_draws_screen <- function(draws){
  
  draws$deltaDraw %>% as_tibble %>%
    set_names(colnames(attributes(draws)$Af)) %>% 
    rowid_to_column(var = 'draw') %>%
    pivot_longer(cols = -any_of('draw'), 
                 names_to = 'attribute_level') 
  
}

#' Generate MU_theta traceplot
#'
#'
#'
#' @param draws ec-draws, output from any ec-model
#' @param burnin burn-in to remove
#' 
#' @seealso [ec_boxplot_MU()] to obtain boxplot
#' 
#' @export
ec_trace_MU <- function(draws, burnin=100){
  
  draws$MUDraw %>% data.frame %>% as_tibble %>%
    set_names(draws$parnames) %>% 
    rowid_to_column(var = 'draw') %>%
    pivot_longer(cols = -any_of('draw'), 
                 names_to = 'par') %>%
    left_join(attributes(draws)$ec_data$attributes %>%
                select(attr_level,attribute,lvl,reference_lvl), 
              by=c('par'='attr_level')) %>%
    ggplot(aes(x=draw, y=value)) +
    geom_line() +guides(color='none')+facet_wrap(~par,scale='free_y')
}


#' Generate Screening probability traceplot
#'
#'
#'
#' @param draws ec-draws, output from any ec-model with screening
#' @param burnin burn-in to remove
#' 
#' @return  draws
#' 
#' @seealso [ec_draws_MU()] to obtain MU_theta draws,
#' [ec_boxplot_screen()] to generate boxplot
#' 
#' @export
ec_trace_screen <- function(draws, burnin=100){
  
  draws$deltaDraw %>% as_tibble %>%
    set_names(colnames(attributes(draws)$Af)) %>% 
    rowid_to_column(var = 'draw') %>%
    pivot_longer(cols = -any_of('draw'), 
                 names_to = 'attribute_level') %>%
    left_join(attributes(draws)$ec_data$attributes %>%
                select(attr_level,attribute,lvl,reference_lvl), 
              by=c('attribute_level'='attr_level')) %>%
    filter(draw>burnin) %>%
    ggplot(aes(x=draw, y=value)) +
    geom_line() +guides(color='none')+facet_wrap(~attribute_level,scale='free_y')
  
}


#' Generate MU_theta boxplot
#'
#'
#'
#' @param draws ec-draws, output from any ec-model
#' @param burnin burn-in to remove
#' 
#' @seealso [ec_trace_MU()] to obtain traceplot
#' 
#' @export
ec_boxplot_MU <- function(draws, burnin=100){
  
  draws$MUDraw %>% as_tibble %>%
    set_names(draws$parnames) %>% 
    rowid_to_column(var = 'draw') %>%
    pivot_longer(cols = -any_of('draw'), 
                 names_to = 'par') %>%
    left_join(attributes(draws)$ec_data$attributes %>%
                select(attr_level,attribute,lvl,reference_lvl), 
              by=c('par'='attr_level')) %>%
    mutate(par=str_replace_all(par,paste0(attribute,":"),"")) %>%
    select(draw,par,value,attribute) %>%
    filter(draw>burnin) %>%
    ggplot(aes(x=par,y=value)) + geom_boxplot() + coord_flip() +
    facet_wrap(~attribute, scales='free')
  
}


#' Generate Screening probability boxplot
#'
#'
#'
#' @param draws ec-draws, output from any ec-model with screening
#' @param burnin burn-in to remove
#' 
#' @return  draws
#' 
#' @seealso [ec_draws_MU()] to obtain MU_theta draws,
#' [ec_trace_screen()] to generate traceplot
#' 
#' @export
ec_boxplot_screen <- function(draws, burnin=100){
  
  draws$deltaDraw %>% as_tibble %>%
    set_names(colnames(attributes(draws)$Af)) %>% 
    rowid_to_column(var = 'draw') %>%
    pivot_longer(cols = -any_of('draw'), 
                 names_to = 'attribute_level') %>%
    left_join(attributes(draws)$ec_data$attributes %>%
                select(attr_level,attribute,lvl,reference_lvl), 
              by=c('attribute_level'='attr_level')) %>%
    filter(draw>burnin) %>%
    ggplot(aes(x=lvl,y=value)) + geom_boxplot() + coord_flip() +
    facet_wrap(~attribute, scales='free_y')
  
}




#' Thin echoice-vd draw objects
#'
#'
#' @param M is an echoice draw object (list)
#' @param burnin_perc how much burn-in to remove
#' @param total_draws how many draws to keep after thinning
#'
#' @return thinned echoice draw object (list)
#'
#' @export
vd_thin_draw=function(M, 
                      burnin_perc=.5, 
                      total_draws=NULL){
  
  R = dim(M$thetaDraw)[3]
  
  #burnin
  first_draw=floor(burnin_perc*R)
  usableR   =length(first_draw:R)
  message("Draws post burnin: ",usableR)
  
  
  if(usableR<total_draws){
    message("More desired draws than original draws")
    total_draws=usableR
    keepdraws=sample(first_draw:R, size = total_draws, replace = TRUE)
  }else{
    all_draws=first_draw:R
    
    #try to thin as much as possible
    chosen_draws = unique(round(seq.int(first_draw, R, length.out=total_draws),0))
    
    #if missing some draws, randomly grab them from remaining draws
    if(length(chosen_draws)<total_draws){
      message("draws in quence: ",length(chosen_draws), " adding more")
      
      added_draws=sample(all_draws[all_draws%in%chosen_draws],total_draws-length(chosen_draws))
      keepdraws=c(chosen_draws,added_draws)
    } else{
      keepdraws=chosen_draws
    }
    
  }
  
  if(is.null(M$tauDraw)){
    M$thetaDraw=  M$thetaDraw[,,keepdraws]
    M$MUDraw=     M$MUDraw[keepdraws,]
    M$SIGMADraw=  M$SIGMADraw[,,keepdraws]
    M$loglike=    M$loglike[keepdraws,,drop=F]
    M$logpost=    M$logpost[keepdraws,,drop=F]
    M$reject=     M$reject[keepdraws,,drop=F]
    
  }else{
    M$thetaDraw=  M$thetaDraw[,,keepdraws]
    M$MUDraw=     M$MUDraw[keepdraws,]
    M$SIGMADraw=  M$SIGMADraw[,,keepdraws]
    M$tauDraw=    M$tauDraw[,,keepdraws]
    M$deltaDraw=  M$deltaDraw[keepdraws,]
    M$loglike=    M$loglike[keepdraws,,drop=F]
    M$logpost=    M$logpost[keepdraws,,drop=F]
    M$reject=     M$reject[keepdraws,,drop=F]
  }
  
  return(M)
}
